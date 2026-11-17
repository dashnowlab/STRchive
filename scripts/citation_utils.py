"""
Citation lookups used by run-manubot.py instead of Manubot.

- PubMed, GeneReviews and OMIM IDs from NCBI E-utilities. GeneReviews chapters are
  PubMed book records tagged with their NBK accession, so they are parsed from the same
  PubMed XML. OMIM titles come from NCBI's OMIM database.
- MONDO IDs from the EBI Ontology Lookup Service (OLS4) API.

Optional environment variables:
- NCBI_API_KEY: raises the NCBI rate limit from 3 to 10 requests/second
- NCBI_EMAIL: contact address NCBI asks automated clients to provide
"""

import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BATCH_SIZE = 200  # IDs per efetch request
SEARCH_BATCH_SIZE = 50  # NBK accessions per esearch query
MAX_TRIES = 5

OLS_API = "https://www.ebi.ac.uk/ols4/api/ontologies/mondo/terms?obo_id={obo_id}"
OLS_LINK = "https://www.ebi.ac.uk/ols4/ontologies/mondo/entities/{iri}"

PUBMED_LINK = "https://www.ncbi.nlm.nih.gov/pubmed/{id}"
GENEREVIEWS_LINK = "https://www.ncbi.nlm.nih.gov/books/{id}"
OMIM_LINK = "https://omim.org/entry/{id}"

MONTHS = {m: i for i, m in enumerate(
    ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"], 1)}

_last_request = 0.0


def _fetch(url, data=None, label="request", min_interval=0.0, check=None):
    """
    GET (or POST if data is given) a URL, retrying transient failures with exponential
    backoff. min_interval: minimum seconds between requests (rate limit).
    check: optional function(body) that raises RuntimeError if the body is an error page.
    Returns the response body.
    """
    global _last_request
    for attempt in range(1, MAX_TRIES + 1):
        wait = _last_request + min_interval - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        _last_request = time.monotonic()
        try:
            with urllib.request.urlopen(url, data=data, timeout=60) as resp:
                body = resp.read()
            if check:
                check(body)
            return body
        except (urllib.error.URLError, TimeoutError, ConnectionError, RuntimeError) as e:
            retryable = not isinstance(e, urllib.error.HTTPError) or e.code == 429 or e.code >= 500
            if not retryable or attempt == MAX_TRIES:
                raise
            delay = 2 ** attempt
            sys.stderr.write(f"WARNING: {label} failed ({e}); retry {attempt}/{MAX_TRIES - 1} in {delay}s\n")
            time.sleep(delay)


def _eutils_error(body):
    # E-utilities sometimes returns HTTP 200 with an error document
    if b"<ERROR>" in body[:500] or b'"error"' in body[:200]:
        raise RuntimeError(body[:200].decode(errors="replace"))


def _request(endpoint, params):
    """POST to an NCBI E-utilities endpoint, respecting NCBI's rate limit"""
    api_key = os.environ.get("NCBI_API_KEY")
    email = os.environ.get("NCBI_EMAIL")
    params = dict(params, tool="STRchive")
    if api_key:
        params["api_key"] = api_key
    if email:
        params["email"] = email
    return _fetch(EUTILS + endpoint, data=urllib.parse.urlencode(params).encode(),
                  label=f"NCBI {endpoint}", min_interval=0.11 if api_key else 0.34,
                  check=_eutils_error)


def _text(elem):
    """Full text of an element including nested markup (<i>, <sup>, ...), whitespace-normalised"""
    if elem is None:
        return ""
    return " ".join("".join(elem.itertext()).split())


def _date(elem):
    """YYYY-MM-DD from a PubMed date element (missing month/day default to 01)"""
    if elem is None:
        return ""
    year = elem.findtext("Year")
    if not year:
        medline = elem.findtext("MedlineDate") or ""  # e.g. "2019 Jul-Aug"
        year = medline[:4] if medline[:4].isdigit() else ""
    if not year:
        return ""
    month = (elem.findtext("Month") or "1").strip()
    month = MONTHS.get(month[:3].lower(), month)
    try:
        month = int(month)
    except ValueError:
        month = 1
    try:
        day = int(elem.findtext("Day") or 1)
    except ValueError:
        day = 1
    return f"{int(year):04}-{month:02}-{day:02}"


def _authors(author_list):
    authors = []
    if author_list is None:
        return authors
    for a in author_list.findall("Author"):
        given = (a.findtext("ForeName") or a.findtext("Initials") or "").strip()
        family = (a.findtext("LastName") or "").strip()
        if not family:
            family = _text(a.find("CollectiveName"))
        if given or family:
            authors.append([given, family])
    return authors


def _abstract(elem):
    if elem is None:
        return ""
    parts = []
    for t in elem.findall("AbstractText"):
        label = t.get("Label")
        txt = _text(t)
        parts.append(f"{label}: {txt}" if label else txt)
    return "\n".join(p for p in parts if p)


def _parse_article(art):
    """Parse a <PubmedArticle> into a citation dict (without id/link)"""
    mc = art.find("MedlineCitation")
    article = mc.find("Article")
    journal = article.find("Journal")
    doi = ""
    for aid in art.findall("PubmedData/ArticleIdList/ArticleId"):
        if aid.get("IdType") == "doi":
            doi = (aid.text or "").strip()
    if not doi:
        for eloc in article.findall("ELocationID"):
            if eloc.get("EIdType") == "doi":
                doi = (eloc.text or "").strip()
    doi = doi.lower()  # DOIs are case-insensitive; lowercase to match previous Manubot output
    date = _date(article.find("ArticleDate")) or _date(journal.find("JournalIssue/PubDate"))
    return mc.findtext("PMID"), {
        "title": _text(article.find("ArticleTitle")) or _text(article.find("VernacularTitle")),
        "type": "article-journal",
        "doi": doi,
        "authors": _authors(article.find("AuthorList")),
        "publisher": _text(journal.find("Title")) or _text(journal.find("ISOAbbreviation")),
        "issn": (journal.findtext("ISSN") or "").strip(),
        "date": date,
        "abstract": _abstract(article.find("Abstract")),
    }


def _parse_book(art):
    """Parse a <PubmedBookArticle> (e.g. a GeneReviews chapter) into a citation dict"""
    doc = art.find("BookDocument")
    book = doc.find("Book")
    nbk = ""
    for aid in doc.findall("ArticleIdList/ArticleId"):
        if aid.get("IdType") == "bookaccession":
            nbk = (aid.text or "").strip()
    # chapter authors, falling back to book editors
    authors = _authors(doc.find("AuthorList")) or _authors(book.find("AuthorList"))
    date = (_date(doc.find("DateRevised")) or _date(doc.find("ContributionDate"))
            or _date(book.find("PubDate")))
    return doc.findtext("PMID"), nbk, {
        "title": _text(doc.find("ArticleTitle")) or _text(book.find("BookTitle")),
        "type": "chapter",
        "doi": "",
        "authors": authors,
        "publisher": _text(book.find("BookTitle")) or _text(book.find("Publisher/PublisherName")),
        "issn": "",
        "date": date,
        "abstract": _abstract(doc.find("Abstract")),
    }


def _efetch(pmids):
    """Fetch PubMed XML for PMIDs in batches. Returns (articles_by_pmid, books_by_nbk)."""
    articles, books = {}, {}
    pmids = list(dict.fromkeys(pmids))
    for i in range(0, len(pmids), BATCH_SIZE):
        batch = pmids[i:i + BATCH_SIZE]
        sys.stderr.write(f"Fetching PubMed records {i + 1}-{i + len(batch)} of {len(pmids)}\n")
        try:
            root = ET.fromstring(_request("efetch.fcgi", {
                "db": "pubmed", "id": ",".join(batch), "retmode": "xml"}))
        except Exception as e:
            sys.stderr.write(f"WARNING: PubMed fetch failed for batch starting {batch[0]}: {e}\n")
            continue
        for art in root.findall("PubmedArticle"):
            try:
                pmid, cite = _parse_article(art)
                articles[pmid] = cite
            except Exception as e:
                sys.stderr.write(f"WARNING: couldn't parse PubMed article: {e}\n")
        for art in root.findall("PubmedBookArticle"):
            try:
                pmid, nbk, cite = _parse_book(art)
                articles[pmid] = cite
                if nbk:
                    books[nbk] = cite
            except Exception as e:
                sys.stderr.write(f"WARNING: couldn't parse PubMed book article: {e}\n")
    return articles, books


def _nbk_to_pmids(nbks):
    """Find the PubMed records for GeneReviews (Bookshelf) NBK accessions"""
    pmids = []
    for i in range(0, len(nbks), SEARCH_BATCH_SIZE):
        batch = nbks[i:i + SEARCH_BATCH_SIZE]
        term = " OR ".join(f"{n}[All Fields]" for n in batch)
        try:
            res = json.loads(_request("esearch.fcgi", {
                "db": "pubmed", "term": term, "retmode": "json", "retmax": 10 * len(batch)}))
            pmids += res["esearchresult"]["idlist"]
        except Exception as e:
            sys.stderr.write(f"WARNING: GeneReviews search failed for batch starting {batch[0]}: {e}\n")
    return pmids


def _bookshelf(nbks):
    """
    Fallback for NBK accessions not indexed in PubMed (e.g. GeneReviews appendices,
    other Bookshelf books): look them up in the Bookshelf (books) database.
    Returns {nbk: citation dict}.
    """
    books = {}
    for nbk in nbks:
        try:
            res = json.loads(_request("esearch.fcgi", {
                "db": "books", "term": f"{nbk}[accn]", "retmode": "json", "retmax": 100}))
            uids = res["esearchresult"]["idlist"]
            if not uids:
                continue
            summ = json.loads(_request("esummary.fcgi", {
                "db": "books", "id": ",".join(uids), "retmode": "json"}))["result"]
        except Exception as e:
            sys.stderr.write(f"WARNING: Bookshelf lookup failed for {nbk}: {e}\n")
            continue
        # the record for the document itself (not a table/figure within it)
        doc = next((summ[u] for u in summ.get("uids", []) if summ[u].get("rid") == nbk), None)
        if doc is None:
            continue
        book_title = ""
        try:
            info = ET.fromstring(doc.get("bookinfo") or "<Info/>")
            src = info.find("Path/Parent[@role='source']/Title")
            book_title = _text(src)
        except ET.ParseError:
            pass
        title = _text(ET.fromstring(f"<t>{doc.get('title', '')}</t>")) if doc.get("title") else ""
        pubdate = (doc.get("pubdate") or "")[:10].replace("/", "-")
        books[nbk] = {
            "title": title,
            "type": "chapter",
            "doi": "",
            "authors": [],
            "publisher": book_title,
            "issn": "",
            "date": pubdate,
            "abstract": "",
        }
    return books


def _omim(mims):
    """
    Look up OMIM entry titles via NCBI's OMIM database (omim.org blocks automated requests).
    Returns {mim: citation dict}.
    """
    entries = {}
    for i in range(0, len(mims), BATCH_SIZE):
        batch = mims[i:i + BATCH_SIZE]
        try:
            res = json.loads(_request("esummary.fcgi", {
                "db": "omim", "id": ",".join(batch), "retmode": "json"}))["result"]
        except Exception as e:
            sys.stderr.write(f"WARNING: OMIM lookup failed for batch starting {batch[0]}: {e}\n")
            continue
        for uid in res.get("uids", []):
            title = " ".join((res[uid].get("title") or "").split())
            if not title:
                continue
            entries[uid] = {
                "title": title,
                "type": "entry",
                "doi": "",
                "authors": [],
                "publisher": "OMIM",
                "issn": "",
                "date": "",
                "abstract": "",
            }
    return entries


def _note(source):
    return f"Generated from NCBI E-utilities ({source})"


def _failed(_id, link, reason):
    return {"id": _id, "manubot_success": False, "link": link,
            "note": f"WARNING: NCBI could not generate citation: {reason}"}


def cite_with_ncbi(ids, bad_title_reason=None):
    """
    Generate citations for pmid:*, genereviews:* and omim:* ids.
    Returns {id: citation} in the same format as run-manubot.py.
    bad_title_reason: optional function(title) -> reason or None; failing titles are marked unsuccessful.
    """
    pmid_ids = [i for i in ids if i.startswith("pmid:")]
    gr_ids = [i for i in ids if i.startswith("genereviews:")]
    omim_ids = [i for i in ids if i.startswith("omim:")]
    sys.stderr.write(f"Processing {len(pmid_ids)} pmid, {len(gr_ids)} genereviews and {len(omim_ids)} omim citations via NCBI E-utilities\n")

    nbks = [i.split(":", 1)[1] for i in gr_ids]
    pmids = [i.split(":", 1)[1] for i in pmid_ids] + _nbk_to_pmids(nbks)
    articles, books = _efetch(pmids)
    missing_nbks = [n for n in nbks if n not in books]
    if missing_nbks:
        sys.stderr.write(f"Looking up {len(missing_nbks)} NBK accessions not in PubMed via Bookshelf\n")
        books.update(_bookshelf(missing_nbks))

    results = {}
    for _id in pmid_ids:
        pmid = _id.split(":", 1)[1]
        link = PUBMED_LINK.format(id=pmid)
        cite = articles.get(pmid)
        if cite is None:
            results[_id] = _failed(_id, link, "PMID not returned by PubMed")
            continue
        results[_id] = {"id": _id, "manubot_success": True, "link": link, **cite,
                        "language": "en", "note": _note(f"pubmed:{pmid}")}
    for _id in gr_ids:
        nbk = _id.split(":", 1)[1]
        link = GENEREVIEWS_LINK.format(id=nbk)
        cite = books.get(nbk)
        if cite is None:
            results[_id] = _failed(_id, link, "NBK accession not found in PubMed")
            continue
        results[_id] = {"id": _id, "manubot_success": True, "link": link, **cite,
                        "language": "en", "note": _note(f"bookaccession:{nbk}")}

    if omim_ids:
        omim = _omim([i.split(":", 1)[1] for i in omim_ids])
        for _id in omim_ids:
            mim = _id.split(":", 1)[1]
            link = OMIM_LINK.format(id=mim)
            cite = omim.get(mim)
            if cite is None:
                results[_id] = _failed(_id, link, "MIM number not found in NCBI OMIM")
                continue
            results[_id] = {"id": _id, "manubot_success": True, "link": link, **cite,
                            "language": "en", "note": _note(f"omim:{mim}")}

    if bad_title_reason:
        for _id, cite in results.items():
            reason = cite["manubot_success"] and bad_title_reason(cite.get("title"))
            if reason:
                results[_id] = _failed(_id, cite["link"], reason)
    return results


def cite_mondo(ids):
    """
    Generate citations for mondo:* ids (e.g. mondo:0025193).
    Returns {id: citation} in the same format as run-manubot.py.
    """
    sys.stderr.write(f"Processing {len(ids)} mondo citations via OLS\n")
    results = {}
    for _id in ids:
        obo_id = "MONDO:" + _id.split(":", 1)[1]
        iri = f"http://purl.obolibrary.org/obo/{obo_id.replace(':', '_')}"
        # same double-encoded link format Manubot produced
        link = OLS_LINK.format(iri=urllib.parse.quote(urllib.parse.quote(iri, safe=""), safe=""))
        try:
            terms = json.loads(_fetch(OLS_API.format(obo_id=obo_id), label="OLS")).get("_embedded", {}).get("terms", [])
            label = " ".join((terms[0].get("label") or "").split()) if terms else ""
        except Exception as e:
            sys.stderr.write(f"WARNING: OLS lookup failed for {_id}: {e}\n")
            label = ""
        if not label:
            results[_id] = {"id": _id, "manubot_success": False, "link": link,
                            "note": f"WARNING: OLS could not generate citation for {obo_id}"}
            continue
        results[_id] = {
            "id": _id,
            "manubot_success": True,
            "link": link,
            "title": label,
            "type": "entry",
            "doi": "",
            "authors": [],
            "publisher": "Mondo Disease Ontology",
            "issn": "",
            "date": "",
            "abstract": "",
            "language": "en",
            "note": f"Generated from the EBI Ontology Lookup Service ({obo_id})",
        }
    return results
