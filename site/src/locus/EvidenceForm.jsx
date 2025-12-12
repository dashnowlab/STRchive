import { useMemo } from "react";
import { FaXmark } from "react-icons/fa6";
import { LuBookCheck, LuSend } from "react-icons/lu";
import { cloneDeep, isEqual, omitBy, startCase, sumBy } from "lodash-es";
import { useLocalStorage } from "@reactuses/core";
import { createPR } from "@/api/pr";
import Alert from "@/components/Alert";
import Button from "@/components/Button";
import { contactSchema } from "@/components/ContactForm";
import Form from "@/components/Form";
import Heading from "@/components/Heading";
import Link from "@/components/Link";
import SchemaForm from "@/components/SchemaForm";
import { repo } from "@/layouts/meta";
import { useQuery } from "@/util/hooks";
import { shortenUrl } from "@/util/string";
import loci from "~/STRchive-loci.json";
import schema from "~/STRtr-kit-evidence.schema.json";

/** add extra fields for evidence metadata */
schema.properties = {
  "evidence-name": {
    section: "Overview",
    ...contactSchema.name,
    type: ["string", "null"],
    default: "",
  },
  "evidence-username": {
    section: "Overview",
    ...contactSchema.username,
    type: ["string", "null"],
    pattern: "^@.+",
    default: "",
  },
  "evidence-email": {
    section: "Overview",
    ...contactSchema.email,
    type: ["string", "null"],
    default: "",
  },
  "evidence-title": {
    section: "Overview",
    title: "Evidence Title",
    description: "Succinct title describing these changes",
    examples: ["Fix mechanism details", "Update disease onset information"],
    type: "string",
    default: null,
  },
  "evidence-description": {
    section: "Overview",
    title: "Evidence Description",
    description:
      "Summary of changes, justification for changes, uncertainty in literature, or anything else we should know for review. Please be detailed. Provide at least 2-3 sentences.",
    examples: [
      "Currently, the disease mechanism details cite a recently retracted paper doi:123456. This evidence corrects the reference and updates...",
    ],
    multiline: true,
    type: ["string", "null"],
    default: null,
  },
  ...schema.properties,
};

schema.required.push("evidence-title");

/** new/evidence locus form */
const EvidenceForm = ({ heading, locus }) => {
  /** confirm with user before leaving page */
  // window.onbeforeunload = () => "";

  /** unique storage key for this page and form */
  const storageKey = `evidence-locus-${locus?.id ?? "new"}`;

  /** form data state */
  let [data, setData] = useLocalStorage(storageKey, {
    ["evidence-title"]: null,
    ["evidence-description"]: null,
    id: locus?.id ?? null,
  });

  /** was data loaded from storage */
  const storageExists = useMemo(() => {
    const fromStorage = window.localStorage.getItem(storageKey);
    /** if saved draft exists and is different from initial data */
    return fromStorage && !isEqual(JSON.parse(fromStorage), data);
  }, []);

  /** submission query */
  const {
    query: submit,
    data: response,
    status,
  } = useQuery(async () => {
    console.info(data);
  });

  /** preview/summary */
  const geneticPoints = sumBy(data.genetic_evidence, "points") || 0;
  const experimentalPoints = sumBy(data.experimental_evidence, "points") || 0;
  const totalPoints = geneticPoints + experimentalPoints;

  return (
    <Form onSubmit={submit}>
      <section>
        <Heading level={1}>
          <LuBookCheck />
          {heading}
        </Heading>

        <Alert type="info">
          Every suggestion is reviewed by our team before inclusion in STRchive.
          Please enter as much accurate information as possible.
        </Alert>

        {storageExists && (
          <div className="row">
            Loaded saved draft
            <Button
              design="plain"
              onClick={() => {
                window.localStorage.removeItem(storageKey);
                window.location.reload();
              }}
            >
              <span>Forget</span>
              <FaXmark />
            </Button>
          </div>
        )}
      </section>

      <SchemaForm
        schema={schema}
        data={data}
        onChange={setData}
        sections={["Overview", "Evidence"]}
        onSubmit={(data) => submit(data)}
      />

      <section>
        <Heading level={2}>Summary</Heading>

        <dl>
          <dt>Genetic Points</dt>
          <dd>{geneticPoints}</dd>

          <dt>Experimental Points</dt>
          <dd>{experimentalPoints}</dd>

          <dt>Total</dt>
          <dd>{totalPoints}</dd>
        </dl>

        <Alert type={status || "info"}>
          {status === "" && (
            <>
              This will make a <strong>public</strong> pull request on{" "}
              <Link to={repo}>our GitHub</Link>. You'll get a link once it's
              created.
            </>
          )}
          {startCase(status)}{" "}
          {response && (
            <Link to={response.link}>{shortenUrl(response.link)}</Link>
          )}
        </Alert>

        <Button type="submit" design="bubble" disabled={status === "success"}>
          <LuSend />
          <span>Submit</span>
        </Button>
      </section>
    </Form>
  );
};

export default EvidenceForm;
