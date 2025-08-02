import { useEffect, useRef, useState } from "react";
import { FaBars, FaXmark } from "react-icons/fa6";
import { useEventListener } from "@reactuses/core";
import { firstInView } from "@/util/dom";
import classes from "./TableOfContents.module.css";

/** all used heading elements */
const headingSelector = "h2, h3, h4";

/**
 * floating table of contents that outlines sections/headings on page. can be
 * turned on/off at route level. singleton.
 */
const TableOfContents = () => {
  const list = useRef(null);

  /** open/closed state */
  const [open, setOpen] = useState(window.innerWidth > 1200);

  /** full heading details */
  const [headings, setHeadings] = useState([]);

  /** active heading id (first in view) */
  const [activeId, setActiveId] = useState("");

  /** on page load */
  useEffect(() => {
    /** read headings from page */
    setHeadings(
      [...document.querySelectorAll(headingSelector)].map((element) => {
        const clone = element.cloneNode(true);
        clone.querySelector("a").remove();
        const { innerHTML, id, tagName } = clone;
        const html = innerHTML.trim();
        const level = (parseInt(tagName.slice(1)) || 0) - 1;
        return { element, html, id, level };
      }),
    );
  }, []);

  /** on window scroll */
  useEventListener("scroll", () => {
    /** get active heading */
    setActiveId(firstInView(headings.map(({ element }) => element))?.id || "");
  });

  /** if not much value in showing toc, hide */
  if (headings.length <= 1) return <></>;

  return (
    <aside className={classes.aside} aria-label="Table of contents">
      <div className={classes.table}>
        {/* toggle button */}
        <button
          type="button"
          className={classes.button}
          aria-expanded={open}
          onClick={() => setOpen(!open)}
          data-tooltip={open ? "Close" : "Open"}
        >
          Table of Contents{open ? <FaXmark /> : <FaBars />}
        </button>

        {/* links */}
        {open && (
          <div ref={list} className={classes.list}>
            {headings.map(({ id, level, html }, index) => (
              <a
                key={index}
                data-active={id === activeId ? "" : undefined}
                className={classes.link}
                href={`#${id}`}
                style={{ paddingLeft: level * 15 }}
                dangerouslySetInnerHTML={{ __html: html }}
              />
            ))}
          </div>
        )}
      </div>
    </aside>
  );
};

export default TableOfContents;
