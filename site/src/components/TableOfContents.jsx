import { useEffect, useRef, useState } from "react";
import { FaBars, FaXmark } from "react-icons/fa6";
import { useEventListener } from "@reactuses/core";
import { firstInView } from "@/util/dom";
import Button from "./Button";

/** all used heading elements */
const headingSelector = "h1[id], h2[id], h3[id], h4[id]";

/**
 * floating table of contents that outlines sections/headings on page. can be
 * turned on/off at route level. singleton.
 */
const TableOfContents = () => {
  const list = useRef(null);

  /** open/closed state */
  const [open, setOpen] = useState(window.innerWidth > 1400);

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
        const level = parseInt(tagName.slice(1)) || 0;
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
    <aside className="sticky top-0 z-1" aria-label="Table of contents">
      <div className="absolute top-0 flex max-w-[300px] flex-col gap-2.5 overflow-hidden rounded-br-md bg-white shadow-md">
        {/* toggle button */}
        <Button
          className="justify-between"
          aria-expanded={open}
          onClick={() => setOpen(!open)}
          data-tooltip={open ? "Close" : "Open"}
        >
          <span>Table of Contents</span>
          {open ? <FaXmark /> : <FaBars />}
        </Button>

        {/* links */}
        {open && (
          <div ref={list} className="flex max-h-[40vh] flex-col overflow-y-auto">
            {headings.map(({ id, level, html }, index) => (
              <a
                key={index}
                data-active={id === activeId ? "" : undefined}
                className="flex items-center gap-2.5 px-[15px] py-[7.5px] text-inherit no-underline transition hover:bg-light-gray hover:text-primary data-[active]:font-bold [&_svg]:opacity-25"
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
