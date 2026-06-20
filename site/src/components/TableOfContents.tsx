import { useEffect, useRef, useState } from "react";
import { firstInView } from "@/util/dom";
import { useEventListener } from "@reactuses/core";
import { IconMenu2, IconX } from "@tabler/icons-react";
import Button from "./Button";

/** all used heading elements */
const headingSelector = "h2[id], h3[id], h4[id]";

/**
 * floating table of contents that outlines sections/headings on page. can be
 * turned on/off at route level. singleton.
 */
export default function TableOfContents() {
  const list = useRef(null);

  /** open/closed state */
  const [open, setOpen] = useState(window.innerWidth > 1600);

  type Heading = {
    element: Element;
    html: string;
    id: string;
    level: number;
  };

  /** full heading details */
  const [headings, setHeadings] = useState<Heading[]>([]);

  /** active heading id (first in view) */
  const [activeId, setActiveId] = useState("");

  /** on page load */
  useEffect(() => {
    /** read headings from page */
    // eslint-disable-next-line
    setHeadings(
      [...document.querySelectorAll(headingSelector)].map((element) => {
        const clone = element.cloneNode(true) as Element;
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
    <aside className="sticky top-0 z-10" aria-label="Table of contents">
      <div className="absolute top-0 flex max-w-80 flex-col gap-2 overflow-hidden rounded-br-md bg-white shadow-md">
        {/* toggle button */}
        <Button
          className="justify-between!"
          aria-expanded={open}
          onClick={() => setOpen(!open)}
        >
          Table of Contents
          {open ? <IconX /> : <IconMenu2 />}
        </Button>

        {/* links */}
        {open && (
          <div
            ref={list}
            className="flex max-h-[40vh] flex-col overflow-y-auto"
          >
            {headings.map(({ id, level, html }, index) => (
              <a
                key={index}
                data-active={id === activeId ? "" : undefined}
                className="flex items-center gap-2 px-4 py-2 leading-normal text-inherit no-underline transition hover:bg-light-gray hover:text-primary data-active:bg-light-gray [&_svg]:opacity-25"
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
}
