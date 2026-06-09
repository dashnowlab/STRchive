import { useEffect, useRef, useState } from "react";
import clsx from "clsx";

/** add "show more/less" control to long lines of content if needed */
const ShowMoreLines = ({ lines = 2, children }) => {
  const ref = useRef(null);
  /** whether to show the control at all */
  const [show, setShow] = useState(false);
  /** whether to show all content */
  const [expanded, setExpanded] = useState(false);

  /** on first render */
  useEffect(() => {
    if (!ref.current) return;

    /** height of full content */
    const { height } = ref.current.getBoundingClientRect();

    /** height of each line */
    const lineHeight = parseFloat(
      window.getComputedStyle(ref.current).lineHeight,
    );

    /** estimate number of lines of full content */
    const lineCount = Math.round(height / lineHeight);

    /** hide control if content can fit within limit */
    setShow(lineCount > lines);
  }, []);

  return show ? (
    <div
      className={clsx(
        "[&:not(:has(:is(span,sup)))]:underline [&:not(:has(:is(span,sup)))]:decoration-dotted [&:not(:has(:is(span,sup)))]:decoration-1 [&:not(:has(:is(span,sup)))]:underline-offset-[2px] [&:has(sup)_span]:underline [&:has(sup)_span]:decoration-dotted [&:has(sup)_span]:decoration-1 [&:has(sup)_span]:underline-offset-[2px] hover:[&:not(:has(:is(span,sup)))]:no-underline hover:[&:has(sup):has(span:hover)_span]:no-underline",
        expanded
          ? "hover:cursor-zoom-out hover:[&:has(:is(sup):hover)]:cursor-auto"
          : "hover:cursor-zoom-in hover:[&:has(:is(sup):hover)]:cursor-auto",
        !expanded &&
          "overflow-hidden [display:-webkit-box] [-webkit-box-orient:vertical] [-webkit-line-clamp:var(--lines,1)]",
      )}
      style={{ "--lines": show && !expanded ? lines : undefined }}
      onClick={() => setExpanded(!expanded)}
      onKeyUp={(event) => event.key === "Enter" && setExpanded(!expanded)}
      role="button"
      tabIndex={0}
      aria-expanded={expanded}
    >
      {children}
    </div>
  ) : (
    <div ref={ref}>{children}</div>
  );
};

export default ShowMoreLines;
