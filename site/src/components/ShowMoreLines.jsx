import { useEffect, useRef, useState } from "react";
import clsx from "clsx";
import classes from "./ShowMoreLines.module.css";

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
        classes.content,
        expanded ? classes.expanded : classes.collapsed,
        !expanded && "truncate-lines",
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
