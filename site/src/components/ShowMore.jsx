import { useEffect, useRef, useState } from "react";
import classes from "./ShowMore.module.css";

/** collapse long lines of text */
const ShowMore = ({ lines = 2, children }) => {
  const ref = useRef();
  const [enabled, setEnabled] = useState(false);
  const [expanded, setExpanded] = useState(false);

  /** on first render */
  useEffect(() => {
    if (!ref.current) return;
    /** select text */
    const range = document.createRange();
    range.setStartBefore(ref.current);
    range.setEndAfter(ref.current);
    /** get number of lines */
    const rects = range.getClientRects();
    /** if line limit exceeded, enable component */
    if (rects.length - 1 > lines) setEnabled(true);
  }, []);

  return (
    <>
      <span
        ref={ref}
        className={enabled && !expanded ? classes.box : undefined}
        style={{ WebkitLineClamp: enabled && !expanded ? lines : undefined }}
      >
        {children}
      </span>
      {enabled && (
        <button
          className={classes.button}
          onClick={() => setExpanded(!expanded)}
        >
          {expanded ? "Show less" : "Show more"}
        </button>
      )}
    </>
  );
};

export default ShowMore;
