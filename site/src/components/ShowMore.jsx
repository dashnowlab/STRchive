import { useEffect, useRef, useState } from "react";
import classes from "./ShowMore.module.css";

/** add "show more/less" control to long lines of content if needed */
const ShowMore = ({ lines = 2, children }) => {
  const ref = useRef();
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

  return (
    <>
      <span
        ref={ref}
        className={show && !expanded ? classes.box : undefined}
        style={{ WebkitLineClamp: show && !expanded ? lines : undefined }}
      >
        {children}
      </span>
      {show && (
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
