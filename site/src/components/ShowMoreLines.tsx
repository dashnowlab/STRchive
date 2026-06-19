import type { ReactNode } from "react";
import { useEffect, useRef, useState } from "react";
import clsx from "clsx";

type Props = {
  lines?: keyof typeof counts;
  children: ReactNode;
};

const counts = {
  1: "line-clamp-1",
  2: "line-clamp-2",
  3: "line-clamp-3",
  4: "line-clamp-4",
  5: "line-clamp-5",
};

/** add "show more/less" control to long lines of content if needed */
export default function ShowMoreLines({ lines = 2, children }: Props) {
  const ref = useRef<HTMLDivElement>(null);
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
    const count = Math.round(height / lineHeight);

    /** hide control if content can fit within limit */
    setShow(count > lines);
  }, [lines]);

  return show ? (
    <div
      className={clsx(
        expanded
          ? "cursor-zoom-out"
          : "cursor-zoom-in decoration-dotted not-hover:underline",
        !expanded && counts[lines],
      )}
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
}
