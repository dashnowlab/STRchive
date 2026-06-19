import type { ComponentProps, ReactNode } from "react";
import { useState } from "react";
import Button from "@/components/Button";

type Props = {
  items: ReactNode[];
  limit?: number;
} & ComponentProps<"button">;

/** add "show more/less" control to group of items if needed */
export default function ShowMoreItems({ items, limit = 9, ...props }: Props) {
  /** whether to show all items */
  const [expanded, setExpanded] = useState(false);

  return (
    <>
      {items.slice(0, expanded ? Infinity : limit)}
      {items.length > limit && (
        <Button
          className="self-center"
          design="plain"
          onClick={() => setExpanded(!expanded)}
          aria-expanded={expanded}
          {...props}
        >
          {expanded
            ? "Show Less"
            : `Show All (+${(items.length - limit).toLocaleString()})`}
        </Button>
      )}
    </>
  );
}
