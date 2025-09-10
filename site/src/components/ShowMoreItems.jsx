import { useState } from "react";
import Button from "@/components/Button";
import classes from "./ShowMoreItems.module.css";

/** add "show more/less" control to group of items if needed */
const ShowMoreItems = ({ items, limit = 9 }) => {
  /** whether to show all items */
  const [expanded, setExpanded] = useState(false);

  return (
    <>
      {items.slice(0, expanded ? Infinity : limit)}
      {items.length > limit && (
        <Button
          design="plain"
          className={classes.button}
          onClick={() => setExpanded(!expanded)}
          aria-expanded={expanded}
        >
          {expanded
            ? "Show Less"
            : `Show All (+${(items.length - limit).toLocaleString()})`}
        </Button>
      )}
    </>
  );
};

export default ShowMoreItems;
