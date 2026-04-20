import { startCase } from "lodash-es";
import Link from "@/components/Link";
import { tagOptions } from "@/data/tags";
import classes from "./Tag.module.css";

/** tag pill link */
const Tag = ({ value }) => {
  const option = tagOptions.find((t) => t.value === value) ?? {};
  const fallback = startCase(value);
  const {
    Icon = () => <></>,
    label = fallback,
    bg = "var(--gray)",
    text = "var(--white)",
    tooltip = fallback,
  } = option;
  return (
    <Link
      className={classes.tag}
      to={`/loci?tag=${value}#table`}
      style={{ backgroundColor: bg, color: text }}
      data-tooltip={`${tooltip}.<br/>See other loci with this tag.`}
    >
      <Icon className={classes.icon} />
      {label}
    </Link>
  );
};

export default Tag;
