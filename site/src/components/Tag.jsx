import { startCase } from "lodash-es";
import Link from "@/components/Link";
import { tagOptions } from "@/data/tags";
import classes from "./Tag.module.css";

/** tag pill link */
const Tag = ({ value, to = "", tooltip = "", small = false }) => {
  const option = tagOptions.find((o) => o.value === value) ?? {};
  const fallback = startCase(value);
  const {
    Icon = () => <></>,
    label = fallback,
    bg = "var(--gray)",
    text = "var(--white)",
    description = fallback,
  } = option;
  return (
    <Link
      className={classes.tag}
      to={to}
      style={
        small ? { color: bg, padding: 0 } : { backgroundColor: bg, color: text }
      }
      data-tooltip={[description, tooltip].flat().filter(Boolean).join("<br/>")}
    >
      <Icon className={classes.icon} />
      {!small && label}
    </Link>
  );
};

export default Tag;
