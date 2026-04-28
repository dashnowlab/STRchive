import { startCase } from "lodash-es";
import Link from "@/components/Link";
import { tagOptions } from "@/data/tags";
import classes from "./Tag.module.css";

/** tag pill link */
const Tag = ({ value, to = "", tooltip = "" }) => {
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
      style={{ backgroundColor: bg, color: text }}
      data-tooltip={[description, tooltip].flat().filter(Boolean).join("<br/>")}
    >
      <Icon className={classes.icon} />
      {label}
    </Link>
  );
};

export default Tag;
