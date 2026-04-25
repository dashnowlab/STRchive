import Link from "@/components/Link";
import { evidenceOptions } from "@/data/evidence";
import classes from "./Tag.module.css";

/** single evidence pill */
const EvidencePill = ({ option, to }) => {
  const {
    Icon = () => <></>,
    label,
    bg = "var(--gray)",
    text = "var(--white)",
    tooltip = label,
  } = option;

  const tooltipText = to
    ? `${tooltip}.<br/>See the criTRia curation for this locus.`
    : `${tooltip}.`;

  const inner = (
    <>
      <Icon className={classes.icon} />
      {label}
    </>
  );

  return to ? (
    <Link
      className={classes.tag}
      to={to}
      style={{ backgroundColor: bg, color: text }}
      data-tooltip={tooltipText}
    >
      {inner}
    </Link>
  ) : (
    <span
      className={classes.tag}
      style={{ backgroundColor: bg, color: text }}
      data-tooltip={tooltipText}
    >
      {inner}
    </span>
  );
};

/**
 * Evidence pills for a locus. `value` is the evidence array from locus data
 * (e.g. ["Moderate"]). Links to criTRia curation page when `critriaId` is set.
 */
const Evidence = ({ value, critriaId = undefined }) => {
  const values = Array.isArray(value) ? value : [value];
  const to = critriaId ? `/critria/${critriaId}` : undefined;

  return (
    <>
      {evidenceOptions
        .filter((opt) => values.includes(opt.value))
        .map((opt) => (
          <EvidencePill key={opt.value} option={opt} to={to} />
        ))}
    </>
  );
};

export default Evidence;
