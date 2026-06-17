import Link from "@/components/Link";
import { tagOptions } from "@/data/tags";
import { startCase } from "lodash-es";

type Props = {
  value: string;
  to?: string;
  tooltip?: string;
  small?: boolean;
};

/** tag pill link */
export default function Tag({
  value,
  to = "",
  tooltip = "",
  small = false,
}: Props) {
  const option = tagOptions.find((o) => o.value === value);
  const fallback = startCase(value);
  const {
    Icon,
    label = fallback,
    bg = "var(--color-gray)",
    text = "var(--color-white)",
    description = fallback,
  } = option || {};
  return (
    <Link
      className="flex items-center gap-2 rounded-full px-2 py-1 no-underline transition hover:opacity-75"
      to={to}
      style={
        small ? { color: bg, padding: 0 } : { backgroundColor: bg, color: text }
      }
      data-tooltip={[description, tooltip].flat().filter(Boolean).join("<br/>")}
    >
      {Icon && <Icon />}
      {!small && label}
    </Link>
  );
}
