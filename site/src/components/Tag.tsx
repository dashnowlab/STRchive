import Link from "@/components/Link";
import { tagOptions } from "@/data/tags";
import clsx from "clsx";
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
  /** look up matching tag options */
  const option = tagOptions.find((option) => option.value === value);

  /** fallback text */
  const fallback = startCase(value);

  /** get tag option props */
  const {
    Icon,
    label = fallback,
    className = "text-gray",
    description = fallback,
  } = option || {};

  /** combine tooltip */
  tooltip = [description, tooltip].flat().filter(Boolean).join("<br/>");

  if (small)
    return (
      <div
        className={clsx("rounded-full p-1", className)}
        data-tooltip={tooltip}
      >
        {Icon && <Icon className="size-3" />}
      </div>
    );
  return (
    <Link
      className={clsx(
        "flex items-center gap-2 rounded-full px-2 py-1 no-underline transition hover:opacity-75",
        className,
      )}
      to={to}
      data-tooltip={tooltip}
    >
      {Icon && <Icon />}
      {label}
    </Link>
  );
}
