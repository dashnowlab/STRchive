import Button from "@/components/Button";
import Link from "@/components/Link";
import Popover from "@/components/Popover";
import { tagOptions } from "@/data/tags";
import clsx from "clsx";
import { startCase } from "lodash-es";

type Props = {
  value: string;
  to?: string;
  tooltip?: string;
  small?: boolean;
  className?: string;
};

/** tag pill link */
export default function Tag({
  value,
  to = "",
  tooltip = "",
  small = false,
  className = "",
}: Props) {
  /** look up matching tag options */
  const option = tagOptions.find((option) => option.value === value);

  /** fallback text */
  const fallback = startCase(value);

  /** get tag option props */
  const {
    Icon,
    label = fallback,
    className: optionClassName = "text-black bg-light-gray",
    description = fallback,
  } = option || {};

  /** combine class */
  className = clsx(optionClassName, className);

  /** combine tooltip */
  const _tooltip = (
    <>
      {description}
      <br />
      {tooltip}
    </>
  );

  if (small)
    return (
      <Popover content={_tooltip}>
        <Button className={clsx("size-5 rounded-full", className)}>
          {Icon && <Icon className="size-4" />}
        </Button>
      </Popover>
    );
  return (
    <Popover content={_tooltip} button={false}>
      <Link
        className={clsx(
          "flex items-center gap-2 rounded-full px-2 py-1 leading-normal no-underline transition hover:opacity-75",
          className,
        )}
        to={to}
      >
        {Icon && <Icon />}
        {label}
      </Link>
    </Popover>
  );
}
