import type { ComponentProps, ReactNode } from "react";
import Loading from "@/assets/loading.svg?react";
import {
  IconAlertCircle,
  IconAlertTriangle,
  IconCircleCheck,
  IconInfoCircle,
} from "@tabler/icons-react";
import clsx from "clsx";

/** available categories of marks and associated styles */
export const types = {
  info: { className: "text-tertiary", icon: <IconInfoCircle /> },
  loading: { className: "text-gray", icon: <Loading /> },
  success: { className: "text-primary", icon: <IconCircleCheck /> },
  warning: { className: "text-secondary", icon: <IconAlertCircle /> },
  error: { className: "text-secondary", icon: <IconAlertTriangle /> },
};

type Type = keyof typeof types;

type Props = {
  type?: string;
  icon?: ReactNode;
} & ComponentProps<"div">;

/** colored box with icon and text */
export default function Alert({
  type = "info",
  icon,
  className = "",
  children,
  ...props
}: Props) {
  return (
    <div
      className={clsx(
        "relative isolate flex shrink-0 items-center gap-4 overflow-hidden rounded-md p-4",
        className,
      )}
      {...props}
    >
      <div
        className={clsx(
          "absolute inset-0 -z-10 bg-current opacity-10",
          types[type as Type]?.className,
        )}
      />
      <div className={types[type as Type]?.className}>
        {icon ?? types[type as Type]?.icon}
      </div>
      <div>{children}</div>
    </div>
  );
}
