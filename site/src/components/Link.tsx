import type { ComponentProps, ReactNode } from "react";
import { LuExternalLink } from "react-icons/lu";

type Props = {
  to: string;
  newTab?: boolean;
  arrow?: boolean;
  children: ReactNode;
} & ComponentProps<"a" | "span">;

export default function Link({
  to,
  newTab = false,
  arrow = false,
  children,
  ...props
}: Props) {
  /** whether link is to external site, or page within this site */
  const external = !!to.match(/^(https|http|ftp|mailto)/);

  const Component = to.trim() ? "a" : "span";

  return (
    <Component
      href={to || ""}
      // whether to open in new tab
      target={(newTab ?? external) ? "_blank" : ""}
      {...(props as ComponentProps<"a"> & ComponentProps<"span">)}
    >
      {children}
      {/* indicate third-party site with icon  */}
      {(arrow ?? external) && !children && (
        <LuExternalLink className="relative ml-1 scale-75" />
      )}
    </Component>
  );
}
