import type { ComponentProps } from "react";
import { LuExternalLink } from "react-icons/lu";

type Props = {
  to: string;
  newTab?: boolean;
  arrow?: boolean;
} & ComponentProps<"a">;

export default function Link({ to, newTab, arrow, children, ...props }: Props) {
  /** whether link is to external site, or page within this site */
  const external = !!to.match(/^(https|http|ftp|mailto)/);

  return (
    <a
      href={to || undefined}
      // whether to open in new tab
      target={(newTab ?? external) ? "_blank" : ""}
      {...props}
    >
      {children}
      {/* indicate third-party site with icon  */}
      {(arrow ?? external) && !!children && (
        <LuExternalLink className="ml-1 inline-block -translate-y-0.5" />
      )}
    </a>
  );
}
