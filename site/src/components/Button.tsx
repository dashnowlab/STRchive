import type { ComponentProps } from "react";
import Link from "@/components/Link";
import clsx from "clsx";

type ButtonProps = ComponentProps<"button">;
type AnchorProps = ComponentProps<typeof Link>;

type Props = {
  design?: "" | "plain" | "bubble";
} & (ButtonProps | AnchorProps);

/** looks like a button, and either does something or goes somewhere */
export default function Button({
  design = "",
  className = "",
  ...props
}: Props) {
  const Component = "to" in props ? Link : "button";

  return (
    <Component
      type="button"
      className={clsx(
        "inline-flex min-h-10 min-w-10 items-center justify-center gap-2 px-[0.75em] py-[0.5em] no-underline",
        design === "" && "hover:text-primary",
        design === "plain" &&
          "rounded-md bg-light-gray text-current hover:text-primary",
        design === "bubble" &&
          "gap-2 rounded-full border-2 border-current text-lg font-medium text-current hover:border-black hover:bg-secondary hover:text-white",
        className,
      )}
      {...(props as ButtonProps & AnchorProps)}
    />
  );
}
