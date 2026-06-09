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
        "inline-flex min-h-10 min-w-10 items-center justify-center gap-2 no-underline not-has-[span]:p-[0.5em] has-[span]:px-[0.75em] has-[span]:py-[0.5em]",
        design === "plain" && "rounded-md bg-light-gray",
        design === "bubble" &&
          "gap-2 rounded-full border-2 border-current font-medium text-lg text-current hover:bg-secondary hover:text-white",
        className,
      )}
      {...(props as ButtonProps & AnchorProps)}
    />
  );
}
