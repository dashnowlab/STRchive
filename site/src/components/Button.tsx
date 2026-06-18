import type { ComponentProps, Ref } from "react";
import Link from "@/components/Link";
import clsx from "clsx";

type ButtonProps = ComponentProps<"button">;
type AnchorProps = ComponentProps<typeof Link>;

type Props = {
  ref?: Ref<HTMLButtonElement>;
  design?: "" | "plain" | "bubble";
} & (ButtonProps | AnchorProps);

/** looks like a button, and either does something or goes somewhere */
export default function Button({
  ref,
  design = "",
  className = "",
  ...props
}: Props) {
  const Component = "to" in props ? Link : "button";

  return (
    <Component
      ref={ref}
      type="button"
      className={clsx(
        "inline-flex min-h-10 max-w-full min-w-10 items-center justify-center gap-2 px-[0.75em] py-[0.5em] no-underline",
        design === "" && "hover:text-primary",
        design === "plain" &&
          "rounded-md bg-light-gray text-current hover:text-primary",
        design === "bubble" &&
          "gap-2 rounded-full text-lg font-medium text-current ring-2 ring-current ring-inset hover:bg-secondary hover:text-white hover:ring-black",
        className,
      )}
      {...(props as ButtonProps & AnchorProps)}
    />
  );
}
