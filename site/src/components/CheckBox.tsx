import type { ComponentProps, ReactNode } from "react";
import { useEffect, useRef } from "react";
import checkUrl from "@/assets/check.svg?inline";
import xUrl from "@/assets/x.svg?inline";
import clsx from "clsx";

type Props = {
  label?: ReactNode;
  checked: boolean | "mixed";
  onChange?: (checked: boolean | "mixed") => void;
  tooltip?: string;
} & Omit<ComponentProps<"input">, "onChange" | "checked">;

/** checkbox with label */
export default function CheckBox({
  label,
  checked,
  onChange,
  className,
  tooltip,
  ...props
}: Props) {
  const ref = useRef<HTMLInputElement>(null);

  useEffect(() => {
    if (ref.current) ref.current.indeterminate = checked === "mixed";
  }, [checked]);

  return (
    <label data-tooltip={tooltip}>
      <input
        ref={ref}
        type="checkbox"
        className={clsx(
          "size-5 appearance-none rounded-md border-2 border-current accent-primary",
          checked === true
            ? "bg-primary"
            : checked === false
              ? "bg-secondary"
              : "bg-white",
          className,
        )}
        style={{
          backgroundImage:
            checked === true
              ? `url("${checkUrl}")`
              : checked === false
                ? `url("${xUrl}")`
                : "none",
        }}
        checked={!!checked}
        onChange={() => {
          if (checked === "mixed") onChange?.(true);
          else if (checked === true) onChange?.(false);
          else onChange?.("mixed");
        }}
        aria-checked={checked}
        aria-label={tooltip}
        {...props}
      />
      {label && <span>{label}</span>}
    </label>
  );
}
