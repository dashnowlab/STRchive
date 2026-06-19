import type { ComponentProps, ReactNode } from "react";
import { useEffect, useRef } from "react";
import { IconCheck, IconX } from "@tabler/icons-react";
import clsx from "clsx";

type Props = {
  label: ReactNode;
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
      <div
        className={clsx(
          "grid size-5 appearance-none place-content-center rounded-sm border-current text-white accent-primary ring-2 ring-black ring-inset",
          checked === true
            ? "bg-primary"
            : checked === false
              ? "bg-secondary"
              : "bg-white",
          className,
        )}
      >
        <input
          ref={ref}
          type="checkbox"
          className="sr-only"
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
        {checked === true && <IconCheck />}
        {checked === false && <IconX />}
      </div>
      {label && <span>{label}</span>}
    </label>
  );
}
