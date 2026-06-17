import type { ComponentProps } from "react";
import { useRef } from "react";
import { FaXmark } from "react-icons/fa6";
import Button from "@/components/Button";
import clsx from "clsx";

type Props = {
  label?: string;
  multi?: boolean;
  onChange?: (value: string) => void;
  className?: string;
  tooltip?: string;
} & Omit<ComponentProps<"input"> | ComponentProps<"textarea">, "onChange">;

/** text input with label */
export default function TextBox({
  label,
  multi,
  onChange,
  className,
  tooltip,
  ...props
}: Props) {
  const Component = multi ? "textarea" : "input";
  const ref = useRef<HTMLInputElement | HTMLTextAreaElement>(null);

  return (
    <label data-tooltip={tooltip}>
      {label && <span>{label}</span>}
      <div className="relative flex w-full items-start has-[textarea]:h-30 has-[textarea]:max-h-full has-[textarea]:min-h-20 has-[textarea]:max-w-full has-[textarea]:resize has-[textarea]:overflow-hidden">
        <Component
          // @ts-expect-error ref type is not properly inferred
          ref={ref}
          type="text"
          className={clsx(
            className,
            "min-h-10 w-full min-w-10 rounded-md bg-white px-4 py-2 pr-10",
            multi && "h-full resize-none",
          )}
          onChange={(event) => onChange?.(event.target.value)}
          {...props}
        />

        <div className="absolute top-0 right-0 flex items-center">
          <Button
            onClick={() => {
              if (ref.current) ref.current.value = "";
              onChange?.("");
            }}
            aria-label="Clear textbox"
          >
            <FaXmark />
          </Button>
        </div>
      </div>
    </label>
  );
}
