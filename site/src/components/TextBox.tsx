import type { ComponentProps, ReactNode } from "react";
import { useRef } from "react";
import Button from "@/components/Button";
import Help from "@/components/Help";
import { IconX } from "@tabler/icons-react";
import clsx from "clsx";

type Props = {
  label?: ReactNode;
  multi?: boolean;
  onChange?: (value: string) => void;
  help?: ReactNode;
} & Omit<ComponentProps<"input"> | ComponentProps<"textarea">, "onChange">;

/** text input with label */
export default function TextBox({
  label,
  multi,
  onChange,
  className,
  help,
  ...props
}: Props) {
  const Component = multi ? "textarea" : "input";
  const ref = useRef<HTMLInputElement | HTMLTextAreaElement>(null);

  return (
    <label className={clsx(!label && !help && "contents")}>
      {(label || help) && (
        <span>
          {label}
          <Help>{help}</Help>
        </span>
      )}
      <div
        className={clsx(
          "relative flex min-h-10 min-w-10 items-start rounded-md bg-white ring-2 ring-inset has-[textarea]:h-30 has-[textarea]:max-h-full has-[textarea]:min-h-20 has-[textarea]:max-w-full has-[textarea]:resize has-[textarea]:overflow-hidden",
          className,
        )}
      >
        <Component
          // @ts-expect-error ref type is not properly inferred
          ref={ref}
          type="text"
          className={clsx(
            "w-full rounded-md px-4 py-2 pr-10 leading-normal",
            multi && "h-full resize-none",
          )}
          onChange={(event) => onChange?.(event.target.value)}
          {...props}
        />

        <div className="absolute top-0 right-0 flex size-10 items-center">
          <Button
            design="hollow"
            onClick={() => {
              if (ref.current) ref.current.value = "";
              onChange?.("");
            }}
            aria-label="Clear textbox"
          >
            <IconX />
          </Button>
        </div>
      </div>
    </label>
  );
}
