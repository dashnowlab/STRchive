import { useRef } from "react";
import { FaXmark } from "react-icons/fa6";
import clsx from "clsx";
import Button from "@/components/Button";

/** text input with label */
const TextBox = ({ label, multi, onChange, className, tooltip, ...props }) => {
  const ref = useRef(null);

  const Component = multi ? "textarea" : "input";

  return (
    <label data-tooltip={tooltip}>
      {label && <span>{label}</span>}
      <div className="relative flex w-full items-start has-[textarea]:max-w-full has-[textarea]:min-h-20 has-[textarea]:h-[120px] has-[textarea]:max-h-full has-[textarea]:resize-both has-[textarea]:overflow-hidden">
        <Component
          ref={ref}
          type="text"
          className={clsx(
            className,
            "min-h-10 min-w-10 w-full rounded-md bg-white px-[15px] py-[7.5px] pr-10 [box-shadow:inset_0_0_0_2px_var(--color-dark-gray)]",
            multi && "h-full resize-none",
          )}
          onChange={(event) => onChange?.(event.target.value)}
          {...props}
        />

        <div className="absolute right-0 top-0 flex items-center">
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
};

export default TextBox;
