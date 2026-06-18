import type { ReactNode } from "react";
import { useRef } from "react";
import { LuCheck, LuChevronDown } from "react-icons/lu";
import Button from "@/components/Button";
import { preserveScroll } from "@/util/dom";
import { Select as _Select } from "@base-ui/react";
import clsx from "clsx";

type Single = {
  multi?: false;
  value?: string;
  defaultValue?: string;
  onChange: (value: string) => void;
};

type Multi = {
  multi: true;
  value?: string[];
  defaultValue?: string[];
  onChange: (value: string[]) => void;
};

type Props = {
  label: ReactNode;
  options: { value: string; label?: ReactNode; tooltip?: string }[];
  tooltip?: string;
} & (Single | Multi);

/** dropdown select with label */
export default function Select({
  label,
  multi,
  options,
  defaultValue,
  value,
  onChange,
  tooltip,
}: Props) {
  const ref = useRef<HTMLButtonElement>(null);

  return (
    <_Select.Root
      value={value}
      onValueChange={(value) => {
        if (value === null) return;
        if (multi) onChange(value as string[]);
        else onChange(value as string);
      }}
      onOpenChange={(open) => {
        if (open) return;
        if (ref.current) preserveScroll(ref.current);
      }}
      multiple={multi}
      defaultValue={defaultValue}
    >
      <label data-tooltip={tooltip}>
        {label}
        <_Select.Trigger render={<Button ref={ref} design="plain" />}>
          <_Select.Value>
            {(value: string | string[]) => {
              let selected: ReactNode = "";
              if (Array.isArray(value)) {
                if (value.length === 0) selected = "";
                else if (value.length === options.length)
                  selected = "All selected";
                else if (value.length === 1)
                  selected =
                    options.find((option) => option.value === value[0])
                      ?.label || value[0];
                else selected = `${value.length} selected`;
              } else {
                selected =
                  options.find((option) => option.value === value)?.label ||
                  value ||
                  "";
              }
              return (
                <span className={clsx("grow", !selected && "text-gray")}>
                  {selected || "No Selection"}
                </span>
              );
            }}
          </_Select.Value>

          <_Select.Icon>
            <LuChevronDown />
          </_Select.Icon>
        </_Select.Trigger>
        <_Select.Portal>
          <_Select.Positioner alignItemWithTrigger={false}>
            <_Select.Popup>
              <_Select.List className="flex min-w-(--button-width) flex-col rounded-md bg-white shadow-md">
                {options.map(({ label, value, tooltip }, index) => (
                  <_Select.Item
                    key={index}
                    value={value}
                    className="flex cursor-pointer items-center gap-2 px-2 py-4 transition *:first:opacity-0 *:first:transition hover:bg-light-gray data-selected:*:first:opacity-100"
                    data-tooltip={tooltip}
                  >
                    <LuCheck className="text-primary" />
                    <_Select.ItemText>{label}</_Select.ItemText>
                  </_Select.Item>
                ))}
              </_Select.List>
            </_Select.Popup>
          </_Select.Positioner>
        </_Select.Portal>
      </label>
    </_Select.Root>
  );
}
