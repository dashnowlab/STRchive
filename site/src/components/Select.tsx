import type { ReactNode } from "react";
import { useRef } from "react";
import Button from "@/components/Button";
import { preserveScroll } from "@/util/dom";
import { Select as _Select } from "@base-ui/react";
import { IconCheck, IconChevronDown } from "@tabler/icons-react";
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
  label?: ReactNode;
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
        {label && <span>{label}</span>}
        <_Select.Trigger render={<Button ref={ref} design="plain" />}>
          <_Select.Value>
            {(value: string | string[]) => {
              let selected: ReactNode = "";
              if (Array.isArray(value)) {
                if (!value.length) selected = "";
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
            <IconChevronDown />
          </_Select.Icon>
        </_Select.Trigger>
      </label>

      <_Select.Portal>
        <_Select.Positioner alignItemWithTrigger={false} collisionPadding={20}>
          <_Select.Popup>
            <_Select.List className="flex flex-col rounded-md bg-white shadow-md">
              {options.map(({ label, value, tooltip }, index) => (
                <_Select.Item
                  key={index}
                  value={value}
                  className="flex cursor-pointer items-center gap-2 px-4 py-2 transition outline-none *:first:opacity-0 *:first:transition hover:bg-light-gray data-highlighted:bg-light-gray data-selected:*:first:opacity-100"
                  data-tooltip={tooltip}
                >
                  <IconCheck className="text-primary" />
                  {label || "–"}
                </_Select.Item>
              ))}
            </_Select.List>
          </_Select.Popup>
        </_Select.Positioner>
      </_Select.Portal>
    </_Select.Root>
  );
}
