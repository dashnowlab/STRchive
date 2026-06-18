import type { ElementType, ReactNode } from "react";
import { useRef } from "react";
import { LuCheck, LuChevronDown } from "react-icons/lu";
import Button from "@/components/Button";
import { preserveScroll } from "@/util/dom";
import {
  Listbox,
  ListboxButton,
  ListboxOption,
  ListboxOptions,
} from "@headlessui/react";
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
  multi = false,
  options,
  defaultValue,
  value,
  onChange,
  tooltip,
}: Props) {
  const labelRef = useRef<HTMLLabelElement>(null);

  /** track open state */
  const isOpen = useRef(false);

  return (
    <label ref={labelRef} data-tooltip={tooltip}>
      {label && <span>{label}</span>}
      <Listbox<ElementType, string | string[]>
        value={value}
        onChange={onChange as (value: string | string[]) => void}
        multiple={multi}
        defaultValue={defaultValue}
      >
        {({ value, open }) => {
          /** on close */
          if (!open && isOpen.current) preserveScroll(labelRef.current);
          isOpen.current = open;

          /** label for selected value */
          let selected: ReactNode = "";
          if (Array.isArray(value)) {
            if (value.length === 0) selected = "";
            else if (value.length === options.length) selected = "All selected";
            else if (value.length === 1)
              selected =
                options.find((option) => option.value === value[0])?.label ||
                value[0];
            else selected = `${value.length} selected`;
          } else {
            selected =
              options.find((option) => option.value === value)?.label ||
              value ||
              "";
          }

          return (
            /** https://github.com/tailwindlabs/headlessui/issues/3597 */
            <div className="contents">
              <ListboxButton
                as={Button}
                design="plain"
                className="relative flex min-h-10 min-w-10 grow items-center gap-2 rounded-md bg-light-gray"
              >
                <span className={clsx("grow", !selected && "text-gray")}>
                  {selected || "No Selection"}
                </span>
                <LuChevronDown />
              </ListboxButton>

              <ListboxOptions
                className="flex w-(--button-width) flex-col rounded-md bg-white shadow-md"
                anchor={{ to: "bottom", gap: 2, padding: 10 }}
              >
                {options.map((option, index) => (
                  <ListboxOption
                    key={index}
                    className="flex cursor-pointer items-center gap-2 px-2 py-4 transition hover:bg-light-gray data-focus:bg-light-gray"
                    value={option.value}
                    data-tooltip={option.tooltip}
                  >
                    <LuCheck
                      className="text-primary"
                      style={{
                        opacity: (
                          multi
                            ? value?.includes(option.value)
                            : value === option.value
                        )
                          ? 1
                          : 0,
                      }}
                    />
                    <span className={clsx(!option.label && "text-gray")}>
                      {option.label || "No Selection"}
                    </span>
                  </ListboxOption>
                ))}
              </ListboxOptions>
            </div>
          );
        }}
      </Listbox>
    </label>
  );
}
