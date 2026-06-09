import { useRef } from "react";
import { FaAngleDown, FaCheck } from "react-icons/fa6";
import clsx from "clsx";
import {
  Listbox,
  ListboxButton,
  ListboxOption,
  ListboxOptions,
} from "@headlessui/react";
import Button from "@/components/Button";
import { preserveScroll } from "@/util/dom";

/** dropdown select with label */
const Select = ({
  label,
  multi = false,
  options,
  defaultValue,
  value,
  onChange,
  tooltip,
}) => {
  const labelRef = useRef(null);

  /** track open state */
  const isOpen = useRef(false);

  return (
    <label ref={labelRef} data-tooltip={tooltip}>
      {label && <span>{label}</span>}
      <Listbox
        value={value}
        onChange={onChange}
        multiple={multi}
        defaultValue={defaultValue}
      >
        {({ value, open }) => {
          /** on close */
          if (!open && isOpen.current) preserveScroll(labelRef.current);
          isOpen.current = open;

          /** label for selected value */
          let selected = "";
          if (multi) {
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
            <div style={{ display: "contents" }}>
              <ListboxButton
                as={Button}
                className="relative flex min-h-10 min-w-10 grow items-center gap-2.5 rounded-md bg-light-gray"
              >
                <span className={clsx("grow", !selected && "text-gray")}>
                  {selected || "No Selection"}
                </span>
                <FaAngleDown />
              </ListboxButton>

              <ListboxOptions
                className="flex w-[var(--button-width)] flex-col rounded-md bg-white shadow-md"
                anchor={{ to: "bottom", gap: 2, padding: 10 }}
              >
                {options.map((option, index) => (
                  <ListboxOption
                    key={index}
                    className="flex cursor-pointer items-center gap-2.5 px-2.5 py-[5px] transition data-[focus]:bg-light-gray hover:bg-light-gray"
                    value={option.value}
                    data-tooltip={option.tooltip}
                  >
                    <FaCheck
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
};

export default Select;
