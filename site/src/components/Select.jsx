import { FaAngleDown, FaCheck } from "react-icons/fa6";
import clsx from "clsx";
import {
  Listbox,
  ListboxButton,
  ListboxOption,
  ListboxOptions,
} from "@headlessui/react";
import Button from "@/components/Button";
import classes from "./Select.module.css";

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
  return (
    <label data-tooltip={tooltip}>
      {label && <span>{label}</span>}
      <Listbox
        value={value}
        onChange={onChange}
        multiple={multi}
        defaultValue={defaultValue}
      >
        {({ value }) => {
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
              <ListboxButton as={Button} className={classes.select}>
                <span className={clsx(!selected && classes.label)}>
                  {selected || "No Selection"}
                </span>
                <FaAngleDown />
              </ListboxButton>

              <ListboxOptions
                className={classes.options}
                anchor={{ to: "bottom", gap: 2, padding: 10 }}
              >
                {options.map((option, index) => (
                  <ListboxOption
                    key={index}
                    className={classes.option}
                    value={option.value}
                    data-tooltip={option.tooltip}
                  >
                    <FaCheck
                      className={classes.icon}
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
                    <span className={clsx(!option.label && classes.label)}>
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
