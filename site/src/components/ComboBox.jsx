import { useState } from "react";
import { FaAngleDown, FaCheck } from "react-icons/fa6";
import clsx from "clsx";
import {
  Combobox,
  ComboboxButton,
  ComboboxInput,
  ComboboxOption,
  ComboboxOptions,
} from "@headlessui/react";
import Button from "@/components/Button";
import classes from "./ComboBox.module.css";

/** textbox with autocomplete */
const ComboBox = ({
  label,
  options = [],
  value,
  onChange,
  placeholder,
  tooltip,
}) => {
  const [query, setQuery] = useState("");

  return (
    <label data-tooltip={tooltip}>
      {label && <span>{label}</span>}
      <Combobox
        value={value}
        onChange={onChange}
        immediate
        onClose={() => setQuery("")}
      >
        {({ value }) => {
          /** find option matching current selected value */
          const matchingOption = options.find(
            (option) => option.value === value,
          );

          if (value && !matchingOption)
            options = [{ value, label: value }, ...options];

          const filtered = options.filter((option) =>
            [option.value, option.label]
              .join(" ")
              .toLowerCase()
              .includes(query.toLowerCase()),
          );

          return (
            /** https://github.com/tailwindlabs/headlessui/issues/3597 */
            <div style={{ display: "contents" }}>
              <div className={classes.select}>
                <span className={clsx("truncate", classes.selected)}>
                  {matchingOption?.label || value}
                </span>
                <ComboboxInput
                  className={classes.input}
                  placeholder={placeholder || "Search"}
                  onChange={(event) => setQuery(event.target.value)}
                />
                <ComboboxButton as={Button} className={classes.button}>
                  <FaAngleDown />
                </ComboboxButton>
              </div>
              <ComboboxOptions
                className={classes.options}
                anchor={{ to: "bottom", gap: 2, padding: 10 }}
              >
                {query.length > 0 && (
                  <ComboboxOption value={query} className={classes.option}>
                    "{query}"
                  </ComboboxOption>
                )}
                {filtered.map((option, index) => (
                  <ComboboxOption
                    key={index}
                    className={classes.option}
                    value={option.value}
                  >
                    {option.label || <>&nbsp;</>}
                  </ComboboxOption>
                ))}
              </ComboboxOptions>
            </div>
          );
        }}
      </Combobox>
    </label>
  );
};

export default ComboBox;
