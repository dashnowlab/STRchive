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

const ComboBox = ({
  multi = false,
  label,
  options = [],
  value = [],
  onChange,
  placeholder,
  tooltip,
}) => {
  const [query, setQuery] = useState("");

  const customOptions = value
    .filter((value) => !options.find((option) => option.value === value))
    .map((value) => ({ value, label: value }));

  options = [...customOptions, ...options];

  const filtered =
    query === ""
      ? options
      : options.filter((option) =>
          [option.value, option.label]
            .join(" ")
            .toLowerCase()
            .includes(query.toLowerCase()),
        );

  let selected = "";
  if (value.length === 0) selected = "None selected";
  else if (value.length === options.length) selected = "All selected";
  else if (value.length === 1)
    selected = options.find((option) => option.value === value[0]).label;
  else if (value.length > 1) selected = `${value.length} selected`;

  return (
    <label data-tooltip={tooltip}>
      {label && <span>{label}</span>}
      <Combobox
        multiple
        value={multi ? value : value.slice(0, 1)}
        onChange={(options) =>
          onChange?.(multi ? options : options.slice(1, 2))
        }
        onClose={() => setQuery("")}
      >
        <div className={classes.select}>
          <span className={clsx("truncate", classes.selected)}>{selected}</span>
          <ComboboxInput
            className={classes.input}
            placeholder={placeholder || "Search"}
            displayValue={(option) => option.label}
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
            <ComboboxOption className={classes.option} value={query}>
              Add "{query}"
            </ComboboxOption>
          )}
          {filtered.map((option, index) => (
            <ComboboxOption
              key={index}
              className={classes.option}
              value={option.value}
            >
              <FaCheck
                className={classes.check}
                style={{ opacity: value.includes(option.value) ? 1 : 0 }}
              />
              <div className={classes["option-label"]}>
                {option.label || <>&nbsp;</>}
              </div>
              {option.extra}
            </ComboboxOption>
          ))}
        </ComboboxOptions>
      </Combobox>
    </label>
  );
};

export default ComboBox;
