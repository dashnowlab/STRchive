import type { ReactNode } from "react";
import { useState } from "react";
import { LuChevronDown } from "react-icons/lu";
import Button from "@/components/Button";
import {
  Combobox,
  ComboboxButton,
  ComboboxInput,
  ComboboxOption,
  ComboboxOptions,
} from "@headlessui/react";
import clsx from "clsx";

type Props = {
  label?: ReactNode;
  options?: { value: string; label?: ReactNode }[];
  value: string;
  onChange: (value: string) => void;
  placeholder?: string;
  tooltip?: string;
};

/** textbox with autocomplete */
export default function ComboBox({
  label,
  options = [],
  value,
  onChange,
  placeholder,
  tooltip,
}: Props) {
  const [query, setQuery] = useState("");

  return (
    <label data-tooltip={tooltip}>
      {label && <span>{label}</span>}
      <Combobox
        value={value}
        onChange={(value) => {
          if (value !== null) onChange(value);
        }}
        immediate
        onClose={() => setQuery("")}
      >
        {({ value }) => {
          /** find option matching current selected value */
          const matchingOption = options.find(
            (option) => option.value === value,
          );

          const _options =
            /** if selected value isn't in options, add entry for it */
            value && !matchingOption
              ? [{ value, label: value }, ...options]
              : [...options];

          /** filter options by search */
          const filtered = _options.filter((option) =>
            [option.value, option.label]
              .join(" ")
              .toLowerCase()
              .includes(query.toLowerCase()),
          );

          return (
            /** https://github.com/tailwindlabs/headlessui/issues/3597 */
            <div className="contents">
              <div className="relative flex min-h-10 w-full min-w-10 items-center gap-2 rounded-md bg-white">
                <span
                  className={clsx(
                    "pointer-events-none absolute inset-0 truncate pr-8",
                  )}
                >
                  {matchingOption?.label || value}
                </span>
                <ComboboxInput
                  className="h-10 w-full px-2 pr-8 text-current"
                  placeholder={placeholder || "Search"}
                  onChange={(event) => setQuery(event.target.value)}
                />
                <ComboboxButton
                  as={Button}
                  className="absolute right-0 grid size-10 place-items-center transition"
                >
                  <LuChevronDown />
                </ComboboxButton>
              </div>
              <ComboboxOptions
                className="flex w-(--input-width) flex-col rounded-md bg-white shadow-md"
                anchor={{ to: "bottom", gap: 2, padding: 10 }}
              >
                {query.length > 0 && (
                  <ComboboxOption
                    value={query}
                    className="flex cursor-pointer items-center px-4 py-2 transition hover:bg-light-gray data-focus:bg-light-gray"
                  >
                    "{query}"
                  </ComboboxOption>
                )}
                {filtered.map((option, index) => (
                  <ComboboxOption
                    key={index}
                    className="flex cursor-pointer items-center px-4 py-2 transition hover:bg-light-gray data-focus:bg-light-gray"
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
}
