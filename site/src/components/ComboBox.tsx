import type { ReactNode } from "react";
import { useState } from "react";
import Help from "@/components/Help";
import { Combobox as _Combobox } from "@base-ui/react";

type Option = {
  value: string;
  label?: ReactNode;
  special?: boolean;
};

type Props = {
  label?: ReactNode;
  options?: Option[];
  value: string;
  onChange: (value: string) => void;
  placeholder?: string;
  help?: ReactNode;
};

/** textbox with autocomplete */
export default function ComboBox({
  label,
  options = [],
  value,
  onChange,
  placeholder,
  help,
}: Props) {
  const [query, setQuery] = useState("");

  /** custom entry */
  if (query.trim())
    options = [...options, { value: query, label: query, special: true }];

  return (
    <_Combobox.Root
      items={options}
      value={value}
      onValueChange={(value) => {
        if (value !== null) onChange(value);
      }}
      inputValue={query}
      onInputValueChange={setQuery}
      onOpenChange={(open) => {
        if (!open) setQuery("");
      }}
    >
      <label>
        <span>
          {label}
          <Help>{help}</Help>
        </span>
        <_Combobox.InputGroup className="flex min-h-10 w-full min-w-10 rounded-md bg-white ring-2 ring-inset">
          <_Combobox.Input
            className="grow px-4 py-2"
            placeholder={placeholder}
          />
        </_Combobox.InputGroup>
      </label>

      <_Combobox.Portal>
        <_Combobox.Positioner collisionPadding={20}>
          <_Combobox.Popup>
            <_Combobox.List className="flex max-h-(--available-height) w-(--anchor-width) flex-col overflow-y-auto rounded-md bg-white shadow-md">
              {(item: Option) => (
                <_Combobox.Item
                  key={item.value}
                  value={item.value}
                  className="flex cursor-pointer items-center gap-2 px-4 py-2 transition hover:bg-light-gray data-highlighted:bg-light-gray"
                >
                  {item.special ? `"${item.label}"` : item.label}
                </_Combobox.Item>
              )}
            </_Combobox.List>
          </_Combobox.Popup>
        </_Combobox.Positioner>
      </_Combobox.Portal>
    </_Combobox.Root>
  );
}
