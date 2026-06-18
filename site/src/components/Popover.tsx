import type { ReactNode } from "react";
import { LuChevronDown } from "react-icons/lu";
import Button from "@/components/Button";
import {
  Popover as _Popover,
  PopoverButton,
  PopoverPanel,
} from "@headlessui/react";

type Props = {
  label: ReactNode;
  button: ReactNode;
  children: ReactNode;
};

/** button that triggers floating panel of arbitrary content */
export default function Popover({ label, button, children }: Props) {
  return (
    <_Popover>
      <label>
        {label}
        <PopoverButton as={Button} design="plain">
          {button}
          <LuChevronDown />
        </PopoverButton>
      </label>

      <PopoverPanel
        className="flex min-w-(--button-width) flex-col gap-1 rounded-md bg-white p-2 shadow-md"
        anchor={{ to: "bottom", gap: 2, padding: 10 }}
      >
        {children}
      </PopoverPanel>
    </_Popover>
  );
}
