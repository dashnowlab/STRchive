import {
  Popover as _Popover,
  PopoverButton,
  PopoverPanel,
} from "@headlessui/react";
import Button from "@/components/Button";

/** button that triggers floating panel of arbitrary content */
const Popover = ({ label, button, children }) => (
  <_Popover>
    <label>
      {label}
      <PopoverButton as={Button} design="plain">
        {button}
      </PopoverButton>
    </label>

    <PopoverPanel
      className="flex min-w-[var(--button-width)] flex-col gap-[5px] rounded-md bg-white p-2.5 shadow-md"
      anchor={{ to: "bottom", gap: 2, padding: 10 }}
    >
      {children}
    </PopoverPanel>
  </_Popover>
);

export default Popover;
