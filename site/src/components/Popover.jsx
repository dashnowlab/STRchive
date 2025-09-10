import {
  Popover as _Popover,
  PopoverButton,
  PopoverPanel,
} from "@headlessui/react";
import Button from "@/components/Button";
import classes from "./Popover.module.css";

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
      className={classes.panel}
      anchor={{ to: "bottom", gap: 2, padding: 10 }}
    >
      {children}
    </PopoverPanel>
  </_Popover>
);

export default Popover;
