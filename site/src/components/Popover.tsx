import type { ReactNode } from "react";
import Button from "@/components/Button";
import { Popover as _Popover } from "@base-ui/react";
import { IconChevronDown } from "@tabler/icons-react";

type Props = {
  label?: ReactNode;
  button: ReactNode;
  children: ReactNode;
  tooltip?: string;
};

/** button that triggers floating panel of arbitrary content */
export default function Popover({ label, button, children, tooltip }: Props) {
  return (
    <_Popover.Root>
      <label data-tooltip={tooltip}>
        {label}
        <_Popover.Trigger render={<Button design="plain" />}>
          {button}
          <IconChevronDown />
        </_Popover.Trigger>
      </label>

      <_Popover.Portal>
        <_Popover.Backdrop />
        <_Popover.Positioner collisionPadding={20}>
          <_Popover.Popup>
            <_Popover.Viewport>
              <div className="flex flex-col gap-1 rounded-md bg-white p-2 shadow-md">
                {children}
              </div>
            </_Popover.Viewport>
          </_Popover.Popup>
        </_Popover.Positioner>
      </_Popover.Portal>
    </_Popover.Root>
  );
}
