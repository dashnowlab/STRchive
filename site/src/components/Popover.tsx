import type { ReactElement, ReactNode } from "react";
import type { PopoverRootActions } from "@base-ui/react";
import { useRef } from "react";
import { sleep } from "@/util/misc";
import { Popover as _Popover } from "@base-ui/react";

type Props = {
  /** popover content */
  content: ReactNode;
  /** popover trigger */
  children?: ReactElement;
  /** whether popover should open on hover as well */
  hover?: boolean;
  /** whether trigger ultimately renders <button> element */
  button?: boolean;
};

/** button that triggers floating panel of arbitrary content */
export default function Popover({
  content,
  children,
  hover = true,
  button = true,
}: Props) {
  const popupRef = useRef<HTMLDivElement>(null);
  const actionsRef = useRef<PopoverRootActions>(null);

  return (
    <_Popover.Root
      actionsRef={actionsRef}
      onOpenChange={async (open) => {
        await sleep();
        if (open && !popupRef.current?.innerText.trim())
          actionsRef.current?.close();
      }}
    >
      <_Popover.Trigger
        render={children}
        openOnHover={hover}
        delay={100}
        nativeButton={button}
      />
      <_Popover.Portal className="z-30">
        <_Popover.Backdrop />
        <_Popover.Positioner side="top" sideOffset={12}>
          <_Popover.Popup ref={popupRef} className="max-w-100">
            <_Popover.Arrow className="z-10 size-0 data-[side=bottom]:rotate-45 data-[side=left]:left-full data-[side=left]:rotate-135 data-[side=right]:rotate-315 data-[side=top]:top-full data-[side=top]:rotate-225">
              <div className="size-3 -translate-1/2 bg-white shadow-md [clip-path:polygon(-100%_-100%,200%_-100%,-100%_200%)]" />
            </_Popover.Arrow>
            <div className="relative flex flex-col gap-1 overflow-hidden rounded-md bg-white px-4 py-2 shadow-md">
              {content}
            </div>
          </_Popover.Popup>
        </_Popover.Positioner>
      </_Popover.Portal>
    </_Popover.Root>
  );
}
