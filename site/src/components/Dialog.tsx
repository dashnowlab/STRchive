import type { ReactElement, ReactNode } from "react";
import { Dialog as _Dialog } from "@base-ui/react";
import { IconX } from "@tabler/icons-react";

type Props = {
  trigger: ReactElement<{ onClick?: () => void }>;
  title: string;
  children: ReactNode;
};

export default function Dialog({ trigger, title, children }: Props) {
  return (
    <_Dialog.Root>
      <_Dialog.Trigger render={trigger} />
      <_Dialog.Portal className="z-10">
        <_Dialog.Backdrop className="fixed inset-0 bg-black/75" />
        <_Dialog.Popup className="pointer-events-none fixed inset-0 grid h-screen w-screen place-items-center p-10">
          <div className="pointer-events-auto flex max-h-full max-w-full flex-col overflow-hidden rounded-md bg-white text-black shadow-md">
            <div className="flex items-center p-4 shadow-md">
              <strong className="grow text-lg">{title}</strong>
              <_Dialog.Close>
                <IconX />
              </_Dialog.Close>
            </div>

            <div className="flex flex-col justify-center-safe gap-8 overflow-y-auto p-8">
              {children}
            </div>
          </div>
        </_Dialog.Popup>
      </_Dialog.Portal>
    </_Dialog.Root>
  );
}
