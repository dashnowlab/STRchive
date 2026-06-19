import type { ReactElement, ReactNode } from "react";
import { cloneElement, createContext, useRef, useState } from "react";
import { LuX } from "react-icons/lu";
import Button from "@/components/Button";
import {
  useClickOutside,
  useEventListener,
  useScrollLock,
} from "@reactuses/core";

export const DialogContext = createContext<{ isOpen: boolean }>({
  isOpen: false,
});

type Props = {
  trigger: ReactElement<{ onClick?: () => void }>;
  title: string;
  children: ReactNode;
};

export default function Dialog({ trigger, title, children }: Props) {
  const dialogRef = useRef<HTMLDialogElement>(null);
  const boxRef = useRef(null);
  const [isOpen, setIsOpen] = useState(false);
  const [, setLocked] = useScrollLock(document.body);

  const open = () => {
    setIsOpen(true);
    dialogRef.current?.showModal();
    setLocked(true);
  };

  const close = () => {
    setIsOpen(false);
    dialogRef.current?.close();
    setLocked(false);
  };

  useEventListener("close", () => setIsOpen(false), dialogRef);
  useClickOutside(boxRef, close);

  return (
    <>
      {/* eslint-disable-next-line -- erroneous, ref is being accessed on click not render */}
      {cloneElement(trigger, { onClick: open })}

      <dialog
        ref={dialogRef}
        className="fixed inset-0 hidden h-dvh max-h-none w-dvw max-w-none place-items-center bg-transparent backdrop:opacity-0 open:grid"
      >
        <div
          ref={boxRef}
          className="flex max-h-[calc(100dvh-(--spacing(20)))] max-w-[calc(100dvw-(--spacing(20)))] flex-col rounded-md bg-white text-black shadow-md"
        >
          <div className="flex items-center p-4 shadow-md">
            <strong className="grow text-lg">{title}</strong>
            <Button onClick={close} autoFocus>
              <LuX />
            </Button>
          </div>

          <div className="flex flex-col justify-center-safe gap-8 overflow-y-auto overscroll-none p-8">
            <DialogContext.Provider value={{ isOpen }}>
              {children}
            </DialogContext.Provider>
          </div>
        </div>
      </dialog>
    </>
  );
}
