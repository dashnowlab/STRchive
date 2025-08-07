import { cloneElement, createContext, useRef, useState } from "react";
import { FaXmark } from "react-icons/fa6";
import clsx from "clsx";
import { useClickOutside, useEventListener } from "@reactuses/core";
import Button from "./Button";
import classes from "./Dialog.module.css";

export const DialogContext = createContext({});

const Dialog = ({ trigger, title, children }) => {
  const dialogRef = useRef(null);
  const boxRef = useRef(null);
  const [isOpen, setIsOpen] = useState(false);

  const open = () => {
    setIsOpen(true);
    dialogRef.current?.showModal();
  };

  const close = () => {
    setIsOpen(false);
    dialogRef.current?.close();
  };

  useEventListener("close", () => setIsOpen(false), dialogRef);
  useClickOutside(boxRef, close);

  return (
    <>
      {cloneElement(trigger, { onClick: open })}

      <dialog ref={dialogRef} className={clsx("col", classes.dialog)}>
        <div ref={boxRef} className={clsx("col", classes.box)}>
          <div className={classes.title}>
            <span>{title}</span>
            <Button onClick={close} autoFocus>
              <FaXmark />
            </Button>
          </div>

          <div className={clsx("col", classes.body)}>
            <DialogContext.Provider value={{ isOpen }}>
              {children}
            </DialogContext.Provider>
          </div>
        </div>
      </dialog>
    </>
  );
};

export default Dialog;
