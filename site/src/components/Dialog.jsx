import { cloneElement, useRef, useState } from "react";
import { FaXmark } from "react-icons/fa6";
import {
  useClickOutside,
  useEventListener,
  useScrollLock,
} from "@reactuses/core";
import classes from "./Dialog.module.css";

const Dialog = ({ trigger, title, children }) => {
  const ref = useRef(null);
  const [, setIsOpen] = useState(false);

  const [, setLock] = useScrollLock(document.body);

  const open = () => {
    setIsOpen(true);
    ref.current?.showModal();
    setLock(true);
  };

  const close = () => {
    setIsOpen(false);
    ref.current?.close();
    setLock(false);
  };

  useEventListener("close", () => setIsOpen(false), ref);
  useClickOutside(ref, close);

  return (
    <>
      {cloneElement(trigger, { onClick: open })}

      <dialog ref={ref} className={classes.dialog}>
        <div className={classes.title}>
          <span>{title}</span>
          <button onClick={close} autoFocus>
            <FaXmark />
          </button>
        </div>

        <div className={classes.content}>{children}</div>
      </dialog>
    </>
  );
};

export default Dialog;
