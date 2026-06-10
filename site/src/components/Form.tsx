import type { SubmitEvent } from "react";
import { useEventListener } from "@reactuses/core";

type Props = {
  onSubmit: (
    data: Record<string, string | number>,
    event: SubmitEvent<HTMLFormElement>,
  ) => void;
};

/** form wrapper around set of fields */
export default function Form({ onSubmit, ...props }: Props) {
  usePreventImplicitSubmit();

  return (
    <form
      className="contents"
      onSubmit={(event) => {
        /** prevent page navigation */
        event.preventDefault();

        /** get data from form */
        const form = event.currentTarget;
        const formData = new FormData(form);

        const data: Record<string, string | number> = {};

        /** transform form data into nicer format */
        for (const [key, value] of formData.entries()) {
          /** if we can parse as number, do it */
          if (
            typeof value === "string" &&
            value.trim() &&
            !Number.isNaN(Number(value))
          )
            data[key] = Number(value);
          else
            /** raw (string) value */
            data[key] = String(value);
        }

        /** call callback with data */
        onSubmit(data, event);
      }}
      {...props}
    />
  );
}

/** prevent implicit form submit */
const usePreventImplicitSubmit = () => {
  useEventListener("keydown", (event) => {
    const { key, target } = event;
    /** only on enter key */
    if (key !== "Enter") return;
    /** only on elements */
    if (!(target instanceof Element)) return;
    /** only on inputs */
    if (!target.matches("input")) return;
    /** allow enter press on submit button to still work */
    if (target.matches("[type='submit']")) return;
    /** prevent submit */
    event.preventDefault();
  });
};
