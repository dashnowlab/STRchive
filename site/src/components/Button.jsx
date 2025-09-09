import clsx from "clsx";
import Link from "@/components/Link";
import classes from "./Button.module.css";

/** looks like a button, and either does something or goes somewhere */
const Button = ({ design, className, ...props }) => {
  const Component = props.to ? Link : "button";
  return (
    <Component
      type="button"
      className={clsx(classes.button, className, classes[design])}
      {...props}
    />
  );
};

export default Button;
