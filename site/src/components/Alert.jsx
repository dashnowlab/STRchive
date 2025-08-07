import {
  FaCircleCheck,
  FaCircleExclamation,
  FaCircleInfo,
  FaTriangleExclamation,
} from "react-icons/fa6";
import clsx from "clsx";
import Loading from "@/assets/loading.svg?react";
import classes from "./Alert.module.css";

/** available categories of marks and associated styles */
export const types = {
  info: { color: "var(--tertiary)", icon: <FaCircleInfo /> },
  loading: { color: "var(--gray)", icon: <Loading /> },
  success: { color: "var(--primary)", icon: <FaCircleCheck /> },
  warning: { color: "var(--secondary)", icon: <FaCircleExclamation /> },
  error: { color: "var(--secondary)", icon: <FaTriangleExclamation /> },
};

/** colored box with icon and text */
const Alert = ({
  type = "info",
  icon,
  className = "",
  style = {},
  children,
  ...props
}) => (
  <div
    className={clsx(className, classes.alert)}
    style={{ "--color": types[type].color, ...style }}
    {...props}
  >
    {icon ?? types[type].icon}
    <div>{children}</div>
  </div>
);

export default Alert;
