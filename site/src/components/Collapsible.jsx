import { FaAngleDown } from "react-icons/fa6";
import classes from "./Collapsible.module.css";

/** collapsible section */
const Collapsible = ({ label, children }) => (
  <details className={classes.details}>
    <summary className={classes.summary}>
      <FaAngleDown className={classes.icon} />
      {label}
    </summary>
    {children}
  </details>
);

export default Collapsible;
