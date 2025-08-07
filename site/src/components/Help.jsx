import { FaRegQuestionCircle } from "react-icons/fa";
import classes from "./Help.module.css";

const Help = ({ children }) => (
  <FaRegQuestionCircle className={classes.help} data-tooltip={children} />
);

export default Help;
