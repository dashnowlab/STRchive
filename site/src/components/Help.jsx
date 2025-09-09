import { FaRegCircleQuestion } from "react-icons/fa6";
import classes from "./Help.module.css";

const Help = ({ children }) => (
  <FaRegCircleQuestion className={classes.help} data-tooltip={children} />
);

export default Help;
