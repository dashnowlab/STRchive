import { FaExternalLinkAlt } from "react-icons/fa";
import classes from "./Link.module.css";

const Link = ({ to, noIcon = false, children, ...props }) => {
  const external = !!to.match(/^(https|http|ftp|mailto)/);

  return (
    <a
      href={to}
      target={to.match(/(https?|ftp|mailto):\/\//) ? "_blank" : ""}
      {...props}
    >
      {children}
      {external && !noIcon && <FaExternalLinkAlt className={classes.icon} />}
    </a>
  );
};

export default Link;
