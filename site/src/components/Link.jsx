import { FaExternalLinkAlt } from "react-icons/fa";
import classes from "./Link.module.css";

const Link = ({
  to,
  newTab = undefined,
  arrow = undefined,
  children,
  ...props
}) => {
  /** whether link is to external site, or page within this site */
  const external = !!to.match(/^(https|http|ftp|mailto)/);

  return (
    <a
      href={to}
      // whether to open in new tab
      target={(newTab ?? external) ? "_blank" : ""}
      {...props}
    >
      {children}
      {/* indicate third-party site with icon  */}
      {(arrow ?? external) && <FaExternalLinkAlt className={classes.icon} />}
    </a>
  );
};

export default Link;
