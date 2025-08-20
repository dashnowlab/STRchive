import { FaLink } from "react-icons/fa6";
import { kebabCase } from "lodash-es";
import Link from "@/components/Link";
import classes from "./Heading.module.css";

const Heading = ({ level, id, children }) => {
  if (!id && typeof children === "string") id = kebabCase(children);
  const Component = `h${level}`;

  return (
    <Component id={id} className={classes.heading}>
      {children}
      {id && (
        <Link
          to={"#" + id}
          className={classes.anchor}
          aria-label="Heading link"
        >
          <FaLink />
        </Link>
      )}
    </Component>
  );
};

export default Heading;
