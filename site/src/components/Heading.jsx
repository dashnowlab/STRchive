import { onlyText } from "react-children-utilities";
import { FaLink } from "react-icons/fa6";
import clsx from "clsx";
import { stripHtml } from "string-strip-html";
import Link from "@/components/Link";
import { slugify } from "@/util/string";
import classes from "./Heading.module.css";

/** heading, with automatic anchors, server or client side */
const Heading = ({ level, id: explicitId, children, className, ...props }) => {
  /** heading level */
  const Component = `h${level}`;

  children =
    children &&
    children.props &&
    children.props.value &&
    children.props.value.toString
      ? /** being rendered by an astro component */
        stripHtml(children.props.value.toString()).result
      : /** being rendered by another react component */
        onlyText(children);

  /** anchor hash */
  const id = slugify(children);

  return (
    <Component
      id={explicitId ?? id}
      className={clsx(className, classes.heading)}
      {...props}
    >
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
