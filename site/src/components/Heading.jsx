import { onlyText } from "react-children-utilities";
import { FaLink } from "react-icons/fa6";
import clsx from "clsx";
import { stripHtml } from "string-strip-html";
import Link from "@/components/Link";
import { slugify } from "@/util/string";

/** heading, with automatic anchors, server or client side */
const Heading = ({ level, id: explicitId = "", children, className = "", ...props }) => {
  /** heading level */
  const Component = `h${level}`;

  /** anchor hash */
  const id = slugify(
    children &&
      children.props &&
      children.props.value &&
      children.props.value.toString
      ? /** being rendered by an astro component */
        stripHtml(children.props.value.toString()).result
      : /** being rendered by another react component */
        onlyText(children),
  );

  return (
    <Component
      id={explicitId ?? id}
      className={clsx(
        className,
        "group flex items-center gap-2 wrap-anywhere [&>svg:first-child]:opacity-50",
      )}
      {...props}
    >
      {children}
      {id && (
        <Link
          to={"#" + id}
          className="w-0 translate-x-0 whitespace-nowrap text-primary opacity-0 transition group-hover:opacity-100 focus:opacity-100"
          aria-label="Heading link"
        >
          <FaLink />
        </Link>
      )}
    </Component>
  );
};

export default Heading;
