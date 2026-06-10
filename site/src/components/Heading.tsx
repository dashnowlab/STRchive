import type { ComponentProps, ReactNode } from "react";
import { onlyText } from "react-children-utilities";
import Link from "@/components/Link";
import { slugify } from "@/util/string";
import clsx from "clsx";
import { stripHtml } from "string-strip-html";

type Level = 1 | 2 | 3 | 4;

type Props = {
  level: Level;
  id?: string;
} & ComponentProps<"h1" | "h2" | "h3" | "h4">;

/** heading, with automatic anchors, server or client side */
function Heading({
  level,
  id = "",
  children,
  className = "",
  ...props
}: Props) {
  /** heading level */
  const Component: `h${Level}` = `h${level}`;

  const _id = id || slugify(getTextContent(children));

  return (
    <Component
      id={_id}
      className={clsx("flex items-center gap-2 wrap-anywhere", className)}
      {...props}
    >
      <Link to={`#${_id}`} className="contents text-current no-underline">
        {children}
      </Link>
    </Component>
  );
}

export function H1(props: Omit<Props, "level">) {
  return <Heading level={1} {...props} />;
}

export function H2(props: Omit<Props, "level">) {
  return <Heading level={2} {...props} />;
}

export function H3(props: Omit<Props, "level">) {
  return <Heading level={3} {...props} />;
}

export function H4(props: Omit<Props, "level">) {
  return <Heading level={4} {...props} />;
}

/** get text content of children */
const getTextContent = (children: ReactNode) => {
  try {
    /** being rendered by astro component */
    // eslint-disable-next-line -- too onerous to do all nested type checks/guards, so just try/catch
    return stripHtml((children as any).props.value.toString()).result;
  } catch {}
  try {
    /** being rendered by another react component */
    return onlyText(children);
  } catch {}
  return "";
};
