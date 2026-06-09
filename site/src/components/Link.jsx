import { FaArrowUpRightFromSquare } from "react-icons/fa6";

const Link = ({
  to,
  newTab = false,
  arrow = false,
  children,
  ...props
}) => {
  /** whether link is to external site, or page within this site */
  const external = !!to.match(/^(https|http|ftp|mailto)/);

  const Component = to.trim() ? "a" : "span";

  return (
    <Component
      href={to || null}
      // whether to open in new tab
      target={(newTab ?? external) ? "_blank" : ""}
      {...props}
    >
      {children}
      {/* indicate third-party site with icon  */}
      {(arrow ?? external) && !children && (
        <FaArrowUpRightFromSquare className="relative top-[0.1em] ml-[0.1em] scale-75" />
      )}
    </Component>
  );
};

export default Link;
