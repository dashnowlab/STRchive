const Link = ({ to, ...props }) => (
  <a
    href={to}
    target={to.match(/(https?|ftp|mailto):\/\//) ? "_blank" : ""}
    {...props}
  />
);

export default Link;
