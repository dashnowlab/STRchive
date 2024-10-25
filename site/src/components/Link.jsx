const Link = ({ to, ...props }) => (
  <a href={to} target={to.startsWith("https://") ? "_blank" : ""} {...props} />
);

export default Link;
