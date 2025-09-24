import { Fragment } from "react/jsx-runtime";

/** make english list from array */
export const makeList = (items, Wrapper = "span", joiner = "or") => {
  items = [items].flat();
  items = items.map((item, index) => <Wrapper key={index}>{item}</Wrapper>);
  if (!items?.length) return <></>;
  if (items.length === 1) return items[0];
  const first = items.slice(0, -1);
  const last = items.slice(-1);
  if (items.length === 2)
    return (
      <>
        {first} {joiner} {last}
      </>
    );
  return (
    <>
      {first.map((item, index) => (
        <Fragment key={index}>{item}, </Fragment>
      ))}
      {joiner} {last}
    </>
  );
};
