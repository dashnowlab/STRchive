import type { JSX } from "react";
import { Fragment } from "react/jsx-runtime";

/** make english list from array */
export const makeList = (
  items: unknown,
  Wrapper: keyof JSX.IntrinsicElements = "span",
  joiner = "or",
) => {
  let list = [items].flat();
  list = list.map((item, index) => (
    <Wrapper key={index}>{String(item)}</Wrapper>
  ));
  if (!list?.length) return <></>;
  if (list.length === 1) return list[0];
  const first = list.slice(0, -1);
  const last = list.slice(-1);
  if (list.length === 2)
    return (
      <>
        {first} {joiner} {last}
      </>
    );
  return (
    <>
      {first.map((item, index) => (
        <Fragment key={index}>{String(item)}, </Fragment>
      ))}
      {joiner} {last}
    </>
  );
};
