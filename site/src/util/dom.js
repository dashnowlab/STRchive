import { sleep } from "@/util/misc";

/** scroll page so that mouse stays at same position in document */
export const preserveScroll = async (element) => {
  const oldY = element.getBoundingClientRect().top;
  await sleep(0);
  const newY = element.getBoundingClientRect().top;
  window.scrollBy({ top: newY - oldY, behavior: "instant" });
};

/** fit svg viewbox to content */
export const fitViewBox = (element, padding = 0) => {
  let { x, y, width, height } = element.getBBox();
  x -= padding;
  y -= padding;
  width += 2 * padding;
  height += 2 * padding;
  element.setAttribute("viewBox", [x, y, width, height].join(" "));
  return { x, y, width, height };
};
