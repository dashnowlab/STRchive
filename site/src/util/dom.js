import { debounce } from "lodash-es";
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

/** find index of first element "in view". model behavior off of wikiwand.com. */
export const firstInView = (elements) => {
  for (const element of elements.reverse())
    if (element.getBoundingClientRect()?.top < 100) return element;
};

/** highlight an element */
export const spotlight = async (selector) => {
  const element = document.querySelector(selector);
  if (!element) return;
  element.focus();
  element.scrollIntoView({ block: "center", behavior: "smooth" });
  await waitForStop(window, "scroll");
  element.animate(
    [
      { boxShadow: "0 0 0 10px transparent, 0 0 0 9999px transparent" },
      { boxShadow: "0 0 0 10px var(--primary), 0 0 0 9999px #000000cc" },
      { boxShadow: "0 0 0 10px transparent, 0 0 0 9999px transparent" },
    ],
    { duration: 2000, easing: "ease-in-out" },
  );
};

/** wait for event to stop being emitted */
const waitForStop = async (target, event, wait = 100) =>
  new Promise(async (resolve) => {
    const listener = debounce(() => {
      target.removeEventListener(event, listener);
      resolve();
    }, wait);
    listener();
    target.addEventListener(event, listener);
  });
