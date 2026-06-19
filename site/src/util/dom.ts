import { sleep } from "@/util/misc";
import { debounce } from "lodash-es";

/** scroll page so that mouse stays at same position in document */
export const preserveScroll = async (element: HTMLElement) => {
  const oldY = element.getBoundingClientRect().top;
  await sleep();
  const newY = element.getBoundingClientRect().top;
  window.scrollBy({ top: newY - oldY, behavior: "instant" });
};

/** fit svg viewbox to content */
export const fitViewBox = (element: SVGSVGElement, padding = 0) => {
  let { x, y, width, height } = element.getBBox();
  x -= padding;
  y -= padding;
  width += 2 * padding;
  height += 2 * padding;
  element.setAttribute("viewBox", [x, y, width, height].join(" "));
  return { x, y, width, height };
};

/** find index of first element "in view". model behavior off of wikiwand.com. */
export const firstInView = (elements: Element[]) => {
  for (const element of elements.reverse())
    if (element.getBoundingClientRect()?.top < 100) return element;
};

/** highlight an element */
export const spotlight = async (element: Element) => {
  if (element instanceof HTMLElement) element.focus();
  element.scrollIntoView({ block: "center", behavior: "smooth" });
  await waitForStop(window, "scroll");
  element.animate(
    [
      { boxShadow: "0 0 0 9999px transparent" },
      { boxShadow: "0 0 0 9999px #000000cc" },
      { boxShadow: "0 0 0 9999px transparent" },
    ],
    { duration: 2000, easing: "ease-in-out" },
  );
};

/** wait for event to stop being emitted */
const waitForStop = async (target: EventTarget, event: string, wait = 100) =>
  new Promise(async (resolve) => {
    const listener = debounce(() => {
      target.removeEventListener(event, listener);
      resolve(true);
    }, wait);
    listener();
    target.addEventListener(event, listener);
  });
