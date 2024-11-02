import { sleep } from "@/util/misc";

/** scroll page so that mouse stays at same position in document */
export const preserveScroll = async (element) => {
  const oldY = element.getBoundingClientRect().top;
  await sleep(0);
  const newY = element.getBoundingClientRect().top;
  window.scrollBy({ top: newY - oldY, behavior: "instant" });
};
