import { UAParser } from "ua-parser-js";

export const userAgent = new UAParser(window.navigator.userAgent).getResult();

export const isFirefox = userAgent.browser.name
  ?.toLowerCase()
  .includes("firefox");

export const isSafari = userAgent.browser.name
  ?.toLowerCase()
  .includes("safari");
