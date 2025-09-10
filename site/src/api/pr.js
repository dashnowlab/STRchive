import { mock, request } from "./";

/** create pr in repo. see /cloud/pr */
export const createPR = async (params) => {
  const headers = new Headers();
  headers.append("Content-Type", "application/json");
  body = JSON.stringify(params);
  const options = { method: "POST", headers, body };
  const url = "https://strchive-pr-467600139623.us-central1.run.app";
  if (mock) {
    console.debug("PR", params);
    return { link: "https://fake-link.com" };
  }
  const created = await request(url, options);
  return { link: created.html_url };
};
