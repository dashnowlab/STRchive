import { mock, request } from "./";

/** create pr in repo. see /cloud/pr */
export const createPR = async (params) => {
  if (mock) {
    console.debug("PR", params);
    return { link: "https://fake-link.com" };
  }
  const headers = new Headers();
  headers.append("Content-Type", "application/json");
  const body = JSON.stringify(params);
  const options = { method: "POST", headers, body };
  const url = "https://strchive-pr-467600139623.us-central1.run.app";
  const created = await request(url, options);
  return { link: created.html_url };
};
