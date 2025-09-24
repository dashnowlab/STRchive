import { mock, request } from "./";

/** create issue in repo. see /cloud/issue */
export const createIssue = async (params) => {
  if (mock) {
    console.debug("Issue", params);
    return { link: "https://fake-link.com" };
  }
  const headers = new Headers();
  headers.append("Content-Type", "application/json");
  body = JSON.stringify(params);
  const options = { method: "POST", headers, body };
  const url = "https://strchive-issue-467600139623.us-central1.run.app";
  const created = await request(url, options);
  return { link: created.html_url };
};
