import { request } from "./";

export const createIssue = async ({ title, body, labels }) => {
  const headers = new Headers();
  headers.append("Content-Type", "application/json");
  body = JSON.stringify({ title, body, labels });
  const options = { method: "POST", headers, body };
  const url = "https://strchive-issue-467600139623.us-central1.run.app";
  const created = await request(url, options);
  return { link: created.html_url };
};
