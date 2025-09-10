import { mock, request } from "./";

export const createPR = async ({ branch, title, body, files, labels }) => {
  const headers = new Headers();
  headers.append("Content-Type", "application/json");
  body = JSON.stringify({ branch, title, body, files, labels });
  const options = { method: "POST", headers, body };
  const url = "https://strchive-pr-467600139623.us-central1.run.app";
  if (mock) {
    console.debug("PR", { branch, title, body, files, labels });
    return { link: "https://fake-link.com" };
  }
  const created = await request(url, options);
  return { link: created.html_url };
};
