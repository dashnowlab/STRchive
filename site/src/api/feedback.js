import { request } from "./";

/** feedback form */
export const submitFeedback = async (title, body) => {
  const headers = new Headers();
  headers.append("Content-Type", "application/json");
  body = JSON.stringify({ title, body });
  const options = { method: "POST", headers, body };
  /** cloud func entry point */
  const url = "https://strchive-feedback-XXXXXXXXXXXX.us-central1.run.app";
  const created = await request(url, options);
  return { link: created.html_url };
};
