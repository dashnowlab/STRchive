const functions = require("@google-cloud/functions-framework");
const { Octokit } = require("octokit");

/** settings */
const owner = "dashnowlab";
const repo = "STRchive";
const labels = ["feedback"];
const auth = process.env.GITHUB_TOKEN;

/** entry point */
functions.http("entrypoint", async (request, response) => {
  response.set("Access-Control-Allow-Origin", "*");
  response.set("Access-Control-Allow-Methods", "OPTIONS, POST");
  response.set("Access-Control-Allow-Headers", "Accept, Content-Type");

  /** handle pre-flight */
  if (request.method == "OPTIONS") return response.status(200).send();

  let octokit;
  try {
    /** auth with github */
    octokit = new Octokit({ auth });
  } catch (error) {
    /** https://github.com/octokit/request-error.js */
    return response.status(error.status).send(error.message);
  }

  /** get params */
  const { title, body } = request.body || {};

  /** check for missing params */
  if (!title) return response.status(400).send("Missing title");
  if (!body) return response.status(400).send("Missing body");

  /** create issue */
  try {
    const { data } = await octokit.rest.issues.create({
      owner,
      repo,
      title,
      body,
      labels,
    });

    /** return created issue */
    return response.status(200).json(data);
  } catch (error) {
    /** https://github.com/octokit/request-error.js */
    return response.status(error.status).send(error.message);
  }
});
