const functions = require("@google-cloud/functions-framework");
const { Octokit } = require("octokit");

/** settings */
const owner = "dashnowlab";
const repo = "STRchive";
const labels = ["locus-edit"];
const auth = process.env.GITHUB_TOKEN;

/** endpoint */
functions.http("entrypoint", async (request, response) => {
  response.set("Access-Control-Allow-Origin", "*");
  response.set("Access-Control-Allow-Methods", "OPTIONS, POST");
  response.set("Access-Control-Allow-Headers", "Accept, Content-Type");

  /** handle pre-flight */
  if (request.method == "OPTIONS") return response.status(200).send();

  let octokit;
  try {
    /** auth with github https://octokit.github.io/rest.js/v22/#authentication */
    octokit = new Octokit({ auth });
  } catch (error) {
    /** https://github.com/octokit/request-error.js */
    return response
      .status(error.status)
      .send(`Couldn't authenticate with GitHub: ${error.message}`);
  }

  /** get params */
  const { branch, title, body, files = [] } = params || {};

  /** check for missing params */
  const missing = [];
  if (!branch) missing.push("branch");
  if (!title) missing.push("title");
  if (!body) missing.push("body");
  if (!files.length) missing.push("files");
  files.forEach(({ path, content }, index) => {
    if (!path) missing.push(`file path ${index}`);
    if (!content) missing.push(`file content ${index}`);
  });
  if (missing.length)
    return response.status(400).send(`Missing params: ${missing.join(", ")}`);

  /** get main branch */
  let main;
  try {
    /** https://octokit.github.io/rest.js/v22/#git-get-ref */
    main = (
      await octokit.rest.git.getRef({
        owner,
        repo,
        ref: "heads/main",
      })
    ).data.object.sha;
  } catch (error) {
    return response
      .status(error.status)
      .send(`Couldn't get main branch: ${error.message}`);
  }

  /** check if branch already exists */
  for (let tries = 1; tries <= 10; tries++) {
    try {
      octokit.rest.repos.getBranch({
        owner,
        repo,
        branch,
      });
      /** branch exists, add suffix to name to avoid duplicate */
      branch = branch.replace(/(-?\d*)$/, `-${tries}`);
    } catch (error) {
      /** branch doesn't exist */
      if (error.status === 404) break;
      return response
        .status(error.status)
        .send(`Couldn't get branch ${branch}: ${error.message}`);
    }
  }

  /** create new branch */
  try {
    /** https://octokit.github.io/rest.js/v22/#git-create-ref */
    await octokit.rest.git.createRef({
      owner,
      repo,
      ref: `refs/heads/${branch}`,
      /** create branch off of main */
      sha: main,
    });
  } catch (error) {
    return response
      .status(error.status)
      .send(`Couldn't create new branch ${branch}: ${error.message}`);
  }

  /** update files */
  for (const { path, content } of files) {
    /** get existing file */
    let existing;
    try {
      /** https://octokit.github.io/rest.js/v22/#repos-get-content */
      existing = (
        await octokit.rest.repos.getContent({
          owner,
          repo,
          branch,
          path,
        })
      ).data.sha;
    } catch (error) {
      console.warn(`Couldn't get existing file ${path}: ${error.message}`);
    }

    /** create/update file */
    try {
      /** https://octokit.github.io/rest.js/v22/#repos-create-or-update-file-contents */
      await octokit.rest.repos.createOrUpdateFileContents({
        owner,
        repo,
        branch,
        path,
        message: title,
        content: Buffer.from(content.replace(/\n?$/m, "\n")).toString("base64"),
        sha: existing,
      });
    } catch (error) {
      return response
        .status(error.status)
        .send(`Couldn't update file ${path}: ${error.message}`);
    }
  }

  let pr;
  /** create pull request */
  try {
    /** https://octokit.github.io/rest.js/v22/#pulls-create */
    pr = await octokit.rest.pulls.update({
      owner,
      repo,
      base: "main",
      head: branch,
      title,
      body,
      maintainer_can_modify: true,
      draft: true,
    });
  } catch (error) {
    return response
      .status(error.status)
      .send(`Couldn't create PR ${branch}: ${error.message}`);
  }

  /** add labels */
  try {
    await octokit.rest.issues.addLabels({
      owner,
      repo,
      issue_number: pr.data.number,
      labels,
    });
  } catch (error) {
    return response
      .status(error.status)
      .send(`Couldn't add label: ${error.message}`);
  }

  return pr;
});
