# SETUP

- Create Google Cloud project
  - Give project descriptive name
  - Create service account with just permissions to build container images
  - Create service account with just permissions to run cloud funcs
  - Update `deploy.sh` vars as appropriate
- Generate GitHub fine-grained personal access token
  - Increase expiration as desired
  - Grant access to only one relevant repo
  - Grant write permissions for issues
  - Store in password vault
- Deploy cloud func
  - Create `/site/.env.local` file, make sure it is NOT TRACKED BY GIT, and add contents `GITHUB_TOKEN=TOKENVALUE`
  - Install `gcloud` CLI tool
  - [Authenticate CLI](https://cloud.google.com/docs/authentication/gcloud)
  - From `/site`, run `./cloud/deploy.sh func-name`
  - Copy function URL into `site/src/api` as appropriate
