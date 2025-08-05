Cloud func setup:

- Generate GitHub fine-grained personal access token
  - Increase expiration as desired
  - Grant access to only one relevant repo
  - Grant write permissions for issues
  - Store in password vault
- Manually deploy cloud func
  - Create `/site/.env.local` file, make sure it is NOT TRACKED BY GIT, and make its contents `GITHUB_TOKEN=TOKENVALUE123`
  - Install `gcloud` CLI tool
  - [Authenticate CLI](https://cloud.google.com/docs/authentication/gcloud)
  - From /site, run `./scripts/deploy.sh script-name`
  - Copy function URL into site/src/api as appropriate
