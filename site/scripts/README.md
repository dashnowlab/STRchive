Cloud func setup:

- Generate GitHub fine-grained personal access token
  - Increase expiration as desired
  - Grant access to only one relevant repo
  - Grant write permissions for issues
  - Store in password vault
- Manually deploy cloud func; no auto-deploy (yet)
  - Go to appropriate project in Google Cloud Console
  - Create a new inline cloud run func
  - Give it descriptive name/URL that includes project name
  - Use Node runtime
  - Allow unauthenticated requests
  - Under container variables & secrets, create a new environment variable with name GITHUB_TOKEN and value of token
  - Set other settings as desired
  - Copy/paste script and package files in this folder into source
  - Rename "function entry point" to match http function name in script file
  - Copy endpoint URL into site/src/api as appropriate
