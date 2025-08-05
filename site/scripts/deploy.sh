#!/bin/bash
# run from /site

script=$1
project=gap-som-dbmi-strchive-app-60
function=strchive-$script
region=us-central1

source .env.local
github=$GITHUB_TOKEN

gcloud run deploy $function \
  --source ./scripts/$script \
  --function entrypoint \
  --base-image nodejs22 \
  --region $region \
  --project $project \
  --build-service-account projects/$project/serviceAccounts/cloudrunfunction-builder@$project.iam.gserviceaccount.com \
  --service-account cloudrunfunction-runtime@$project.iam.gserviceaccount.com \
  --max-instances 2 \
  --no-invoker-iam-check \
  --set-env-vars GITHUB_TOKEN=$github
