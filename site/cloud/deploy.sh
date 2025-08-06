#!/bin/bash
# run from /site

name=$1
project=gap-som-dbmi-strchive-app-60
function=strchive-$name
build=cloudrunfunction-builder
run=cloudrunfunction-runtime
region=us-central1

source .env.local
github=$GITHUB_TOKEN

gcloud run deploy $function \
  --source ./cloud/$name \
  --function entrypoint \
  --base-image nodejs22 \
  --region $region \
  --project $project \
  --build-service-account projects/$project/serviceAccounts/$build@$project.iam.gserviceaccount.com \
  --service-account $run@$project.iam.gserviceaccount.com \
  --max-instances 2 \
  --no-invoker-iam-check \
  --set-env-vars GITHUB_TOKEN=$github
