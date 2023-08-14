#!/bin/bash

# check if an github artifact exists

# Name of the artifact
ARTIFACT_NAME=$1
# GitHub repository, like user/repo
REPO_NAME=canoe
# GitHub API token for authentication
TOKEN=none

# Fetch artifact details using GitHub API
ARTIFACT_JSON=$(curl -s -H "Authorization: token $TOKEN" \
  "https://api.github.com/repos/$REPO_NAME/actions/artifacts")

# Check if the artifacts list is null
if [[ $(echo "$ARTIFACT_JSON" | jq '.artifacts') == "null" ]]; then
  echo "false"
else
  ARTIFACT_URL=$(echo "$ARTIFACT_JSON" \
  | jq -r --arg ARTIFACT_NAME "$ARTIFACT_NAME" \
  '.artifacts[] | select(.name==$ARTIFACT_NAME) | .url')

  # Check and set environment variable
  if [[ -z "$ARTIFACT_URL" ]]; then
    echo "false"
  else
    echo "true"
  fi
fi
