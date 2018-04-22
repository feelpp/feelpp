#!/bin/bash

TOKEN=$GITHUB_OAUTH
# we require jq for parsing json files
LATEST_VERSION=`curl https://api.github.com/repos/feelpp/feelpp/releases/latest | jq '.tag_name' | awk 'BEGIN{FS="\""}{print $2}'`
echo "computing ChangeLog since ${LATEST_VERSION}..."
github_changelog_generator -p feelpp -t $TOKEN --release-branch develop \
  --since-tag $LATEST_VERSION \
  --no-issues-wo-labels \
  --future-release v0.104.0 \
  --enhancement-labels type:feature,type:optimisation,type:question,type:refactoring,type:cleanup,type:update-3rdparty\
  --bug-labels type:bug,type:ftbs,type:fte

