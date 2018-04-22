#!/bin/bash

TOKEN=$GITHUB_OAUTH
# we require jq for parsing json files
LATEST_VERSION=`curl https://api.github.com/repos/feelpp/feelpp/releases/latest | jq '.tag_name' | awk 'BEGIN{FS="\""}{print $2}'`
echo "computing ChangeLog since ${LATEST_VERSION}..."
github_changelog_generator -p feelpp -t $TOKEN --release-branch develop \
  --between-tags v0.100.0,v0.101.0,v0.102.0,v0.103.0,v0.103.1,v0.103.2,v0.104.0 \
  --no-issues-wo-labels \
  --future-release v0.104.0 \
  --enhancement-labels type:feature,type:optimisation,type:question,type:refactoring,type:cleanup,type:update-3rdparty\
  --bug-labels type:bug,type:ftbs,type:fte

