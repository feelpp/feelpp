#!/bin/bash

TOKEN=$GITHUB_OAUTH
VERSION_MAJOR=`cat cmake/modules/feelpp.version.cmake | grep "FEELPP_VERSION_MAJOR " | awk 'BEGIN{FS="\""}{print $2}'`
VERSION_MINOR=`cat cmake/modules/feelpp.version.cmake | grep "FEELPP_VERSION_MINOR " | awk 'BEGIN{FS="\""}{print $2}'`
VERSION_MICRO=`cat cmake/modules/feelpp.version.cmake | grep "FEELPP_VERSION_MICRO " | awk 'BEGIN{FS="\""}{print $2}'`
VERSION="${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_MICRO}"

# we require jq for parsing json files
LATEST_VERSION=`curl https://api.github.com/repos/feelpp/feelpp/releases/latest | jq '.tag_name' | awk 'BEGIN{FS="\""}{print $2}'`
echo "computing ChangeLog since ${LATEST_VERSION}..."
# we require github_changelog_generator for changelog generation
github_changelog_generator -p feelpp -t $TOKEN --release-branch develop \
                           --since-tag $LATEST_VERSION \
                           --no-issues-wo-labels \
                           --enhancement-labels type:feature,type:optimisation,type:question,type:refactoring,type:cleanup,type:update-3rdparty\
                           --bug-labels type:bug,type:ftbs,type:fte

CHANGELOG=`cat CHANGELOG.md`

API_JSON=$(printf '{"tag_name": "v%s","target_commitish": "develop","name": "v%s","body": "${CHANGELOG}","draft": false,"prerelease": false}' $VERSION $VERSION $VERSION)
echo $API_JSON
curl --data "$API_JSON" https://api.github.com/repos/feelpp/feelpp/releases?access_token=$TOKEN

# merge develop into master
git checkout master
git merge develop
git checkout develop

# update version now but only minor version
NEW_VERSION_MINOR=`expr ${VERSION_MINOR} + 1`
echo $NEW_VERSION_MINOR
sed -i "s/FEELPP_VERSION_MINOR \"[0-9]*\"/FEELPP_VERSION_MINOR \"${NEW_VERSION_MINOR}\"/g" cmake/modules/feelpp.version.cmake | grep "FEELPP_VERSION_MINOR "
NEW_VERSION="${VERSION_MAJOR}.${NEW_VERSION_MINOR}.${VERSION_MICRO}"
git add cmake/modules/feelpp.version.cmake
git commit -m"bump up version from v${VERSION} to v${NEW_VERSION} [ci skip]" cmake/modules/feelpp.version.cmake
git push https://$TOKEN@github.com/feelpp/feelpp.git

