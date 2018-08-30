#!/bin/bash

set -x
TOKEN=$GITHUB_OAUTH
ASSET=$2
echo $(basename $ASSET)
VERSION=$1
VERSIONID=`curl -s https://$TOKEN@api.github.com/repos/feelpp/feelpp/releases/tags/$VERSION | jq '.id'`
echo $VERSIONID
curl -X POST --header "Content-Type:$(file -b --mime-type $ASSET)" --data-binary @$ASSET "https://$TOKEN@uploads.github.com/repos/feelpp/feelpp/releases/$VERSIONID/assets?name=$(basename $ASSET)"

#&access_token=$TOKEN"


#curl  --data $ASSET https://$TOKEN@api.github.com/repos/feelpp/feelpp/releases/$VERSION/assets?name=$ASSET
