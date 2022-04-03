#!/bin/bash

set -eo pipefail
set -x
source $(dirname $0)/feelpp_pkg_common.sh

echo $BUILDKITE_PASSPHRASE > pp
echo "aptly publish update -passphrase-file=pp -force-overwrite ${DIST} s3:apt.feelpp.org:${FLAVOR}"
rm pp
