#!/bin/bash

set -eo pipefail
set -x
source $(dirname $0)/feelpp_pkg_common.sh

echo $BUILDKITE_PASSPHRASE > $(pwd)/pp
cat $(pwd)/pp
aptly publish update -batch -passphrase-file=pp -force-overwrite ${DIST} s3:apt.feelpp.org:${FLAVOR}
rm $(pwd)/pp
