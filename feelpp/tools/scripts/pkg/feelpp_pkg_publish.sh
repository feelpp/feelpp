#!/bin/bash

set -eo pipefail
set -x
source $(dirname $0)/feelpp_pkg_common.sh

aptly publish update -batch -passphrase="${GPG_PASSPHRASE}" -force-overwrite  -gpg-key=${GPG_KEY} ${DIST} s3:apt.feelpp.org:${FLAVOR}/${DIST}
