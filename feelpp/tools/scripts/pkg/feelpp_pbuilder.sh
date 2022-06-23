set -eo pipefail
#set -x

# this script must be executed at the top level of the Feel++ directories

scriptdir=$PWD/$(dirname $0)
source $(dirname $0)/feelpp_pkg_common.sh

builddeps=$(cat $DIST | tr "\n" " ")


pbuilder-dist $DIST login --save-after-login << EOF
echo "--- apt update"
apt-get update
apt-get -y install apt-transport-https ca-certificates gnupg software-properties-common wget

echo "--- get repo signaures"
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc  | apt-key add -
wget -O - http://apt.feelpp.org/apt.gpg | apt-key add -

echo "--- add-apt-repository  "
if test ! "$DIST" = "bullseye"; then
    add-apt-repository  'deb [trusted=yes] http://apt.feelpp.org/$FLAVOR/$DIST $DIST $CHANNEL'
    add-apt-repository  'deb https://apt.kitware.com/$FLAVOR/ $DIST main'
fi
if [ "$DIST" = "bullseye" ]; then
add-apt-repository  'deb http://deb.debian.org/debian $DIST-backports main'
fi

echo "--- apt update"
apt-get update

echo $builddeps

echo "--- apt install"
apt-get -y install $builddeps
if [ "$DIST" = "bullseye" ]; then
    apt-get -y install -t bullseye-backports  cmake
fi
EOF