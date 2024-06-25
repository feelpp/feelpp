set -eo pipefail
#set -x

# this script must be executed at the top level of the Feel++ directories

scriptdir=$PWD/$(dirname $0)
source $(dirname $0)/feelpp_pkg_common.sh

builddeps=$(cat $DIST | tr "\n" " ")


feelpp-pbuilder-dist $DIST login --save-after-login << EOF
echo "--- apt update"
apt-get update
apt-get -y install apt-transport-https ca-certificates gnupg software-properties-common wget


echo "--- add-apt-repository  "
add-apt-repository  'deb [trusted=yes] http://apt.feelpp.org/$FLAVOR/$DIST $DIST $CHANNEL'
echo "--- get repo signatures"
wget -O - http://apt.feelpp.org/apt.gpg | apt-key add -

if test "$DIST" = "focal"; then    
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc  | apt-key add -
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