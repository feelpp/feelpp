#! /bin/bash

set -eo pipefail
#set -x

# this script must be executed at the top level of the Feel++ directories

scriptdir=$PWD/$(dirname $0)
BUILDKITE_AGENT_NAME=${BUILDKITE_AGENT_NAME:-default}
# default values
CHANNEL=latest
if [ "$BUILDKITE_BRANCH" = "develop" -o  "$BRANCH" = "develop" ]; then
    CHANNEL=latest
fi
if [ "$BUILDKITE_BRANCH" = "master" -o  "$BRANCH" = "master" ]; then
    CHANNEL=stable
fi 
DIST=${DIST:-focal}
if [ "$DIST" = "bionic" -o "$DIST" = "eoan" -o "$DIST" = "focal"  -o "$DIST" = "groovy"  -o "$DIST" = "hirsute"  ]; then
   FLAVOR=ubuntu
elif [ "$DIST" = "buster" -o "$DIST" = "bullseye" -o "$DIST" = "sid" ]; then
    FLAVOR=debian
fi



COMPONENT=${COMPONENT:-feelpp}
#OTHERMIRROR=
#if [ "$COMPONENT" = "feelpp-toolboxes" ]; then
#    OTHERMIRROR="deb https://dl.bintray.com/feelpp/ubuntu $DIST $CHANNEL"
#fi
echo "Building $FLAVOR/$DIST packages for channel $CHANNEL Feel++ component $COMPONENT"

# update apt tools
echo "--- update apt tools"
sudo apt update
sudo apt install -y ubuntu-dev-tools
if [ "$COMPONENT" = "feelpp-toolboxes" -o "$COMPONENT" = "feelpp-mor" ]; then
    echo "--- update feelpp packages"
    # this is necessary for toolboxes and mor to retrieve the proper feelpp version
    sudo apt install -y --reinstall libfeelpp-dev feelpp-tools 
fi

#PBUILDER_RESULTS=/var/lib/buildkite-agent/pbuilder/${DIST}_result_${BUILDKITE_AGENT_NAME}
# local debug build
PBUILDER_RESULTS=$HOME/pbuilder/${DIST}_result_${BUILDKITE_AGENT_NAME}
if [ ! -f $HOME/pbuilder/${DIST}_base.tgz ]; then
    echo "--- creating distribution $DIST"
    pbuilder-dist $DIST create
fi
echo "--- start from clean slate"
if [ -d build-$DIST ]; then rm -rf build-$DIST; fi

if [ -n "$(ls -A ${PBUILDER_RESULTS}/ 2>/dev/null)" ];
then
    echo "removing previous builds in $PBUILDER_RESULTS";
    ls -1 ${PBUILDER_RESULTS}/
    rm -f "${PBUILDER_RESULTS}/*";
else
    echo "no files in ${PBUILDER_RESULTS}/";
fi

echo "--- update for pbuilder $DIST"
pbuilder-dist $DIST update
if [ "$DIST" = "bionic" ]; then
echo "--- fixes for pbuilder $DIST"
    pbuilder-dist $DIST login --save-after-login << EOF
apt-get update
apt-get install apt-transport-https ca-certificates gnupg software-properties-common wget

echo "deb https://apt.kitware.com/ubuntu/ bionic main" >> /etc/apt/sources.list
echo "deb http://archive.ubuntu.com/ubuntu bionic universe"  >> /etc/apt/sources.list
echo "deb http://archive.ubuntu.com/ubuntu bionic-updates main" >> /etc/apt/sources.list
echo "deb http://archive.ubuntu.com/ubuntu bionic-updates universe" >> /etc/apt/sources.list
echo "deb http://archive.ubuntu.com/ubuntu bionic-backports main" >> /etc/apt/sources.list
echo "deb http://archive.ubuntu.com/ubuntu bionic-security main" >> /etc/apt/sources.list
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc  | apt-key add -
apt-get update
EOF
fi
if [ "$DIST" = "buster" ]; then
echo "--- fixes for pbuilder $DIST"
    pbuilder-dist $DIST login --save-after-login << EOF
apt-get update
apt-get install -y apt-transport-https ca-certificates gnupg software-properties-common wget

echo "deb http://deb.debian.org/debian buster-backports main" >> /etc/apt/sources.list
apt-get update
EOF
fi
# pbuilder-dist $DIST update
if [ "$COMPONENT" = "feelpp-toolboxes" -o "$COMPONENT" = "feelpp-mor" ]; then
    echo "--- add bintray key for $COMPONENT"
    export DIST
    export CHANNEL
    export FLAVOR
    pbuilder-dist $DIST login --save-after-login << EOF
apt-get install -y apt-transport-https ca-certificates gnupg software-properties-common wget
if [ "$DIST" = "buster" ]; then
    echo "deb http://deb.debian.org/debian buster-backports main" >> /etc/apt/sources.list
fi
echo "deb http://apt.feelpp.org/$FLAVOR $DIST $CHANNEL" | tee -a /etc/apt/sources.list
wget -qO - http://apt.feelpp.org/apt.gpg | apt-key add
# wget -qO  - https://feelpp.jfrog.io/artifactory/api/security/keypair/gpg-debian/public | apt-key add -
apt update
EOF
fi

echo "--- setting directory build-$DIST to build source tarball"
#git clone https://github.com/feelpp/feelpp /tmp/feelpp
mkdir build-$DIST
if [ "$COMPONENT" = "feelpp" ]; then
    cd build-$DIST && ../configure -r --enable-toolboxes --enable-mor --cmakeflags="-DFEELPP_ENABLE_GIT=OFF -DLIBBSON_DIR=/usr -DLIBMONGOC_DIR=/usr"
elif [ "$COMPONENT" = "feelpp-toolboxes" ]; then
    cd build-$DIST && ../configure -r --root=../toolboxes
elif [ "$COMPONENT" = "feelpp-mor" ]; then
    cd build-$DIST && ../configure -r --root=../mor    
fi
make package_source
echo "--- cloning feelpp.pkg"
git clone -q https://github.com/feelpp/feelpp.pkg.git
# local debug build
#ln -s ../../Debian/feelpp.pkg



main_version=$(echo ${COMPONENT}-*.tar.gz | sed  "s/${COMPONENT}-\([0-9.]*\)-*\([a-z.0-9]*\).tar.gz/\1/g")
extra_version=$(echo ${COMPONENT}-*.tar.gz | sed "s/${COMPONENT}-\([0-9.]*\)-*\([a-z.0-9]*\).tar.gz/\2/g")

if [ -z $extra_version ]; then
    version=$(echo ${COMPONENT}-*.tar.gz | sed  "s/${COMPONENT}-\([0-9.]*\)-*\([a-z.0-9]*\).tar.gz/\1/g" )
    rename_archive=$(echo ${COMPONENT}-*.tar.gz | sed "s/${COMPONENT}-\([0-9.]*\)-*\([a-z.0-9]*\).tar.gz/${COMPONENT}_\1.orig.tar.gz/g" )
else
    rename_archive=$(echo ${COMPONENT}-*.tar.gz | sed "s/${COMPONENT}-\([0-9.]*\)-*\([a-z.0-9]*\).tar.gz/${COMPONENT}_\1~\2.orig.tar.gz/g")
    version=$(echo ${COMPONENT}-*.tar.gz | sed  "s/${COMPONENT}-\([0-9.]*\)-*\([a-z.0-9]*\).tar.gz/\1~\2/g" )
fi
echo "--- building archive $rename_archive for debian"
cp ${COMPONENT}-*.tar.gz feelpp.pkg/${COMPONENT}/$rename_archive
cd feelpp.pkg/${COMPONENT}/$DIST && tar xzf ../$rename_archive --strip 1

echo "--- update changelog ${COMPONENT}  $version-1"
dch -v "$version-1" --distribution "unstable" -b "New upstream commits"

echo "--- add source ${COMPONENT}  $version-1"
dpkg-source -b .

echo "--- building ${COMPONENT} debian version $version-1"
pbuilder-dist $DIST build --buildresult ${PBUILDER_RESULTS}  ../${COMPONENT}_${version}-1.dsc

echo "+++ uploading ${PBUILDER_RESULTS} to bintray $COMPONENT $FLAVOR/$DIST"
ls  -1 ${PBUILDER_RESULTS}
aptly repo add feelpp-$DIST-$CHANNEL ${PBUILDER_RESULTS}
aptly publish update ${DIST} s3:apt.feelpp.org:${FLAVOR}
#echo "$scriptdir/publish.sh $main_version ${PBUILDER_RESULTS} $HOME/debian $DIST $CHANNEL $COMPONENT"
#$scriptdir/publish.sh  $main_version ${PBUILDER_RESULTS} $HOME/debian $DIST $CHANNEL $COMPONENT
#repreprocmd=reprepro -Vb $HOME/debian/$DIST -C $COMPONENT 

## echo "../upload_bintray.sh $main_version ${PBUILDER_RESULTS} $FLAVOR $DIST $CHANNEL $COMPONENT"
## ../upload_bintray.sh $main_version ${PBUILDER_RESULTS} $FLAVOR $DIST $CHANNEL $COMPONENT

