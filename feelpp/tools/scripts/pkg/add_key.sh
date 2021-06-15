#!/bin/sh

pbuilder-dist $DIST login --save-after-login << EOF
apt install -y wget gnupg ca-certificates
wget -qO - https://feelpp.jfrog.io/artifactory/api/security/keypair/gpg-debian/public | gpg --dearmor > feelpp.gpg
mv feelpp.gpg /etc/apt/trusted.gpg.d
echo "deb [signed-by=/etc/apt/trusted.gpg.d/feelpp.gpg] https://feelpp.jfrog.io/artifactory/$FLAVOR $DIST $CHANNEL" | tee -a feelpp.list
mv feelpp.list /etc/apt/sources.list.d/
apt update
EOF
