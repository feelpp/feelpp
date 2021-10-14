#!/bin/bash

version=$1
incoming=${2:-.}
repo=${3:-$HOME/debian}
distribution=${4:-focal}
channel=${5:-latest}
component=${6:-feelpp}
arch=${7:-amd64}
repository=$repo/$distribution

for i in ${incoming}/*.deb; do
    newpackage=`dpkg-deb -I $i | grep Package | awk '{print $2}'`
    newversion=`dpkg-deb -I $i | grep Version | awk '{print $2}'`
    currentpackage=`reprepro -Vb $repository -C $channel list $distribution | grep "\b$newpackage\b"`
    if ! test  -z  "$currentpackage"; then 
        #check version
        currentversion=`reprepro -Vb $repository -C $channel list $distribution | grep "\b$newpackage\b" | awk '{print $3}'`
        if test $currentversion = $newversion; then
            # remove package before adding the new version
            reprepro -Vb $repository -C $channel remove $distribution $newpackage
        fi
    fi
    reprepro -Vb $repository -C $channel includedeb $distribution $i
done