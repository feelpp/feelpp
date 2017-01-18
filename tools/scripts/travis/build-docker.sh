#! /bin/bash

docker --version
sudo apt-get update
sudo apt-get install -y -q -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" docker-engine
git clone --depth=1 https://github.com/feelpp/docker
cd docker/feelpp-libs && bash mkimg.sh -f ${TARGET_IMAGE} --jobs 1 --branch ${TARGET_BRANCH} --cmakeflags "${CMAKE_FLAGS}" --cxxflags "${CXXFLAGS}" --cxx "${FEELPP_CXX}" --cc "${FEELPP_CC}" --travis true
docker images
