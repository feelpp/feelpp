# -*- mode: Dockerfile -*-

ARG TAG
FROM ghcr.io/feelpp/feelpp:$TAG

USER root
RUN apt-get update;\
        apt-get install -y --no-install-recommends git git-lfs jq curl pybind11-dev python3-pybind11; \
        apt-get install -y 'libfeelpp*dev';\
        rm -rf /var/lib/apt/lists/* ;

RUN json=$(curl -s "https://api.github.com/repos/Kitware/CMake/releases"); \
        latest_tag=$(echo $json | jq -r '.[] | select(.prerelease==false).tag_name' | sort --unique --version-sort | grep -v "rc" | tail -1);\
	sh_url=$(echo $json | jq -r ".[] | select(.tag_name==\"${latest_tag}\").assets[].browser_download_url | select(endswith(\"inux-x86_64.sh\"))");\
	curl -sL ${sh_url} -o cmakeinstall.sh \
	&& chmod +x cmakeinstall.sh \
	&& ./cmakeinstall.sh --prefix=/usr/local --exclude-subdir \
	&& rm cmakeinstall.sh

USER feelpp
