MAINTAINER Feel++ Support <support@feelpp.org>

ARG PROJECT="project"
ARG URL=""
ARG BRANCH="master"
ARG BUILD_JOBS=1
ARG CMAKE_FLAGS=""
ARG MAKE_FLAGS=""
ARG CXX=clang++-3.9
ARG CC=clang-3.9
ARG GITHUB_OAUTH=

USER feelpp
ENV HOME /home/feelpp

COPY feelpp_build.sh $HOME/feelpp_build.sh 
# Install Feel++
RUN bash -c "echo" ${PROJECT} "-" ${URL} "-"  ${BRANCH} "-"  ${BUILD_JOBS} "-" ${CMAKE_FLAGS}
RUN bash -c "source $HOME/feelpp_build.sh && feelpp_build ${PROJECT} ${URL} ${BRANCH} ${BUILD_JOBS} \"${CMAKE_FLAGS}\" \"${MAKE_FLAGS}\""

# COPY WELCOME $HOME/WELCOME
USER root
#ENTRYPOINT ["/sbin/my_init","--quiet","--","sudo","-u","feelpp","/bin/sh","-c"]
CMD ["/usr/local/bin/start.sh"]
