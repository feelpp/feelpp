#!/usr/bin/env bash
SOURCE_FILE_NAME=$1
TARGET_FILE_NAME=${1/times/time_per_core}
if [ -f ${SOURCE_FILE_NAME} ];
then
    printf "processes, hsize, points, elements, time\n" > ${TARGET_FILE_NAME}
    sort -n ${SOURCE_FILE_NAME} >> ${TARGET_FILE_NAME}
else
    printf "Could not find ${SOURCE_FILE_NAME}.\n"
fi

