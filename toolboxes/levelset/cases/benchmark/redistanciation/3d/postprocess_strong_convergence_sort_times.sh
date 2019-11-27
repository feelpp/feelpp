ROOT_DIR=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/redistanciation/3d
H_SIZES="0.1024 0.0512 0.0256 0.0128 0.0064 0.0032 0.0016 0.0008 0.0004 0.0002 0.0001"
H_SIZES="0.1024 0.0512 0.0256 0.0128 0.01 0.0064 0.0032 0.0016 0.0008 0.0004 0.0002 0.0001"
SOURCE_FILE_NAME_BASE=redistanciation_3d_timer_h
TARGET_FILE_NAME_BASE=results_3d_strong_convergence_time_per_core_h
for HSIZE in ${H_SIZES}
do
    SOURCE_FILE_NAME=${ROOT_DIR}/${SOURCE_FILE_NAME_BASE}_$HSIZE.csv
    TARGET_FILE_NAME=${ROOT_DIR}/${TARGET_FILE_NAME_BASE}_$HSIZE.csv
    if [ -f ${SOURCE_FILE_NAME} ];
    then
        printf "processes, hsize, points, elements, time\n" > ${TARGET_FILE_NAME}
        sort -n ${SOURCE_FILE_NAME} >> ${TARGET_FILE_NAME}
    else
        printf "Could not find ${SOURCE_FILE_NAME}.\n"
    fi
done
