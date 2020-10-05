TABLE_FILE_NAME="hsize_elements_2d.csv"
printf "hsize, elements\n" > ${TABLE_FILE_NAME}
grep elements meshes/*/partitioner.log --color | sed s,meshes/,, | sed -e s#/partitioner.log#,# -e 's|:      number of elements in memory : | |' | sort -r >> ${TABLE_FILE_NAME}
