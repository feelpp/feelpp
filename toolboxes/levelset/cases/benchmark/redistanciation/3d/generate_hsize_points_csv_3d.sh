TABLE_FILE_NAME="hsize_points_3d.csv"
printf "hsize, points\n" > ${TABLE_FILE_NAME}
grep points meshes/*/partitioner.log --color | sed s,meshes/,, | sed -e s#/partitioner.log#,# -e 's|:      number of points  in memory : | |' | sort -r >> ${TABLE_FILE_NAME}
