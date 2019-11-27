TABLE_FILE_NAME="hsize_points_elements_2d.csv"
printf "hsize, points, elements\n" > ${TABLE_FILE_NAME}

H_SIZES=`grep points meshes/*/partitioner.log --color | sed s,meshes/,, | sed -e 's|/partitioner.log||' -e 's|:      number of points  in memory||' -e 's|:.*||' | sort -r | tr -d "\n" `

#printf "H_SIZES = ${H_SIZES}"

for HSIZE in ${H_SIZES}
do
    NPOINTS=`grep points meshes/${HSIZE}/partitioner.log --color | sed 's#      number of points  in memory : ##'`    
    NELEMENTS=`grep elements meshes/${HSIZE}/partitioner.log --color | sed 's#      number of elements in memory : ##'`
    printf "${HSIZE}, ${NPOINTS}, ${NELEMENTS}\n" >> ${TABLE_FILE_NAME}
done
