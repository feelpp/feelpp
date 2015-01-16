#$/bin/bash

if [ -a tmpFile ]
then
  rm tmpFile
fi
touch tmpFile
today=$(date +"%F")
echo "/**" >> tmpFile
echo "\\page  DocCmakeLists List of CMake Options" >> tmpFile
echo "\\\\author Feel++ Consortium" >> tmpFile
echo "\\date $today" >> tmpFile
echo "" >> tmpFile
echo "In this section are presented the different cmake option named 'FEELPP_ENABLE_*'." >> tmpFile
echo "" >> tmpFile
echo "You are encouraged to directly read the CMakeLists.txt if some option have ambiguous comportement." >> tmpFile
echo "" >> tmpFile
echo "Many option (linked to Boost, Petsc or Gmsh) has to be defined by the user, or set in the environment variable, or defined at cmake time." >> tmpFile
echo "" >> tmpFile
echo "cmake \$FEELPP_SRC_DIR -DGMSH_DIR=/my/custom/gmsh/install" >> tmpFile
echo "" >> tmpFile

echo "<table class=\"manual\">" >> tmpFile
echo "<tr><th>Option Name</th><th>Description</th><th>Default value</th></tr>" >> tmpFile
find $1 \( -name CMakeLists.txt -o -name "*.cmake" \) -exec grep -i 'option(' {} \; | grep -i feelpp_enable | grep -v \# | sed "s/ [oO][pP][tT][iI][oO][nN]*(//" | sed "s/\t//g" | sed "s/^[ \t]*//g" | sed "s/)$//g" | awk '{printf "<tr><th>"$1"</th><th>" ; for(i=2;i<NF;i++) printf $i" "; print "</th><th>"$NF"</th></tr>"}' | sort >> tmpFile
echo "</table>" >> tmpFile

echo "Here are now given the list of used Env Variable" >> tmpFile

echo "<table class=\"manual\">" >> tmpFile
echo "<tr><th>Option Name</th></tr>" >> tmpFile
find $1 \( -name CMakeLists.txt -o -name "*.cmake" \) -exec grep -io 'ENV{.*}' {} \; | sed "s/}/} /g" |awk '{print $1}' | sort | uniq | awk '{print "<tr><th>"$1"</th></tr>"}' | sort >> tmpFile
echo "</table>" >> tmpFile

echo "">> tmpFile
echo "*/" >> tmpFile

mv tmpFile $1/doc/api/GettingStarted/DocCmakeLists.doc
