#$/bin/bash

rm tmpFile
touch tmpFile
echo "/**" >> tmpFile
echo "\\page  DocCmakeLists List of CMake Options" >> tmpFile
echo "\\author Feel++ Consortium" >> tmpFile
echo "\\date 2013-05-06" >> tmpFile
echo "" >> tmpFile
echo "In this section are presented the different cmake option available." >> tmpFile

echo "<table class=\"manual\">" >> tmpFile
echo "<tr><th>Option Name</th><th>Description</th><th>Default value</th></tr>" >> tmpFile

find $1 \( -name CMakeLists.txt -o -name "*.cmake" \) -exec grep -i 'option(' {} \; | grep -i feelpp_enable | grep -v \# | sed "s/ *option(//I" | sed "s/\t//g" | sed "s/^[ \t]*//g" | sed "s/)$//g" | awk '{printf "<tr><th>"$1"</th><th>" ; for(i=2;i<NF;i++) printf $i" "; print "</th><th>"$NF"</th></tr>"}' >> tmpFile
echo "</table>" >> tmpFile

echo "">> tmpFile
echo "*/" >> tmpFile

mv tmpFile $1/doc/api/GettingStarted/DocCmakeLists.doc
