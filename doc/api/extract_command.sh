#$/bin/bash

rm tmpFile
touch tmpFile
echo "/**" >> tmpFile
echo "\\page  DocCmakeLists List of CMake Options" >> tmpFile
echo "\\author Feel++ Consortium" >> tmpFile
echo "\\date 2013-05-06" >> tmpFile
echo "" >> tmpFile
echo "In this section are presented the different cmake option available." >> tmpFile

find $1 \( -name CMakeLists.txt -o -name "*.cmake" \) -exec grep -i 'option(' {} \; | grep -i feelpp_enable | grep -v \# | sed "s/ *option//I" | sed "s/\t//g" >> tmpFile

echo "">> tmpFile
echo "*/" >> tmpFile

mv tmpFile $1/doc/api/GettingStarted/DocCmakeLists.doc
