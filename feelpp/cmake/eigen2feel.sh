find . -name CMakeLists.txt  | xargs -istr sed -ibak "s/DESTINATION include /DESTINATION include\/feel /g" str
find . -name CMakeLists.txt  | xargs -istr sed -ibak "s/include\/unsupported/include\/feel\/unsupported/g" str
find . -name CMakeLists.txt  | xargs -istr sed -ibak "s/include\/Eigen/include\/feel\/Eigen/g" str
