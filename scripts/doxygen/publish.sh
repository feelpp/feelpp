#! /opt/local/bin/bash

function builddox
{
    branch=$1
    cd feelpp.git && git pull && git checkout $branch && cd ..
    [ -d doxygen-$branch ] && rm -rf doxygen-$branch
    mkdir doxygen-$branch
    cd doxygen-$branch
    cmake -DCMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.6 -DCMAKE_C_COMPILER=/opt/local/bin/gcc-mp-4.6 ../feelpp.git
    make doxygen
    cd ../

    # now work in feelpp.docs to push the newly created doxygen files
    mkdir -p gh-pages/$branch
    cd gh-pages  && git pull
    
    rsync -avz ../doxygen-$branch/doc/api/html/ $branch/
    git add -f $branch/* 
    git commit -m"update Feel++ online documentation of branch $branch" -a
    git push origin gh-pages
    cd ..
}


if [ ! -d gh-pages ]; then mkdir gh-pages; cd gh-pages; git init; git remote add -t gh-pages -f origin https://github.com/feelpp/feelpp.git; git checkout gh-pages; cd ..; fi
cd gh-pages && git pull && cd ..
#if [ ! -d feelpp.git ]; then git clone  https://github.com/feelpp/feelpp.git feelpp.git; fi
cd feelpp.git && git pull && cd ..

# checkout in master branch
# builddox master
builddox develop
