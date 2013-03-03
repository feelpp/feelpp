#! /opt/local/bin/bash

function builddox
{
    branch=$1
    cd feelpp.git && git pull && git checkout $branch && cd ..
    [ -d doxygen-$branch ] && rm -rf doxygen-$branch
    mkdir doxygen-$branch
    cd doxygen-$branch
    cmake ../feelpp.git
    make doxygen
    cd doc/manual/ 
#    make feelpp-manual_pdf
    cd ../../..

    # now work in feelpp.docs to push the newly created doxygen files
    mkdir -p feelpp.docs/api/$branch
    cd feelpp.docs && git pull
    
    rsync -avz ../doxygen-$branch/doc/api/html/ api/$branch/html/
    mkdir api/$branch/pdfs
#    cp ../doxygen-$branch/doc/manual/feelpp-manual.pdf  api/$branch/pdfs
    git add api/$branch/html/* 
#    git add api/$branch/pdfs/* 
    git commit -m"update $branch doxygen and user manual documentation" -a
    git push
    cd ..
}


if [ ! -d feelpp.docs ]; then git clone  https://code.google.com/p/feelpp.docs/ feelpp.docs; fi
cd feelpp.docs && git pull && cd ..
if [ ! -d feelpp.git ]; then git clone  https://github.com/feelpp/feelpp.git feelpp.git; fi
cd feelpp.git && git pull && cd ..

# checkout in master branch
# builddox master
builddox develop
