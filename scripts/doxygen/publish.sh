#/bin/sh

function builddox
{
    branch = $1
    cd feelpp.git && git pull && git checkout $branch && cd ..
    [ -d doxygen-$branch ] && rm -rf doxygen-$branch
    mkdir doxygen-$branch
    cd doxygen-$branch
    cmake ../feelpp.git
    make doxygen
    cd ..

    # now work in feelpp.docs to push the newly created doxygen files
    mkdir -p feelpp.docs/api/$branch
    cd feelpp.docs && git pull
    # cleanup html directory
    if [ -d api/$branch/html]; then
        rm -rf api/$branch/html;
    fi
    cp -r ../doxygen-$branch/doc/api/html api/$branch/
    git add api/$branch/html
    git push
}

if [! -d feelpp.docs ]; then git clone  https://code.google.com/p/feelpp.docs/ feelpp.docs; fi
cd feelpp.docs && git pull && cd ..
if [! -d feelpp.git ]; then git clone  https://github.com/feelpp/feelpp.git feelpp.git; fi
cd feelpp.git && git pull && cd ..

# checkout in master branch
builddox master
builddox develop
