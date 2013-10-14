#! /opt/local/bin/bash

function builddox
{
  doxygen_dir=$HOME/doxygen-$1
  feelpp_source=$2
  gh_pages=$3

  cd $feelpp_source
  git checkout $1

  if [ ! -d ${doxygen_dir} ]; 
  then 
    mkdir ${doxygen_dir};
  fi 
  cd $doxygen_dir

  cmake $feelpp_source -DFEELPP_ENABLE_DOXYGEN=ON
  make doxygen #generate the doc associated to the branch $1 in ${doxygen_dir}/doc/api/html

  # now work in feelpp.docs to push the newly created doxygen files
  if [ ! -d  $gh_pages/$branch ];
  then 
    mkdir -p $gh_pages/$branch;
  fi
  cd $gh_pages;
  rsync -avz $doxygen_dir/doc/api/html/ $branch/
  git commit -m "update Feel++ online documentation of branch $branch" -a
}

base_dir=${1:-$HOME}
#Where the sources are stored
feelpp_source=$base_dir/${2:-feelpp}
#Where the gh-pages copy is
gh_pages=$base_dir/${3:-gh-pages}

#Create and/or update the ${gh-pages}/feelpp clone's repo
if [ ! -d ${gh_pages} ]; 
then 
  mkdir ${gh_pages}; 
  cd ${gh_pages}; 
  git clone -b gh-pages --single-branch https://github.com/feelpp/feelpp.git
else
  cd ${gh_pages}/feelpp;
  git pull
fi

#Create and/or update the feelpp's copy
if [ ! -d ${feelpp_source} ]; 
then 
  mkdir ${feelpp_source}; 
  cd ${feelpp_source};
  git clone https://github.com/feelpp/feelpp.git
else
  cd $feelpp_source;
  git pull
fi

#Create in ${gh_pages}/feelpp the associated doc of the ${branch}
builddox develop $feelpp_source $gh_pages
builddox release/version-0.92 $feelpp_source $gh_pages
builddox release/v0.95.0 $feelpp_source $gh_pages

cd $feelpp_source
git checkout develop
cd ${gh_pages}
git push origin gh-pages


