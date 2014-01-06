#!/bin/bash

function builddox
{
  branch=${1/\//_}
  doxygen_dir=$HOME/doxygen-$branch
  feelpp_source=$2
  gh_pages=$3

  cd $feelpp_source
  git checkout $1

  if [ ! -d $doxygen_dir ]; 
  then 
    mkdir $doxygen_dir;
  fi
  echo $doxygen_dir 
  cd $doxygen_dir
  pwd

  cmake $feelpp_source -DFEELPP_ENABLE_DOXYGEN=ON
  make doxygen #generate the doc associated to the branch $1 in ${doxygen_dir}/doc/api/html

  # now work in feelpp.docs to push the newly created doxygen files
  if [ ! -d  $gh_pages/$branch ];
  then 
    mkdir -p $gh_pages/$branch;
  fi
  cd $gh_pages
  rsync -avz $doxygen_dir/doc/api/html/ $branch/
  git add -A $branch
  git commit -m "update Feel++ online documentation of branch $branch"
}

base_dir=${1:-$HOME}
#Where the sources are stored
feelpp_source=$base_dir/${2:-feelpp.git}
#Where the gh-pages copy is
gh_pages=$base_dir/${3:-gh-pages}

#Create and/or update the ${gh-pages}/feelpp clone's repo
if [ ! -d $gh_pages ]; 
then 
  mkdir $gh_pages
  cd $gh_pages
  git clone -b gh-pages --single-branch https://github.com/feelpp/feelpp.git $gh_pages
else
  cd $gh_pages
  git pull
fi

#Create and/or update the feelpp's copy
if [ ! -d ${feelpp_source} ]; 
then 
  mkdir ${feelpp_source}; 
  cd ${feelpp_source};
  git clone https://github.com/feelpp/feelpp.git $feelpp_source
else
  cd $feelpp_source
  git pull
fi

#Update hifimagnet
if [ -d $feelpp_source/research/hifimagnet ]; 
then 
  cd $feelpp_source/research/hifimagnet
  git pull
  cd 
fi

#Create in ${gh_pages}/feelpp the associated doc of the ${branch}
builddox develop $feelpp_source $gh_pages
#builddox release/version-0.92 $feelpp_source $gh_pages
#builddox release/v0.95.0 $feelpp_source $gh_pages
#builddox release/v0.96.0 $feelpp_source $gh_pages

cd $feelpp_source
git checkout develop
cd ${gh_pages}
git push origin gh-pages

