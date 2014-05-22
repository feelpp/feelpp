#!/bin/bash

function builddox
{
  branch=${1/\//_}
  doxygen_dir=$HOME/doxygen-$branch
  feelpp_source=$2
  gh_pages=$3

  cpt=0
  STR=""
  # allow to generate doc for project in feelpp/research directory
  if [ "$#" -ge 3 ]
  then
    for i in "$@"
    do
      if [ $cpt -ge 3 ]
      then
        project=$(echo $i | awk '{print toupper($0)}')
        STR="$STR -DFEELPP_ENABLE_"$project"_DOCUMENTATION=ON"
        if [ ! -d $feelpp_source/research/$i ]; 
        then
          cd $feelpp_source/research/
          git clone https://github.com/feelpp/$i\.git
        else
          cd $feelpp_source/research/$i
          git pull
        fi
      fi
      cpt=$(($cpt+1))
    done
  fi

  # Update feelpp source
  cd $feelpp_source
  git checkout $1

  # create doxygen_dir
  if [ ! -d $doxygen_dir ]; 
  then 
    mkdir $doxygen_dir;
  fi
  echo $doxygen_dir 
  cd $doxygen_dir

  # generate doc for branch $1 and specified project
  # the results are in ${doxygen_dir}/doc/api/html
  cmake $feelpp_source -DFEELPP_ENABLE_DOXYGEN=ON $STR
  make doxygen #generate the doc associated to the branch $1 in ${doxygen_dir}/doc/api/html

  # now work in feelpp.docs to push the newly created doxygen files
  if [ ! -d  $gh_pages/$branch ];
  then 
    mkdir -p $gh_pages/$branch;
  fi
  cd $gh_pages
  rsync -avz $doxygen_dir/doc/api/html/ $branch/
  # git add -A $branch
  # git commit -m "update Feel++ online documentation of branch $branch"
}

base_dir=$1
#Where the sources are stored
feelpp_source=$2
#Where the gh-pages copy is
gh_pages=$3

#Create and/or update the ${gh-pages}/feelpp clone's repo
if [ ! -d $gh_pages ]; 
then 
  mkdir $gh_pages
  #cd $gh_pages
  #git clone -b gh-pages --single-branch https://github.com/feelpp/feelpp.git $gh_pages
fi
cd $gh_pages

#Create and/or update the feelpp's copy
if [ ! -d ${feelpp_source} ]; 
then 
  mkdir ${feelpp_source}; 
  cd ${feelpp_source};
  git clone https://github.com/feelpp/feelpp.git $feelpp_source
else
  cd $feelpp_source
  git checkout .
  git pull
fi

#Update hifimagnet
if [ -d $feelpp_source/research/hifimagnet ]; 
then 
  cd $feelpp_source/research/hifimagnet
  git pull
  cd 
fi

#Create in ${gh_pages}/feelpp the associated doc for ${branch}
cpt=0
if [ "$#" -ge 3 ]
then
  for i in "$@"
  do
    if [ $cpt -ge 3 ]
    then
      builddox $i $feelpp_source $gh_pages # cemosis bubble
    fi
      cpt=$(($cpt+1))
  done
fi

#cd $feelpp_source
#git checkout develop
#cd ${gh_pages}
#git push origin gh-pages

