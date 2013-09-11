#! /opt/local/bin/bash

#usage: $1 = repo, $2 = branch_name
function check_repo
{
  if [ ! -d $1 ]; then mkdir $1; cd $1; git init; git remote add -t $2 -f origin https://github.com/feelpp/feelpp.git; git checkout $2; fi
  cd $1 && git pull
}

function builddox
{
  branch=$1
  feel_git=${2:-feelpp.git}
  cxx_compiler=${3:-/opt/local/bin/g++-mp-4.6}
  c_compiler=${4:-/opt/local/bin/gcc-mp-4.6}
  doxygen_dir=$HOME/doxygen-$branch

  cd $feel_git && git checkout $branch
  [ -d $doxygen_dir ] && rm -rf $doxygen_dir
  mkdir $doxygen_dir
  cd $doxygen_dir
 #cmake -DCMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.6 -DCMAKE_C_COMPILER=/opt/local/bin/gcc-mp-4.6 ../$feel_git
  cmake $feel_git -DFEELPP_ENABLE_DOXYGEN=ON
  make doxygen

  # now work in feelpp.docs to push the newly created doxygen files
  mkdir -p $gh_pages/$branch
  cd $gh_pages  && git pull

  rsync -avz $doxygen_dir/doc/api/html/ $branch/
  git add -f $branch/* 
  git commit -m "update Feel++ online documentation of branch $branch" -a
  git push origin gh-pages
}

branch=${1:-develop}
base_dir=${2:-$HOME}
#Where the sources are stored
feelpp_source=$base_dir/${3:-feelpp.git}
#Where the gh-pages copy is
gh_pages=$base_dir/${4:-gh-pages}

cd $feelpp_source
behind=$(git rev-list --left-right --count $branch...HEAD | awk '{print $1}')
ahead=$(git rev-list --left-right --count $branch...HEAD | awk '{print $2}')

#TODO -> Check the differences at the rsync time to check if there are differencies.
#if [ $behind -gt 0 ]; then

  #Create the gh-pages and feelpp_sources copies if they does not exist.
  #echo "if [ ! -d $gh_pages ]; then mkdir $gh_pages; cd $gh_pages; git init; git remote add -t gh-pages -f origin https://github.com/feelpp/feelpp.git; git checkout gh-pages; fi"
  check_repo $gh_pages gh-pages
  #echo "if [ ! -d $fellpp_source ]; then mkdir $fellpp_source; cd $fellpp_source; git init; git remote add -t $branch -f origin https://github.com/feelpp/feelpp.git; git checkout $branch; fi"
  check_repo $feelpp_source $branch

  # checkout in master branch
  # builddox master
  builddox $branch $feelpp_source $gh_pages
#fi
