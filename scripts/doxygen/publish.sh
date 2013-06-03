#! /opt/local/bin/bash

#usage: $1 = repo, $2 = branch_name
function check_repo
{
  echo "if [ ! -d $1 ]; then mkdir $1; cd $1; git init; git remote add -t $2 -f origin https://github.com/feelpp/feelpp.git; git checkout $2; fi"
  echo "cd $1 && git pull"
}

function builddox
{
  branch=$1
  feel_git=${2:-feelpp.git}
  cxx_compiler=${3:-/opt/local/bin/g++-mp-4.6}
  c_compiler=${4:-/opt/local/bin/gcc-mp-4.6}
  doxygen_dir=$HOME/doxygen-$branch

  echo "cd $feel_git && git checkout $branch"
  echo "[ -d $doxygen_dir ] && rm -rf $doxygen_dir"
  echo "mkdir $doxygen_dir"
  echo "cd $doxygen_dir"
  #echo "cmake -DCMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.6 -DCMAKE_C_COMPILER=/opt/local/bin/gcc-mp-4.6 ../$feel_git"
  echo "cmake $feel_git"
  echo "make doxygen"

  # now work in feelpp.docs to push the newly created doxygen files
  echo "mkdir -p $gh_pages/$branch"
  echo "cd $gh_pages  && git pull"

  echo "rsync -avz $doxygen_dir/doc/api/html/ $branch/"
  echo "git add -f $branch/* "
  echo "git commit -m "update Feel++ online documentation of branch $branch" -a"
  echo "git push origin gh-pages"
}

branch=${1:-develop}
base_dir=${2:-$HOME}
#Where the sources are stored
feelpp_source=${3:-$base_dir/feelpp.git}
#Where the gh-pages copy is
gh_pages=${4:-$base_dir/gh-pages}

#Create the gh-pages and feelpp_sources copies if they does not exist.
#echo "if [ ! -d $gh_pages ]; then mkdir $gh_pages; cd $gh_pages; git init; git remote add -t gh-pages -f origin https://github.com/feelpp/feelpp.git; git checkout gh-pages; fi"
check_repo $gh_pages gh-pages
#echo "if [ ! -d $fellpp_source ]; then mkdir $fellpp_source; cd $fellpp_source; git init; git remote add -t $branch -f origin https://github.com/feelpp/feelpp.git; git checkout $branch; fi"
check_repo $feelpp_source $branch

# checkout in master branch
# builddox master
builddox $branch $feelpp_source $gh_pages
