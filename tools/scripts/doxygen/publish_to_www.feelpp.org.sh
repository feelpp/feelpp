#!/bin/bash

#the base dir
base_dir=${1:-$HOME}
#Where the sources are stored
feelpp_source=$base_dir/${2:-feelpp.git}
#Where the gh-pages copy is
gh_pages=$base_dir/${3:-gh-pages}

# First: generate the doc extracted from feelpp/feelpp
/home/vhuber/feelpp/tools/scripts/doxygen/publish.sh $base_dir $feelpp_source $gh_pages master develop 

cd $base_dir

if [ ! -d "$base_dir/www.feelpp.org" ]; 
then
  git clone http://github.com/feelpp/www.feelpp.org $base_dir/www.feelpp.org
else
  cd $base_dir/www.feelpp.org
  git pull
fi
cd $base_dir/www.feelpp.org

git checkout master
# the resulting static pages are in _site that is not tracked by github
jekyll build
git checkout gh-pages

# Copy the generated doc from doxygen
rsync -avz $gh_pages/ docs
git add -A
git commit -m "Update documentation" -a

rsync -avz _site ./
git add -A
git commit -m "Update web site" -a

git push

git checkout master
