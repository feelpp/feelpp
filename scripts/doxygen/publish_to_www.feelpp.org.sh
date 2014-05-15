#!/bin/bash

#the base dir
base_dir=${1:-$HOME}
#Where the sources are stored
feelpp_source=$base_dir/${2:-feelpp.git}
#Where the gh-pages copy is
gh_pages=$base_dir/${3:-gh-pages}

# First: generate the doc extracted from feelpp/feelpp
# ./publish.sh $base_dir $feelpp_source $gh_pages master develop release/version-0.92 release/v0.95.0 release/v0.96.0

cd $base_dir

if [ ! -d "$base_dir/www.feelpp.org" ]; 
then
  git clone http://github.com/feelpp/www.feelpp.org $base_dir/www.feelpp.org
else
  cd $base_dir/www.feelpp.org
  git pull
fi
cd www.feelpp.org

git checkout gh-pages

# Copy the generated doc
rsync -avz $gh_pages/ docs

# Generate the static web site
jekyll build

# the resulting static pages are in _site
rsync -avz _site ./

git add -A
git commit -m "deploy web site" -a
git pull origin gh-pages
git push origin gh-pages
git checkout master
