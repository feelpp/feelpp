#!/bin/bash

base_dir=$HOME

# First: generate the doc extracted from feelpp/feelpp
./publish.sh

if [ ! -d www.feelpp.org ]; 
then
  git clone http://github.com/feelpp/www.feelpp.org
else
  cd www.feelpp.org
  git pull
fi
cd www.feelpp.org

# Copy the generated doc
rsync -avz $HOME/gh-pages docs

# Generate the static web site
jekyll build

# the resulting static pages are in _site
git checkout gh-pages
rsync -avz _site ./
git add -A
git commit -m "deploy web site" -a
git pull origin gh-pages
git push origin gh-pages
git checkout master
