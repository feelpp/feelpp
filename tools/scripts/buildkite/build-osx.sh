#! /bin/bash

set -eo pipefail

echo '--- tapping in homebrew-science'
brew tap homebrew/homebrew-science
echo '--- tapping in homebrew-feelpp'
brew tap feelpp/homebrew-feelpp
echo '--- updating formulas'
brew update
brew upgrade
echo '--- install feelpp HEAD'
brew install feelpp --HEAD
echo '--- cleaning up'
brew remove feelpp

