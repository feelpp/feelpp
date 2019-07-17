#! /bin/bash

set -eo pipefail

echo '--- tapping in homebrew-feelpp'
brew tap feelpp/homebrew-feelpp
echo '--- updating formulas'
# brew update
echo '--- install feelpp HEAD'
brew install feelpp --HEAD
echo '--- cleaning up'
brew remove feelpp

