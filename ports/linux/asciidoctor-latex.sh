#!/bin/bash

# First setup and use ruby2.0
# On Ubuntu 14.04.2:
# There seem to be a bug with ruby2.0 in 14.04.2, see:
# https://bugs.launchpad.net/ubuntu/+source/ruby2.0/+bug/1310292
# sudo ln -sf /usr/bin/ruby2.0 /usr/bin/ruby
# sudo ln -sf /usr/bin/gem2.0 /usr/bin/gem

git clone https://github.com/otavio/asciidoctor-latex.git
gem build asciidoctor-latex.gemspec
gem install *.gem
sudo gem install *.gem
