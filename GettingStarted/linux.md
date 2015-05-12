Feel++ on Linux 
===============

# Debian

Debian is the platform of choice for Feel++, it was developed mainly
on it. The commands to install Feel++ on Debian are

```
  sudo apt-get update
  sudo apt-get install feel++-apps libfeel++-dev feel++-doc
```

The interested user is encouraged to follow the Feel++ PTS page
* Feel++ <a href="http://packages.qa.debian.org/f/feel%2B%2B.html">Debian Packages Tracking System</a>

At the moment Feel++ compiles and is available on the following Debian
platforms:
* Feel++ <a href="https://buildd.debian.org/status/package.php?p=feel%2b%2b">Build results</a>

#  Ubuntu
Feel++ has its own PPA, checkout <a href="https://launchpad.net/~feelpp/+archive/ppa">feelpp ppa for Trusty, Saucy and Precise</a>.

The interested user might want to look at the following Ubuntu Launchpad Feel++ page <a href="https://launchpad.net/ubuntu/+source/feel++">Ubuntu Source
  Page for all Ubuntu versions</a>

## precise Ubuntu Precise - 12.04
```
	sudo add-apt-repository -y ppa:feelpp/ppa
	sudo add-apt-repository -y ppa:feelpp/petsc
	sudo add-apt-repository -y ppa:mapnik/boost-backports-1-54
	sudo add-apt-repository -y ppa:kalakris/eigen
	sudo apt-get -qq update
	sudo apt-get install -qq libboost1.54-all-dev mpi-default-dev mpi-default-bin libeigen3-dev libpetsc3.4.2-dev libcln-dev petsc-dev libxml2-dev gmsh bison flex doxygen doxygen-latex transfig imagemagick libtbb-dev libann-dev libglpk-dev automake libtool
	sudo apt-get install feel++-apps libfeel++-dev
```

\remark 
It allows us also to provide a recent version to compile the Feel++ projects on <a href="https://travis-ci.org/feelpp/feelpp">Travis-ci</a>, which is a continuous integration tool. 
The installation procedure is currently <a href="https://github.com/feelpp/feelpp/blob/develop/.travis.yml">as follows</a>.

## Ubuntu Trusty - 14.04
```
	sudo add-apt-repository ppa:feelpp/ppa
	sudo apt-get -qq update
	sudo apt-get install feel++-apps libfeel++-dev
```


