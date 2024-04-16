Prerequisites
=============

* C++11-compliant C++ compiler. GNU C++ compiler (version >= 5.1) is recommended.

* [Recommended] GNU multiprecision library (http://gmplib.org), version > 4.0.

* CMake, version 3.10 or newer.

* Ninja (https://ninja-build.org), version 1.8 or newer

* [Optional] To build CLN documentation texinfo and TeX are necessary.

* If you install CLN from git, you also need git (http://git-scm.com),
  version >= 1.8.

Installation from a source .tar.bz2 distribution
================================================

1) Unpack the archive

   $ tar xaf cln-x.y.z-tar.bz2

2) Create a build directory

   $ mkdir cln_build

3) Run CMake to generate build files

   $ cd cln_build
   $ cmake -GNinja ../cln-x.y.z

4) Actually build CLN

   $ cmake --build .

5) Run the test suite (not mandatory, but strongly recommended)

   $ cmake --build . --target check

6) Install CLN

   [become root if necessary]
   # cmake --build . -t install


Installation from git
=====================

The steps are essentially the same as compiling from the tarball, the only
difference is using git to get the code.

1) Download the code:

   $ git clone git://ginac.de/cln.git cln

   or, if you have already cloned the repository,

   $ git pull

Subsequent steps are the same as compiling from a source distribution.

