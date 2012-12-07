PREREQUISITES
=============

1. A decent ANSI-compliant C++ compiler. GCC (version >= 4.4) is recommended.
2. CLN library (http://www.ginac.de/CLN), version >= 1.2.2
3. CMake, version >= 2.8 (version 2.6.x might work too).
4. Python, version >= 2.6
5. (optional) pkg-config utility (http://pkg-config.freedesktop.org)
6. (optional) bison (>= 2.3), flex (2.5.33), GNU readline (>= 4.3). These
   are necessary to build ginsh, the GiNaC interactive shell, and can be
   skipped if you don't intend to use ginsh (i.e. you need the GiNaC library
   for compiling another piece of a software).
7. (optional) To build the GiNaC tutorial and reference manual the doxygen
   utility (it can be downloaded from http://www.stack.nl/~dimitri/doxygen)
   and TeX are necessary.

INSTALLATION
============

To install from a source distribution:
--------------------------------------

1) Unpack the tarball

 $ tar xaf GiNaC-x.y.z.tar.bz2

2) Create a build directory

 $ mkdir ginac_build

3) Run CMake to generate Makefile's

 $ cd ginac_build
 $ cmake ../GiNaC-x.y.z

4) Actually build GiNaC

 $ make

5) Run the test and benchmark suite (not mandatory, but strongly recommended)

 $ make check

6) Install GiNaC

 [become root if necessary]
 # make install

To install from git: 
--------------------

The steps are essentially the same as compiling from the tarball, the only
difference is using git to get the code.

1) Download the code:

  $ git clone git://www.ginac.de/ginac.git ginac
  
  or, if you have already cloned the repository,

  $ git pull

Subsequent steps are the same as compiling from a source tarball.

