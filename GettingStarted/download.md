Download Feel++ sources
=======================

# Getting the source via an archive

Feel++ is distributed as a tarball once in a while. The tarballs are available :
at

[https://github.com/feelpp/feelpp/releases](https://github.com/feelpp/feelpp/releases)

Download the latest tarball, then :
```
tar -xzf feelpp-X.YY.0.tar.gz
cd feelpp-X.YY.0
```
Once downloaded, you can compile and install it via the cmake procedure (see above).

# Getting the source via Git

In order to download the sources of Feel++, you can download it
directly from [the source repository](https://github.com/feelpp/feelpp)
thanks to Git.

To make it possible, you can download them anonymously or with an
account in Github that you have created. As an open-source project, we
strongly suggest you to create an account and take part of the project
with sharing your ideas, developments or suggests. For now, if you
want to get the sources without an account, open a command-line and
type

```
git clone --depth=xx https://github.com/feelpp/feelpp.git
```
with `xx` the number of last commits you want to save.
`xx=1` reduce the time to clone.
Then you can go to the Feel++ top directory with
```
  cd feel
```
You should obtain further directories such as:
```
applications/   # functional applications
benchmarks/  # applications under test
cmake/   # do not touch, used for compilation
contrib/
data/   # geometric data
doc/   # tutorial and examples
feel/   # Feel++ library
ports/   # used for Mac OS X installation
quickstart/   # basic examples
research/   # research projects using Feel++
scripts/ # various scripts
testsuite/ # Feel++ unit tests testsuite
CMakeLists.txt   # the file for cmake to build, do not modify
...
```
