<!-- -*- mode: markdown  -->

Feel++ Installation
===================

This file describes some basic installation procedures for Feel++. A more complete documentation is available in
the [Feel++ online documentation](http://www.feelpp.org/docs/develop/BuildingP.html).

## Requirements

 - G++ (>=4.6) or Clang (>= 3.1)
 - Boost C++ libraries (>=1.46)
 - PETSc (>=3.1)
 - Gmsh (>= 2.8)

## Optional packages

 - SLEPc (>=3.1), Library to solve standard and generalized eigenvalue problems.
   Note that the version of SLEPc installed depends on the version of PETSc
   (they have should the same release number),
 - ANN (>=1.1.2) Approximate Nearest Neighbor
 - GLPK (>= 4.45) GNU Linear Programming ToolKit

## Building Feel++

Note that in-source builds of Feel++ are not possible.

### Out source build

A `configure` script is provided to simplify the setup of the programming environment: a small set of options
is available :
 - `-r`: Release mode
 - `-d`: Debug mode
By default `configure` looks for the default C++ compiler available (e.g. typically `/usr/bin/c++`),
to change the default C++ compiler, define the environment variable CXX (e.g. `export CXX=/usr/bin/clang++`).

Here is an example:
```
cd <feelpp toplevel source directory>
mkdir opt
cd opt
../configure -r
make
make install
```

At the end of the process, you have

 - the Feel++ libraries and its dependencies compiled
 - a quickstart example solving $- \Delta u = f$ in $[0,1]^2$ with $u=g$ on $\partial \Omega$   located in the quickstart/ directory

More details are available in Feel++ manual and Feel++ website (http://www.feelpp.org)

# Running Feel++ examples

## Environment variables

If Feel++ was installed in a directory other than
```
/opt/local
/usr/local
/usr
```
then you need to set the environment variable **FEELPP_DIR** to ensure that next
**cmake** commands find the required components to compile the Feel++ examples

For example, if Feel++ was installed in `/opt/feelpp`, then set
```
# using sh type shells
export FEELPP_DIR=/opt/feelpp
# using csh type shells
setenv FEELPP_DIR /opt/feelpp
```


## Compiling Feel++ examples

The Feel++ examples are in `FEELPP_DIR/share/doc/feel/examples` organized in
sub-directories following the layout in the Feel++ source code in `doc/manual`.
A `CMakeLists.txt` is available to compile the examples.

```
mkdir <feelpp examples build dir>
cd <feelpp examples build dir>
cmake $FEELPP_DIR/share/doc/feel/examples
```

Once the `cmake` process is done, type for example
```
make tutorial
```
to build all the Feel++ tutorial examples from the manual.
