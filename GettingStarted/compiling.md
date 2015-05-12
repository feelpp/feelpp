Compiling Feel++ from sources
=============================

# Quick start

```bash
	cmake ../feel -DCMAKE_CXX_COMPILER=/usr/bin/clang++ // Generate the makefile
	sudo make install // actually compile and install the library
	make -j 4 feelpp_qs_laplacian // Compile the library and an example
	./quickstart/feelpp_qs_laplacian // execute one example
```

# General information

Feel++ uses <a href="http://www.cmake.org">`cmake`</a> as its build system. Check that \cmake is using `gcc4.7` (or a higher version) or `clang++` as C++ compiler (you can use the option `CMAKE_CXX_COMPILER=<path>/g++-4.7` where the `path` depends on your OS, it's probably `/usr/bin` or
`/opt/local/bin` but you can also change it with the command `ccmake`
and press `t` for advanced options).
Feel++, using `cmake`, can be built either in source and out of source and different
build type:
* minsizerel: minimal size release
* release: release
* debug: debug
* none(default)

In order to install Feel++ on Unix systems (other than Mac OS X, if
you have a Macintosh, please go to \ref Mac), you have to install many
dependencies before. Those libraries and programs are necessary for
the compilation and installation of the Feel++ libraries.

This is the list of all the libraries you must have installed on your
computer, and the `*-dev` packages for some of them.

Required packages:
* g++ (4.7, 4.8, 4.9, 5.1) OR Clang (3.3, 3.4, 3.5, 3.6, 3.7)
* MPI: openmpi (preferred) or mpich
* Boost >= 1.55
* Petsc >= 3.2
* CMake >= 2.8.7
* Gmsh >= 2.8.3 (can be found at <a href="http://www.geuz.org/gmsh">http://www.geuz.org/gmsh</a>)
* Libxml2

Optional packages:
* Superlu
* Suitesparse(umfpack)
* Metis: scotch with the metis interface (preferred), metis (non-free)
* Google perftools
* Paraview, this is not stricly required to run Feel++ programs but it is somehow
  necessary for visualization. (Can be found at <a href="http://www.paraview.org">http://www.paraview.org</a>)
* Python (2.5 or later) for the validation tools

Note that all these packages are available under Debian GNU/Linux and
Ubuntu. Once you have installed those dependencies, you can jump to \ref Compiling.



## Compiling out the source

**It is not allowed to build the library in the source directory **: The best way is to have a directory (`FEEL` for example) in which you have:
\code
ls FEEL
feel/ // Sources
feel.opt/ // Build directory
\endcode
where `feel` is the top directory where the source have been downloaded, using git or trackballs.
`cmake` is used as it:
```
	cd FEEL/feel.opt
	cmake ../feel -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_BUILD_TYPE=RelWithDebInfo
```

You can customize the build type in any way.
Please have a look to
```
	FEEL/feel/CMakeLists.txt
```
or at the `List of Cmake Options` section.

## Per-system specific configuration

If you need to set specific CMake variables to fit an exotic system
(e.g. supercomputer), you can create a machine file fitting this
system.  To do so, simply create a file named
`feelpp.machines.&lt;name&gt;.cmake` in the
`cmake/machines` directory, where `&lt;name&gt;` is the
machine name returned by the execution of `uname -n`.  The file
will automatically be recognized and loaded, when you will be
compiling Feel++ on this system.  </p>

Once CMake has done its work, you are now able to compile the library with:
```bash
  make
```

You can speed up the compilation process, if you have a multicore processor. To
do so, you have to specify the number of parallel jobs `make` will be
allowed to spawn by using the `-j` flag:

```bash
  make -j <nbjobs>
```

**Important**: from now on, all commands should be typed in
  `feel.opt` or its subdirectories.


## Unit tests

Feel++ comes with built-in unit tests to ensure that the library is
working correctly. You can launch basic tests based on the quickstart
application by typing:


```bash
  make check
```

Or you can execute the whole testsuite, by typing the following
command: (This might take a while)

```bash
  make check-full
```
</p>

<p>

<a href="#" class="top">top</a>
<hr>
