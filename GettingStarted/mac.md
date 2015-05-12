Feel++ on Mac OS X
==================

Feel++ is also supported on Mac operating systems, starting from OS X 10.9 Mavericks. The way to make it work is quite different.

# Compilers
In order to make Feel++ and `cmake` work properly, you first have to install *a Xcode for Gcc/clang. <br />
If your computer is recent, you can install it with your DVD that came with your machine (not the
OS DVD, but the applications one). You don't have to install the complete Xcode (you can uncheck iOS SDK for example,
it's not necessary here and requires a lot of memory). You can also install Xcode through the App Store.
Xcode will provide your computer all the basic tools to compile, such as gcc and clang.
It's the first step, you'll see later how to easily install the necessary tools to build Feel++.

## Package management systems

Currently, building Feel++ is only available through <a
href="http://brew.sh">Homebrew</a>. <a
href="http://www.macports.org/install.php">MacPorts</a> versions are
also present in the source tree of Feel++. You can find a version using
gcc in &lt;path to top feel directory&gt;/ports/macosx/macports, and
an experimental version using clang/libc++ in &lt;path to top feel
directory&gt;/ports/macosx/macports-xc5.

## Homebrew

**Introduction**: Homebrew is a free/open source software
introduced to simplify the installation of other free/open source
software on the MacOS X ecosystem. It is distributed under the <a
href="https://github.com/mxcl/homebrew/blob/master/Library/Homebrew/LICENSE">BSD
2 Clause (NetBSD) license</a>. For more information, visit their <a
href="http://brew.sh">website</a>.

**Installation**: To install the latest version of Homebrew, simply
visit their <a href="http://brew.sh">website</a> and follow the
instructions. Each new package Homebrew installs is built into an
intermediate place called the Cellar (usually /usr/local/Cellar) and
then the packages are symlinked into /usr/local (default).

**Key commands**: Homebrew base command is `brew`. Here is a list of useful commands to use Homebrew:
* `brew doctor`: Check if the system has any problem with the current installation of brew,
* `brew install mypackage`: This action installs the package mypackage,
* `brew install [--devel|--HEAD] mypackage`: These action respectively installs either the development version or the HEAD version of the package mypackage, if such versions are specified in the Formula file,
* `brew uninstall mypackage`: This action uninstalls the package mypackage.

**Formula**: A Formula is a <a href="https://www.ruby-lang.org">Ruby</a> script describing to Homebrew how to install a package. Feel++ uses specific Formulae that you can get in the Feel++ github repository: <a href="https://github.com/feelpp/homebrew-science">feelpp/homebrew-science</a>.

###  First time installation: Homebrew and Feel++
<p>
This section is aimed at users that do not have Homebrew already installed. If it is not your case, keep in mind the repository to tap and see the next section so you do not miss any important dependency. </p>
In order to build Feel++ from Homebrew, you have to do the following steps:
```
# Install Homebrew:
ruby -e "$(curl -fsSL https://raw.github.com/feelpp/homebrew/go/install)"
# Check that everything is ok or fix warnings/errors if necessary:
brew doctor
# Get the Formulae specific to Feel++ on the github of feelpp:
brew tap feelpp/science

# Finally either:
# - Install Feel++ and all its dependencies:
brew install feelpp
# - Or install only the dependencies:
# (useful if you plan to use a different version from the master one, e.g. develop)
brew install --only-dependencies feelpp
```

<!--
or in a more detailed way:
```
# Install Homebrew:
ruby -e "$(curl -fsSL https://raw.github.com/feelpp/homebrew/go)"
# Check that everything is ok or fix warnings/errors if necessary:
brew doctor
# Get the Formulae specific to Feel++ on the github of feelpp:
brew tap feelpp/science
# Install openmpi with c++11 support:
brew install open-mpi --c++11
# Install boost:
brew install boost --without-python --without-single --without-static --with-mpi --c++11
# Install Ann, Petsc, Gmsh and HDF5:
brew install ann && brew install petsc && brew install gmsh && brew install hdf5
# Install Feel++:
brew install feelpp
```
-->

And you should be set to use Feel++ on MacOS X.

Be careful: you have to specifically use clang++ (and not clang).

<!--
cmake ../where/is/my/code -DCMAKE_CXX_COMPILER=`which clang++`
-->

###  Advanced usage: Homebrew and Feel++ dependencies

If Homebrew is already installed on your system, you might want to customize your installation for the correct dependencies to be met for Feel++.

You can browse Feel++ dependencies using the following command:
```
brew deps feelpp
```

If you want to customize the compilation process for a dependency (Set debug mode, Remove checking steps, Remove the link with certain libraries, etc.),
you can access to the building options with the `info` flag. For exemple, with open-mpi:
```
user@server: brew info open-mpi

open-mpi: stable 1.7.5
http://www.open-mpi.org/
Conflicts with: lcdf-typetools, mpich2
/usr/local/Cellar/open-mpi/1.7.5 (761 files, 23M) *
  Built from source
From: https://github.com/Homebrew/homebrew/commits/master/Library/Formula/open-mpi.rb
==> Dependencies
Required: libevent
==> Options
--c++11
    Build using C++11 mode
--disable-fortran
    Do not build the Fortran bindings
--enable-mpi-thread-multiple
    Enable MPI_THREAD_MULTIPLE
```

You then just have to pass the needed flags, when installing the dependency.<br />
When you install feelpp with `brew install feelpp`, the default flags for each dependency will be used.
If you require specific flags, first install the dependency with the correct flags before installing Feel++.

**Important:** Needed dependencies for Feel++: <br />
<ul>
<li>`boost` needs to be installed with the following flags: `--without-python --without-single --without-static --with-mpi --c++11`</li>
<li>`mumps` needs to be installed with the following flag: `--with-scotch5`</li>
</ul>

**Tips:** Reducing the compilation time: <br />
<ul>
<li>`scalapack` can be installed with the following flag: `--without-check`</li>
</ul>

<a href="#" class="top">top</a>
<hr>

## MacPorts

**Introduction**: MacPorts is an open-source community projet which
  aims to design an easy-to-use system for compiling, installing and
  upgrading open-source software on Mac OS X operating system. It is
  distributed under <a
  href="http://opensource.org/licenses/bsd-license.php">BSD
  License</a> and facilitate the access to thousands of ports
  (software) without installing or compiling open-source software.
  MacPorts provides a single software tree which includes the latest
  stable releases of approximately 17700 ports targeting the current
  Mac OS X release (10.9). If you want more information, please visit
  their <a href="http://www.macports.org/">website</a>.

**Installation**: To install the latest version of MacPorts, please go
to <a href="http://www.macports.org/install.php">Installing
MacPorts</a> page and follow the instructions. The simplest way is to
install it with the Mac OS X Installer using the `pkg` file
provided on their website. It is recommended that you install X11 (X
Window System) which is normally used to display X11 applications.<br>
If you have installed with the package installer
(`MacPorts-2.x.x.pkg`) that means MacPorts will be installed in
`/opt/local`. From now on, we will suppose that macports has
been installed in `/opt/local` which is the default MacPorts
location. Note that from now on, all tools installed by MacPorts will
be installed in `/opt/local/bin` or `/opt/local/sbin`
for example (that's here you'll find gcc4.7 or later e.g
`/opt/local/bin/g++-mp-4.7` once being installed).

**Key commands**: In your command-line, the software MacPorts is
  called by the command `port`. Here is a list of key commands
  for using MacPorts, if you want more informations please go to <a
  href="http://guide.macports.org/#using.port">MacPorts Commands</a>.

 * `sudo port -v selfupdate`: This action should be used regularly to update the local tree with the global MacPorts ports. The option `-v` enables verbose which generates verbose messages.
 * `port info mypackage`: This action is used to get information about a port (description, license, maintainer, etc.)
 * `sudo port install mypackage`: This action install the port mypackage
 * `sudo port uninstall mypackage`: This action uninstall the port mypackage
 * `port installed`: This action displays all ports installed and their versions, variants and activation status. You can also use the `-v` option to also display the platform and CPU  architecture(s) for which the ports were built, and any variants which were explicitly negated.
 * `sudo port upgrade mypackage`: This action updgrades installed ports and their dependencies when a `Portfile` in the repository has been updated. To avoid the upgrade of a port's dependencies, use the option `-n`.

**Portfile**: A Portfile is a TCL script which usually contains simple
keyword values and TCL expressions. Each package/port has a
corresponding Portfile but it's only a part of a port description.
Feel++ provides some mandatory Portfiles for its compilation which are
either not available in MacPorts or are buggy but Feel++ also provides
some Portfiles which are already available in MacPorts such as gmsh or
petsc. They usually provide either some fixes to ensure Feel++ works
properly or new version not yet available in MacPorts.  These
Portfiles are installed in `ports/macosx/macports`.


### MacPorts and Feel++

To be able to install Feel++, add the following line in
`/opt/local/etc/macports/source.conf` at the top of the file
before any other sources:

```
file:///<path to feel top directory>/ports/macosx/macports
```

Once it's done, type in a command-line:

```
  cd <your path to feel top directory>/ports/macosx/macports
  sudo portindex -f
```

You should have an output like this:

```
Reading port index in $<$your path to feel top directory$>$/ports/macosx/macports
Adding port science/feel++
Adding port science/gmsh
Adding port science/petsc

Total number of ports parsed:   3
Ports successfully parsed:      3
Ports failed:                   0
Up-to-date ports skipped:       0
```

Your are now able to type

```bash
  sudo port install feel++
```

It might take some time (possibly an entire day) to compile all the
requirements for Feel++ to compile properly. If you have several cores
on your MacBook Pro, iMac or MacBook we suggest that you configure
macports to use all or some of them.

To do that uncomment the following line in the file
`/opt/local/etc/macports/macports.conf`


```bash
buildmakejobs	0 $\#$ all the cores
```

At the end of the `sudo port install feel++`, you have all
dependencies installed. To build all the Makefile, \cmake is
automatically launched but can have some libraries may not be found
but they are not mandatory for build Feel++, only the features related
to the missing libraries will be missing.

### Missing ports

`cmake` can build Makefiles even if some packages are missing
(latex2html, VTK ...). It's not necessary to install them but you can
complete the installation with MacPorts, `cmake` will find them
by itself once they have been installed.

### MacPorts and XCode 5

There is an experimental version of ports for Feel++ in &lt;path to top
feel directory&gt;/ports/macosx/macports-xc5. Using these ports will
set up the compilation using clang and libc++. The process is similar
to the one previously described for MacPorts, except for one point:
Before starting to install packages, you must switch to the llvm c++
standard library by adding the following line to your macports.conf
file:

```bash
cxx_stdlib  libc++
```

This requires MacPorts to be at least on version 2.2.1 for the flag to
be recognized and will normally cause all the packages you will
install to be recompiled using libc++ instead of libstdc++.
