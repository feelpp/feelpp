Feel++ on Mac OS X
==================

\feel is also supported on Mac operating systems, starting from OS X 10.9 Mavericks. The way to make it work is quite different.

# Compilers
In order to make \feel and <tt>cmake</tt> work properly, you first have to install *a Xcode for Gcc/clang. <br />
If your computer is recent, you can install it with your DVD that came with your machine (not the
OS DVD, but the applications one). You don't have to install the complete Xcode (you can uncheck iOS SDK for example,
it's not necessary here and requires a lot of memory). You can also install Xcode through the App Store.
Xcode will provide your computer all the basic tools to compile, such as gcc and clang.
It's the first step, you'll see later how to easily install the necessary tools to build \feel.

## Package management systems

Currently, building \feel is only available through <a
href="http://brew.sh">Homebrew</a>. <a
href="http://www.macports.org/install.php">MacPorts</a> versions are
also present in the source tree of \feel. You can find a version using
gcc in &lt;path to top feel directory&gt;/ports/macosx/macports, and
an experimental version using clang/libc++ in &lt;path to top feel
directory&gt;/ports/macosx/macports-xc5.

## Homebrew

<b>Introduction**: Homebrew is a free/open source software
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

**Key commands**: Homebrew base command is <tt>brew</tt>. Here is a list of useful commands to use Homebrew:
* <tt>brew doctor</tt>: Check if the system has any problem with the current installation of brew,
* <tt>brew install mypackage</tt>: This action installs the package mypackage,
* <tt>brew install [--devel|--HEAD] mypackage</tt>: These action respectively installs either the development version or the HEAD version of the package mypackage, if such versions are specified in the Formula file,
* <tt>brew uninstall mypackage</tt>: This action uninstalls the package mypackage.

**Formula**: A Formula is a <a href="https://www.ruby-lang.org">Ruby</a> script describing to Homebrew how to install a package. \feel uses specific Formulae that you can get in the \feel github repository: <a href="https://github.com/feelpp/homebrew-science">feelpp/homebrew-science</a>.

###  First time installation: Homebrew and Feel++
<p>
This section is aimed at users that do not have Homebrew already installed. If it is not your case, keep in mind the repository to tap and see the next section so you do not miss any important dependency. </p>
In order to build \feel from Homebrew, you have to do the following steps:
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

And you should be set to use \feel on MacOS X.

Be careful: you have to specifically use clang++ (and not clang).

<!--
cmake ../where/is/my/code -DCMAKE_CXX_COMPILER=`which clang++`
-->

###  Advanced usage: Homebrew and Feel++ dependencies

If Homebrew is already installed on your system, you might want to customize your installation for the correct dependencies to be met for \feel.

You can browse \feel dependencies using the following command:
```
brew deps feelpp
```

If you want to customize the compilation process for a dependency (Set debug mode, Remove checking steps, Remove the link with certain libraries, etc.),
you can access to the building options with the <tt>info</tt> flag. For exemple, with open-mpi:
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
When you install feelpp with <tt>brew install feelpp</tt>, the default flags for each dependency will be used.
If you require specific flags, first install the dependency with the correct flags before installing \feel.

**Important:** Needed dependencies for \feel: <br />
<ul>
<li><tt>boost</tt> needs to be installed with the following flags: <tt>--without-python --without-single --without-static --with-mpi --c++11</tt></li>
<li><tt>mumps</tt> needs to be installed with the following flag: <tt>--with-scotch5</tt></li>
</ul>

**Tips:** Reducing the compilation time: <br />
<ul>
<li><tt>scalapack</tt> can be installed with the following flag: <tt>--without-check</tt></li>
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
install it with the Mac OS X Installer using the <tt>pkg</tt> file
provided on their website. It is recommended that you install X11 (X
Window System) which is normally used to display X11 applications.<br>
If you have installed with the package installer
(<tt>MacPorts-2.x.x.pkg</tt>) that means MacPorts will be installed in
<tt>/opt/local</tt>. From now on, we will suppose that macports has
been installed in <tt>/opt/local</tt> which is the default MacPorts
location. Note that from now on, all tools installed by MacPorts will
be installed in <tt>/opt/local/bin</tt> or <tt>/opt/local/sbin</tt>
for example (that's here you'll find gcc4.7 or later e.g
<tt>/opt/local/bin/g++-mp-4.7</tt> once being installed).

**Key commands**: In your command-line, the software MacPorts is
  called by the command <tt>port</tt>. Here is a list of key commands
  for using MacPorts, if you want more informations please go to <a
  href="http://guide.macports.org/#using.port">MacPorts Commands</a>.

 * <tt>sudo port -v selfupdate</tt>: This action should be used regularly to update the local tree with the global MacPorts ports. The option <tt>-v</tt> enables verbose which generates verbose messages.
 * <tt>port info mypackage</tt>: This action is used to get information about a port (description, license, maintainer, etc.)
 * <tt>sudo port install mypackage</tt>: This action install the port mypackage
 * <tt>sudo port uninstall mypackage</tt>: This action uninstall the port mypackage
 * <tt>port installed</tt>: This action displays all ports installed and their versions, variants and activation status. You can also use the <tt>-v</tt> option to also display the platform and CPU  architecture(s) for which the ports were built, and any variants which were explicitly negated.
 * <tt>sudo port upgrade mypackage</tt>: This action updgrades installed ports and their dependencies when a <tt>Portfile</tt> in the repository has been updated. To avoid the upgrade of a port's dependencies, use the option <tt>-n</tt>.

**Portfile**: A Portfile is a TCL script which usually contains simple
keyword values and TCL expressions. Each package/port has a
corresponding Portfile but it's only a part of a port description.
\feel provides some mandatory Portfiles for its compilation which are
either not available in MacPorts or are buggy but \feel also provides
some Portfiles which are already available in MacPorts such as gmsh or
petsc. They usually provide either some fixes to ensure \feel works
properly or new version not yet available in MacPorts.  These
Portfiles are installed in <tt>ports/macosx/macports</tt>.


### MacPorts and Feel++

To be able to install \feel, add the following line in
<tt>/opt/local/etc/macports/source.conf</tt> at the top of the file
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
requirements for \feel to compile properly. If you have several cores
on your MacBook Pro, iMac or MacBook we suggest that you configure
macports to use all or some of them.

To do that uncomment the following line in the file
<tt>/opt/local/etc/macports/macports.conf</tt>


```bash
buildmakejobs	0 $\#$ all the cores
```

At the end of the <tt>sudo port install feel++</tt>, you have all
dependencies installed. To build all the Makefile, \cmake is
automatically launched but can have some libraries may not be found
but they are not mandatory for build Feel++, only the features related
to the missing libraries will be missing.

### Missing ports

<tt>cmake</tt> can build Makefiles even if some packages are missing
(latex2html, VTK ...). It's not necessary to install them but you can
complete the installation with MacPorts, <tt>cmake</tt> will find them
by itself once they have been installed.

### MacPorts and XCode 5

There is an experimental version of ports for \feel in &lt;path to top
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
