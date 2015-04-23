/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/*! \page TutorialEnv Setting up the Feel++ Environment

\tableofcontents


See section \ref ProgEnv for more information about \feel installation.

\section Minimal Minimal Example

Let's begin with our first program using the \feel framework (source:
<tt>doc/manual/tutorial/myapp.cpp</tt>).  Before all, you have to include the
\feel headers.

We use the C++ <tt>namespace</tt> to avoid <tt>Feel::</tt> prefix before
\feel objects.
\snippet myapp.cpp marker

We initialize the environment variables through the \feel <tt>Environment</tt> class.

<a href="#" class="top">top</a>
<hr>
\section MyAppWithOption Adding options

\snippet myappwithoptions.cpp marker

\li We pass command line options using the <a href="http://www.boost.org/doc/libs/1_53_0/doc/html/program_options.html">Boost Program Options</a>, library using the prefix <tt>po::</tt> which is a \feel alias for the Boost::program_options namespace. To add a new \feel option, we must create a new

\feel <tt>options_description</tt>. You must add the default \feel options
and the new one that we choose here as a double value. Note that the default
value will be assigned if not specified by the user.

<a href="#" class="top">top</a>
<hr>
\section Compilation Compilation, execution, logs
To compile a tutorial, just use the GNU make command.
\verbatim
  make feelpp_doc_<appname>
\endverbatim

where <tt><appname></tt> is the name of the application you wish to compile (here, <tt>myapp</tt>). Go to the execution directory as specified in the program, and execute it. You can change your option value.
\verbatim
  ./feelpp_doc_myapp [--value 6.6]
\endverbatim

You can list the log files created.
\verbatim
  ls /tmp/<your login>/feelpp_doc_myapp/
\endverbatim

If you open one of these log, you should be able to see your value and the processor number used to compute. You can run your application on several processors using MPI.
\verbatim
  mpirun -np 2 feelpp_doc_myapp
\endverbatim

Note that there will be one log for each processor in that case.

<a href="#" class="top">top</a>
<hr>
\section Config Config files

A config file can be parsed to the program to profile your options. The default config paths are,
    \li current dir
    \li <tt>$HOME/feel/config/</tt>
    \li <tt>$INSTALL_PREFIX/share/feel/config/</tt>

then you have to write inside one of these folders a file called
<tt><app_name>.cfg</tt> or <tt>feelpp_<app_name>.cfg</tt>. For example, our
<tt>myapp.cfg</tt> would looks like,

\verbatim
value=0.53
\endverbatim

Note that you can specify the config file through the option <tt>--config-file=<path></tt>

It's also possible to give several configuration files with the option <tt>--config-files <path1> <path2> <path3></tt>
\verbatim
 ./feelpp_doc_myapp --config-files ex1.cfg ex2.cfg ex3.cfg
\endverbatim

In the case where some options are duplicated in the files, the priority is given at the end :
\li <tt>ex3.cfg</tt> can overwrite options in <tt>ex2.cfg</tt> and <tt>ex3.cfg</tt>
\li <tt>ex2.cfg</tt> can overwrite options in <tt>ex1.cfg</tt>

All files in <tt> --config-files </tt> can overwrite options given by <tt> --config-file </tt>.
And all options in the command line can overwrite all options given in cfg files.

<a href="#" class="top">top</a>
<hr>
\section Initializing Initializing PETSc and Trilinos

PETSc is a suite of data structures and routines for the scalable (parallel)
solution of scientific applications modeled by partial differential
equations. It employs the MPI standard for parallelism.

\feel supports the PETSc framework, the <tt>Environment</tt> takes care of initializing the associated PETSc environment.

<a href="#" class="top">top</a>
<hr>



*/
}
