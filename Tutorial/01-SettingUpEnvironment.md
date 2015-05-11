Setting up the Feel++++ Environment 
=================================


# Minimal Example

Let's begin with our first program using the \Feel++ framework (source:
<tt>doc/manual/tutorial/myapp.cpp</tt>).  Before all, you have to include the
\Feel++ headers.

We use the C++ <tt>namespace</tt> to avoid <tt>Feel++::</tt> prefix before
\Feel++ objects.
\snippet myapp.cpp marker

We initialize the environment variables through the \Feel++ <tt>Environment</tt> class.


# Adding options

\snippet myappwithoptions.cpp marker

* We pass command line options using the <a href="http://www.boost.org/doc/libs/1_53_0/doc/html/program_options.html">Boost Program Options</a>, library using the prefix <tt>po::</tt> which is a \Feel++ alias for the Boost::program_options namespace. To add a new \Feel++ option, we must create a new

\Feel++ <tt>options_description</tt>. You must add the default \Feel++ options
and the new one that we choose here as a double value. Note that the default
value will be assigned if not specified by the user.


# Compilation execution and logs

To compile a tutorial, just use the GNU make command.
```
  make Feel++pp_doc_<appname>
```

where <tt><appname></tt> is the name of the application you wish to compile (here, <tt>myapp</tt>). Go to the execution directory as specified in the program, and execute it. You can change your option value.
```
  ./Feel++pp_doc_myapp [--value 6.6]
```

You can list the log files created.
```
  ls /tmp/<your login>/Feel++pp_doc_myapp/
```

If you open one of these log, you should be able to see your value and the processor number used to compute. You can run your application on several processors using MPI.
```
  mpirun -np 2 Feel++pp_doc_myapp
```

Note that there will be one log for each processor in that case.

<a href="#" class="top">top</a>
<hr>
# Config files

A config file can be parsed to the program to profile your options. The default config paths are,
    * current dir
    * <tt>$HOME/Feel++/config/</tt>
    * <tt>$INSTALL_PREFIX/share/Feel++/config/</tt>

then you have to write inside one of these folders a file called
<tt><app_name>.cfg</tt> or <tt>Feel++pp_<app_name>.cfg</tt>. For example, our
<tt>myapp.cfg</tt> would looks like,

```
value=0.53
```

Note that you can specify the config file through the option <tt>--config-file=<path></tt>

It's also possible to give several configuration files with the option <tt>--config-files <path1> <path2> <path3></tt>
```
 ./Feel++pp_doc_myapp --config-files ex1.cfg ex2.cfg ex3.cfg
```

In the case where some options are duplicated in the files, the priority is given at the end :
  * <tt>ex3.cfg</tt> can overwrite options in <tt>ex2.cfg</tt> and <tt>ex3.cfg</tt>
  * <tt>ex2.cfg</tt> can overwrite options in <tt>ex1.cfg</tt>

All files in <tt> --config-files </tt> can overwrite options given by <tt> --config-file </tt>.
And all options in the command line can overwrite all options given in cfg files.

<a href="#" class="top">top</a>
<hr>
# Initializing PETSc, SLEPc and other third party libraries

PETSc is a suite of data structures and routines for the scalable (parallel)
solution of scientific applications modeled by partial differential
equations. It employs the MPI standard for parallelism.

\Feel++ supports the PETSc framework, the <tt>Environment</tt> takes care of initializing the associated PETSc environment.



