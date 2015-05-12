Setting up the Feel++++ Environment
=================================


# Minimal Example

Let's begin with our first program using the Feel++ framework.  To
start, you include the Feel++ headers.

We use the C++ `namespace` to avoid `Feel++::` prefix before
Feel++ objects.

We initialize the environment variables through the Feel++ `Environment` class.

!CODEFILE "code/environment.cpp"


# Adding options

 We pass command line options using the
 [Boost Program Options](http://www.boost.org/doc/libs/1_53_0/doc/html/program_options.html),
 library using the prefix `po::` which is a Feel++ alias for the
 Boost::program_options namespace. To add a new Feel++ option, we must
 create a new

Feel++ `options_description`. You must add the default Feel++ options
and the new one that we choose here as a double value. Note that the default
value will be assigned if not specified by the user.


# Compilation execution and logs

To compile a tutorial, just use the GNU make command.
```bash
  make Feel++pp_doc_<appname>
```

where `<appname>` is the name of the application you wish to compile (here, `myapp`). Go to the execution directory as specified in the program, and execute it. You can change your option value.
```bash
  ./Feel++pp_doc_myapp [--value 6.6]
```

You can list the log files created.
```
  ls /tmp/<your login>/Feel++pp_doc_myapp/
```

If you open one of these log, you should be able to see your value and the processor number used to compute. You can run your application on several processors using MPI.
```bash
  mpirun -np 2 Feel++pp_doc_myapp
```

Note that there will be one log for each processor in that case.



# Config files

A config file can be parsed to the program to profile your options. The default config paths are,
    * current dir
    * `$HOME/Feel++/config/`
    * `$INSTALL_PREFIX/share/Feel++/config/`

then you have to write inside one of these folders a file called
`<app_name>.cfg` or `Feel++pp_<app_name>.cfg`. For example, our
`myapp.cfg` would looks like,

```
value=0.53
```

Note that you can specify the config file through the option `--config-file=<path>`

It's also possible to give several configuration files with the option `--config-files <path1> <path2> <path3>`
```bash
 ./Feel++pp_doc_myapp --config-files ex1.cfg ex2.cfg ex3.cfg
```

In the case where some options are duplicated in the files, the priority is given at the end :
  * `ex3.cfg` can overwrite options in `ex2.cfg` and `ex3.cfg`
  * `ex2.cfg` can overwrite options in `ex1.cfg`

All files in ` --config-files ` can overwrite options given by ` --config-file `.
And all options in the command line can overwrite all options given in cfg files.



# Initializing PETSc, SLEPc and other third party libraries

PETSc is a suite of data structures and routines for the scalable (parallel)
solution of scientific applications modeled by partial differential
equations. It employs the MPI standard for parallelism.

Feel++ supports the PETSc framework, the `Environment` takes care of initializing the associated PETSc environment.
