# Onefeel

This part of the repository contains the script interface a Feel++ application with Gmsh through Onelab.

## Initial setup

Make sure you compiled the patched Gmsh version, specifically provided for Feel++, if you intend to use remote execution.
The patched version is located in the [feelpp/feelpp.tps](https://github.com/feelpp/feelpp.tps/tree/master) respository (see [link](https://github.com/feelpp/feelpp.tps/tree/master/gmsh)).

## Execution syntax

The base script for using the Onelab interface is onefeel.py.
The syntax of the script is the following:
```python onefeel.py [-d n] [-n np] [-r host] [-c chroot] [-g gmsh] -- /app_path/app_name [app_option ...]```  
With the following options:
- d (optional): Set the debug level for debugging the script (defaults to 0 for off)
- n (optional): Set the number of MPI processes to uses (Can be modified through the interface, defaults to 1)
- r (optional): Set the remote host to use (Can be modified through the interface, defaults to localhost)
- c (optional): Set the chroot to use (defaults to empty)
- g (optional): Set the path to the patched gmsh executable (defaults to /usr/bin/gmsh)

The `--` special flag is here to separate the script options on its left and the application executable and options on its right.
