#!/usr/bin/python
# Script checking out a Gmsh version working with Onefeel
# Authors: Alexandre Ancel

import os
import subprocess

def main():
    # Initial checks
    """
    if(os.path.exists("./gmsh")):
        print("Gmsh directory already exists")
        return 1
    """

    if(not os.path.exists("./onefeel.patch")):
        print("Onefeel patch does not exist (onefeel.patch)")
        return 1

    if(not os.path.exists("./gmsh")):
        print "Checking out Gmsh ..."
        cmd = ["svn", "checkout", "-r", "17790", "https://geuz.org/svn/gmsh/trunk", "gmsh"]
        retval = subprocess.call(cmd)
        if(retval != 0):
            exit(1)

    os.chdir("./gmsh");

    cmd = ["cp", "../onefeel.patch", "."]
    retval = subprocess.call(cmd)
    if(retval != 0):
        exit(1)

    print "Patching Gmsh ..."
    cmd =  ["patch", "-p0", "-i", "onefeel.patch" ]
    retval = subprocess.call(cmd)
    if(retval != 0):
        exit(1)

    print "Building Gmsh ..."
    if(not os.path.exists("build")):
        os.mkdir("build")
    os.chdir("build")
    cmd = ["cmake", "..", "-DENABLE_MPI=ON", "-DENABLE_FLTK=ON"]
    retval = subprocess.call(cmd)
    if(retval != 0):
        exit(1)

    cmd = ["make"]
    retval = subprocess.call(cmd)
    if(retval != 0):
        exit(1)


main()
