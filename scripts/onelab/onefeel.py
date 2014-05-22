#!/usr/bin/python
# Script executing a Feel++ app with Gmsh through OneLab
# Initial contributors: Carolina Diaz, Jerome Boeglin, Sebastien Landre
# Authors: Alexandre Ancel

import os
import sys
import argparse
import subprocess

# Check if we have the executable accessible via the PATH variable
def checkExecutable(executable):
    ret = 0
    try:
        with open(os.devnull, "w") as fnull:
            ret = subprocess.call(executable, stdout = fnull)
    except OSError as e:
        return False

    # test return code
    if(ret != 0):
        return False

    return True

def main():

    sargs=[]
    pargs=[]
    # path to the script for synchronizing data
    syncData = os.path.dirname(os.path.abspath(__file__)) + "/syncData.py"

    # Check that we have the option separator
    hasDDash = False
    for i in range(0, len(sys.argv)):
        if(sys.argv[i] == "--"):
            hasDDash = True

    if(not hasDDash):
        print "Missing -- separator"
        print "usage: " + sys.argv[0] + " [script options] -- [application executable] [application options]"
        return 1

    for i in range(0, len(sys.argv)):
        if(sys.argv[i] == "--"):
            sargs=sys.argv[1:i]
            pargs=sys.argv[i+1:]

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--debug", help="Debug mode", type=int, default=0)
    parser.add_argument("-n", "--np", help="Number of MPI processes", type=int, default=1)
    parser.add_argument("-r", "--remote", help="Remote host on which the program will be launched", default="")
    parser.add_argument("-c", "--chroot", help="The chroot in which the program has to be launched", default="")
    parser.add_argument("-g", "--gmsh", help="Gmsh executable to use", default="/usr/bin/gmsh")
    args = parser.parse_args(sargs)

    if(args.debug > 0):
        print sys.argv
        print sargs
        print pargs


    sys.stdout.write("Checking for Gmsh ... ")
    if(checkExecutable([args.gmsh, "-"])):
        sys.stdout.write("OK\n")
    else:
        sys.stdout.write("Not found\n")
        return 1

    sys.stdout.write("Checking for data synchronization script ... ")
    cmd = ["ls", syncData]
    if(checkExecutable(cmd)):
        sys.stdout.write("OK\n")
    else:
        sys.stdout.write("Not found\n")
        return 1

    # check for the executable to use
    if(args.remote != ""):
        cmd = ["ssh", args.remote, "ls", pargs[0]];
        if(args.debug > 0):
            print cmd
        sys.stdout.write("Checking for executable ... ")
        if(checkExecutable(cmd)):
            sys.stdout.write("OK\n")
        else:
            sys.stdout.write("Not found\n")
            return 1
    else:
        cmd = ["ls", pargs[0]];
        if(args.debug > 0):
            print cmd
        sys.stdout.write("Checking for executable ... ")
        if(checkExecutable(cmd)):
            sys.stdout.write("OK\n")
        else:
            sys.stdout.write("Not found\n")
            return 1

    maincommand=""
    feelcommand=""
    for a in pargs:
        feelcommand = feelcommand + " " + a
    feelcommand = feelcommand + " --onelab.enable=1 --onelab.sync.script=" + syncData + " "

    # are we connecting on a remote computer ?
    if(args.remote != ""):
        maincommand = maincommand + "ssh " + str(args.remote) + " "
        feelcommand = feelcommand + "--onelab.remote=" + str(args.remote) + " "

    # do we want chroots ?
    if(args.chroot != ""):
        maincommand = maincommand + "schroot -c " + str(args.chroot) + " -- "
        feelcommand = feelcommand + "--onelab.chroot=" + str(args.chroot) + " "

    # do we want MPI ?
    if(args.np <= 0):
        print "Invalid number of MPI processes"
        exit(1)
    elif(args.np > 1):
        maincommand = maincommand + "mpiexec -np " + str(args.np) + " "

    maincommand = maincommand + feelcommand

    print "Creating Onelab files ..."
    if(args.debug > 0):
        print maincommand.split()
    retval = subprocess.call(maincommand.split())
    if(retval != 0):
        exit(1)

    # Checking for remote executable
    if(args.remote != ""):
        print "Copying remote config file ..."
        cmd = ["scp", args.remote + ":" + pargs[0] + ".ol", "."]
        if(args.debug > 0):
            print cmd
        retval = subprocess.call(cmd)
        if(retval != 0):
            exit(1)
        cmd = ["scp", args.remote + ":" + pargs[0] + ".onelab.cfg.ol", "."]
        if(args.debug > 0):
            print cmd
        retval = subprocess.call(cmd)
        if(retval != 0):
            exit(1)
    # we copy the config file for the local config only if we are not already in the same path
    elif(os.path.abspath(os.path.dirname(pargs[0])) != os.getcwd()):
        print "Copying config file ..."
        cmd = ["cp", pargs[0] + ".ol", "."]
        if(args.debug > 0):
            print cmd
        retval = subprocess.call(cmd)
        if(retval != 0):
            exit(1)
        cmd = ["cp", pargs[0] + ".onelab.cfg.ol", "."]
        if(args.debug > 0):
            print cmd
        retval = subprocess.call(cmd)
        if(retval != 0):
            exit(1)


    print "Executing Gmsh ..."
    cmd = [args.gmsh, "./" + os.path.basename(pargs[0]) + ".ol"]
    if(args.debug > 0):
        print cmd
    retval = subprocess.call(cmd)
    if(retval != 0):
        exit(1)

main()
