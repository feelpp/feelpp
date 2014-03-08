#!/usr/bin/python
# Script that syncs data after the execution of the application
# this is to compensate the fact that we cannot use several .out commands for an app
# it will be executed, before the actual merges of data
# Authors: Alexandre Ancel

import sys
import subprocess

def main():
    # TODO argument check

    # read first line of file with files to sync
    f = open(sys.argv[1], 'r')
    line = f.readline()

    print line

    # delete comment character and remove blanks before and after
    line = line.translate(None, '#').lstrip().rstrip()

    # get filenames
    files = line.split()

    print "External call: Syncing data ..."
    for f in files:
        cmd = ["scp", f, "."]
        retval = subprocess.call(cmd)
        if(retval != 0):
            exit(1)

main()
