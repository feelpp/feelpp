#!/usr/bin/python
# Script that syncs data after the execution of the application
# this is to compensate the fact that we cannot use several .out commands for an app
# it will be executed, before the actual merges of data
# Authors: Alexandre Ancel

import os
import sys
import subprocess

def main():
    # TODO argument check
    if(len(sys.argv) != 2):
        print "Invalid number of arguments"
        return 1

    # check that the file exists
    if(os.path.exists(sys.argv[1])):
        # read first line of file with files to sync
        f = open(sys.argv[1], 'r')

        # read all the lines
        lines = f.readlines()
        for line in lines:
            # if the line is a comment, we load the corresponding files
            if line.startswith("#"):
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
    else:
        print "The file named " + sys.argv[1] + " does not exist"

    return 0

main()
