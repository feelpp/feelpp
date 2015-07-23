#!/usr/bin/python

import argparse
import json
import os

params = {
    "dry_run": False,
    "verbose": False
}

def check(mesh):
    # Check that the geo file exists
    if(not os.path.exists(mesh["geo"])):
        return 1

    # Check that we have the required keywords
    if(("geo" not in mesh) or ("backend" not in mesh)):
        return 1

    return 0;

def generateGMSH(mesh):
    # get filename
    out = os.path.basename(mesh["geo"])
    name, ext = os.path.splitext(out)

    # Start building commands
    lstb = []
    lstb.append("gmsh " + mesh["geo"] + " -3")
    lstc = lstb

    # Check if we input the hmin/hmax
    #if "h" in mesh:
        #lstb = []
        ## we iterate over the previously built commands
        #for cmd in lstc:
            #c = ""
            ## we iterate over the different hmin, hmax elements
            #for i in range(len(mesh["h"])):
                ## if we specified only one h component, consider it as hmin
                #if(len(mesh["h"][i]) == 1):
                    #c = cmd + " -clmin " + str(mesh["h"][i][0])
                #if(len(mesh["h"][i]) == 2):
                    #c = cmd + " -clmin " + str(mesh["h"][i][0])
                    #c = c + " -clmax " + str(mesh["h"][i][1])
                #lstb.append(c)

        #lstc = lstb

    if "partitions" in mesh:
        lstb = []
        #cmd = cmd + " -part " + mesh["partitions"]
        # we iterate over the previously built commands
        for cmd in lstc:
            c = ""
            # we iterate over the different partition numbers
            for i in range(len(mesh["partitions"])):
                c = cmd + " -part " + str(mesh["partitions"][i])
                lstb.append(c)

        lstc = lstb

    # set output filename
    for cmd in lstc:
        if(params["verbose"] or params["dry_run"]):
            print cmd
        if(not params["dry_run"]):
            os.system(cmd)

def generateMG(mesh):
    # Get filename
    out = os.path.basename(mesh["geo"])
    name, ext = os.path.splitext(out)

    # Mesh in 2d and convert to mesh format
    cmd = "gmsh " + mesh["geo"]
    cmd = cmd + " -2"

    # set output filename
    cmd = cmd + " -format mesh -o " + name + ".mesh"
    if(params["verbose"]):
        print cmd

    # Generate tetrahedralization with meshgems
    cmd = "run_mg-tetra.sh" + " -i " + name + ".mesh"
    if(params["verbose"]):
        print cmd

    cmd = "gmsh" + name + "_tetra.mesh" + " -3 -format msh -o " + name + ".msh"
    if(params["verbose"]):
        print cmd

def configureGEO(mesh):
    # get file info
    dirname = os.path.dirname(mesh["geo"])
    basename = os.path.basename(mesh["geo"])
    name, ext = os.path.splitext(basename)

    # build filename
    filename = name
    if("name" in mesh):
        filename = filename + "_" + str(mesh["name"]) + "_"
    if("h" in mesh):
        filename = filename + "_hmin" + str(mesh["h"][0]) + "_hmax" + str(mesh["h"][1])
    if("reynolds" in mesh):
        filename = filename + "_reynolds" + str(mesh["reynolds"])
    if("lo" in mesh):
        filename = filename + "_lo" + str(mesh["lo"][0])
    filename = filename

    # Check if we did not already create this file before
    # ensure that a new file is created
    i = 0
    while(os.path.exists(filename + "." + str(i) + ext)):
        i = i + 1

    filename = filename + "." + str(i) + ext

    if(params["verbose"]):
        print "Creating " + filename

    fin = open(mesh["geo"], "r")
    fout = open(filename, "w")

    for line in fin.readlines():
        cleanline = line.translate(None, "\n").lstrip().rstrip()
        l = cleanline.split("=")

        newln = ""

        # Set the hmin/hmax value of points
        # if we find a point
        if(l[0].find("Point(") == 0):
            # Get the id of the point
            pid = int(l[0].translate(None, "Point()"))
            # The following only need one replace in cas we are in the wrong case
            # if this is a point that we have to set to hmin, we do it
            if(pid in mesh["hminPoints"]):
                newln = l[0] + "=" + l[1].replace("hmax", "hmin")
            # otherwise we set it to hmax
            else:
                newln = l[0] + "=" + l[1].replace("hmin", "hmax")
            # output line to file
            if(cleanline != newln and params["verbose"]):
                print "\"" + cleanline + "\""
                print "\"" + newln + "\""
                print "Changed line " + l[0]
        # change hmin
        elif(l[0].find("hmin") == 0):
            newln = "hmin = " + str(mesh["h"][0]) + ";"
        # change hmax
        elif(l[0].find("hmax") == 0):
            newln = "hmax = " + str(mesh["h"][1]) + ";"
        # change hmax
        elif(l[0].find("lo") == 0):
            newln = "lo = " + str(mesh["lo"][0]) + ";"
        # by default we output the line to the file
        else:
            newln = cleanline

        fout.write(newln + "\n")

    fin.close()
    fout.close()

    return filename

def partition(mesh):
    print "test"

def main():
    # generate 3d mesh
    # gmsh -3 ../../../../benchmarks/navierstokes/nozzle/fda-3d.geo -o fda-3d.msh

    # argument parser
    parser = argparse.ArgumentParser(description='Generate FDA meshes with json configuration')
    parser.add_argument('-n', '--dry-run', action="store_true", help='dry run')
    parser.add_argument('-v', '--verbose', action="store_true", help='verbose mode')

    args = parser.parse_args()

    # configure app params
    params["dry_run"] = args.dry_run
    params["verbose"] = args.verbose

    f = open("mesh.json", "r")
    j = json.load(f)
    meshes = j["meshes"]

    print meshes
    for i in range(len(meshes)):
        if(check(meshes[i]) == 0):
            filename = configureGEO(meshes[i])
            print "Created new geo file " + filename

            # Create a copy of the subdictionary and fill it with the new geo file
            newMesh = meshes[i]
            newMesh["geo"] = filename

            if("gmsh" in newMesh["backend"]):
                if(params["verbose"]):
                    print "Generating using GMSH"
                generateGMSH(newMesh)
            if("meshgems" in newMesh["backend"]):
                if(params["verbose"]):
                    print "Generating using Meshgems"
                generateMG(newMesh)
        else:
            print "Entry " + str(i) + " is invalid"

    f.close()

# Check that we are not using the script as a module
if __name__ == '__main__':
    main()
