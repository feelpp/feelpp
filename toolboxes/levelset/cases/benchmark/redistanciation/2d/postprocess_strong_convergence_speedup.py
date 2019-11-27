#!/usr/bin/env python3

if __name__ == '__main__':

    import csv
    import sys
    
    dictrows=[]
    arglength = len(sys.argv)
    #print("length = " + str(arglength))
    sourcefilename = sys.argv[1]
    targetfilename = sourcefilename.replace("time_per_core", "speedup_per_core")
    
    #print("source filename = " + sourcefilename)
    #print("target filename = " + targetfilename)
    
    with open(sourcefilename, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            #print(row)
            dictrows.append(row)

    # We only compute speedup if we find sequential time
    if dictrows[0]["processes"] == "1":
        sequential_time = float(dictrows[0][" time"])
        with open(targetfilename, mode='w') as targetcsv:
            targetcsv.write("processes, speedup\n")
            for d in dictrows:
                targetcsv.write(d["processes"] + ", " + str( sequential_time / float(d[" time"])) + "\n" )

    else:
        print("Could not find sequential time in " + sourcefilename)
    #print("sequential time is : " + str(sequential_time))

