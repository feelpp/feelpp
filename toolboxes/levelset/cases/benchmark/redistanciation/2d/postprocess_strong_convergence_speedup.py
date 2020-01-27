#!/usr/bin/env python3

if __name__ == '__main__':

    import csv
    import sys
    
    dictrows=[]
    arglength = len(sys.argv)
    sourcefilename = sys.argv[1]
    targetfilename = sourcefilename.replace("time_per_core", "speedup_per_core")
    
    with open(sourcefilename, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            dictrows.append(row)

    # We only compute speedup if we find sequential time
    if dictrows[0]["processes"] == "1":
        sequentialtime = float(dictrows[0][" time"])
        with open(targetfilename, mode='w') as targetcsv:
            targetcsv.write("processes, speedup, incrementalspeedup, linearincrementalspeedup\n")
            targetcsv.write("1, 1.0, , \n" )# first line
            lasttime = sequentialtime
            lastnproc = 1
            dictrows.pop(0)
            for d in dictrows:
                currentnproc = int(d["processes"])
                currenttime = float(d[" time"])
                speedup = sequentialtime / currenttime
                incrementalspeedup = lasttime / currenttime
                linearincrementalspeedup = currentnproc / lastnproc
                targetcsv.write(str(currentnproc) + ", " + str(speedup) + ", " + str(incrementalspeedup) + ", " + str(linearincrementalspeedup) + "\n" )
                lasttime = currenttime
                lastnproc = currentnproc
    else:
        print(__file__ + ": could not find sequential time in " + sourcefilename)
    
