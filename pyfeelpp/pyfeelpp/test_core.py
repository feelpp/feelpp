import core
import sys,time

e=core.Environment(sys.argv)
print("pid:",e.worldComm().localRank() )

