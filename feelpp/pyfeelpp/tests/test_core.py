import feelpp 
import sys

print(sys.argv)
e=feelpp.Environment(sys.argv)

print("pid:",e.worldComm().localRank() )
print("isMasterRank:",e.isMasterRank() )

