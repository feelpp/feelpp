import pyfeelpp.core as core
import sys

e=core.Environment(sys.argv)

print("pid:",e.worldComm().localRank() )
print("isMasterRank:",e.isMasterRank() )

