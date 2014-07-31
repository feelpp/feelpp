## to move into the folders where Python libraries are created

#!/usr/bin/python

from mpi4py import MPI
import libPyFunctSpace
import sys


libPyFunctSpace.wrap(sys.argv)
