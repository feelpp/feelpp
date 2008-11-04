BOOST = -I$(HOME)/local_disk/ext/boost
BOOSTSANDBOX = -I$(HOME)/local_disk/svn/boost-sandbox
BLAS = -L$(HOME)/local_disk/ext/lib -lf77blas -latlas /usr/lib/libf2c.a
LAPACK = -L$(HOME)/local_disk/ext/lapack -llapack-gfortran411 -llapack_tools
CPP = $(HOME)/local_disk/ext/bin/g++ -g
MUMPS = -I$(HOME)/local_disk/ext/mumps/MUMPS_4.6.4/include -L$(HOME)/local_disk/ext/mumps/MUMPS_4.6.4/lib -ldmumps -lsmumps -lcmumps -lzmumps -lpord -L$(HOME)/local_disk/ext/mumps/MUMPS_4.6.4/libseq -lmpiseq

mumps_ublas: mumps_ublas.cpp
	$(CPP) -o mumps_ublas mumps_ublas.cpp $(BOOST) $(BOOSTSANDBOX) $(MUMPS) $(BLAS) -lgfortran -lpthread
