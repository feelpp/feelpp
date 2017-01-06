<!-- -*- mode: markdown  -->

Feel++ TestSuite
===================

# Machine IRMA-ATLAS

## Feel++ Version
Wed Jan 4 09:01:32 2017 +0100 (commmit d204f6503b4483e455249067fbf6e7b3ac604ba7)

### Environment
Loaded Modulefiles:  
```gcc/4.9.0```  
```openmpi/1.10```  
```boost/1.59```  
```petsc/3.6.3```  
```slepc/3.6.3```  
```gmsh/2.10.1```  
```hdf5/1.8.15-patch1```  
```VTK/5.10.1```  
```ParaView/5.0.1```  
```fftw/3.3.4```  
```CMake/3.5.2```  

### CMake options
```
-DCMAKE_BUILD_TYPE=Release
-DCMAKE_CXX_COMPILER=clang++-3.8 -DCMAKE_C_COMPILER=clang-3.8 -DFEELPP_ENABLE_GSL=OFF
```

### Results
96% tests passed, 14 tests failed out of 333  

Label Time Summary:  
testalg              =  38.68 sec (22 tests)  
testcore             =  56.22 sec (36 tests)  
testcrb              =  61.31 sec (11 tests)  
testdiscr            = 431.53 sec (94 tests)  
testfilters          =  32.81 sec (10 tests)  
testfit              =   6.67 sec (4 tests)  
testintegration      =  30.11 sec (14 tests)  
testinterpolation    =  96.34 sec (25 tests)  
testleaks            =  30.83 sec (16 tests)  
testls               =   6.17 sec (3 tests)  
testmaterial         =   2.72 sec (2 tests)  
testmath             =   2.70 sec (2 tests)  
testmesh             =  63.41 sec (24 tests)  
testmodels           =   3.92 sec (2 tests)  
testopt              =   2.70 sec (2 tests)  
testpde              =  13.65 sec (4 tests)  
testpoly             =  45.15 sec (18 tests)  
testts               =  10.80 sec (5 tests)  
testvf               = 434.79 sec (38 tests)  

Total Test time (real) = 1372.59 sec  

The following tests FAILED:  
	249 - feelpp\_test\_element0D-np-6 (Failed)  
	250 - feelpp\_test\_element0D-np-1 (Failed)  
	264 - feelpp\_test\_bdf2-np-6 (Failed)  
	265 - feelpp\_test\_bdf2-np-1 (Failed)  
	266 - feelpp\_test\_bdf3-np-1 (Failed)  
	273 - feelpp\_test\_laplacianv-np-6 (Failed)  
	274 - feelpp\_test\_laplacianv-np-1 (Failed)  
	275 - feelpp\_test\_laplacian-np-6 (Failed)  
	276 - feelpp\_test\_laplacian-np-1 (Failed)  
	277 - feelpp\_test\_laplaciant-np-6 (Failed)  
	278 - feelpp\_test\_laplaciant-np-1 (Failed)  
	289 - feelpp\_test\_convolve-np-6 (Timeout)  
	330 - feelpp\_test\_nofeel\_interpolator-np-6 (Failed)  
	331 - feelpp\_test\_nofeel\_interpolator-np-1 (Failed)  