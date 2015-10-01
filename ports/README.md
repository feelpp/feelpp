L'ordre de compilation doit être respecté

Gmsh ne dépend ni de Petsc ni d'openMPI !

## OpenMPI 
````
./configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 --prefix=/data/software/install/openmpi-1.10.0
make -j all install
```

