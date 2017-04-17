# Building

## Full install.

To compile all librairies just type

```
./install.sh
```

You can modify the profile name.

It will create 3 directories
- src/          compile libraries here (clean after install)
- install/      install libraries here
- archive/      keep an archive of install here


## Install per library

To install by hand a library do for example for petsc

```
./petsc.sh gcc630
```
where gcc630 is the profile name.
To remove files

```
./petsc.sh -r gcc630
```



