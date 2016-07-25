# Two implementation for the magnetostatic problem

For the moment, only strong bc are available

## AMS
We solve here the regularized problem: curl(1/mu curl (u)) + regul u = j
The command line to make it run is:
```sh
./feelpp_test_ams --config-files cube.cfg backend_stab.cfg
```

## Saddle
We solve here the regularized problem: curl(1/mu curl (u)) + grad(p) = j; div(u) = 0
The command line to make it run is:

```sh
./feelpp_test_saddle --config-files cube.cfg backend_stab.cfg
```
