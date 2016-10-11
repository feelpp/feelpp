## Tests

The test load `data.txt`.

`feelpp_test_nofeel_interpolator`: reads `data.txt` and generates `data_interp.txt`.
Vizualization with gnuplot
```sh
gnuplot
plot `data_interp.txt` u 1:2 w l t "P0", `data_interp.txt` u 1:3 w l t "P1", `data_interp.txt` u 1:4 w l t "Spline", `data_interp.txt` u 1:5 w l t "Akima", `data_interp.txt` u 1:6 w l t "GSL"
```

`feelpp_test_interpolator`: if the data file stands for the measured relation `T - k(T)`, it evaluates `k` and `kd` (the derivative) given a `--functions.f` option.
