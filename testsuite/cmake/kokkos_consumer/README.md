# Installed Kokkos integration consumer

After installing Feel++, activate its dependency environment and configure this
standalone project against the installation:

```bash
cmake -S testsuite/cmake/kokkos_consumer -B build/kokkos-consumer \
  -DCMAKE_PREFIX_PATH=/path/to/install \
  -DFeel++_DIR=/path/to/install/share/feelpp/feel/cmake/modules
cmake --build build/kokkos-consumer
ctest --test-dir build/kokkos-consumer --output-on-failure
```

The host executable links the ordinary CPU target. With Kokkos enabled, the
second executable compiles and links the explicit integration target and checks
that linking has not initialized Kokkos. A dependency sentinel verifies that
Kokkos compile definitions reach the integration consumer and stay out of the
ordinary CPU consumer. No parallel kernel is launched.
