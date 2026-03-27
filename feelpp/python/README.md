# feelpp-core

Python bindings for [Feel++](https://www.feelpp.org), a C++ library for solving partial differential equations using finite element methods. Feel++ provides a domain-specific embedded language (DSEL) in C++ that closely mirrors the mathematical formulation of variational problems, making it accessible for scientists and engineers working in computational mechanics, fluid dynamics, heat transfer, and more.

## Installation

```bash
pip install feelpp-core
```

**Requirements:** Linux x86_64 with glibc >= 2.39 (Ubuntu 24.04+, Debian 13+, Fedora 40+).

## Quick Start

```python
import feelpp.core as fppc
import sys

# Initialize Feel++ environment
e = fppc.Environment(sys.argv, opts=fppc.backend_options("laplacian"))

# Create a 2D mesh
m = fppc.mesh(dim=2, geo=1, realdim=2)

# Create a P1 continuous function space
Vh = fppc.functionSpace(mesh=m, space="Pch", order=1)

print(f"Feel++ {fppc.__version__}")
print(f"Mesh elements: {m.numGlobalElements()}")
print(f"DOFs: {Vh.nDof()}")
```

## Features

- Meshes: simplex and hypercube in 1D, 2D, 3D with geometric orders 1 and 2
- Function spaces: Pch (continuous), Pdh (discontinuous), Pchv (vector continuous)
- Variational formulations, integration, and interpolation
- Parallel computing via MPI
- Linear and nonlinear solvers (PETSc/SLEPc backend)
- Data export (Ensight, VTK, XDMF/HDF5)
- Remote data management

## Optional Dependencies

```bash
pip install feelpp-core[gmsh]   # For mesh generation via Gmsh
pip install feelpp-core[test]   # For running tests
```

## Documentation

Full documentation: [docs.feelpp.org](https://docs.feelpp.org)

## License

LGPL-3.0-or-later. See [LICENSE](https://github.com/feelpp/feelpp/blob/develop/LICENSE).
