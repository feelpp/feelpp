Spaces and elements {#TutorialSpaces}
============================



You've learned how to discretize the space you want to compute on.
You now have to learn how to define and use function spaces and elements of functions spaces.

The source code is available in `myfunctionspace.cpp`
(The listing is given at the end).

# Loading a Mesh in 2D {#load}

We recall how to load a mesh :
\snippet myfunctionspace.cpp mesh


# Constructing a function space {#fs}

For basic function spaces, we have predetermined constructors:
\snippet myfunctionspace.cpp space

One can also use :
- `Pdh<ORDER>(mesh)` : Polynomial Discontinuous
- `Pvh<ORDER>(mesh)` : Polynomial Continuous Vectorial
- `Pdhv<ORDER>(mesh)` : Polynomial Discontinuous Vectorial
- `Pchm<ORDER>(mesh)` : Polynomial Continuous Matrix
- `Ned1h<ORDER>(mesh)` : Nedelec function spaces

# Defining an element {#elem}

Elements are basically defined and created like that :
\snippet myfunctionspace.cpp element

# Code with other features {#code}

\snippet myfunctionspace.cpp all
