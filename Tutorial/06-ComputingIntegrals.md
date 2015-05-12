Computing integrals over mesh {#TutorialIntegrals}
==================================

[TOC]

The next step is to compute integrals over the mesh. The source code is available in `myintegrals.cpp.`

# Step by step explanations {#TutorialIntegralsSteps}

We start by loading a Mesh in 2D
\snippet myintegrals.cpp mesh

then we define the Feel++ expression that we are going to integrate using the \c
soption function that retrieves the command line option string `functions.g.`  We then transform this string into a Feel++ expression using `expr().`

\snippet myintegrals.cpp expression

then We compute two integrals over the domain and its boundary respectively

 * $$\int_\Omega g$$

 * $$\int_{\partial \Omega} g$$

and we print the results to the screen.

\snippet myintegrals.cpp integrals

Only the rank 0 process (thanks to `Environment` isMasterRank() prints
to the screen as the result is the same over all mpi processes if the
application was run in parallel. Note also that the code actually prints the
expression passed by the user through the command line option `functions.g.`

# Some results {#TutorialIntegralsResults}

We start with the following function $$g=1$$. Recall that by default the
domain is the unit square, hence the $$\int_\Omega g$$ and $$\int_{\partial
\Omega} g$$ should be equal to 1 and 4 respectively.

```
shell> ./feelpp_doc_myintegrals --functions.g=1
int_Omega 1 = 1
int_{boundary of Omega} 1 = 4
```

Now we try $$g=x$$. We need to tell Feel++ what are the symbols associated
with the expression: here the symbol is `x`  and it works as follows

```
shell> ./feelpp_doc_myintegrals --functions.g=x:x
int_Omega x = 0.5
int_{boundary of Omega} x = 2
```

Recall that \c : is a separator between the expression and each symbol composing it.

Now we try $$g=x y$$. We need to tell Feel++ what are the symbols associated
with the expression: here the symbol is `x`  and `y`  and it works as follows

```
shell> ./feelpp_doc_myintegrals --functions.g=x*y:x:y
int_Omega x*y = 0.25
int_{boundary of Omega} x*y = 1
```

More complicated functions are of course doable, such as $$g=\sin( x y )$$.

```
./feelpp_doc_myintegrals --functions.g="sin(x*y):x:y"
int_Omega sin(x*y) = 0.239812
int_{boundary of Omega} sin(x*y) = 0.919395
```

Here is the last example in parallel over 4 processors which returns, of
course, the exact same results as in sequential

```
mpirun -np 4 ./feelpp_doc_myintegrals --functions.g="sin(x*y):x:y"
int_Omega sin(x*y) = 0.239812
int_{boundary of Omega} sin(x*y) = 0.919395
```

Finally we can change the type of domain and compute the area and perimeter of the unit disk as follows
```
./feelpp_doc_myintegrals --functions.g="1:x:y" --gmsh.domain.shape=ellipsoid --gmsh.hsize=0.05
int_Omega 1 = 0.784137
int_{boundary of Omega} 1 = 3.14033
```

Note that we don't get the exact results due to the fact that
$$\Omega_h = \cup_{K \in \mathcal{T}_h} K$$ which we use for the numerical
integration is different from the exact domain $$\Omega = \{ (x,y)\in
\mathbb{R}^2 | x^2+y^2 < 1\}$$.

#  Complete code {#TutorialIntegralsCode}

The complete code reads as follows

\snippet myintegrals.cpp all

to compile just type
```
make feelpp_doc_myintegrals
```
