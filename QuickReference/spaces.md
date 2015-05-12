Function Spaces
===============





# QRFunctionSpace Function Space

We now turn to the next crucial mathematical ingredient: the function space,
whose definition depends on $$\Omega_h$$ --- or more precisely its partitioning
$$\calT_h$$ --- and the choice of basis function. Function spaces in Feel++
follow the same definition, see listing~\ref fspace1, and Feel++ provides
support for continuous and discontinuous Galerkin methods and in particular
approximations in $$L^2$$, $$H^1$$-conforming and $$H^1$$-nonconforming, $$H^2$$,
$$H(\mathrm{div})$$ and $$H(\mathrm{curl})$$\footnote{At the time of writing, $$H^2$$,
$$H(\mathrm{div})$$ and $$H(\mathrm{curl})$$ approximations are in experimental
support.

\anchor fspace1
```cpp
 // space of continuous piecewise
 // $$\P_3$$ functions defined on a mesh
 // of order 2 triangles in 3D
 FunctionSpace<Mesh<Simplex<2,2,3>,
               bases<Lagrange<3> > > Xh;
```


The `FunctionSpace`  class

 -  constructs the table of degrees of freedom which maps local (elementwise) degrees of
  freedom to the global ones with respect to the geometrical entities,
 -  embeds the definition of the elements of the function space allowing for a
  tight coupling between the elements and their function spaces,
 -  stores an interpolation data structure (\eg region tree) for rapid
  localisation of point sets (determining in which element they reside).


We introduce the following spaces
$$
  \begin{aligned}
    \mathbb{W}_h &= \{v_h \in L^2(\Omega_h): \ \forall K \in \mathcal{T}_h, v_h|_K
    \in \mathbb{P}_K\},\\
    \mathbb{V}_h &= \mathbb{W}_h \cap C^0(\Omega_h)= \{ v_h \in \mathbb{W}_h: \ \forall F \in
    \mathcal{F}^i_h\ \jump{v_h}_F = 0\}\\
    \mathbb{H}_h &= \mathbb{W}_h \cap C^1(\Omega_h)= \{ v_h \in \mathbb{W}_h: \ \forall F \in
    \mathcal{F}^i_h\ \jump{v_h}_F = \jump{\nabla v_h}_F = 0\}\\
    \CR_h &= \{ v_h \in L^2(\Omega_h):\ \forall K \in \calT_h, v_h|_K \in
    \P_1; \forall F \in \calF^i_h\ \int_F \jump{v_h} = 0 \}\\
    \RaTu_h &= \{ v_h \in L^2(\Omega_h):\ \forall K \in \calT_h, v_h|_K \in
    \Span{1,x,y,x^2-y^2}; \forall F \in \calF^i_h\ \int_F \jump{v_h} = 0 \}\\
    \RT_h&=\{\bm{v}_h \in [L^2(\Omega_h)]^d:\ \forall K \in \calT_h, v_h|_K \in
    \RT_k; \forall F \in \calF^i_h\ \jump{\bm{v}_h \cdot \normal}_F = 0 \}\\
    \N_h&=\{\bm{v}_h \in [L^2(\Omega_h)]^d:\ \forall K \in \calT_h, v_h|_K \in
    \N_k; \forall F \in \calF^i_h\ \jump{\bm{v}_h \times \normal}_F = 0 \}
  \end{aligned}
$$
where $$\RT_k$$ and $$\N_k$$ are respectively the Raviart-Thomas and N&eacute;d&eacute;lec finite
elements of degree $$k$$.

The Legrendre and Dubiner basis yield implicitely discontinuous
approximations, the Legendre and Dubiner boundary adapted basis,
see~\cite MR1696933, were designed to handle continuous approximations
whereas the Lagrange basis can yield either discontinuous or continuous
(default behavior) approximations.  $$\RT_h$$ and $$\N_h$$ are implicitely spaces
of vectorial functions $$\bm{f}$$ \st $$\bm{f}: \Omega_h \subset \R^d \mapsto
\R^d$$. As to the other basis functions, \ie Lagrange, Legrendre, Dubiner,
etc., they are parametrized by their values namely `Scalar`  ,
`Vectorial`  or `Matricial.`   Note that
`FunctionSpace`  handles also products of function spaces.  This is
very powerful to describe complex multiphysics problems when coupled with
operators, functionals and forms described in the next section. Extracting
subspaces or component spaces are part of the interface.

```cpp
// continuous piecewise P3
// approximations
FunctionSpace<Mesh<Simplex<2> >,
  bases<Lagrange<3,Scalar,
                 Continuous> > > P3ch;
// discontinuous piecewise P3
// approximations
FunctionSpace<Mesh<Simplex<2> >,
  bases<Lagrange<3,Scalar,
                 Discontinuous> > > P3dh;
// mixed (P2 vectorial, P1 scalar,
// P1 Scalar) approximation
FunctionSpace<Mesh<Simplex<2>>,
 bases<Lagrange<2,Vectorial>,
       Lagrange<1,Scalar>,
       Lagrange<1,Scalar>>> P2P1P1;
```

The most important feature in `FunctionSpace`  is that it embeds the
definition of element which allows for the strict definition of an \c
Element of a `FunctionSpace`  and thus ensures the correctness of the
code.  An element has its representation as a vector --- also in the
case of product of multiple spaces. --- The vector representation is
parametrized by one of the linear algebra backends. Other supported
operations are interpolation and extraction of components --- be it a
product of function spaces element or a vectorial/matricial element:

```cpp
FunctionSpace<Mesh<Simplex<2> >,
  bases<Lagrange<3,Scalar, Continuous> > > P3ch;
// get an element from P3ch
auto u = P3ch.element();
FunctionSpace<Mesh<Simplex<2> >,
 bases<Lagrange<2,Vectorial>, Lagrange<1,Scalar>,
       Lagrange<1,Scalar> > > P2P1P1;
auto U = P2P1P1.element();
// Views: changing a view changes U and vice versa
// view on element associated to P2
auto u = U.element<0>();
// extract view of first component
auto ux = u.comp(X);
// view on element associated to 1st P1
auto p = U.element<1>();
// view on element associated to 2nd P1
auto q = U.element<2>();
```

Finally Feel++ provides the Lagrange, $$\Iclag, \Idlag$$, Crouzeix-Raviart, $$\Icr$$,
Raviart-Thomas, $$\Irt$$ and N&eacute;d&eacute;lec, $$\In$$ global interpolation operators.
In abstract form, they read
$$
  \calI : \X \ni v \mapsto \sum_{i=1}^{\opdim\X} \ell_i(v) \phi_i
$$
where $$\X$$ is the infinite dimensional space, $$(\ell_i)_{i=1,...,\opdim\X}$$ are
the linear forms and $$(\phi_i)_{i=1...\opdim\X}$$ the basis function associated
with the various approximations.

# Function Space helper functions

Function        | Description
----------------|------------------------------
`Pch<N>(mesh)`  | generates $$P^N_{c,h}$$
`Pchv<N>(mesh)` | generates $$[P^N_{c,h}]^d$$
`THch<N>(mesh)` | generates $$[P^{N+1}_{c,h}]^d \times P^N_{c,h}$$

