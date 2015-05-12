Defining a Model {#TutorialModel}
================



# Introduction {#Intro}

It is well known an equation can describe a huge range of physical
problems.  Each of theses problems will have a particular environment,
but the equation to solve will be the same.  To make our program
applicable to theses range of problem, we have defined a model.

# What is a model {#What}

A model is defined by :
- a Name
- a Description
- a Model
- [Parameters](@ref Parameters)
- [Materials](@ref Materials)
- [Boundary Conditions](@ref BoundaryConditions)
- [Post Processing](@ref PostPro)

## Parameters {#Parameters}
A parameter is a non physical property for a model.

## Materials {#Materials}
To retrieve the materials properties, we use :
\snippet aniso_laplacian.cpp get_mat
And to apply them :
\snippet aniso_laplacian.cpp materials

## BoundaryConditions {#BoundaryConditions}
Thanks to GiNaC, we handle boundary conditions (Dirichlet, Neumann, Robin) as expression.
You have to indicate in the json file the quantity to handle (velocity, pressure...) and the associated expression.
\snippet aniso_laplacian.cpp get_bc
We can apply theses boundary condition this way
\snippet aniso_laplacian.cpp boundary

## PostProcessing {#PostPro}

# Example {#Example}
We have set up an example : an anisotropic laplacian.

\snippet aniso_laplacian.cpp global
