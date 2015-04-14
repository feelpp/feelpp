
Defining a Model {#TutorialModel}
================

[TOC]

# Introduction {#Intro}

It is well known an equation can describe a huge range of physical problems.
Each of theses problems will have a particular environment, but the equation to solve will be the same.
To make our program applicable to theses range of problem, we have defined a model.

# What is a model {#What}

A model is defined by :
- a Name
- a Description
- a Modèle
- Parameters
- Materials
- Boundary Conditions
- Post Processing

## Parameters {#Parameters}
A parameter is 

## Materials {#Materials}

## BoundaryConditions {#BoundaryConditions}
Thanks to GiNaC, we handle boundary conditions (Dirichlet, Neumann, Robin) as expression.
You have to indicate in the json file the quantity to handle (velocity, pressure...) and the associated expression.

## PostProcessing {#PostPro}

# Example {#Example}
We have set up an example : an anisotropic laplacian.

\snippet aniso_laplacian.cpp global

