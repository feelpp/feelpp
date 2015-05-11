Using function {#TutorialFunctions}
=====================

[TOC]

Once you  have created an element, you may want to give it a value, that can depends on a lot of parameters (mainly spaces, but other may apply).

To do so, Feel++ relies on expression.
We may use tree kind of expression :

- [Parsed](#parsed)
- [Build-in](#build-in)
- [Hard coded](#hc)


# Parsed {#parsed}

Thanks to [GiNaC](http://www.ginac.de), we can parse expression like that :
```sh
./feelpp_myexpression --functions.f="2*x*y+cos(x+y):x:y"
```

Step by step explanations
------------ 

We start by loading a Mesh in 2D
\snippet myintegrals.cpp mesh

then we define some expression through the command line or config file: \c g is a scalar field and \c f is a vector field
\snippet myexpression.cpp expr

here is an example how to enter them, more are available below
```c++
feelpp_doc_myexpression --a=3 --functions.g="a*x*y:x:y:a" --functions.f="{sin(pi*x),cos(pi*y)}:x:y"
```

\remark You can print back the expression to the screen to check that everything is ok.
\remark You want to use as expression `a*x+b*y`, you have to define `a` and `b` as option (either in your code, either in the library).

then we compute the gradient of \c g and \c f
\snippet myexpression.cpp grad

Notice that template argument are given to \c grad to specify the shape of the
gradient: in the case of $$\nabla g$$ it is $$1\times2$$ and $$\nabla f$$
$$2\times 2$$ since we are in 2D.

then we compute the laplacian of \c g and \c f
\snippet myexpression.cpp laplacian

then we compute the divergence of \c f
\snippet myexpression.cpp div

and the curl of \c f
\snippet myexpression.cpp curl

Finally we evaluate these expression at one point given by the option \c x and \c y
\snippet myexpression.cpp eval

# Build-in {#build-in}

Instead of defining an expression from a string, you can use
\snippet myexporter.cpp expr
The list of the Feel++ Keyword is [here](Keywords.html).

# Hard Coded {#hc}

One other method to define function is described here.

\snippet myfunctor.cpp all

