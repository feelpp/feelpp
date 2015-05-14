Defining and using expressions {#TutorialExpr}
================================



The next step is to construct a function space over the mesh. The source code is
available in `myexpression.cpp.`

# Step by step explanations {#ex}

We start by loading a Mesh in 2D
!CODEFILE "code/myintegrals.cpp" mesh

then we define some expression through the command line or config file: `g`  is a scalar field and `f`  is a vector field
!CODEFILE "code/myexpression.cpp" expr

here is an example how to enter them, more are available below
```c++
feelpp_doc_myexpression --a=3 --functions.g="a*x*y:x:y:a" --functions.f="{sin(pi*x),cos(pi*y)}:x:y"
```

You can print back the expression to the screen to check that everything is ok.
You want to use as expression `a*x+b*y`, you have to define `a` and `b` as option (either in your code, either in the library).

then we compute the gradient of `g`  and `f`
!CODEFILE "code/myexpression.cpp" grad

Notice that template argument are given to `grad`  to specify the shape of the
gradient: in the case of $$\nabla g$$ it is $$1\times2$$ and $$\nabla f$$
$$2\times 2$$ since we are in 2D.

then we compute the laplacian of `g`  and `f`
!CODEFILE "code/myexpression.cpp" laplacian

then we compute the divergence of `f`
!CODEFILE "code/myexpression.cpp" div

and the curl of `f`
!CODEFILE "code/myexpression.cpp" curl

Finally we evaluate these expression at one point given by the option `x`  and `y`
!CODEFILE "code/myexpression.cpp" eval

# Some results {#res}

We start with the following function $$g=1$$ and $$f=(1,1)$$.

```bash
shell> ./feelpp_doc_myexpression --functions.g=1:x:y --functions.f="{1,1}:x:y
g=1
f={x,-y}
grad(g)=[[0]]
grad(f)=[[1,0],[0,-1]]
laplacian(g)=[[0]]
laplacian(f)=[[0],[0]]
div(f)=[[0]]
curl(f)=[[0]]
Evaluation  at  (0,0):
           g(x,y)=1
           f(x,y)= 0
-0
Gradient:
     grad(g)(x,y)= 0 -0
     grad(f)(x,y)= 1  0
 0 -1
Divergence:
      div(f)(x,y)=0
Curl:
     curl(f)(x,y)=-3.14159
Laplacian:
laplacian(g)(x,y)=0
laplacian(f)(x,y)=0
0
```

The symbolic calculus system worked as expected.


# Complete code {#code}

The complete code reads as follows

!CODEFILE "code/myexpression.cpp" all

to compile just use the `make` command in your compilation directory
```bash
make feelpp_doc_myexpression
```
