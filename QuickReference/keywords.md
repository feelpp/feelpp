tr<Keywords
========

One of Feel++ assets is it finite element embedded language. The language follows the C++ grammar, and provides keywords as well as operations between objects which are, mathematically, tensors of rank 0, 1 or 2.

In all following tables we use the notations:<br>
$$f: \mathbb{R}^n \mapsto \mathbb{R}^{m\times p}$$  with $$n=1,2,3$$, $$m=1,2,3$$, $$p=1,2,3$$<br>
$$\Omega^e$$ current mesh element.

The genesis of the language can be found in \cite prudhomme05:_domain_specif_embed_languag_c and an update on Feel++ is available in \cite PRUDHOMME:2012:HAL-00662868:3.

# Keywords_Points Points

## Current Point:

|Feel++ Keyword | Math Object | Description | Dimension|
|---------------|---------------------------------------
|```cpp{.cpp} P()```|$$\overrightarrow{P}$$| $$(P_x, P_y, P_z)^T$$|$$d \times 1$$|
|```cpp{.cpp} Px()```|$$P_x$$|$$x$$ coordinate of $$\overrightarrow{P}$$|$$1 \times 1$$|
|```cpp{.cpp} Py()```|$$P_y$$|$$y$$ coordinate of $$\overrightarrow{P}$$<br>(value is 0 in 1D)|$$1 \times 1$$|
|```cpp{.cpp} Pz()```|$$P_z$$|$$z$$ coordinate of $$\overrightarrow{P}$$<br>(value is 0 in 1D and 2D)|$$1 \times 1$$|


## Element Barycenter Point:

|Feel++ Keyword | Math Object | Description | Dimension | 
|```cpp{.cpp} C()```|$$\overrightarrow{C}$$| $$(C_x, C_y, C_z)^T$$|$$d \times 1$$|
|```cpp{.cpp} Cx()```|$$C_x$$|$$x$$ coordinate of $$\overrightarrow{C}$$|$$1 \times 1$$|
|```cpp{.cpp} Cy()```|$$C_y$$|$$y$$ coordinate of $$\overrightarrow{C}$$<br>(value is 0 in 1D)|$$1 \times 1$$|
|```cpp{.cpp} Cz()```|$$C_z$$|$$z$$ coordinate of $$\overrightarrow{C}$$<br>(value is 0 in 1D and 2D)|$$1 \times 1$$|


Normal at Current Point:
<table class="manual">
|Feel++ Keyword | Math Object | Description | Dimension | 
|```cpp{.cpp} N()```|$$\overrightarrow{N}$$| $$(N_x, N_y, N_z)^T$$|$$d \times 1$$|
|```cpp{.cpp} Nx()```|$$N_x$$|$$x$$ coordinate of $$\overrightarrow{N}$$|$$1 \times 1$$|
|```cpp{.cpp} Ny()```|$$N_y$$|$$y$$ coordinate of $$\overrightarrow{N}$$<br>(value is 0 in 1D)|$$1 \times 1$$|
|```cpp{.cpp} Nz()```|$$N_z$$|$$z$$ coordinate of $$\overrightarrow{N}$$<br>(value is 0 in 1D and 2D)|$$1 \times 1$$|



<a href="#" class="top">top</a>
<hr>
<br>
# Keywords_Array Vectors and Matrix
## BuildingVectors Building Vectors
Usual syntax to create vectors:
<table class="manual">
|Feel++ Keyword | Math Object | Description | Dimension | 
|```cpp{.cpp} vec<n>(v_1,v_2,...,v_n)```|$$\begin{pmatrix} v_1\\v_2\\ \vdots \\v_n \end{pmatrix}$$|Column Vector with $$n$$ rows<br>entries being expressions|$$n \times 1$$|

You can also use expressions and the unit base vectors:
<table class="example">
|Feel++ Keyword | Math Object | Description | 
|\c oneX() | $$\begin{pmatrix} 1\\0\\0 \end{pmatrix}$$|Unit vector $$\overrightarrow{i}$$|
|\c oneY() | $$\begin{pmatrix} 0\\1\\0 \end{pmatrix}$$|Unit vector $$\overrightarrow{j}$$|
|\c oneZ() | $$\begin{pmatrix} 0\\0\\1 \end{pmatrix}$$|Unit vector $$\overrightarrow{k}$$|



## Building Matrix
<table class="manual">
|Feel++ Keyword | Math Object | Description | Dimension | 
|```cpp{.cpp} mat<m,n>(m_11,m_12,...,m_mn)```|$$\begin{pmatrix} m_{11} & m_{12} & ...\\ m_{21} & m_{22} & ...\\ \vdots & & \end{pmatrix}$$|$$m\times n$$ Matrix<br> entries beeing expressions |$$m \times n$$|
|```cpp{.cpp} ones<m,n>()```|$$\begin{pmatrix} 1 & 1 & ...\\ 1 & 1 & ...\\ \vdots & & \end{pmatrix}$$|$$m\times n$$ Matrix <br>Filled with 1 |$$m \times n$$|
|```cpp{.cpp} zero<m,n>()```|$$\begin{pmatrix} 0 & 0 & ...\\ 0 & 0 & ...\\ \vdots & & \end{pmatrix}$$|$$m\times n$$ Matrix <br>Filled with 0 |$$m \times n$$|
|```cpp{.cpp} constant<m,n>(c)```|$$\begin{pmatrix} c & c & ...\\ c & c & ...\\ \vdots & & \end{pmatrix}$$|$$m\times n$$ Matrix <br>Filled with a constant c |$$m \times n$$|
|```cpp{.cpp} eye<n>()```|$$\begin{pmatrix} 1 & 0 & ...\\ 0 & 1 & ...\\ \vdots & & \end{pmatrix}$$|Unit diagonal Matrix <br> of size$$n\times n$$ |$$n \times n$$|
|```cpp{.cpp} Id<n>()```|$$\begin{pmatrix} 1 & 0 & ...\\ 0 & 1 & ...\\ \vdots & & \end{pmatrix}$$|Unit diagonal Matrix <br> of size$$n\times n$$ |$$n \times n$$|


## Manipulating Vectors and Matrix
Let $$A$$ be a square matrix of size $$n$$.
<table class="manual">
|Feel++ Keyword | Math Object | Description | Dimension | 
|```cpp inv(A)```|$$A^{-1}$$|Inverse of matrix $$A$$ |$$n \times n$$|
|```cpp det(A)```|$$\det (A)$$|Determinant of matrix $$A$$ |$$1 \times 1$$|
|```cpp sym(A)\eco|$$\text{Sym}(A)$$|Symmetric part of matrix $$A$$: $$\frac{1}{2}(A+A^T)$$<br> |$$n \times n$$|
|```cpp antisym(A)\eco|$$ \text{Asym}(A)$$|Antisymmetric part of  $$A$$: $$\frac{1}{2}(A-A^T)$$<br> |$$n \times n$$|


Let A and B be two matrix (or two vectors) of same dimension $$m \times n$$.
<table class="manual">
|Feel++ Keyword | Math Object | Description | Dimension | 
|```cpp trace(A)\eco|$$\text{tr}(A)$$|Trace of matrix $$A$$<br>Generalized on non-squared Matrix<br>Generalized on Vectors |$$1 \times 1$$|
|```cpp trans(B)\eco|$$B^T$$|Transpose of matrix $$B$$<br>Can be used on non-squared Matrix<br>Can be used on Vectors |$$n \times m$$|
|```cpp inner(A,B)\eco|$$ A.B \\ A:B = \text{tr}(A*B^T)$$|Scalar product of two vectors<br>Generalized scalar product of two matrix |$$1 \times 1$$|
|```cpp cross(A,B)\eco|$$ A\times B$$|Cross product of two vectors|$$n \times 1$$|





<a href="#" class="top">top</a>
<hr>
<br>
# Keywords_Expr Expressions
Following tables present tools to declare and manipulate expressions.
<table class="example">
|Feel++ Keyword | Description | 
|```cpp{.cpp}
Px()
Py()
Pz()
cst( c )
```|
Variable $$x$$<br>
Variable $$y$$<br>
Variable $$z$$<br>
Constant function equal to $$c$$
|


You can of course use all current operators ( + - / * ) and the usual following functions:
<table class="example">
|Feel++ Keyword | Math Object | Description</th>
|```cpp{.cpp} abs(expr) ```|$$|f(\overrightarrow{x})|$$|element wise absolute value of $$f$$|
|```cpp{.cpp} cos(expr)```|$$\cos(f(\overrightarrow{x}))$$|element wise cos value of $$f$$|
|```cpp{.cpp} sin(expr)```|$$\sin(f(\overrightarrow{x}))$$|element wise sin value of $$f$$|
|```cpp{.cpp} tan(expr)```|$$\tan(f(\overrightarrow{x}))$$|element wise tan value of $$f$$|
|```cpp{.cpp} acos(expr)```|$$\acos(f(\overrightarrow{x}))$$|element wise acos value of $$f$$|
|```cpp{.cpp} asin(expr)```|$$\asin(f(\overrightarrow{x}))$$|element wise asin value of $$f$$|
|```cpp{.cpp} atan(expr)```|$$\atan(f(\overrightarrow{x}))$$|element wise atan value of $$f$$|
|```cpp{.cpp} cosh(expr)```|$$\cosh(f(\overrightarrow{x}))$$|element wise cosh value of $$f$$|
|```cpp{.cpp} sinh(expr)```|$$\sinh(f(\overrightarrow{x}))$$|element wise sinh value of $$f$$|
|```cpp{.cpp} tanh(expr)```|$$\tanh(f(\overrightarrow{x}))$$|element wise tanh value of $$f$$|
|```cpp{.cpp} exp(expr)```|$$\exp(f(\overrightarrow{x}))$$|element wise exp value of $$f$$|
|```cpp{.cpp} log(expr)```|$$\log(f(\overrightarrow{x}))$$|element wise log value of $$f$$|
|```cpp{.cpp} sqrt(expr)```|$$\sqrt{f(\overrightarrow{x})}$$|element wise sqrt value of $$f$$|
|```cpp{.cpp} ceil(expr)```|$$\lceil{f(\overrightarrow{x})}\rceil$$|element wise ceil of $$f$$|
|```cpp{.cpp} floor(expr)```|$$\lfloor{f(\overrightarrow{x})}\rfloor$$|element wise floor of $$f$$|
|```cpp{.cpp} sign(expr)```|$$\begin{cases} 1 & \text{if}\ f(\overrightarrow{x}) \geq 0\\-1 & \text{if}\ f(\overrightarrow{x}) < 0\end{cases}$$|element wise sign value of $$f$$|
|```cpp chi(expr)```|$$\chi(f(\overrightarrow{x}))=\begin{cases}0 & \text{if}\ f(\overrightarrow{x}) = 0\\1 & \text{if}\ f(\overrightarrow{x}) \neq 0\\\end{cases}$$|element wise boolean test of $$f$$|


<a href="#" class="top">top</a>
<hr>
<br>
# Keywords_Operators Operators
## Operators_Operations Operations
You can use the usual operations and logical operators.
<table class="example">
|Feel++ Keyword | Math Object | Description | 
|```cpp{.cpp} + ``` |$$ f+g$$|tensor sum|
|```cpp{.cpp} - ``` |$$ f-g$$|tensor substraction|
|```cpp{.cpp} * ``` |$$ f*g$$|tensor product|
|```cpp{.cpp} / ``` |$$ f/g$$|tensor tensor division <br>($$g$$ scalar field)|
|```cpp{.cpp} < ``` |$$ f<g$$|element wise less|
|```cpp{.cpp} <= ``` |$$ f<=g$$|element wise less or equal|
|```cpp{.cpp} > ``` |$$ f>g$$|element wise greater|
|```cpp{.cpp} >= ``` |$$ f>=g$$|element wise greater or equal|
|```cpp{.cpp} == ``` |$$ f==g$$|element wise equal|
|```cpp{.cpp} != ``` |$$ f!=g$$|element wise not equal|
|```cpp{.cpp} - ``` |$$ -g$$|element wise unary minus|
|```cpp{.cpp} && ``` |$$ f$$ and $$g$$|element wise logical and |
|```cpp{.cpp} || ``` |$$ f$$ or $$g$$|element wise logical or|
|```cpp{.cpp} ! ``` |$$ !g$$|element wise logical not|


## Operators_Differential Differential Operators
Feel++ finit element language use <em>test</em> and <em>trial</em> functions. Keywords are different according to the kind of the manipulated function.<br>
<strong>Usual operators</strong> are for <strong>test</strong> functions.<br>
<strong>t-operators</strong> for <strong>trial</strong> functions.<br>
<strong>v-operators</strong> to get an <strong>evaluation</strong>.
<table class="example">
|Feel++ Keyword | Math Object | Description | Rank | Dimension | 
|```cpp{.cpp} id(f)``` | $$f$$ | test function | $$\mathrm{rank}(f(\overrightarrow{x}))$$ | $$m \times p $$|
|```cpp{.cpp} idt(f)```| $$f$$ | trial function | $$\mathrm{rank}(f(\overrightarrow{x}))$$ | $$m \times p $$|
|```cpp{.cpp} idv(f)```| $$f$$ | evaluation function   | $$\mathrm{rank}(f(\overrightarrow{x}))$$ | $$m \times p $$|
|```cpp{.cpp} grad(f)``` | $$\nabla f$$ | gradient of test function | $$\mathrm{rank}(f(\overrightarrow{x}))+1$$ | $$m \times n $$ <br> $$p=1$$|
|```cpp{.cpp} gradt(f)```| $$\nabla f$$ | grdient of trial function | $$\mathrm{rank}(f(\overrightarrow{x}))+1$$ |$$m \times n $$<br> $$p=1$$|
|```cpp{.cpp} gradv(f)```| $$\nabla f$$ | evaluation function gradient  | $$\mathrm{rank}(f(\overrightarrow{x}))+1$$ |$$m \times n $$<br> $$p=1$$|
|```cpp{.cpp} div(f)``` | $$\nabla\cdot f$$ | divergence of test function | $$\mathrm{rank}(f(\overrightarrow{x}))-1$$ | $$1 \times 1 $$|
|```cpp{.cpp} divt(f)```| $$\nabla\cdot f$$ | divergence of trial function | $$\mathrm{rank}(f(\overrightarrow{x}))-1$$ |$$1 \times 1 $$|
|```cpp{.cpp} divv(f)```| $$\nabla\cdot f$$ | evaluation function divergence  | $$\mathrm{rank}(f(\overrightarrow{x}))-1$$ |$$1 \times 1 $$|
|```cpp{.cpp} curl(f)``` | $$\nabla\times f$$ | curl of test function |1| $$n \times 1 $$<br>$$m=n$$|
|```cpp{.cpp} curlt(f)```| $$\nabla\times f$$ | curl of trial function |1 |$$n \times 1 $$<br>$$m=n$$|
|```cpp{.cpp} curlv(f)```| $$\nabla\times f$$ | evaluation function curl  |1 |$$n \times 1 $$<br>$$m=n$$|
|```cpp{.cpp} hess(f)```| $$\nabla^2 f$$ | hessian of test function  |2 |$$n \times n $$<br>$$m=p=1$$|


## Operators_TwoValued Two Valued Operators
<table class="example">
|Feel++ Keyword | Math Object | Description | Rank | Dimension | 
|```cpp{.cpp} jump(f)``` |  $$[f]=f_0\overrightarrow{N_0}+f_1\overrightarrow{N_1}$$ | jump of test function |0| $$n \times 1 $$<br>$$m=1$$|
|```cpp{.cpp} jump(f)``` |  $$[\overrightarrow{f}]=\overrightarrow{f_0}\cdot\overrightarrow{N_0}+\overrightarrow{f_1}\cdot\overrightarrow{N_1}$$ | jump of test function |0| $$1 \times 1 $$<br>$$m=2$$|
|```cpp{.cpp} jumpt(f)``` |  $$[f]=f_0\overrightarrow{N_0}+f_1\overrightarrow{N_1}$$ | jump of trial function |0| $$n \times 1 $$<br>$$m=1$$|
|```cpp{.cpp} jumpt(f)``` |  $$[\overrightarrow{f}]=\overrightarrow{f_0}\cdot\overrightarrow{N_0}+\overrightarrow{f_1}\cdot\overrightarrow{N_1}$$ | jump of trial function |0| $$1 \times 1 $$<br>$$m=2$$|
|```cpp{.cpp} jumpv(f)``` |  $$[f]=f_0\overrightarrow{N_0}+f_1\overrightarrow{N_1}$$ | jump of function evaluation |0| $$n \times 1 $$<br>$$m=1$$|
|```cpp{.cpp} jumpv(f)``` |  $$[\overrightarrow{f}]=\overrightarrow{f_0}\cdot\overrightarrow{N_0}+\overrightarrow{f_1}\cdot\overrightarrow{N_1}$$ | jump of function evaluation|0| $$1 \times 1 $$<br>$$m=2$$|
|```cpp{.cpp} average(f)``` |  $${f}=\frac{1}{2}(f_0+f_1)$$ | average of test function|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times n $$<br>$$m=n$$|
|```cpp{.cpp} averaget(f)``` |  $${f}=\frac{1}{2}(f_0+f_1)$$ | average of trial function|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times n $$<br>$$m=n$$|
|```cpp{.cpp} averagev(f)``` |  $${f}=\frac{1}{2}(f_0+f_1)$$ | average of function evaluation|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times n $$<br>$$m=n$$|
|```cpp{.cpp} leftface(f)``` |  $$f_0$$ |left test function|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times n $$<br>$$m=n$$|
|```cpp{.cpp} leftfacet(f)``` |  $$f_0$$ |left trial function|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times n $$<br>$$m=n$$|
|```cpp{.cpp} leftfacev(f)``` |  $$f_0$$ |left function evaluation|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times n $$<br>$$m=n$$|
|```cpp{.cpp} rightface(f)``` |  $$f_1$$ |right test function|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times n $$<br>$$m=n$$|
|```cpp{.cpp} rightfacet(f)``` |  $$f_1$$ |right trial function|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times n $$<br>$$m=n$$|
|```cpp{.cpp} rightfacev(f)``` |  $$f_1$$ |right function evaluation|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times n $$<br>$$m=n$$|
|```cpp{.cpp} maxface(f)``` |  $$\max(f_0,f_1)$$ |maximum of right and left<br>test function|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times p $$|
|```cpp{.cpp} maxfacet(f)``` |  $$\max(f_0,f_1)$$ |maximum of right and left<br>trial function|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times p $$|
|```cpp{.cpp} maxfacev(f)``` |  $$\max(f_0,f_1)$$ |maximum of right and left<br>function evaluation|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times p $$|
|```cpp{.cpp} minface(f)``` |  $$\min(f_0,f_1)$$ |minimum of right and left<br>test function|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times p $$|
|```cpp{.cpp} minfacet(f)``` |  $$\min(f_0,f_1)$$ |minimum of right and left<br>trial function|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times p $$|
|```cpp{.cpp} minfacev(f)``` |  $$\min(f_0,f_1)$$ |minimum of right and left<br>function evaluation|$$\mathrm{rank}( f(\overrightarrow{x}))$$| $$n \times p $$|



<a href="#" class="top">top</a>
<hr>
<br>
# Kyewords_Geometric Geometric Transformations
## Matrix_Jacobian Jacobian Matrix
You can access to the jacobian matrix, $$J$$, of the geometric transformation, using the keyword:
```cpp J() \eco
There are some tools to manipulate this jacobian.
<table class="manual">
|Feel++ Keyword | Math Object | Description | 
|```cpp{.cpp} detJ()```|$$\det(J)$$|Determinant of jacobian matrix |
|```cpp{.cpp} invJT()```|$$(J^{-1})^T$$|Transposed inverse of jacobian matrix |


