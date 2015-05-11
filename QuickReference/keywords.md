Keywords
========

One of \feel assets is it finite element embedded language. The language follows the C++ grammar, and provides keywords as well as operations between objects which are, mathematically, tensors of rank 0, 1 or 2.

In all following tables we use the notations:<br>
$$f: \mathbb{R}^n \mapsto \mathbb{R}^{m\times p}$$  with $$n=1,2,3$$, $$m=1,2,3$$, $$p=1,2,3$$<br>
$$\Omega^e$$ current mesh element.

The genesis of the language can be found in \cite prudhomme05:_domain_specif_embed_languag_c and an update on Feel++ is available in \cite PRUDHOMME:2012:HAL-00662868:3.

# Keywords_Points Points

Current Point:
<table class="manual">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th><th>Dimension</th></tr>
<tr><td>```cpp{.cpp} P()```</td><td>$$\overrightarrow{P}$$</td><td> $$(P_x, P_y, P_z)^T$$</td><td>$$d \times 1$$</td></tr>
<tr><td>```cpp{.cpp} Px()```</td><td>$$P_x$$</td><td>$$x$$ coordinate of $$\overrightarrow{P}$$</td><td>$$1 \times 1$$</td></tr>
<tr><td>```cpp{.cpp} Py()```</td><td>$$P_y$$</td><td>$$y$$ coordinate of $$\overrightarrow{P}$$<br>(value is 0 in 1D)</td><td>$$1 \times 1$$</td></tr>
<tr><td>```cpp{.cpp} Pz()```</td><td>$$P_z$$</td><td>$$z$$ coordinate of $$\overrightarrow{P}$$<br>(value is 0 in 1D and 2D)</td><td>$$1 \times 1$$</td></tr>
</table>

Element Barycenter Point:
<table class="manual">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th><th>Dimension</th></tr>
<tr><td>```cpp{.cpp} C()```</td><td>$$\overrightarrow{C}$$</td><td> $$(C_x, C_y, C_z)^T$$</td><td>$$d \times 1$$</td></tr>
<tr><td>```cpp{.cpp} Cx()```</td><td>$$C_x$$</td><td>$$x$$ coordinate of $$\overrightarrow{C}$$</td><td>$$1 \times 1$$</td></tr>
<tr><td>```cpp{.cpp} Cy()```</td><td>$$C_y$$</td><td>$$y$$ coordinate of $$\overrightarrow{C}$$<br>(value is 0 in 1D)</td><td>$$1 \times 1$$</td></tr>
<tr><td>```cpp{.cpp} Cz()```</td><td>$$C_z$$</td><td>$$z$$ coordinate of $$\overrightarrow{C}$$<br>(value is 0 in 1D and 2D)</td><td>$$1 \times 1$$</td></tr>
</table>

Normal at Current Point:
<table class="manual">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th><th>Dimension</th></tr>
<tr><td>```cpp{.cpp} N()```</td><td>$$\overrightarrow{N}$$</td><td> $$(N_x, N_y, N_z)^T$$</td><td>$$d \times 1$$</td></tr>
<tr><td>```cpp{.cpp} Nx()```</td><td>$$N_x$$</td><td>$$x$$ coordinate of $$\overrightarrow{N}$$</td><td>$$1 \times 1$$</td></tr>
<tr><td>```cpp{.cpp} Ny()```</td><td>$$N_y$$</td><td>$$y$$ coordinate of $$\overrightarrow{N}$$<br>(value is 0 in 1D)</td><td>$$1 \times 1$$</td></tr>
<tr><td>```cpp{.cpp} Nz()```</td><td>$$N_z$$</td><td>$$z$$ coordinate of $$\overrightarrow{N}$$<br>(value is 0 in 1D and 2D)</td><td>$$1 \times 1$$</td></tr>
</table>


<a href="#" class="top">top</a>
<hr>
<br>
# Keywords_Array Vectors and Matrix
## BuildingVectors Building Vectors
Usual syntax to create vectors:
<table class="manual">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th><th>Dimension</th></tr>
<tr><td>```cpp{.cpp} vec<n>(v_1,v_2,...,v_n)```</td><td>$$\begin{pmatrix} v_1\\v_2\\ \vdots \\v_n \end{pmatrix}$$</td><td>Column Vector with $$n$$ rows<br>entries being expressions</td><td>$$n \times 1$$</td></tr>
</table>
You can also use expressions and the unit base vectors:
<table class="example">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th></tr>
<tr><td>\c oneX() </td><td> $$\begin{pmatrix} 1\\0\\0 \end{pmatrix}$$</td><td>Unit vector $$\overrightarrow{i}$$</td></tr>
<tr><td>\c oneY() </td><td> $$\begin{pmatrix} 0\\1\\0 \end{pmatrix}$$</td><td>Unit vector $$\overrightarrow{j}$$</td></tr>
<tr><td>\c oneZ() </td><td> $$\begin{pmatrix} 0\\0\\1 \end{pmatrix}$$</td><td>Unit vector $$\overrightarrow{k}$$</td></tr>
</table>


## BuildingMatrix Building Matrix
<table class="manual">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th><th>Dimension</th></tr>
<tr><td>```cpp{.cpp} mat<m,n>(m_11,m_12,...,m_mn)```</td><td>$$\begin{pmatrix} m_{11} & m_{12} & ...\\ m_{21} & m_{22} & ...\\ \vdots & & \end{pmatrix}$$</td><td>$$m\times n$$ Matrix<br> entries beeing expressions </td><td>$$m \times n$$</td></tr>
<tr><td>```cpp{.cpp} ones<m,n>()```</td><td>$$\begin{pmatrix} 1 & 1 & ...\\ 1 & 1 & ...\\ \vdots & & \end{pmatrix}$$</td><td>$$m\times n$$ Matrix <br>Filled with 1 </td><td>$$m \times n$$</td></tr>
<tr><td>```cpp{.cpp} zero<m,n>()```</td><td>$$\begin{pmatrix} 0 & 0 & ...\\ 0 & 0 & ...\\ \vdots & & \end{pmatrix}$$</td><td>$$m\times n$$ Matrix <br>Filled with 0 </td><td>$$m \times n$$</td></tr>
<tr><td>```cpp{.cpp} constant<m,n>(c)```</td><td>$$\begin{pmatrix} c & c & ...\\ c & c & ...\\ \vdots & & \end{pmatrix}$$</td><td>$$m\times n$$ Matrix <br>Filled with a constant c </td><td>$$m \times n$$</td></tr>
<tr><td>```cpp{.cpp} eye<n>()```</td><td>$$\begin{pmatrix} 1 & 0 & ...\\ 0 & 1 & ...\\ \vdots & & \end{pmatrix}$$</td><td>Unit diagonal Matrix <br> of size$$n\times n$$ </td><td>$$n \times n$$</td></tr>
<tr><td>```cpp{.cpp} Id<n>()```</td><td>$$\begin{pmatrix} 1 & 0 & ...\\ 0 & 1 & ...\\ \vdots & & \end{pmatrix}$$</td><td>Unit diagonal Matrix <br> of size$$n\times n$$ </td><td>$$n \times n$$</td></tr>
</table>

## Matrix_Manipulating Manipulating Vectors and Matrix
Let $$A$$ be a square matrix of size $$n$$.
<table class="manual">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th><th>Dimension</th></tr>
<tr><td>```cpp{.cpp} inv(A)```</td><td>$$A^{-1}$$</td><td>Inverse of matrix $$A$$ </td><td>$$n \times n$$</td></tr>
<tr><td>```cpp{.cpp} det(A)```</td><td>$$\det (A)$$</td><td>Determinant of matrix $$A$$ </td><td>$$1 \times 1$$</td></tr>
<tr><td>\co sym(A)\eco</td><td>$$\text{Sym}(A)$$</td><td>Symmetric part of matrix $$A$$: $$\frac{1}{2}(A+A^T)$$<br> </td><td>$$n \times n$$</td></tr>
<tr><td>\co antisym(A)\eco</td><td>$$ \text{Asym}(A)$$</td><td>Antisymmetric part of  $$A$$: $$\frac{1}{2}(A-A^T)$$<br> </td><td>$$n \times n$$</td></tr>
</table>

Let A and B be two matrix (or two vectors) of same dimension $$m \times n$$.
<table class="manual">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th><th>Dimension</th></tr>
<tr><td>\co trace(A)\eco</td><td>$$\text{tr}(A)$$</td><td>Trace of matrix $$A$$<br>Generalized on non-squared Matrix<br>Generalized on Vectors </td><td>$$1 \times 1$$</td></tr>
<tr><td>\co trans(B)\eco</td><td>$$B^T$$</td><td>Transpose of matrix $$B$$<br>Can be used on non-squared Matrix<br>Can be used on Vectors </td><td>$$n \times m$$</td></tr>
<tr><td>\co inner(A,B)\eco</td><td>$$ A.B \\ A:B = \text{tr}(A*B^T)$$</td><td>Scalar product of two vectors<br>Generalized scalar product of two matrix </td><td>$$1 \times 1$$</td></tr>
<tr><td>\co cross(A,B)\eco</td><td>$$ A\times B$$</td><td>Cross product of two vectors</td><td>$$n \times 1$$</td></tr>
</table>




<a href="#" class="top">top</a>
<hr>
<br>
# Keywords_Expr Expressions
Following tables present tools to declare and manipulate expressions.
<table class="example">
<tr><th>\feel Keyword</th><th>Description</th></tr>
<tr><td>```cpp{.cpp}
Px()
Py()
Pz()
cst( c )
```</td><td>
Variable $$x$$<br>
Variable $$y$$<br>
Variable $$z$$<br>
Constant function equal to $$c$$
</td></tr>
</table>

You can of course use all current operators ( + - / * ) and the usual following functions:
<table class="example">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th>
<tr><td>```cpp{.cpp} abs(expr) ```</td><td>$$|f(\overrightarrow{x})|$$</td><td>element wise absolute value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} cos(expr)```</td><td>$$\cos(f(\overrightarrow{x}))$$</td><td>element wise cos value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} sin(expr)```</td><td>$$\sin(f(\overrightarrow{x}))$$</td><td>element wise sin value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} tan(expr)```</td><td>$$\tan(f(\overrightarrow{x}))$$</td><td>element wise tan value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} acos(expr)```</td><td>$$\acos(f(\overrightarrow{x}))$$</td><td>element wise acos value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} asin(expr)```</td><td>$$\asin(f(\overrightarrow{x}))$$</td><td>element wise asin value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} atan(expr)```</td><td>$$\atan(f(\overrightarrow{x}))$$</td><td>element wise atan value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} cosh(expr)```</td><td>$$\cosh(f(\overrightarrow{x}))$$</td><td>element wise cosh value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} sinh(expr)```</td><td>$$\sinh(f(\overrightarrow{x}))$$</td><td>element wise sinh value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} tanh(expr)```</td><td>$$\tanh(f(\overrightarrow{x}))$$</td><td>element wise tanh value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} exp(expr)```</td><td>$$\exp(f(\overrightarrow{x}))$$</td><td>element wise exp value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} log(expr)```</td><td>$$\log(f(\overrightarrow{x}))$$</td><td>element wise log value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} sqrt(expr)```</td><td>$$\sqrt{f(\overrightarrow{x})}$$</td><td>element wise sqrt value of $$f$$</td></tr>
<tr><td>```cpp{.cpp} ceil(expr)```</td><td>$$\lceil{f(\overrightarrow{x})}\rceil$$</td><td>element wise ceil of $$f$$</td></tr>
<tr><td>```cpp{.cpp} floor(expr)```</td><td>$$\lfloor{f(\overrightarrow{x})}\rfloor$$</td><td>element wise floor of $$f$$</td></tr>
<tr><td>```cpp{.cpp} sign(expr)```</td><td>$$\begin{cases} 1 & \text{if}\ f(\overrightarrow{x}) \geq 0\\-1 & \text{if}\ f(\overrightarrow{x}) < 0\end{cases}$$</td><td>element wise sign value of $$f$$</td></tr>
<tr><td>```cpp chi(expr)```</td><td>$$\chi(f(\overrightarrow{x}))=\begin{cases}0 & \text{if}\ f(\overrightarrow{x}) = 0\\1 & \text{if}\ f(\overrightarrow{x}) \neq 0\\\end{cases}$$</td><td>element wise boolean test of $$f$$</td></tr>
</table>

<a href="#" class="top">top</a>
<hr>
<br>
# Keywords_Operators Operators
## Operators_Operations Operations
You can use the usual operations and logical operators.
<table class="example">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th></tr>
<tr><td>```cpp{.cpp} + ``` </td><td>$$ f+g$$</td><td>tensor sum</td></tr>
<tr><td>```cpp{.cpp} - ``` </td><td>$$ f-g$$</td><td>tensor substraction</td></tr>
<tr><td>```cpp{.cpp} * ``` </td><td>$$ f*g$$</td><td>tensor product</td></tr>
<tr><td>```cpp{.cpp} / ``` </td><td>$$ f/g$$</td><td>tensor tensor division <br>($$g$$ scalar field)</td></tr>
<tr><td>```cpp{.cpp} < ``` </td><td>$$ f<g$$</td><td>element wise less</td></tr>
<tr><td>```cpp{.cpp} <= ``` </td><td>$$ f<=g$$</td><td>element wise less or equal</td></tr>
<tr><td>```cpp{.cpp} > ``` </td><td>$$ f>g$$</td><td>element wise greater</td></tr>
<tr><td>```cpp{.cpp} >= ``` </td><td>$$ f>=g$$</td><td>element wise greater or equal</td></tr>
<tr><td>```cpp{.cpp} == ``` </td><td>$$ f==g$$</td><td>element wise equal</td></tr>
<tr><td>```cpp{.cpp} != ``` </td><td>$$ f!=g$$</td><td>element wise not equal</td></tr>
<tr><td>```cpp{.cpp} - ``` </td><td>$$ -g$$</td><td>element wise unary minus</td></tr>
<tr><td>```cpp{.cpp} && ``` </td><td>$$ f$$ and $$g$$</td><td>element wise logical and </td></tr>
<tr><td>```cpp{.cpp} || ``` </td><td>$$ f$$ or $$g$$</td><td>element wise logical or</td></tr>
<tr><td>```cpp{.cpp} ! ``` </td><td>$$ !g$$</td><td>element wise logical not</td></tr>
</table>

## Operators_Differential Differential Operators
\feel finit element language use <em>test</em> and <em>trial</em> functions. Keywords are different according to the kind of the manipulated function.<br>
<strong>Usual operators</strong> are for <strong>test</strong> functions.<br>
<strong>t-operators</strong> for <strong>trial</strong> functions.<br>
<strong>v-operators</strong> to get an <strong>evaluation</strong>.
<table class="example">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th><th>Rank</th><th>Dimension</th></tr>
<tr><td>```cpp{.cpp} id(f)``` </td><td> $$f$$ </td><td> test function </td><td> $$\mathrm{rank}(f(\overrightarrow{x}))$$ </td><td> $$m \times p $$</td></tr>
<tr><td>```cpp{.cpp} idt(f)```</td><td> $$f$$ </td><td> trial function </td><td> $$\mathrm{rank}(f(\overrightarrow{x}))$$ </td><td> $$m \times p $$</td></tr>
<tr><td>```cpp{.cpp} idv(f)```</td><td> $$f$$ </td><td> evaluation function   </td><td> $$\mathrm{rank}(f(\overrightarrow{x}))$$ </td><td> $$m \times p $$</td></tr>
<tr><td>```cpp{.cpp} grad(f)``` </td><td> $$\nabla f$$ </td><td> gradient of test function </td><td> $$\mathrm{rank}(f(\overrightarrow{x}))+1$$ </td><td> $$m \times n $$ <br> $$p=1$$</td></tr>
<tr><td>```cpp{.cpp} gradt(f)```</td><td> $$\nabla f$$ </td><td> grdient of trial function </td><td> $$\mathrm{rank}(f(\overrightarrow{x}))+1$$ </td><td>$$m \times n $$<br> $$p=1$$</td></tr>
<tr><td>```cpp{.cpp} gradv(f)```</td><td> $$\nabla f$$ </td><td> evaluation function gradient  </td><td> $$\mathrm{rank}(f(\overrightarrow{x}))+1$$ </td><td>$$m \times n $$<br> $$p=1$$</td></tr>
<tr><td>```cpp{.cpp} div(f)``` </td><td> $$\nabla\cdot f$$ </td><td> divergence of test function </td><td> $$\mathrm{rank}(f(\overrightarrow{x}))-1$$ </td><td> $$1 \times 1 $$</td></tr>
<tr><td>```cpp{.cpp} divt(f)```</td><td> $$\nabla\cdot f$$ </td><td> divergence of trial function </td><td> $$\mathrm{rank}(f(\overrightarrow{x}))-1$$ </td><td>$$1 \times 1 $$</td></tr>
<tr><td>```cpp{.cpp} divv(f)```</td><td> $$\nabla\cdot f$$ </td><td> evaluation function divergence  </td><td> $$\mathrm{rank}(f(\overrightarrow{x}))-1$$ </td><td>$$1 \times 1 $$</td></tr>
<tr><td>```cpp{.cpp} curl(f)``` </td><td> $$\nabla\times f$$ </td><td> curl of test function </td><td>1</td><td> $$n \times 1 $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} curlt(f)```</td><td> $$\nabla\times f$$ </td><td> curl of trial function </td><td>1 </td><td>$$n \times 1 $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} curlv(f)```</td><td> $$\nabla\times f$$ </td><td> evaluation function curl  </td><td>1 </td><td>$$n \times 1 $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} hess(f)```</td><td> $$\nabla^2 f$$ </td><td> hessian of test function  </td><td>2 </td><td>$$n \times n $$<br>$$m=p=1$$</td></tr>
</table>

## Operators_TwoValued Two Valued Operators
<table class="example">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th><th>Rank</th><th>Dimension</th></tr>
<tr><td>```cpp{.cpp} jump(f)``` </td><td>  $$[f]=f_0\overrightarrow{N_0}+f_1\overrightarrow{N_1}$$ </td><td> jump of test function </td><td>0</td><td> $$n \times 1 $$<br>$$m=1$$</td></tr>
<tr><td>```cpp{.cpp} jump(f)``` </td><td>  $$[\overrightarrow{f}]=\overrightarrow{f_0}\cdot\overrightarrow{N_0}+\overrightarrow{f_1}\cdot\overrightarrow{N_1}$$ </td><td> jump of test function </td><td>0</td><td> $$1 \times 1 $$<br>$$m=2$$</td></tr>
<tr><td>```cpp{.cpp} jumpt(f)``` </td><td>  $$[f]=f_0\overrightarrow{N_0}+f_1\overrightarrow{N_1}$$ </td><td> jump of trial function </td><td>0</td><td> $$n \times 1 $$<br>$$m=1$$</td></tr>
<tr><td>```cpp{.cpp} jumpt(f)``` </td><td>  $$[\overrightarrow{f}]=\overrightarrow{f_0}\cdot\overrightarrow{N_0}+\overrightarrow{f_1}\cdot\overrightarrow{N_1}$$ </td><td> jump of trial function </td><td>0</td><td> $$1 \times 1 $$<br>$$m=2$$</td></tr>
<tr><td>```cpp{.cpp} jumpv(f)``` </td><td>  $$[f]=f_0\overrightarrow{N_0}+f_1\overrightarrow{N_1}$$ </td><td> jump of function evaluation </td><td>0</td><td> $$n \times 1 $$<br>$$m=1$$</td></tr>
<tr><td>```cpp{.cpp} jumpv(f)``` </td><td>  $$[\overrightarrow{f}]=\overrightarrow{f_0}\cdot\overrightarrow{N_0}+\overrightarrow{f_1}\cdot\overrightarrow{N_1}$$ </td><td> jump of function evaluation</td><td>0</td><td> $$1 \times 1 $$<br>$$m=2$$</td></tr>
<tr><td>```cpp{.cpp} average(f)``` </td><td>  $${f}=\frac{1}{2}(f_0+f_1)$$ </td><td> average of test function</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times n $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} averaget(f)``` </td><td>  $${f}=\frac{1}{2}(f_0+f_1)$$ </td><td> average of trial function</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times n $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} averagev(f)``` </td><td>  $${f}=\frac{1}{2}(f_0+f_1)$$ </td><td> average of function evaluation</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times n $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} leftface(f)``` </td><td>  $$f_0$$ </td><td>left test function</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times n $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} leftfacet(f)``` </td><td>  $$f_0$$ </td><td>left trial function</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times n $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} leftfacev(f)``` </td><td>  $$f_0$$ </td><td>left function evaluation</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times n $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} rightface(f)``` </td><td>  $$f_1$$ </td><td>right test function</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times n $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} rightfacet(f)``` </td><td>  $$f_1$$ </td><td>right trial function</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times n $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} rightfacev(f)``` </td><td>  $$f_1$$ </td><td>right function evaluation</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times n $$<br>$$m=n$$</td></tr>
<tr><td>```cpp{.cpp} maxface(f)``` </td><td>  $$\max(f_0,f_1)$$ </td><td>maximum of right and left<br>test function</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times p $$</td></tr>
<tr><td>```cpp{.cpp} maxfacet(f)``` </td><td>  $$\max(f_0,f_1)$$ </td><td>maximum of right and left<br>trial function</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times p $$</td></tr>
<tr><td>```cpp{.cpp} maxfacev(f)``` </td><td>  $$\max(f_0,f_1)$$ </td><td>maximum of right and left<br>function evaluation</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times p $$</td></tr>
<tr><td>```cpp{.cpp} minface(f)``` </td><td>  $$\min(f_0,f_1)$$ </td><td>minimum of right and left<br>test function</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times p $$</td></tr>
<tr><td>```cpp{.cpp} minfacet(f)``` </td><td>  $$\min(f_0,f_1)$$ </td><td>minimum of right and left<br>trial function</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times p $$</td></tr>
<tr><td>```cpp{.cpp} minfacev(f)``` </td><td>  $$\min(f_0,f_1)$$ </td><td>minimum of right and left<br>function evaluation</td><td>$$\mathrm{rank}( f(\overrightarrow{x}))$$</td><td> $$n \times p $$</td></tr>
</table>


<a href="#" class="top">top</a>
<hr>
<br>
# Kyewords_Geometric Geometric Transformations
## Matrix_Jacobian Jacobian Matrix
You can access to the jacobian matrix, $$J$$, of the geometric transformation, using the keyword:
\co J() \eco
There are some tools to manipulate this jacobian.
<table class="manual">
<tr><th>\feel Keyword</th><th>Math Object</th><th>Description</th></tr>
<tr><td>```cpp{.cpp} detJ()```</td><td>$$\det(J)$$</td><td>Determinant of jacobian matrix </td></tr>
<tr><td>```cpp{.cpp} invJT()```</td><td>$$(J^{-1})^T$$</td><td>Transposed inverse of jacobian matrix </td></tr>
</table>

<a href="#" class="top">top</a>
<hr>

* \b Next: \ref Environment



*/
}
