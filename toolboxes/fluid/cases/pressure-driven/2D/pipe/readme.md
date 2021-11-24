# Pressure driven 2D pipe test case

## Geometry

The domain is the rectangle $\Omega=[0;L]\times[0;H]\subset\mathbb R^2$. Default is $L=1,H=4$. Top, bottom, left and right boundaries are labeled accordingly : $\Gamma_{top}=[0;L]\times\{H\}$, $\Gamma_{bottom}=[0;L]\times\{H\}$, $\Gamma_{left}=\{0\}\times[0;H]$ and $\Gamma_{right}=\{L\}\times[0;H]$.

The parameter are:

* `pin` : inlet pressure on $\Gamma_{in}
* `pout` : outlet pressure on $\Gamma_{out}$
* `L` : pipe width
* `H` : pipe length
* `nu` : dynamic viscosity

## Equations

Line 6 of `pipe-2d.json`, switching between 

- `"equations":"Stokes"`
- `"equations":"Navier-Stokes"`

allows to solve for either problem.

## Analytical solution

The exact solution for Stokes' problem writes :

$$\mathbf u_{ex} = \left(\begin{matrix}0\\\frac{p_{in}-p_{out}}{2H\nu}(L-x)x\end{matrix}\right)$$
$$p_{ex} = p_{in} + \frac{p_{out}-p_{in}}Hy$$

An easy computation shows the convective term $(\mathbf u\cdot\nabla)\mathbf u$ is zero. Thus, the analytical solution is the same for the Navier-Stokes' problem.