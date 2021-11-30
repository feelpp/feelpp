# Pressure driven 2D curved pipe test case

## Geometry

The domain in polar coordinates $(r,\theta)$ is $\Omega=[r_1;r_2]\times[\frac\pi2;\frac\pi6]\subset\mathbb R^2$. Default is $r_1=1.9,r_2=2.1$. The curved boundaries are labeled `wall` and straight boundaries are labeled as follows : $\Gamma_{in}=[r_1;r_2]\times\{\frac\pi2+\frac\pi6\}$, $\Gamma_{out}=[r_1;r_2]\times\{\frac\pi2\}$.

The parameter are:

* `pin` : inlet pressure on $\Gamma_{in}
* `pout` : outlet pressure on $\Gamma_{out}$
* `L` : pipe width
* `H` : pipe length
* `nu` : dynamic viscosity

## Equations

Line 6 of `curved-pipe-2d.json`, switching between 

- `"equations":"Stokes"`
- `"equations":"Navier-Stokes"`

allows to solve for either problem.

## Analytical solution

In polar coordinates, the exact solution for Stokes' problem writes

$$\mathbf u_{ex} = \left(\frac{p_{in}-p_{out}}{\frac\pi6}\left(\frac12r\text{ln}(r)+\frac Cr+Dr\right)\right)\mathbf e_\theta$$
$$p_{ex} = \frac{p_{in}(\theta-\frac\pi6)-p_{out}\theta}{\frac\pi6}$$

where

$$C=\frac{r_1^2r_2^2}2\frac{\text{ln}(r_2)-\text{ln}(r_1)}{r_2^2-r_1^2}$$
$$D=-\frac12\frac{r_2^2\text{ln}(r_2)-r_1^2\text{ln}(r_1)}{r_2^2-r_1^2}$$

This yields, for default values of $r_1,r_2$, the following solution in cartesian coordinates.

$$\mathbf u_{ex} = \left(\begin{matrix} -\sin(\text{acos}(\frac x{\sqrt{x^2+y^2}}))(19.0985931710274-1.90985931710274)(0.5\sqrt{x^2+y^2}log(\sqrt{x^2+y^2})+\frac{0.995836667858137}{\sqrt{x^2+y^2}}-0.596781975733881\sqrt{x^2+y^2}) \\ \cos(\text{acos}(\frac x{\sqrt{x^2+y^2}}))(19.0985931710274-1.90985931710274)(0.5\sqrt{x^2+y^2}log(\sqrt{x^2+y^2})+\frac{0.995836667858137}{\sqrt{x^2+y^2}}-0.596781975733881\sqrt{x^2+y^2}) \end{matrix}\right)$$

$$p_{ex} = (19.0985931710274(\text{acos}(\frac x{\sqrt{x^2+y^2}})-1.5707963267949)+1.90985931710274(2.09439510239320-\text{acos}(\frac x{\sqrt{x^2+y^2}})))$$


A computation with sympy shows the convective term $(\mathbf u\cdot\nabla)\mathbf u=(C_x,C_y,C_z)^T$ is as follows.

```
Cx = 

8.59436692696235*x*(-8.59436692696235*x**2*y*(0.5*sqrt(x**2 + y**2)*log(sqrt(x**2 + y**2)) - 0.596781975733881*sqrt(x**2 + y**2) + 0.995836667858137/sqrt(x**2 + y**2))/((x**2 + y**2)**2*sqrt(-x**2/(x**2 + y**2) + 1)) - 8.59436692696235*sqrt(-x**2/(x**2 + y**2) + 1)*(0.5*y*log(sqrt(x**2 + y**2))/sqrt(x**2 + y**2) - 0.096781975733881*y/sqrt(x**2 + y**2) - 0.995836667858137*y/(x**2 + y**2)**(3/2)))*(0.5*sqrt(x**2 + y**2)*log(sqrt(x**2 + y**2)) - 0.596781975733881*sqrt(x**2 + y**2) + 0.995836667858137/sqrt(x**2 + y**2))/sqrt(x**2 + y**2) - 8.59436692696235*sqrt(-x**2/(x**2 + y**2) + 1)*(-8.59436692696235*sqrt(-x**2/(x**2 + y**2) + 1)*(0.5*x*log(sqrt(x**2 + y**2))/sqrt(x**2 + y**2) - 0.096781975733881*x/sqrt(x**2 + y**2) - 0.995836667858137*x/(x**2 + y**2)**(3/2)) - 8.59436692696235*(x**3/(x**2 + y**2)**2 - x/(x**2 + y**2))*(0.5*sqrt(x**2 + y**2)*log(sqrt(x**2 + y**2)) - 0.596781975733881*sqrt(x**2 + y**2) + 0.995836667858137/sqrt(x**2 + y**2))/sqrt(-x**2/(x**2 + y**2) + 1))*(0.5*sqrt(x**2 + y**2)*log(sqrt(x**2 + y**2)) - 0.596781975733881*sqrt(x**2 + y**2) + 0.995836667858137/sqrt(x**2 + y**2)) 
```

```
Cy = 

 8.59436692696235*x*(-8.59436692696235*x*y*(0.5*sqrt(x**2 + y**2)*log(sqrt(x**2 + y**2)) - 0.596781975733881*sqrt(x**2 + y**2) + 0.995836667858137/sqrt(x**2 + y**2))/(x**2 + y**2)**(3/2) + 8.59436692696235*x*(0.5*y*log(sqrt(x**2 + y**2))/sqrt(x**2 + y**2) - 0.096781975733881*y/sqrt(x**2 + y**2) - 0.995836667858137*y/(x**2 + y**2)**(3/2))/sqrt(x**2 + y**2))*(0.5*sqrt(x**2 + y**2)*log(sqrt(x**2 + y**2)) - 0.596781975733881*sqrt(x**2 + y**2) + 0.995836667858137/sqrt(x**2 + y**2))/sqrt(x**2 + y**2) - 8.59436692696235*sqrt(-x**2/(x**2 + y**2) + 1)*(0.5*sqrt(x**2 + y**2)*log(sqrt(x**2 + y**2)) - 0.596781975733881*sqrt(x**2 + y**2) + 0.995836667858137/sqrt(x**2 + y**2))*(-8.59436692696235*x**2*(0.5*sqrt(x**2 + y**2)*log(sqrt(x**2 + y**2)) - 0.596781975733881*sqrt(x**2 + y**2) + 0.995836667858137/sqrt(x**2 + y**2))/(x**2 + y**2)**(3/2) + 8.59436692696235*x*(0.5*x*log(sqrt(x**2 + y**2))/sqrt(x**2 + y**2) - 0.096781975733881*x/sqrt(x**2 + y**2) - 0.995836667858137*x/(x**2 + y**2)**(3/2))/sqrt(x**2 + y**2) + 8.59436692696235*(0.5*sqrt(x**2 + y**2)*log(sqrt(x**2 + y**2)) - 0.596781975733881*sqrt(x**2 + y**2) + 0.995836667858137/sqrt(x**2 + y**2))/sqrt(x**2 + y**2)) 
```