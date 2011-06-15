# -*- mode: octave -*-
x=[0.1   # hsize
   1     # shape (0:simplex, 1:hypercube, 2:ellipsoid)
   2     # function (1, 2)
   3     # parameter alpha
   10    # parameter beta
   0     # weakdir (0: weak dirichlet condition, 1 strong)
   50    # penaldir (value for weak dirichlet parameter)
   ]
L2=[];
H1=[];
eL2=[];
eH1=[];
#h=[0.1 0.05 0.025 0.0125 0.005]
h=[0.1 0.05 0.025]
for H=h
  x(1)=H;
  y=residualestimator_3_1(x);
  L2=[L2 y(1)];
  H1=[H1 y(2)];
  eL2=[eL2 y(3)];
  eH1=[eH1 y(4)];
end

loglog(h,L2,'-*;L_2;',h,h.*h,'-+;h^2;',
       h,H1,'-o;H_1;',h,h,'-x;h;',
       h,eL2,'-*;estimated L_2;',
       h,eH1,'-o;estimated H_1;');
title('-Delta u =f in 1D');

