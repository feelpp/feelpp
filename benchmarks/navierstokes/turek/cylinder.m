function cylinder (Ti,Tf,dt, file,color)
  %%  Copyright (C) 2007 Université Joseph Fourier
  %% 
  %%  This file is part of feel_cylinder.
  %% 
  %%  feel_cylinder is free software; you can redistribute it and/or modify
  %%  it under the terms of the GNU General Public License as published by
  %%  the Free Software Foundation; either version 2, or (at your option)
  %%  any later version.
  %% 
  %%  feel_cylinder is distributed in the hope that it will be useful, but
  %%  WITHOUT ANY WARRANTY; without even the implied warranty of
  %%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  %%  General Public License for more details.
  %% 
  %%  You should have received a copy of the GNU General Public License
  %%  along with Octave; see the file COPYING.  If not, write to the Free
  %%  Software Foundation, 59 Temple Place - Suite 330, Boston, MA
  %%  02111-1307, USA.
  
  %% usage:  cylinder (dt)
  %%
  %%

  %%  Author: Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
  %%  Keywords: CD, CL, DP and Frequency calculation
  %%  Maintainer: Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>

  figure(1);
  
  if ( file == "" )
    file='quantities.dat';
  endif 
  if ( color == "" )
    color='r*-';
  endif 
  
  x=load(file);
  
  Nsteps=(Tf-Ti)/dt;
  NTi=Ti/dt;
  NT=Tf/dt;
  
  nrows=3;
  ncols=2;
  
  %%top_title("Cylinder Test Case Re=100, dt=0.0125");
  
  %% CD
  %%minCD=5.57*ones(size(x(NTi:NT,1),1),1);
  %%maxCD=5.59*ones(size(x(NTi:NT,1),1),1);
  minCD=3.22*ones(size(x(NTi:NT,1),1),1);
  maxCD=3.24*ones(size(x(NTi:NT,1),1),1);

  subplot(nrows,ncols,1);
  plot(x(NTi:NT,1),x(NTi:NT,2),[ color, ';CD(t);'],
       x(NTi:NT,1),minCD,'@11;min CD(t);',
       x(NTi:NT,1),maxCD,'@12;max CD(t);');
  disp(['CD plotted']);
  
  %% CL
  %minCL=0.0104*ones(size(x(NTi:NT,1),1),1);
  %maxCL=0.0110*ones(size(x(NTi:NT,1),1),1);
  minCL=0.99*ones(size(x(NTi:NT,1),1),1);
  maxCL=1.01*ones(size(x(NTi:NT,1),1),1);

  subplot(nrows,ncols,2);
  plot(x(NTi:NT,1),x(NTi:NT,3), [color, ';CL(t);'],
       x(NTi:NT,1),minCL,'@11;min CL(t);',
       x(NTi:NT,1),maxCL,'@12;max CL(t);');
  
  disp(['CL plotted']);
  
  %% D P
  subplot(nrows,ncols,3);
  plot(x(NTi:NT,1),x(NTi:NT,5), [color, ';DP(t);']);
  
  disp(['DP plotted']);
  
  %%
  F=fft(x(NTi:NT,3));
  N=length(F);
  P=F(1:N/2).*conj(F(1:N/2))/(N/2);
  freq=(1:N/2)/(N/2)*0.5/dt;
  subplot(nrows,ncols,4);
  plot(freq,P,[color, ';P=f(CL(t));']);
  
  disp(['Frequencies plotted']);
  
  %% strouhal
  [maxf,imaxf]=max(P)
  max_CL_frequency=freq(imaxf)
  DT=1/freq(imaxf)
  strouhal=0.1*freq(imaxf)/(2*1.5/3)
  
  %%P=polyfit(
  
