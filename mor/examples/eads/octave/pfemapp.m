# -*- octave -*-

inP(1) = 10   # kIC : thermal conductivity (default: 2)
inP(2) = 7e-3 # D : fluid flow rate (default: 5e-3)
inP(3) = 1e6  # Q : heat flux (default: 1e6)
inP(4) = 100  # r : conductance (default: 100)
inP(5) = 4e-3 # ea : length air flow channel (default: 4e-3)
inP(6) = 1 # meshsize times 1e-3 (default: 1)
inP(7) = 2 # integer : polynomial degree for the temperature in {1,2,3,4}

inP;

for D=linspace(1e-4,1e-2,10)
  inP(2)=D;
  y= [y pfemapp( inP )];

end

D=linspace(1e-4,1e-2,10);
plot(D,y(1:2:20));
y= pfemapp( inP )



