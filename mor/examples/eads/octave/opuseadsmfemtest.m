# -*- octave -*- 

  inP(1) = 1.0e+1;  # kIC : thermal conductivity (default: 2)
  inP(2) = 7.0e-3;  # D : fluid flow rate (default: 5e-3)
  inP(3) = 1.0e+6;  # Q : heat flux (default: 1e6)
  inP(4) = 1.0e+2;  # r : conductance (default: 100)
  inP(5) = 4.0e-3;  # ea : length air flow channel (default: 4e-3)
  inP(6) = 1;       # meshsize times 1e-3 (default: 1)
  inP(7) = 2;       # integer : polynomial degree for the temperature in {1,2,3,4}

inP;

y=[];
D=[];
N=20;
for i=1:N 
	inP(1)= 0.2+(i-1)*(150-0.2)/N;
      #  exp(log(0.2)+i*log(150)/4);
      #  exp(log(2e-3)+i*(log(1e-2)-log(2e-3))/(N-1));
      #  D=linspace(1e-4,1e-2,10)
      #  inP(2)=D;
D=[D inP(1)];
inP;
y= [y opuseadsmfem( inP )];

end

#D=linspace(1e-4,1e-2,10);


figure(1)
axis([0 150 300 350]);
plot(D,y(1:2:40),'--og', D,y(2:2:40),'--rs');
xlabel('K_{ic}');
ylabel('T(K)');  
saveas(figure(1),'Temperature-Conduc');

