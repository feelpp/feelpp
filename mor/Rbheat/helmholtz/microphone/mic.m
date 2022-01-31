# -*- octave -*-

inP(1) = 1.5;  # mu1 \in [0.8;1.5]
inP(2) = 10;  # mu2 \in [10;50]
inP(3) = 0.1;
inP(4) = 2;
inP(5) = 1e-2;
inP;

yfem=[];ycrb=[];
D=[];
N=20;
z=[];y=[];x=[];
resize(x,N);
resize(y,N);
resize(z,N,N);
for i=1:N
  inP(1)= .8+(i-1)*(1.5-0.8)/N;
  x(i)=inP(1);
  for j=1:N
    inP(2)= 10+(j-1)*(50-10)/N;
    y(j)=inP(2);
    ##D=[D inP(2)];
    ##inP;
    ##yfem= [yfem opusawpfem( inP )];
    z(i,j)= opusawcrb( inP )(1);
    e(i,j)= opusawcrb( inP )(2);
  end
end
[xx,yy]=meshgrid(x,y);
#plot(D,yfem(1:2:2*N),'-xb',D,ycrb(1:2:2*N),'-xr');
#plot(D,ycrb(1:2:2*N),'-xr');
figure(1)
surf(x,y,z)
figure(2)
surf(x,y,e)
#xlabel('\mu_2');
#ylabel('s');



