# -*- octave -*-


inP(1) = 380;
inP(2) = 380;
inP(3) = 0.001;

inP;

yfem=[]; ycrb=[];
D=[];
N=10;
z=[];y=[];x=[];
resize(x,N);
resize(y,N);
resize(z,N,N);
for i=1:N 
  inP(1)= 0.2+(i-1)*(50-0.2)/N;
  x(i)=inP(1);
  for j=1:N
    intP(2)=0.8+(i-1)*(50-0.8)/N;
    y(j) = inP(2);
    z(i,j) = opusheat1dcrb( inP )(1);
    e(i,j) = opusheat1dcrb( inP )(2);
  end
end
[xx,yy]=meshgrid(x,y);
figure(1)
surf(x,y,z)
figure(2)
surf(x,y,e)


