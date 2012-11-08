n=100
lm=9
m=10^lm
\precision=lm+5
root3=sqrt(3)
root5=sqrt(5)
p1 = sum(0,j=0,2*n,divres(floor(root5*m*j),m)[2]*(-1)^j*x^j)
p2 = sum(0,j=0,n,divres(floor(root3*m*j),m)[2]*x^j)
#
p1*p2
p1*p2
divres(p1,p2)
divres(p1,p2)
gcd(p1,p2)
gcd(p1,p2)
