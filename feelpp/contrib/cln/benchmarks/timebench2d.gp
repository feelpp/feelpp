n=500
\precision=floor(n*log(2)/log(10)+5)
p=nextprime(ceil(pi*2^500))
#
factor(mod(x^n+x+1,p))
