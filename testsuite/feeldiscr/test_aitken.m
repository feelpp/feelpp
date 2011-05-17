M=load('history.dat')
semilogy(M(:,1),M(:,2),';abs;',M(:,1),M(:,3),';rel;')
polyfit(M(:,1),log(M(:,2)),1)
