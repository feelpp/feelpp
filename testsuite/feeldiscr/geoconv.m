M=load('conv.dat');
%%ex=1+1./(2*pi);
%ex=exp(1)-1;
%%ex=1+pi*0.5*0.5/2;
%%ex=(%e-1)*(8*%pi^2-%pi+2)/(2*(4*%pi^2+1))
%ex=(e-1)*(8*pi^2-pi+2)/(2*(4*pi^2+1))
% ex=1+1/2-1/3;
%ex=31/(2^2*3);
ex=0;
loglog(M(:,1),abs(M(:,2)-ex),';P_1;',
       M(:,1), abs(M(:,3)-ex),';P_2;',
       M(:,1),abs(M(:,4)-ex),';P_3;',
       M(:,1),abs(M(:,5)-ex),';P_4;',
       M(:,1),abs(M(:,6)-ex),';P_5;');

polyfit(log(M(:,1)),log(abs(M(:,3)-ex)),1)
polyfit(log(M(:,1)),log(abs(M(:,4)-ex)),1)
polyfit(log(M(:,1)),log(abs(M(:,5)-ex)),1)
polyfit(log(M(:,1)),log(abs(M(:,6)-ex)),1)