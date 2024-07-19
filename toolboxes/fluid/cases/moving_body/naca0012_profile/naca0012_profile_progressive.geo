// Using Chebyshev nodes to refine better in the front and back of the wing
//+
L=1.009;
Ls =1/128;//1/128
//+
Point(1) = {-4, -2, 0, 80*Ls};
//+
Point(2) = {-4, 2, 0, 80*Ls};
//+
Point(3) = {16, 2, 0, 80*Ls};
//+
Point(4) = {16, -2, 0, 80*Ls};
Point(5) = {0, 0, 0, Ls};
Point(6) = {1.009, 0, 0, Ls};
nPoints=100;
Xi=0;
CM_x=0;
CM_y=0;
pList[0] = 6;
For i In {2 : nPoints}
  X = 0.5*1.009 + 0.5*1.009*Cos((2*i-1)/(2*nPoints)*Pi);//Xi + L*i/(nPoints + 1);
  pList[i-1]=i+10;
  Point(pList[i-1]) = {X,
	             0.6*(0.2969*Sqrt(X)-0.1260*X-0.3516*X^2+0.2843*X^3-0.1015*X^4),
	             0,
	             Ls};
  CM_x=CM_x+X;
EndFor	
pList[nPoints] = 5;
CM_x=CM_x+1.009;
For i In {1 : nPoints-1}
  X = L - (0.5*1.009 + 0.5*1.009*Cos((2*i-1)/(2*nPoints)*Pi));
  pList[nPoints+i]= nPoints+1+i+10;
  Point(pList[nPoints+i]) = {X,
	             -0.6*(0.2969*Sqrt(X)-0.1260*X-0.3516*X^2+0.2843*X^3-0.1015*X^4),
	             0,
	             Ls};
  CM_x=CM_x+X;
EndFor
pList[2*nPoints]=6;
Line(100) = pList[];
Point(10)={CM_x/(2*nPoints+1),0,0,Ls};
//+
Line(101) = {2, 1};
//+
Line(102) = {1, 4};
//+
Line(103) = {4, 3};
//+
Line(104) = {3, 2};
//+
Curve Loop(1) = {104, 101, 102, 103};
//+
Curve Loop(2) = {100};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(2) = {2};
//+
Physical Surface("AirfoilVolume") = {2};
//+
Physical Surface("Fluid") = {1};
//+
Physical Curve("Inlet") = {101};
//+
Physical Curve("Walls") = {104, 102};
//+
Physical Curve("Outlet") = {103};
//+
Physical Curve("AirfoilSurface") = {100};
