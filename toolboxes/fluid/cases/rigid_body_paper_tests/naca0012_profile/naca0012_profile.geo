// Gmsh project created on Wed Oct  6 14:10:28 2021
//+
L=1.009;
Ls = 1/64;
//+
Point(1) = {-4, -2, 0, Ls};
//+
Point(2) = {-4, 2, 0, Ls};
//+
Point(3) = {16, 2, 0, Ls};
//+
Point(4) = {16, -2, 0, Ls};
Point(5) = {0, 0, 0, Ls};
Point(6) = {1.009, 0, 0, Ls};
nPoints=500;
Xi=0;
CM_x=0;
CM_y=0;
pList[0] = 5;
For i In {1 : nPoints}
  X = Xi + L*i/(nPoints + 1);
  pList[i]=i+10;
  Point(pList[i]) = {X,
	             0.6*(0.2969*Sqrt(X)-0.1260*X-0.3516*X^2+0.2843*X^3-0.1015*X^4),
	             0,
	             Ls};
  CM_x=CM_x+X;
EndFor	
pList[nPoints+1] = 6;
CM_x=CM_x+1.009;
For i In {1 : nPoints}
  X = L - L*i/(nPoints + 1);
  pList[nPoints+1+i]= nPoints+1+i+10;
  Point(pList[nPoints+1+i]) = {X,
	             -0.6*(0.2969*Sqrt(X)-0.1260*X-0.3516*X^2+0.2843*X^3-0.1015*X^4),
	             0,
	             Ls};
  CM_x=CM_x+X;
EndFor
pList[2*nPoints+2]=5;
Spline(100) = pList[];
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
