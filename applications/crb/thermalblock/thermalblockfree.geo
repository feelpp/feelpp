h = 0.1;
nx=3 ;
ny=3 ;
hsize=0.05;
dx = 1/nx;
dy = 1/ny;
t1=0;
x=0;
y=0;
For j In {0:ny}
  For i In {0:nx}
  y = j*dy ;
  x = i*dx ;
  t1+=1;
  Point(t1) = {x, y, 0,  hsize} ;
  EndFor
EndFor
t2=0;
For j In {0:ny}
  For i In {0:nx-1}
  t2 +=1;
  Line(t2)={t2+j,t2+j+1};
//  Physical Line(t2+1)={t2};
  EndFor
EndFor
Physical Line("south_domain-1")={1};
Physical Line("south_domain-2")={2};
Physical Line("south_domain-3")={3};
  Physical Line("north_domain-7") = {10};
  Physical Line("north_domain-8") = {11};
  Physical Line("north_domain-9") = {12};
t3 = (ny+1)*nx;
t4 = 0;
For i In {0:nx}
  For j In {0:ny-1}
  t3 +=1;
  t4 +=1;
 Line(t3)={t4,t4+nx+1};
// Physical Line(t3+1)={t3};
 EndFor
EndFor
Physical Line("west_domain-1")={13};
Physical Line("west_domain-4")={17};
Physical Line("west_domain-7")={21};
Physical Line("east_domain-3")={16};
Physical Line("east_domain-6")={20};
Physical Line("east_domain-9")={24};
t5 = 0;
ne = (ny+1)*nx+1;
For j In {0:ny-1}
  For i In {0:nx-1}
  t5 +=1;
  Line Loop(t5)={t5,(ne+t5+j),-(nx+t5),-(ne+t5+j-1)};
  Plane Surface(t5)={t5};
 EndFor
EndFor
  Physical Surface("domain-1")={1};
  Physical Surface("domain-2")={2};
  Physical Surface("domain-3")={3};
  Physical Surface("domain-4")={4};
  Physical Surface("domain-5")={5};
  Physical Surface("domain-6")={6};
  Physical Surface("domain-7")={7};
  Physical Surface("domain-8")={8};
  Physical Surface("domain-9")={9};
