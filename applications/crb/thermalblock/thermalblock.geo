/*********************************************************************
 *
 *  Gmsh code
 *
 *
 *
 *********************************************************************/

hsize = .03;
nx = 5;
ny = 5;
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

Physical Line("south")={1:nx};
Physical Line("north")={nx*(ny+1)-nx+1:nx*(ny+1)};

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

u=0;

For i In {0:ny-1}

tab[i]=nx*(ny+1)+u+1;

tabl[i]=nx*(ny+1)+u+nx+1;

u += nx+1;

EndFor

Physical Line("east")=tab[];

Physical Line("west")=tabl[];

t5 = 0;
ne = (ny+1)*nx+1;

For j In {0:ny-1}

  For i In {0:nx-1}

t5 +=1;

 Line Loop(t5)={t5,(ne+t5+j),-(nx+t5),-(ne+t5+j-1)};

 Plane Surface(t5)={t5};

 Physical Surface(t5)={t5};

 EndFor

EndFor














