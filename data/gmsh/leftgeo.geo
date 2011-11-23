h =0.1;
// generation of circle
Point ( 1 ) = { 0 , 0 , 0 , h } ;
Point ( 2 ) = { 1 , 0 , 0 , h } ;
Point ( 3 ) = { 0 , 1 , 0 , h } ;
Point ( 4 ) = { -1 ,0 , 0 , h/5 } ;
Point ( 5 ) = { 0 ,-1 , 0 , h/5 } ;
Circle (1) ={2 ,1 ,3 };
//Line (1) ={ 2 ,3 };
Circle (2) ={3 ,1 ,4 };
Circle (3) ={4 ,1 ,5 };
Circle (4) ={5 ,1 ,2 };
// generation of shape
Line Loop ( 5 ) = { 1 , 2 , 3 , 4 } ;
// generation of surface
Plane Surface ( 6 ) ={5};
//Physical Line("inside") = { 1 } ;
//Physical Line("outside")= { 2, 3, 4 } ;
Physical Line(1) = { 2 };
Physical Line(2) = { 3 };
Physical Line(3) = { 1 };
Physical Line(4) = { 4 };
Physical Surface("Mat1") = {6};
