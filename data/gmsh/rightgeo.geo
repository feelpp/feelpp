h =0.1;
// generation of rectangle
Point ( 1 ) ={ 1 ,0 ,0 , h } ;
Point ( 2 ) ={ 2 ,0 ,0 , h } ;
Point ( 3 ) ={ 2 ,1 ,0 , h/2 } ;
Point ( 4 ) ={ 0 ,1 ,0 , h/2 } ;
Line( 1 ) ={1 ,2};
Line( 2 ) ={2 ,3};
Line( 3 ) ={3 ,4};
Line( 4 ) ={4 ,1};
// generation of shape
Line Loop ( 5 ) = { 1 , 2 , 3 , 4 } ;
// generation de la surface
Plane Surface ( 6 ) ={5};
//Physical Line("inside") = { 4 } ;
//Physical Line("outside")= { 1, 2, 3 } ;
Physical Line(1) = { 4 };
Physical Line(2) = { 1 };
Physical Line(3) = { 2 };
Physical Line(4) = { 3 };
Physical Surface("Mat1") = {6};
