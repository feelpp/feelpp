h=0.1; 

re=2.;
ri=1.;
l=1.;

Point(0) = {ri,0,0,h};
Point(1) = {re,0,0,h};
Line(0) = {0,1};
out[] = Extrude{{0,0,1},{0,0,0}, Pi/2} { Line{0};};
outV[] = Extrude{0,0,l} {Surface{out[1]};};

Physical Surface("top") = {outV[0]};
Physical Surface("bottom") = {out[1]};
Physical Surface("in") = {outV[2]};
Physical Surface("ext") = {outV[3]};
Physical Surface("out") = {outV[4]};
Physical Surface("int") = {outV[5]};


// Physical Surface("all") = {outV[0], out[1],outV[2],outV[3],outV[4],outV[5]};
Physical Volume("omega") = {outV[1]};




