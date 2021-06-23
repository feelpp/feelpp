h=0.1;

re=2.;
ri=1.;
angle=Pi/2.;

Point(0) = {ri,0,0,h};
Point(1) = {re,0,0,h};

Line(1) = {0,1};

out[] = Extrude{{0,0,1},{0,0,0}, angle} { Line{1};};

Physical Line("in") = {1};
Physical Line("out") = {out[0]};
Physical Line("ext") = {out[2]};
Physical Line("int") = {out[3]};
Physical Surface("omega") = {out[1]};

Mesh.ElementOrder = 2;

