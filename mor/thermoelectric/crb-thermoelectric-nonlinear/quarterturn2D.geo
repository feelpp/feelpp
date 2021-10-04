h = 0.1;

Point(1) = {1,0,0,h};
Point(2) = {2,0,0,h};

Line(1) = {1,2};
out[] = Extrude{{0,0,1},{0,0,0},Pi/2} {Line{1}; };

Physical Surface("omega") = {out[1]};
Physical Line("V0") = {1};
Physical Line("V1") = {out[0]};
Physical Line("R") = {out[2], out[3]};
