h = 0.1;
Point(1) = {0,0,0,h};
Point(2) = {1,0,0,h};
Line(1) = {1,2};
out[] = Extrude{0,1,0} {Line{1};};

Physical Surface("omega") = {out[1]};
Physical Line("gamma1") = {1};
Physical Line("gamma2") = {out[0]};
Physical Line("gamma3") = {out[2]};
Physical Line("gammb1") = {out[3]};
