h=0.1;

Point(0) = {0,-0.5,0,h};
Point(1) = {0,1.5,0,h};

Line(0) = {0,1};

out[] = Extrude{2,0,0} {Line{0}; Layers{2/h};};

Physical Surface("fluid") = {out[1]};
Physical Line("left") = {0};
Physical Line("top") = {out[2]};
Physical Line("right") = {out[0]};
Physical Line("bottom") = {out[3]};
