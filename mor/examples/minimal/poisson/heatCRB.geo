h = 0.1;

center=newp; Point(center) = {0,0,0,h};
leftCircle = newp; Point(leftCircle) = {-0.5,0,0,h};
rightCircle = newp; Point(rightCircle) = {0.5,0,0,h};
bottomLeft = newp; Point(bottomLeft) = {-1,-1,0,h};
bottomRight = newp; Point(bottomRight) = {1,-1,0,h};
topRight = newp; Point(topRight) = {1,1,0,h};
topLeft = newp; Point(topLeft) = {-1,1,0,h};

circ1 = newl; Circle(circ1) = {leftCircle,center,rightCircle};
circ2 = newl; Circle(circ2) = {rightCircle,center,leftCircle};
baseLine = newl; Line(baseLine) = {bottomLeft,bottomRight};
sideLine1 = newl; Line(sideLine1) = {bottomRight,topRight};
topLine = newl; Line(topLine) = {topRight,topLeft};
sideLine2 = newl; Line(sideLine2) = {topLeft,bottomLeft};

circLoop = newll; Line Loop(circLoop) = {circ1,circ2};
squareLoop = newll; Line Loop(squareLoop) = {baseLine,sideLine1,topLine,sideLine2};

circ = news; Plane Surface(circ) = {circLoop};
square = news; Plane Surface(square) = {squareLoop,circLoop};

Physical Line("base") = {baseLine};
Physical Line("side") = {sideLine1,sideLine2};
Physical Line("top") = {topLine};
Physical Surface("omega0") = {circ};
Physical Surface("omega1") = {square};
