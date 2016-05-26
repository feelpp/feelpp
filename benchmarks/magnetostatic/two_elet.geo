h=10;
p1 = newp; Point(p1) = {0,0,0,h};
p2 = newp; Point(p2) = {1,0,0,h};
p3 = newp; Point(p3) = {0,1,0,h};
p4 = newp; Point(p4) = {0,0,1,h};
p5 = newp; Point(p5) = {-1,0,0,h};

l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p1};
l4 = newl; Line(l4) = {p2,p4};
l5 = newl; Line(l5) = {p4,p1};
l6 = newl; Line(l6) = {p3,p4};

l7 = newl; Line(l7) = {p4,p5};
l8 = newl; Line(l8) = {p3,p5};
l9 = newl; Line(l9) = {p1,p5};

ll1 = newll; Line Loop(ll1) = {l1,l2,l3};
ll2 = newll; Line Loop(ll2) = {l1,l4,l5};
ll3 = newll; Line Loop(ll3) = {l2,l6,-l4};
ll4 = newll; Line Loop(ll4) = {l3,-l5,-l6};

ll5 = newll; Line Loop(ll5) = {l9,-l7,l5};
ll6 = newll; Line Loop(ll6) = {-l3,l8,-l9};
ll7 = newll; Line Loop(ll7) = {l6,l7,-l8};

s1 = news; Plane Surface(s1) = {ll1};
s2 = news; Plane Surface(s2) = {ll2};
s3 = news; Plane Surface(s3) = {ll3};

s4 = news; Plane Surface(s4) = {ll4};

s5 = news; Plane Surface(s5) = {ll5};
s6 = news; Plane Surface(s6) = {ll6};
s7 = news; Plane Surface(s7) = {ll7};

Physical Surface("s1")={s1};
Physical Surface("s2")={s2};
Physical Surface("s3")={s3};
Physical Surface("s4")={s4};
Physical Surface("s5")={s5};
Physical Surface("s6")={s6};
Physical Surface("s7")={s7};

sl1 = newsl; Surface Loop(sl1) = {s1,s2,s3,s4};
sl2 = newsl; Surface Loop(sl2) = {s5,s6,s7,s4};
v1 = newv; Volume(v1) = {sl1};
v2 = newv; Volume(v2) = {sl2};

Physical Volume("v1") = {v1};
Physical Volume("v2") = {v2};
