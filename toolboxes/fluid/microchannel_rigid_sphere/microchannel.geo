// Full 3D model of a cell in a RT-DC experiment
// Here, the cell is simply a rigid sphere

// /!\ The length unit is the millimetre

// We will need OpenCASCADE to define the sphere 
SetFactory("OpenCASCADE");

// mesh sizes adjusted to millimetre (reference paper values)
//hbulk = 0.001;
//hmembrane = 0.0002;

// for testing
hbulk = 0.005;
hmembrane = 0.005;

hsolvent = hbulk;
hcytosol = hbulk;


// Microchannel dimensions
// (thin) channel length : 300µm
// square cross-section side : 20µm
halfside=0.010;
channellength = 0.300;

// rectangle bounds
xmin = 0.0;
xmax = channellength;
ymin = -halfside;
ymax = halfside;
zmin = -halfside;
zmax = halfside;

// first square (inlet) vertices
Point(1) = {xmin,ymin,zmin,hsolvent};
Point(2) = {xmin,ymin,zmax,hsolvent};
Point(3) = {xmin,ymax,zmax,hsolvent};
Point(4) = {xmin,ymax,zmin,hsolvent};

// second square (outlet) vertices
Point(5) = {xmax,ymin,zmin,hsolvent};
Point(6) = {xmax,ymin,zmax,hsolvent};
Point(7) = {xmax,ymax,zmax,hsolvent};
Point(8) = {xmax,ymax,zmin,hsolvent};

// first square rectangle edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// second square rectangle edges
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// square to square edges (channel walls edges
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// first square perimeter
Curve Loop(14) = {1, 2, 3, 4};
Plane Surface(1) = {14};
// second square perimeter
Curve Loop(15) = {5, 6, 7, 8};
Plane Surface(2) = {15};
// bottom channel wall
Curve Loop(16) = {1, 10, -5, -9};
Plane Surface(3) = {16};
// right side channel wall
Curve Loop(17) = {2, 11, 6, 10};
Plane Surface(4) = {17};
// top channel wall
Curve Loop(18) = {3, 12, 7, 11};
Plane Surface(5) = {18};
// left side channel wall
Curve Loop(19) = {4, 9, -8, -12};
Plane Surface(6) = {19};

// inlet shrinking part length : (arbitrary) 10µm
shrinkinglength = 0.020;

// inlet large section square (larger)
// length : (arbitrary) 40µm
// square side : (arbitrary) 40µm
inletlength = 0.040;
inlethalfside = 0.5*inletlength;

// inlet square cuboid bounds
inxmin = xmin - shrinkinglength - inletlength;
inxmax = xmin - shrinkinglength;
inymin = -inlethalfside;
inymax = inlethalfside;
inzmin = -inlethalfside;
inzmax = inlethalfside;

// inlet square cuboid points
Point(9) = {inxmin,inymin,inzmin,hsolvent};
Point(10) = {inxmin,inymin,inzmax,hsolvent};
Point(11) = {inxmin,inymax,inzmax,hsolvent};
Point(12) = {inxmin,inymax,inzmin,hsolvent};

Point(13) = {inxmax,inymin,inzmin,hsolvent};
Point(14) = {inxmax,inymin,inzmax,hsolvent};
Point(15) = {inxmax,inymax,inzmax,hsolvent};
Point(16) = {inxmax,inymax,inzmin,hsolvent};

// inlet square cuboid edges
Line(20) = {9,10};
Line(21) = {10,11};
Line(22) = {11,12};
Line(23) = {12,9};

Line(24) = {13,14};
Line(25) = {14,15};
Line(26) = {15,16};
Line(27) = {16,13};

Line(28) = {9,13};
Line(29) = {10,14};
Line(30) = {11,15};
Line(31) = {12,16};

// shrinking part edges
Line(32) = {13,1};
Line(33) = {14,2};
Line(34) = {15,3};
Line(35) = {16,4};

// inlet cuboid and shrinking part faces perimeters
// Squares
Curve Loop(36) = {20,21,22,23};
Curve Loop(37) = {24,25,26,27};
Plane Surface(7) = {36};
Plane Surface(8) = {37};
// Cuboid walls
Curve Loop(38) = {20,29,-24,-28};
Curve Loop(39) = {21,30,-25,-29};
Curve Loop(40) = {22,31,-26,-30};
Curve Loop(41) = {23,28,-27,-31};
Plane Surface(9) = {38};
Plane Surface(10) = {39};
Plane Surface(11) = {40};
Plane Surface(12) = {41};
// Shrinking part walls
Curve Loop(42) = {24,33,-1,-32};
Curve Loop(43) = {25,34,-2,-33};
Curve Loop(44) = {26,35,-3,-34};
Curve Loop(45) = {27,32,-4,-35};
Plane Surface(13) = {42};
Plane Surface(14) = {43};
Plane Surface(15) = {44};
Plane Surface(16) = {45};



// =================== outlet ===================

// outlet shrinking part length : (arbitrary) 10µm
expandinglength = 0.020;
outletlentgh = inletlength;
outlethalfside = 0.5*outletlentgh;

// outlet square cuboid bounds
outxmin = xmax + expandinglength;
outxmax = xmax + expandinglength + outletlentgh;
outymin = -outlethalfside;
outymax = outlethalfside;
outzmin = -outlethalfside;
outzmax = outlethalfside;

// outlet square cuboid points
Point(17) = {outxmin,outymin,outzmin,hsolvent};
Point(18) = {outxmin,outymin,outzmax,hsolvent};
Point(19) = {outxmin,outymax,outzmax,hsolvent};
Point(20) = {outxmin,outymax,outzmin,hsolvent};

Point(21) = {outxmax,outymin,outzmin,hsolvent};
Point(22) = {outxmax,outymin,outzmax,hsolvent};
Point(23) = {outxmax,outymax,outzmax,hsolvent};
Point(24) = {outxmax,outymax,outzmin,hsolvent};

// outlet square cuboid edges
Line(46) = {17,18};
Line(47) = {18,19};
Line(48) = {19,20};
Line(49) = {20,17};

Line(50) = {21,22};
Line(51) = {22,23};
Line(52) = {23,24};
Line(53) = {24,21};

Line(54) = {17,21};
Line(55) = {18,22};
Line(56) = {19,23};
Line(57) = {20,24};

// expanding part edges
Line(58) = {5,17};
Line(59) = {6,18};
Line(60) = {7,19};
Line(61) = {8,20};

// outlet cuboid and expanding part faces perimeters
// Squares
Curve Loop(62) = {46,47,48,49};
Curve Loop(63) = {50,51,52,53};
Plane Surface(17) = {62};
Plane Surface(18) = {63};
// Cuboid walls
Curve Loop(64) = {46,55,-50,-54};
Curve Loop(65) = {47,56,-51,-55};
Curve Loop(66) = {48,57,-52,-56};
Curve Loop(67) = {49,54,-53,-57};
Plane Surface(19) = {64};
Plane Surface(20) = {65};
Plane Surface(21) = {66};
Plane Surface(22) = {67};

// Expanding part walls
Curve Loop(68) = {46,-59,-5,58};
Curve Loop(69) = {47,-60,-6,59};
Curve Loop(70) = {48,-61,-7,60};
Curve Loop(71) = {49,-58,-8,61};
Plane Surface(23) = {68};
Plane Surface(24) = {69};
Plane Surface(25) = {70};
Plane Surface(26) = {71};


// sphere parameters
// radius = 7.66µm
r = 0.00766;
sphereinitialoffset = xmin - shrinkinglength - inletlength  + 0.001 ;
centerx = sphereinitialoffset + r;
centery = ymin + (ymax-ymin)/2.0;
centerz = zmin + (zmax-zmin)/2.0;
Point(25) = {centerx, centery, centerz, hcytosol};

Sphere(1) = {centerx, centery, centerz, r, -Pi/2, Pi/2, 2*Pi};
// Now point 26 and 27 are created with the sphere
Characteristic Length {25} = hcytosol;
Characteristic Length {26} = hmembrane;
Characteristic Length {27} = hmembrane;

Physical Surface("ChannelInlet") = {7};
Physical Surface("ChannelOutlet") = {18};

// Channel wall 
Physical Surface("ChannelWall") = {3,4,5,6,9,10,11,12,13,14,15,16,19,20,21,22,23,24,25,26};

// Sphere surface
Physical Surface("Membrane") = {27};

// Solvent outer boundary = channel wall U inlet U outlet 
Surface Loop(3) = {7, 9, 10, 11, 12, 13, 14, 15, 16, 13, 4, 5, 6, 3, 25, 26, 23, 24, 20, 21, 22, 19, 18};

// Solvent inner boundary = fluid-structure interface = sphere = cell membrane outer boundary
Surface Loop(4) = {27};

// Cell volume
Physical Volume("Cell") = {1};

// Solvent volume
Volume(2) = {3, 4};
Physical Volume("Solvent") = {2};

