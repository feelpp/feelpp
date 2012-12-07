/********************************************************************* 
 *
 *  Gmsh 
 * 
 *  Variables, elementary entities (points, lines, surfaces), physical
 *  entities (points, lines, surfaces)
 *
 *********************************************************************/

// The simplest construction in Gmsh's scripting language is the
// `affectation'. The following command defines a new variable `lc':

lc1 = 0.2;
lc1d = 0.15;
lc1b = 0.03;
lc1c= 0.1;
lc2 = 0.2;

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is defined by a list of
// four numbers: three coordinates (X, Y and Z), and a characteristic
// length (lc) that sets the target element size at the point:

// The distribution of the mesh element sizes is then obtained by
// interpolation of these characteristic lengths throughout the
// geometry. Another method to specify characteristic lengths is to
// use a background mesh (see `t7.geo' and `bgmesh.pos').

// We can then define some additional points as well as our first
// curve.  Curves are Gmsh's second type of elementery entities, and,
// amongst curves, straight lines are the simplest. A straight line is
// defined by a list of point numbers. In the commands below, for
// example, the line 1 starts at point 1 and ends at point 2:

Point(1) = {-1.9, -0.05, 0, lc1};
Point(2) = {1.9, -0.05,  0, lc1} ;
Point(3) = {1.9, 0.05, 0, lc1} ;
Point(4) = {1.9993749, 0.05,0, lc1} ;
Point(5) = {1.812616, 0.845236,  0, lc1} ;	
Point(6) = {1.3, 1.519868, 0, lc1} ;
Point(7) = {0.46915, 2, 0, lc1} ;	// point on major axis of right side ellipse
Point(8) = {0.125, 1.996089, 0, lc1b};
Point(9) = {0.05, 1.9993749, 0, lc1};
Point(10) = {0, 2, 0, lc1b};
Point(11) = {-1.9, 0.05, 0, lc1};
Point(12) = {-1.9993749,0.05,0, lc1} ;
Point(13) = {-1.812616, 0.845236,  0, lc1} ;
Point(14) = {-1.3, 1.519868, 0, lc1} ;
Point(15) = {-0.46915, 2, 0, lc1} ;	// point on major axis of left side ellipse
Point(16) = {-0.125, 1.996089, 0, lc1b} ;
Point(17) = {-0.05, 1.9993749, 0, lc1} ;
Point(18) = {-1.3, 2, 0, lc1};		// left side ellipse center
Point(19) = {1.3, 2, 0, lc1};		// right side ellipse center
Point(20) = {0, 0, 0, lc1};		// center of circle
Point(21) = {0.05, 1.875, 0, lc1};
Point(22) = {-0.05, 1.875, 0, lc1};

Point (23) = {-1.55, 0.50, 0, lc1b};	//center of P1
Point (24) = {-1.55, 0.80, 0, lc1b};	// point on major axis of P1
Point (25) = {-1.55, 1.00, 0, lc1b};
Point (26) = {-1.55, 0.00, 0, lc1b};
Point (27) = {-1.35, 0.50, 0, lc1c};
Point (28) = {-1.75, 0.50, 0, lc1c};

Point (29) = {-1.05, 0.50, 0, lc1b};	//center of P2
Point (30) = {-1.05, 0.80, 0, lc1b};	// point on major axis of P2
Point (31) = {-1.05, 1.00, 0, lc1b};
Point (32) = {-1.05, 0.00, 0, lc1b};
Point (33) = {-0.85, 0.50, 0, lc1c};
Point (34) = {-1.25, 0.50, 0, lc1c};

Point (35) = {-0.55, 0.50, 0, lc1b};	//center of P3
Point (36) = {-0.55, 0.80, 0, lc1b};	// point on major axis of P3
Point (37) = {-0.55, 1.00, 0, lc1b};
Point (38) = {-0.55, 0.00, 0, lc1b};
Point (39) = {-0.35, 0.50, 0, lc1d};
Point (40) = {-0.75, 0.50, 0, lc1c};

Point (41) = { 0.55, 0.50, 0, lc1b};	//center of P4
Point (42) = { 0.55, 0.80, 0, lc1b};	// point on major axis of P4
Point (43) = { 0.55, 1.00, 0, lc1b};
Point (44) = { 0.55, 0.00, 0, lc1b};
Point (45) = { 0.75, 0.50, 0, lc1c};
Point (46) = { 0.35, 0.50, 0, lc1d};

Point (47) = { 1.05, 0.50, 0, lc1b};	//center of P5
Point (48) = { 1.05, 0.80, 0, lc1b};	// point on major axis of P5
Point (49) = { 1.05, 1.00, 0, lc1b};
Point (50) = { 1.05, 0.00, 0, lc1b};
Point (51) = { 1.25, 0.50, 0, lc1c};
Point (52) = { 0.85, 0.50, 0, lc1c};

Point (53) = { 1.55, 0.50, 0, lc1b};	//center of P6
Point (54) = { 1.55, 0.80, 0, lc1b};	// point on major axis of P6
Point (55) = { 1.55, 1.00, 0, lc1b};
Point (56) = { 1.55, 0.00, 0, lc1b};
Point (57) = { 1.75, 0.50, 0, lc1c};
Point (58) = { 1.35, 0.50, 0, lc1c};

Point(59)={-1.9993749,-0.05, 0, lc1};		//Point(59)={-1.9993749,-0.05,0, lc1};
Point(60)={-1.130265, -1.65, 0, lc2};
Point(61)={ 1.130265, -1.65, 0, lc2};
Point(62)={ 1.9993749,-0.05, 0, lc1};		//Point(62)={ 1.9993749,-0.05,0, lc1};

Point(63)={-0.15,-1.65,0, lc2};          // bay inlet
Point(64)={ 0.20,-1.65,0, lc2};
Point(65)={ 0.20,-1.55,0, lc2};
Point(66)={-0.15,-1.55,0, lc2};

Point(67)={-1.0,-1.5,0, lc1d};            // Core
Point(68)={-0.5,-1.5,0, lc1d};
Point(69)={-0.5,-0.8,0, lc1d};
Point(70)={-1.0,-0.8,0, lc1d};

 
Point(71)={ 0.5,-1.5,0, lc1};            // Bay
Point(72)={ 0.9,-1.5,0, lc1};           // Bay
Point(772)={ 0.7,-1.5,0, lc1c};
Point(73)={ 0.9,-0.4,0, lc1};
Point(74)={ 0.5,-0.4,0, lc1};
Point(774)={ 0.7,-0.4,0, lc1c};

Point(75)={-1.0,-0.2,0, lc2};            // left bay outlet
Point(76)={-0.55,-0.2,0, lc2};         // left bay outlet
Point(776)={-0.775,-0.2,0, lc1c};
Point(77)={-0.55,-0.1,0, lc2};
Point(78)={-1.0,-0.1,0, lc2};
Point(778)={-0.775,-0.1,0, lc1c};


Point(79)={0.5,-0.2,0, lc2};             // right bay outlet
Point(880)={0.725,-0.2,0, lc2};           // right bay outlet
Point(80)={0.95,-0.2,0, lc2};
Point(81)={0.95,-0.1,0, lc2};
Point(882)={0.725,-0.1,0, lc2};
Point(82)={0.5,-0.1,0, lc2};


//Creation des points pour les formes circulaires

Point(83)={-1.5,-0.2,0, lc1};     //creation centre cercle	  // EH1  
Point(84)={-1.5,-0.25,0, lc1};	  //point debut premier arc de cercle	
Point(85)={-1.45,-0.2,0, lc1};	  //point debut deuxieme arc de cercle	
Point(86)={-1.5,-0.15,0, lc1};	  //point debut troisieme arc de cercle 
Point(87)={-1.55,-0.2,0, lc1};	  //point quatrieme arc de cercle	

Point(89)={0,-0.2,0, lc1};               // EH2
Point(90)={0,-0.25,0, lc1};
Point(91)={0.05,-0.2,0, lc1};
Point(92)={0,-0.15,0, lc1};
Point(93)={-0.05,-0.2,0, lc1};

Point(94)={1.6,-0.2,0, lc1};             // EH3
Point(95)={1.6,-0.25,0, lc1};
Point(96)={1.65,-0.2,0, lc1};
Point(97)={1.6,-0.15,0, lc1};
Point(98)={1.55,-0.2,0, lc1};

Point(99)={-1.8,-0.2,0, lc1};            // BD1
Point(100)={-1.8,-0.25,0, lc1};
Point(101)={-1.75,-0.2,0, lc1};
Point(102)={-1.8,-0.15,0, lc1};
Point(103)={-1.85,-0.2,0, lc1};

Point(104)={-0.2,-0.2,0, lc1};           // BD2
Point(105)={-0.2,-0.25,0, lc1};
Point(106)={-0.15,-0.2,0, lc1};
Point(107)={-0.2,-0.15,0, lc1};
Point(108)={-0.25,-0.2,0, lc1};

Point(109)={ 1.3,-0.2,0, lc1};            // BD3
Point(110)={ 1.3,-0.25,0, lc1};
Point(111)={ 1.35,-0.2,0, lc1};
Point(112)={ 1.3,-0.15,0, lc1};
Point(113)={ 1.25,-0.2,0, lc1};


//creation contour cabine+plancher
Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {11,1} ;
Line(5) = {12,11} ;
Circle (6) = {4, 20, 6};
Ellipse (7) = {6, 19, 7, 8};	     // Right HL ellipse
Circle (8) = {8, 20, 9};
Circle (9) = {14, 20, 12};
Ellipse (10) = {16, 18, 15, 14};	 // Left HL ellipse
Circle (11) = {17, 20, 16};
Line(12) = {9,21} ;
Line(13) = {21,22} ;
Line(14) = {22,17} ;
Circle (15) = {6, 20, 8};
Circle (16) = {9, 20, 17};
Circle (17) = {16, 20, 14};
//creation passagers
Ellipse (18) = {26, 23, 24, 27};     // P1
Ellipse (19) = {27, 23, 24, 25};
Ellipse (20) = {25, 23, 24, 28};
Ellipse (21) = {28, 23, 24, 26};
//
Ellipse (22) = {32, 29, 30, 33};     // p2
Ellipse (23) = {33, 29, 30, 31};
Ellipse (24) = {31, 29, 30, 34};
Ellipse (25) = {34, 29, 30, 32};
//
Ellipse (26) = {38, 35, 36, 39};     // p3
Ellipse (27) = {39, 35, 36, 37};
Ellipse (28) = {37, 35, 36, 40};
Ellipse (29) = {40, 35, 36, 38};
//
Ellipse (30) = {44, 41, 42, 45};     // p4
Ellipse (31) = {45, 41, 42, 43};
Ellipse (32) = {43, 41, 42, 46};
Ellipse (33) = {46, 41, 42, 44};
//
Ellipse (34) = {50, 47, 48, 51};     // p5
Ellipse (35) = {51, 47, 48, 49};
Ellipse (36) = {49, 47, 48, 52};
Ellipse (37) = {52, 47, 48, 50};
//
Ellipse (38) = {56, 53, 54, 57};      // p6
Ellipse (39) = {57, 53, 54, 55};
Ellipse (40) = {55, 53, 54, 58};
Ellipse (41) = {58, 53, 54, 56};

//creation contour bay+entrees
Line(42)={59,1};
Line(43)={2,62};

Circle(44)={12,20,59};     
Circle(45)={59,20,60};		          //Circle(45)={59,20,60};

Line(46)={60,63};
Line(47)={63,64};
Line(48)={64,61};

Circle(49)={61,20,62};   	          //Circle(49)={61,20,62};  
Circle(50)={62,20,4};      

Line(51)={64,65};                     // bay inlet
Line(52)={65,66};
Line(53)={66,63};

//creation element internes bay
Line(54)={67,68};                     // Core
Line(55)={68,69};
Line(56)={69,70};
Line(57)={70,67};

//Line(58)={71,72};   
Line(58)={71,772};    
Line(558)={772,72};                     // bay (element interne rectangle)
Line(59)={72,73};
//Line(60)={73,74};
Line(60)={73,774};
Line(660)={774,74};
Line(61)={74,71};

//Line(62)={75,76};                     // left bay outlet //LBO
Line(62)={75,776};                     // left bay outlet //LBO
Line(662)={776,76};                     // left bay outlet //LBO
Line(63)={76,77};
//Line(64)={77,78};
Line(64)={77,778};
Line(664)={778,78};
Line(65)={78,75};

//Line(66)={79,80}; 
Line(66)={79,880}; 
Line(666)={880,80};                    // right bay outlet //RBO
Line(67)={80,81};
//Line(68)={81,82};
Line(68)={81,882};
Line(668)={882,82};
Line(69)={82,79};



Circle(70)={84,83,85};               // EH1
Circle(71)={85,83,86};		//generation des arcs de cercles
Circle(72)={86,83,87};
Circle(73)={87,83,84};

Circle(74)={90,89,91};               // EH2
Circle(75)={91,89,92};
Circle(76)={92,89,93};
Circle(77)={93,89,90};

Circle(78)={95,94,96};               // EH3
Circle(79)={96,94,97};
Circle(80)={97,94,98};
Circle(81)={98,94,95};

Circle(82)={100,99,101};             // BD1
Circle(83)={101,99,102};
Circle(84)={102,99,103};
Circle(85)={103,99,100};

Circle(86)={105,104,106};            // BD2
Circle(87)={106,104,107};
Circle(88)={107,104,108};
Circle(89)={108,104,105};

Circle(90)={110,109,111};            // BD3
Circle(91)={111,109,112};
Circle(92)={112,109,113};
Circle(93)={113,109,110};


Line Loop(1) = {1,2,3,6,7,8,12,13,14,11,10,9,5,4} ;     // Cabin Fluid
Line Loop(2) = {18,19,20,21} ;                          // P1
Line Loop(3) = {22,23,24,25} ;                          // P2
Line Loop(4) = {26,27,28,29};                           // P3
Line Loop(5) = {30,31,32,33};                           // P4
Line Loop(6) = {34,35,36,37};                           // P5
Line Loop(7) = {38,39,40,41};                           // P6

Plane Surface(8) = {1,2,3,4,5,6,7};		                // cabin fluid

Line Loop(15) = {45,46,-53,-52,-51,48,49,-43,-1,-42};   // eBay   
Line Loop(16) = {54,55,56,57};                          // Core
//Line Loop(17) = {58,59,60,61};                          // Bay
Line Loop(17) = {58,558,59,60,660,61};                          // Bay
Line Loop(18) = {62,662,63,64,664,65};                          // left outlet
Line Loop(19) = {66,666,67,68,668,69};                          // right outlet



Line Loop(20) = {70,71,72,73};      //raccord des arcs de cercle                    // EH1
Line Loop(21) = {74,75,76,77};                          // EH2
Line Loop(22) = {78,79,80,81};                          // EH3
Line Loop(23) = {82,83,84,85};                          // BD1
Line Loop(24) = {86,87,88,89};                          // BD2
Line Loop(25) = {90,91,92,93};                          // BD3


//Plane Surface(26) = {15,16,17,18,19,20,21,22,23,24,25};  // eBay avec cercles
Plane Surface(26) = {15,16,17,18,19}; // bay sans les cercles




Physical Surface("cabin") = {8};
 
Physical Surface("bay_Surface") = {26};
Physical Line("passenger1") = {20, 19, 21, 18};
Physical Line("passenger2") = {24, 23, 25, 22};
Physical Line("passenger3") = {28, 27, 29, 26};
Physical Line("passenger4") = {32, 31, 33, 30};
Physical Line("passenger5") = {36, 35, 37, 34};
Physical Line("passenger6") = {40, 39, 41, 38};
Physical Line("cabin_wall") = {9, 6, 10, 7, 8, 11};
Physical Line("plancher") = {1};
Physical Line("bay_wall") = {45, 49, 48, 46};
Physical Line("BD1") = {84, 83, 82, 85};
Physical Line("EH1") = {72, 71, 70, 73};
//Physical Line("LBO") = {65, 64, 63, 62};
Physical Line("LBO") = {65,664, 64, 63, 62,662};
Physical Line("BD2") = {88, 87, 86, 89};
Physical Line("EH2") = {76, 75, 74, 77};
Physical Line("RBO") = {69, 668,68, 67, 666,66};
//Physical Line("RBO") = {69, 68, 67, 66};
Physical Line("BD3") = {92, 91, 90, 93};
Physical Line("EH3") = {80, 79, 78, 81};
Physical Line("Core") = {56, 55, 54, 57};
Physical Line("bay") = {60, 59, 558,58,61,660};
//Physical Line("bay") = {60, 59, 58, 61};
Physical Line("bay_inlet") = {53, 52, 51};
Physical Line("left_outlet") = {5, 4};//Physical Line("left_outlet") = {5, 4, 42};
Physical Line("right_outlet") = {3, 2};//Physical Line("right_outlet") = {3, 2, 43};
Physical Line("cabine_inlet") = {14, 13, 12};
