// -*- mode: c++ -*-
h=1;
Point(0) ={0    , 0   , 0 , h};
Point(1) ={1/3  , 0   , 0 , h};
Point(2) ={2/3  , 0   , 0 , h};
Point(3) ={3/3  , 0   , 0 , h};
Point(4) ={0    , 1/3 , 0 , h};
Point(5) ={1/3  , 1/3 , 0 , h};
Point(6) ={2/3  , 1/3 , 0 , h};
Point(7) ={3/3  , 1/3 , 0 , h};
Point(8) ={0    , 2/3 , 0 , h};
Point(9) ={1/3  , 2/3 , 0 , h};
Point(10)={2/3  , 2/3 , 0 , h};
Point(11)={3/3  , 2/3 , 0 , h};
Point(12)={0    , 3/3 , 0 , h};
Point(13)={1/3  , 3/3 , 0 , h};
Point(14)={2/3  , 3/3 , 0 , h};
Point(15)={3/3  , 3/3 , 0 , h};
Line(101)={0,1};//south_domain-1
Line(102)={1,5};
Line(103)={5,4};
Line(104)={4,0};
Line(105)={1,2};//south_domain-2
Line(106)={2,6};
Line(107)={6,5};
Line(108)={2,3};//south_domain-3
Line(109)={3,7};
Line(110)={7,6};
Line(111)={5,9};
Line(112)={9,8};
Line(113)={8,4};
Line(114)={6,10};
Line(115)={10,9};
Line(116)={7,11};
Line(117)={11,10};
Line(118)={9,13};
Line(119)={13,12};//north_domain-7
Line(120)={12,8};
Line(121)={10,14};
Line(122)={14,13};//north_domain-8
Line(123)={11,15};
Line(124)={15,14};//north_domain-9
Line Loop (201) = {101,102,103,104};//domain-1
Line Loop (202) = {105,106,107,-102};//domain-2
Line Loop (203) = {108,109,110,-106};//domain-3
Line Loop (204) = {-103,111,112,113};//domain-4
Line Loop (205) = {-107,114,115,-111};//domain-5
Line Loop (206) = {-110,116,117,-114};//domain-6
Line Loop (207) = {-112,118,119,120};//domain-7
Line Loop (208) = {-115,121,122,-118};//domain-8
Line Loop (209) = {-117,123,124,-121};//domain-9
Plane Surface (301) = {201};
Plane Surface (302) = {202};
Plane Surface (303) = {203};
Plane Surface (304) = {204};
Plane Surface (305) = {205};
Plane Surface (306) = {206};
Plane Surface (307) = {207};
Plane Surface (308) = {208};
Plane Surface (309) = {209};
Physical Line ("south_domain-1")={101};
Physical Line ("south_domain-2")={105};
Physical Line ("south_domain-3")={108};
Physical Line ("north_domain-7")={119};
Physical Line ("north_domain-8")={122};
Physical Line ("north_domain-9")={124};
Physical Surface ("domain-1") = {301};
Physical Surface ("domain-2") = {302};
Physical Surface ("domain-3") = {303};
Physical Surface ("domain-4") = {304};
Physical Surface ("domain-5") = {305};
Physical Surface ("domain-6") = {306};
Physical Surface ("domain-7") = {307};
Physical Surface ("domain-8") = {308};
Physical Surface ("domain-9") = {309};








