Mesh.Algorithm3D = 4; // frontal
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFactor = 0.15;

Merge "aneurysm.stl";
CreateTopology;

//Extrude outward and inward with 4 layers of thickness 0.5
out1[] = Extrude{Surface{1}; Layers{4, 0.5}; Using Index[0]; };
out2[] = Extrude{Surface{1}; Layers{4, -0.5}; Using Index[1]; };

Printf("inward extrusion: top_surf=%g volume=%g lateral_surf=%g,%g,%g", 
       out2[0], out2[1], out2[2], out2[3], out2[4]);
b1[] = Boundary{ Surface{out2[2]}; };
b2[] = Boundary{ Surface{out2[3]}; };
b3[] = Boundary{ Surface{out2[4]}; };
Printf("interior curves of lateral surfaces = %g %g %g",
       b1[2], b2[2], b3[2]);

//create inlet faces
Line Loop(200)={b1[2]};
Plane Surface(201)={200};

Line Loop(300)={b2[2]};
Plane Surface(301)={300};

//create outlet face
Line Loop(400)={b3[2]};
Plane Surface(401)={400};

//create inside volume
Surface Loop(999)={out2[0], 201,301,401};
Volume(1000) = {999};

//save only physicals
Physical Surface("inlet")={out2[2], 201, out2[3], 301};
Physical Surface("outlet")={out2[4], 401};
Physical Volume("fluid BL")={out2[1]};
Physical Volume("wall")={out1[1]};
Physical Volume("inside")={1000};
