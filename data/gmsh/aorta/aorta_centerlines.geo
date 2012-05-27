Mesh.Algorithm = 6; //(1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
Mesh.Algorithm3D = 1; //(1=tetgen, 4=netgen, 5=FrontalDel, 6=FrontalHex, 7=MMG3D, 9=R-tree

//Mesh.RecombineAll = 1;
//Mesh.Bunin = 150;

Merge "aorta2.stl";

Field[1] = Centerline;
Field[1].FileName = "centerlinesAORTA.vtk";
Field[1].nbPoints = 25;
Field[1].closeVolume =1;

Field[1].extrudeWall =1;
Field[1].nbElemLayer = 4;
Field[1].hLayer = 0.2;//percent of vessel radius

Field[1].reMesh =1;

Field[1].run;

Background Field = 1;

