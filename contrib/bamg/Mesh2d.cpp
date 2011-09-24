#include <cassert>
#include <fstream>
#include <iostream>
#include "ufunction.hpp"
#include "Mesh2d.hpp"

Mesh2::Mesh2(const char * filename)
  : nv(0),nt(0),nbe(0),
    area(0),peri(0),
    vertices(0),triangles(0),borderelements(0)
{ // read the mesh
  int i,iv[3],ir;
  ifstream f(filename);
  if(!f) {cerr << "Mesh2::Mesh2 Erreur openning " << filename<<endl;exit(1);}
  cout << " Read On file \"" <<filename<<"\""<<  endl;
  f >> nv >> nt >> nbe ;
  cout << " Nb of Vertex " << nv << " Nb of Triangle2s " << nt 
       << " Nb of Border Seg : " << nbe << endl;
  assert(f.good() && (nt>0) && (nv>0) &&( nbe>0) ) ;
  triangles = new Triangle2 [nt];
  vertices = new Vertex[nv];
  borderelements = new Seg[nbe];
  area=0;
  assert(triangles && vertices && borderelements);
  for (i=0;i<nv;i++)    
    {
      f >> vertices[i];
      assert(f.good());
    }
  for (i=0;i<nt;i++) 
    { 
      f >> iv[0] >> iv[1] >> iv[2] >> ir;
      assert(f.good() && iv[0]>0 && iv[0]<=nv && iv[1]>0 && iv[1]<=nv && iv[2]>0 && iv[2]<=nv);
      for (int v=0;v<3;++v) iv[v]--;
      triangles[i].init(vertices,iv,ir); 
      area += triangles[i].area;
    }
  for (i=0;i<nbe;i++) 
    { 
      f >> iv[0] >> iv[1]  >> ir;
      assert(f.good() && iv[0]>0 && iv[0]<=nv && iv[1]>0 && iv[1]<=nv);
      for (int v=0;v<2;++v) iv[v]--;
      borderelements[i].init(vertices,iv,ir); 
      peri += borderelements[i].l;
    }

  cout << " End of read: area = " << area << "  perimeter: " << peri << endl;  
}
