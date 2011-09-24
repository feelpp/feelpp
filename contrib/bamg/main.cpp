#include <iostream>
using namespace std;
#include "Mesh2.h"
#include "RNM.hpp"
#include "Mesh2d.hpp"
#include <set>
#include "bamg-gmsh.hpp"
long verbosity=1;
int main(int argc,const char ** argv)
{
  assert(argc>1);
  Mesh2 Th(argv[1]);
  KN<double> data(30);
  data = -1e200;
  int nv = Th.nv;
  KN<double> m11(nv),m12(nv),m22(nv);
  double h=0.01;
  m11=1/(h*h);
  m22=1/(h*h);
  m12=0.;
  
}
