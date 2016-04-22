#include <feel/feel.hpp>
#include <feel/feelfit/fit.hpp>
#include <feel/feelfit/fitdiff.hpp>

using namespace Feel;

int main(int argc, char **argv)
{
  Environment env( _argc=argc, _argv=argv,
      _about=about( _name="test_fit" ,
        _author="Feel++ Consortium",
        _email="feelpp-devel@feelpp.org" ) );

  auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
  auto e = exporter(_mesh=mesh);
  auto Xh = Pch<1>(mesh);
  auto T = Xh->element(); // the temperature (say)
  auto K = Xh->element(); // K(T) - the dependant of the temperature conductivity
  auto Kd= Xh->element(); // K'(T)
  T.on(_range=elements(mesh), _expr=expr(soption("functions.f")));

  for(int i = 0; i < 4; i++)
  {
    // evaluate K(T) with the interpolation from the datafile
    K.on(_range=elements(mesh), _expr=fit(idv(T),soption("fit.datafile"),i));
    Kd.on(_range=elements(mesh), _expr=fitDiff(idv(T),soption("fit.datafile"), i) );

    e->step(i)->add("T",T);
    e->step(i)->add("K",K);
    e->step(i)->add("Kd",Kd);
    e->save();
  }
  return 0;
}

