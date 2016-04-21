#include <feel/feel.hpp>
#include <feel/feelfit/fit.hpp>
#include <feel/feelfit/fitdiff.hpp>

using namespace Feel;

int main(int argc, char **argv)
{
  std::vector<std::pair<double, double>> data;
  std::ifstream infile("data.txt");
  double a, b;
  while (infile >> a >> b)
  {
    data.push_back({a,b});
  }

  /*
   * Launch the Feel++ env
   */
  Environment env( _argc=argc, _argv=argv,
      _about=about( _name="test_fit" ,
        _author="Feel++ Consortium",
        _email="feelpp-devel@feelpp.org" ) );

  auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
  auto Xh = Pch<1>(mesh);
  auto T = Xh->element();
  auto K = Xh->element();
  auto Kd= Xh->element();
  auto M_ak = InterpolatorAkima(data);
  T.on(_range=elements(mesh), _expr=expr(soption("functions.f")));
  // evaluate K(T) with the interpolation from the datafile
  K.on(_range=elements(mesh), _expr=fit(idv(T),M_ak) );
  Kd.on(_range=elements(mesh), _expr=fitDiff(idv(T),M_ak) );

  auto e = exporter(_mesh=mesh);
  e->add("T",T);
  e->add("K",K);
  e->add("Kd",Kd);
  e->save();
  return 0;
}

