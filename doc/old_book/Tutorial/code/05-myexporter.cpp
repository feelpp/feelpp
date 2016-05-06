#include <feel/feel.hpp>

int main(int argc, char**argv )
{
  using namespace Feel;
  Environment env( _argc=argc, _argv=argv,
      _about=about(_name="myexporter",
        _author="Christophe Prud'homme",
        _email="christophe.prudhomme@feelpp.org"));

  // circle - geometrical order: 2
  auto mesh = unitCircle<2>(); 
  
  // circle - geometrical order: 1
  auto meshp1 = unitCircle<1>(); 
  
  // \( \mathbb{p}_2 \) space
  auto Xh = Pch<2>( mesh ); 

  auto myExpr = sin(pi*Px());

  auto v = project( _space=Xh, _range=elements(mesh),
      _expr=myExpr);

  auto exhi = exporter( _mesh=mesh, _name="exhi" );
  auto exlo = exporter( _mesh=meshp1, _name="exlo" );
  auto exhilo = exporter( _mesh=lagrangeP1(_space=Xh)->mesh(),_name="exhilo");


  int max = 10; double dt = 0.1;
  double time = 0;
  for (int i = 0; i<max; i++)
  {
    exhilo->step( time )->add( "vhilo", v );
    exlo->step( time )->add( "vlo", v );
    exhi->step( time )->add( "vhi", v );
    time += dt;
    //! [save]	
    exhi->save();
    exlo->save();
    exhilo->save();
    //! [save]	
  }



}

