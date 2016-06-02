#include <feel/feel.hpp>

int main( int argc, char* argv[] )
{
  using namespace Feel;

  po::options_description advectionoptions( "advection options" );
  advectionoptions.add_options()
      ("epsilon", Feel::po::value<double>()->default_value(0.), "diffusion coefficient")
      ("mu", Feel::po::value<double>()->default_value(0.), "reaction coefficient")
      ("betax", Feel::po::value<double>()->default_value(1.), "x advection velocity")
      ("betay", Feel::po::value<double>()->default_value(0.), "y advection velocity")
      ;

  Environment env( 
          _argc=argc, _argv=argv,
          _desc=advectionoptions,
          _about=about( _name="feelpp_advection",
              _author="Feel++ Consortium",
              _email="feelpp-devel@feelpp.org") );
  // create mesh
  auto mesh = unitSquare();
  // function space
  auto Xh = Pch<1>( mesh );
  auto u = Xh->element( "u" );
  auto v = Xh->element( "v" );
  // diffusion coeff.
  double epsilon = doption(_name="epsilon");
  // reaction coeff.
  double mu = doption(_name="mu");
  auto beta = vec( cst(doption(_name="betax")),
		   cst(doption(_name="betay")) );
  auto f = cst(1.);
  // left hand side
  auto a = form2( _test=Xh, _trial=Xh );
  a += integrate( _range=elements( mesh ),
		  _expr=( epsilon*gradt( u )*trans( grad( v ) )
			  + ( gradt( u )*beta )*id(v)
			  + mu*idt( u )*id( v ) ) );
  // right hand side
  auto l = form1( _test=Xh );
  l+= integrate( _range=elements( mesh ), _expr=f*id( v ) );
    
  // boundary condition
  a +=  on( _range=boundaryfaces( mesh ), _rhs=l, _element=u,
	    _expr=cst(0.) );
           
  // solve the system
  a.solve( _rhs=l, _solution=u );
  // export results
  auto e = exporter( _mesh=mesh );
  e->add("u",u);
  e->save(); 
} // end main
