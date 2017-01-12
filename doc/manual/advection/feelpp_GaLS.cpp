#include <feel/feel.hpp>

int main( int argc, char* argv[] )
{
  using namespace Feel;

  Environment env( _argc=argc, _argv=argv,
		   _about=about( _name="env",
				 _author="Feel++ Consortium",
				 _email="feelpp-devel@feelpp.org") );

  auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
  auto Vh = Pch<1>( mesh );
  auto u = Vh->element();
  auto v = Vh->element();
  //# endmarker2 #
   
  auto f = expr( soption(_name="f",_prefix="functions") );
  auto beta_x = expr( soption(_name="beta_x",_prefix="functions") );
  auto beta_y = expr( soption(_name="beta_y",_prefix="functions") );
  auto beta = vec( beta_x, beta_y );
  auto epsilon = expr( soption(_name="epsilon",_prefix="functions") );
  auto mu = expr( soption(_name="mu",_prefix="functions") );
  auto stable = expr( soption(_name="delta",_prefix="functions") );
  auto delta = stable*constant(1.0)/(1.0/h() + epsilon/(h()*h()));

  auto  cal = -epsilon*trace(hess(v))+ grad(v)*beta + mu*id(v);
  auto  calc = -epsilon*trace(hesst(u))+ gradt(u)*beta + mu*idt(u);

  //# marker3 #
  auto l = form1( _test=Vh );
  l = integrate(_range=elements(mesh),
		_expr=f*(id(v) + delta*cal));

  auto a = form2( _trial=Vh, _test=Vh);
  a = integrate(_range=elements(mesh),
		_expr=(gradt(u)*beta)*id(v)+epsilon*gradt(u)*trans(grad(v))+mu*idt(u)*id(v) );
  a+= integrate(_range=elements(mesh),
		_expr=delta*cal*calc );
  a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
	_expr=constant(0.) );
  a.solve(_rhs=l,_solution=u);
  //# endmarker3 #

  //# marker4 #
  auto e = exporter( _mesh=mesh );
  e->add( "u", u );
  auto f_elt = Vh->element();
  f_elt.on(_range=elements(mesh),_expr=f);
  e->add( "f", f_elt );
  e->save();
  return 0;
  //# endmarker4 # 

}
