//.......
using namespace Feel;

int
main( int argc, char** argv )
{
  Environment env( _argc=argc, _argv=argv,
		   _desc=opts,
		   _about=about(_name="myadvection",
				_author="kyoshe winstone",
				_email="wkyoshe@gmail.com") );

  auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
  auto Vh = Pch<3>( mesh );
  auto u = Vh->element();
  auto v = Vh->element();
  //# endmarker2 #

  auto vars = symbols<2>();
  auto f = expr( option(_name="f",_prefix="functions").as<std::string>(), vars );
  auto beta_x = expr( option(_name="beta_x",_prefix="functions").as<std::string>(), vars );
  auto beta_y = expr( option(_name="beta_y",_prefix="functions").as<std::string>(), vars );
  auto beta = vec( beta_x, beta_y );
  auto epsilon = expr( option(_name="epsilon",_prefix="functions").as<std::string>(), vars );
  auto mu = expr( option(_name="mu",_prefix="functions").as<std::string>(), vars );
  auto stable = expr( option(_name="delta",_prefix="functions").as<std::string>(), vars );
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
  e->save();
  return 0;
  //# endmarker4 #
}
