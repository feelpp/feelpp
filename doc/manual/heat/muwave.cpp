#include <feel/feel.hpp>

using namespace Feel;
	

int main(int argc, char** argv)
{
	Environment env(_argc=argc, _argv=argv,
			_about=about(_name="muwave",
				_author="Feel++ consortium",
				_email="vincent.huber@cemosis.fr"));

	static const int Dim = 2;
	static const int Order = 2;

	typedef double value_type;
	typedef typename std::complex<value_type> complex_type;

	typedef Simplex<Dim> convex_type;
	
	typedef Mesh<convex_type> mesh_type;
	typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

	typedef bases<Lagrange<Order,Scalar,Continuous>> basis_type;

	typedef FunctionSpace<mesh_type,basis_type,complex_type> space_type;
	typedef boost::shared_ptr<space_type> space_ptrtype;

	auto f = expr( soption(_name="functions.f") );
	auto g = expr( soption(_name="functions.g") );

	mesh_ptrtype mesh = loadMesh(_mesh=new mesh_type);
	space_ptrtype Vh = space_type::New(mesh);

	auto u = Vh->element("u");
	auto theta = Vh->element("theta");

	auto a_complex_expr=expr(cst(complex_type(1.,0)));
	auto xx = vf::project(Vh,elements(mesh),a_complex_expr);

	// Helmhotz
	auto f11 = form1(Vh);
	auto f12 = form2(Vh,Vh);
	// heat
	auto f21 = form1(Vh);
	auto f22 = form2(Vh,Vh);

	f12  = integrate(_range=markedelements(mesh,"oven") ,_expr=id(u)*idt(u) );
	f12 += integrate(_range=markedelements(mesh,"steak"),_expr=2.*id(u)*idt(u) );
 	f12 += integrate(_range=elements(mesh),_expr=cst(complex_type(1,-0.5))*gradt( u )*trans( grad( u ) ) );
	
	f12 += on (markedfaces(mesh,"border"),_rhs=f11,_element=u,_expr=f);
	f12 += on (markedfaces(mesh,"source"),_rhs=f11,_element=u,_expr=g);
	f12.solve(_rhs=f11,_solution=u);

#if 0
	auto ff = vf::project(Vh,element(mesh),real(u)*real(u)) + vf::project(elements(mesh),imag(u)*imag(u));
#else
	auto ff = vf::project(Vh,elements(mesh),cst(1.));
#endif
	f22  = integrate(_range=elements(mesh),_expr=grad(theta)*trans(gradt(theta)));
	f21  = integrate(_range=markedelements(mesh,"steak"),_expr=idv(ff)*id(theta));
	f22 += on(_range=markedfaces(mesh,"border"),_rhs=f21,_element=theta,_expr=cst(0.));
	f22 += on(_range=markedfaces(mesh,"source"),_rhs=f21,_element=theta,_expr=cst(0.));

	auto e = exporter(mesh);
	e->add("u",u);
	e->add("theta",theta);
	e->save();

	return 0;
}
