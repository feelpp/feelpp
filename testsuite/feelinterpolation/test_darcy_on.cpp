#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_darcy_on

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;
using namespace Feel::vf;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( darcy_on_suite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    typedef Mesh<Simplex<FEELPP_DIM> > mesh_type;
    typedef FunctionSpace<mesh_type, bases<RaviartThomas<1>, Lagrange<2, Scalar, Continuous> > > space_type;

    auto mesh = loadMesh(_mesh=new mesh_type);

    auto Vh = space_type::New( _mesh=mesh );

    auto e = expr(soption("functions.e")); // error with respect to h
    BOOST_CHECK_EQUAL( e.dynamicContext(), vm::POINT|vm::DYNAMIC );
    auto g = expr(soption("functions.g")); // p exact
    BOOST_CHECK_EQUAL( g.dynamicContext(), vm::POINT|vm::DYNAMIC );

    auto h = grad<FEELPP_DIM>(g);          // u exact
    auto lhs = laplacian(g);               // f
    if( Environment::isMasterRank() )
        std::cout << g << std::endl
                  << h << std::endl
                  << lhs << std::endl;

    auto U = Vh->element();
    auto u = U.template element<0>();
    auto p = U.template element<1>();
    auto V = Vh->element();
    auto v = V.template element<0>();
    auto q = V.template element<1>();

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=-trans(idt(u))*id(v)
                  -divt(u)*id(q)
                  -idt(p)*div(v)
                  // CGLS stabilization terms
                  -1./2*(trans(idt(u)) - gradt(p))*(-id(v)+trans(grad(q)))
                  -1./2*divt(u)*div(v)
                  -1./2*trans(curlt(u))*curl(v));

    auto l = form1( _test=Vh );
    l = integrate( _range=elements(mesh), _expr= -lhs*id(q)
                   // CGLS stabilization terms
                   -1./2*lhs*div(v));

    // boundary conditions u.n
    a += on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=trans(h));

    a.solve(_rhs=l, _solution=U);

    auto err1 = normL2( _range=boundaryfaces(mesh), _expr=trans(idv(u))*N() );
    auto err2 = normL2( _range=boundaryfaces(mesh), _expr=h*N());
    auto err = normL2( _range=elements(mesh), _expr=idv(u)- trans(h) );
    if ( Environment::isMasterRank() )
    {
        auto est = e.evaluate({{"x",doption("gmsh.hsize")}})(0,0);
        BOOST_CHECK_SMALL( err, est );
        std::cout << "||u-h||_L2 =" << err << " I1=" << err1 << " I2=" << err2 << "\n";
    }

    auto UU = Vh->element();
    auto uu =UU.template element<0>();
    uu.on(_range=elements(mesh), _expr=trans(h));
    auto ex = exporter(_mesh=mesh);
    ex->add("u", u);
    ex->add("p", p);
    ex->add("h", uu);
    ex->save();
}

BOOST_AUTO_TEST_SUITE_END()
