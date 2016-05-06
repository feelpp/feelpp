#include <feel/feel.hpp>

int main(int argc, char**argv )
{
     using namespace Feel;

     Environment env( _argc=argc, _argv=argv,
                      _about=about(_name="qs_laplacian",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=id(v));

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=cst(0.) );
    a.solve(_rhs=l,_solution=u);

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
}