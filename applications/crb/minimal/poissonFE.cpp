#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

int main(int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv );

    auto mesh = loadMesh( new Mesh<Simplex<2> > );
    auto Xh = Pch<1>(mesh);
    auto u = Xh->element();
    auto v = Xh->element();
    auto kappa = doption("parameters.kappa");
    auto gamma = doption("parameters.gamma");
    auto flux = doption("parameters.f");

    auto a = form2(_test=Xh, _trial=Xh);
    a = integrate( markedelements(mesh, "omega1"),
                   inner(gradt(u),grad(v)) );
    a+= integrate( markedelements(mesh, "omega0"),
                   kappa*inner(gradt(u),grad(v)) );

    auto f = form1(_test=Xh);
    f = integrate( markedfaces(mesh, "base"),
                   flux*id(v) );

    a+= on( _range=markedfaces(mesh, "top"), _element=u, _rhs=f, _expr=cst(0.) );

    a.solve( _rhs=f, _solution=u );

    auto e = exporter(_mesh=mesh);
    e->add("u", u);
    e->save();

    return 0;
}
