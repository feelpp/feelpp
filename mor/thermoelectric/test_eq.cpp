#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelcrb/empiricalquadrature.hpp>

#include <algorithm>

#if defined(FEELPP_HAS_GLPK_H)
#include <glpk.h>
#endif /* FEELPP_HAS_GLPK_H */

using namespace Feel;
using namespace Feel::vf;

int main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="biotsavart-nonlinear"),
                     _desc=feel_options().add(eq_options()) );

    auto mesh = loadMesh( new Mesh<Simplex<3> > );
    int numGlobalElts = mesh->numGlobalElements();
    int numElts = mesh->numElements();
    Feel::cout << "#elts = " << numGlobalElts << " (" << numElts << ")" << std::endl;

    auto Dmu = ParameterSpace<>::New(0);
    Dmu->setDimension(2);
    auto mu_min = Dmu->element();
    mu_min << 4e3, 5;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 8e3, 20;
    Dmu->setMax( mu_max );

    auto mu = Dmu->element();
    mu << 4.3e3, 9;

    auto r = markedelements(mesh,"Cu");

    // electric problem
    auto Xh = Pch<1>(mesh, r);
    auto u = Xh->element();
    auto lhs = form2(_test=Xh,_trial=Xh);
    auto rhs = form1(_test=Xh);
    lhs = integrate(_range=r, _expr=inner(gradt(u),grad(u)) );
    lhs += on(_rhs=rhs, _element=u, _range=markedfaces(mesh,"V0"), _expr=cst(0.));
    lhs += on(_rhs=rhs, _element=u, _range=markedfaces(mesh,"V1"), _expr=cst_ref(mu(1)));
    lhs.solve(_rhs=rhs, _solution=u);
    using fn_t = std::function<void(decltype(mu) const&, decltype(Xh->element())&)>;
    fn_t lambda = [&Xh,&mesh,&r](decltype(mu) const& mu, decltype(u)& u) mutable {
                      auto lhs = form2(_test=Xh,_trial=Xh);
                      auto rhs = form1(_test=Xh);
                      lhs = integrate(_range=r, _expr=inner(gradt(u),grad(u)) );
                      lhs += on(_rhs=rhs, _element=u, _range=markedfaces(mesh,"V0"), _expr=cst(0.));
                      lhs += on(_rhs=rhs, _element=u, _range=markedfaces(mesh,"V1"), _expr=cst_ref(mu(1)));
                      lhs.solve(_rhs=rhs, _solution=u);
                  };

    // biotsavart
    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI; //SI unit : H.m-1 = m.kg.s-2.A-2
    auto dist = inner( -P(), -P(),
                       mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
    auto ex = -mu0*coeff*cst_ref(mu(0))*cross(trans(gradv(u)),P())/(dist*dist*dist);

    auto eq = EmpiricalQuadrature(r, Dmu, 2);
    eq.addExpression(ex, mu, u, lambda );

    auto ex2 = expr(soption("functions.f"))*cst_ref(mu(0));
    eq.addExpression(ex2, mu);

    tic();
    int err = eq.offline();
    toc("offline");

    mu << 5.5e3, 15;

    tic();
    auto a = integrate(_range=r, _expr=ex).evaluate()(2,0);
    Feel::cout << "a[0] = " << a << std::endl;
    toc("integrate");
    tic();
    auto a2 = integrate(_range=r, _expr=ex2).evaluate()(0,0);
    Feel::cout << "a[1] = " << a2 << std::endl;
    toc("integrate2");
    double d;
    if( !err )
    {
        for( int m = 0; m < 2; ++m )
        {
            tic();
            d = eq.evaluate(mu, m);
            toc("online");
            Feel::cout << "d[" << m << "] = " << d << std::endl;
        }
        lhs = integrate(_range=r, _expr=inner(gradt(u),grad(u)) );
        lhs += on(_rhs=rhs, _element=u, _range=markedfaces(mesh,"V0"), _expr=cst(0.));
        lhs += on(_rhs=rhs, _element=u, _range=markedfaces(mesh,"V1"), _expr=cst_ref(mu(1)));
        lhs.solve(_rhs=rhs, _solution=u);
        auto ex3 = -mu0*coeff*cst_ref(mu(0))*cross(trans(gradv(u)),P())/(dist*dist*dist);
        tic();
        auto a3 = integrate(_range=r, _expr=ex3).evaluate()(2,0);
        Feel::cout << "a[2] = " << a3 << std::endl;
        toc("integrate3");
    }
    else
        Feel::cout << "error " << err << std::endl;

    return 0;
}
