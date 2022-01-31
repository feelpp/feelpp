#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelcrb/empiricalquadrature.hpp>

#include <algorithm>
#include <boost/math/special_functions/binomial.hpp>

using namespace Feel;
using namespace Feel::vf;

using parameter_type = typename ParameterSpace<>::element_type;

std::string alphaRef()
{
    using boost::math::binomial_coefficient;

    double zb = -0.08651;
    double zt = 0.08651;
    auto paramNames = std::vector<std::string>({"p1","p2"});

    int n = paramNames.size() + 1;

    std::stringstream alphaStream;
    alphaStream << "0 + (z<" << zt << ")*(z>" << zb << ")*(";
    for( int k = 1; k < n; ++k )
    {
        double ckn = binomial_coefficient<double>(n,k);
        alphaStream << " + " << ckn
                    << "*((z-" << zb << ")/" << zt-zb << ")^" << k
                    << "*(1-(z-" << zb << ")/" << zt-zb << ")^" << n-k
                    << "*" << paramNames[k-1];
    }
    alphaStream << "):x:y:z";
    for( auto const& p : paramNames )
        alphaStream << ":" << p;

    return alphaStream.str();
}

std::string alphaPrimeRef()
{
    using boost::math::binomial_coefficient;

    double zb = -0.08651;
    double zt = 0.08651;
    auto paramNames = std::vector<std::string>({"p1","p2"});
    paramNames.insert(paramNames.begin(), "0");
    paramNames.push_back("0");

    int n = paramNames.size() - 1;

    std::stringstream alphaPrimeStream;
    alphaPrimeStream << "0 + (z<" << zt << ")*(z>" << zb << ")*(";
    alphaPrimeStream << n << "/" << zt-zb << "*(";
    for( int k = 0; k < n; ++k )
    {
        double ckn1 = binomial_coefficient<double>(n-1,k);
        alphaPrimeStream << " + (" << paramNames[k+1] << "-" << paramNames[k] << ")*" << ckn1
                         << "*((z-" << zb << ")/" << zt-zb << ")^" << k
                         << "*(1-(z-" << zb << ")/" << zt-zb << ")^" << n-1-k;
    }
    alphaPrimeStream << ")):x:y:z";
    for( int k = 1; k < n; ++k )
        alphaPrimeStream << ":" << paramNames[k];

    return alphaPrimeStream.str();
}

int main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="test_eq"),
                     _desc=feel_options().add(eq_options()) );

    auto Dmu = ParameterSpace<>::New(0);
    Dmu->setDimension(2);
    auto mu_min = Dmu->element();
    mu_min << -4, 4;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 4, 10;
    Dmu->setMax( mu_max );
    Dmu->setParameterName(0, "p1");
    Dmu->setParameterName(1, "p2");

    auto mu = Dmu->element();

    auto mesh = loadMesh( new Mesh<Simplex<3> > );
    auto r = markedelements(mesh,"Cu");
    // auto r = elements(mesh);

    auto eq = EmpiricalQuadrature(r,Dmu);

    auto Xh = Pch<1>(mesh, r);
    auto u = Xh->element();
    auto v = Xh->element();

    auto paramNames = std::vector<std::string>({"p1","p2"});
    std::vector<Expr<Feel::vf::Cst<boost::reference_wrapper<double> > >> paramRefs;
    for(int i = 0; i < paramNames.size(); ++i )
        paramRefs.push_back(cst_ref(mu.parameterNamed(paramNames[i])));
    auto alphaStr = alphaRef();
    auto alphaExpr = expr(alphaStr, paramNames, paramRefs);
    auto alphaPrimeStr = alphaPrimeRef();
    auto alphaPrimeExpr = expr(alphaPrimeStr, paramNames, paramRefs);

    auto Jinv = mat<3,3>( cos(alphaExpr), sin(alphaExpr), alphaPrimeExpr*(cos(alphaExpr)*Py()-sin(alphaExpr)*Px()),
                          -sin(alphaExpr), cos(alphaExpr), -alphaPrimeExpr*(cos(alphaExpr)*Px()+sin(alphaExpr)*Py()),
                          cst(0.), cst(0.), cst(1.) );

    auto psi = vec( cos(alphaExpr)*Px() + sin(alphaExpr)*Py(),
                    -sin(alphaExpr)*Px() + cos(alphaExpr)*Py(),
                    Pz() );

    using fn_t = std::function<bool(decltype(mu) const&, decltype(Xh->element())&)>;
    fn_t lambda = [&Xh,&mesh,&r](decltype(mu) const& mu, decltype(u)& v) mutable -> bool {
                      auto lhs = form2(_test=Xh,_trial=Xh);
                      auto rhs = form1(_test=Xh);
                      auto paramNames = std::vector<std::string>({"p1","p2"});
                      std::vector<Expr<Feel::vf::Cst<boost::reference_wrapper<const double> > >> paramRefs;
                      for(int i = 0; i < paramNames.size(); ++i )
                          paramRefs.push_back(cst_ref(mu.parameterNamed(paramNames[i])));
                      auto alphaStr = alphaRef();
                      auto alphaExpr = expr(alphaStr, paramNames, paramRefs);
                      auto alphaPrimeStr = alphaPrimeRef();
                      auto alphaPrimeExpr = expr(alphaPrimeStr, paramNames, paramRefs);
                      auto Jinv = mat<3,3>( cos(alphaExpr), sin(alphaExpr), alphaPrimeExpr*(cos(alphaExpr)*Py()-sin(alphaExpr)*Px()),
                                            -sin(alphaExpr), cos(alphaExpr), -alphaPrimeExpr*(cos(alphaExpr)*Px()+sin(alphaExpr)*Py()),
                                            cst(0.), cst(0.), cst(1.) );
                      lhs = integrate(_range=r, _expr=inner(gradt(v)*Jinv,grad(v)*Jinv) );
                      lhs += on(_rhs=rhs, _element=v, _range=markedfaces(mesh,"V0"), _expr=cst(0));
                      lhs += on(_rhs=rhs, _element=v, _range=markedfaces(mesh,"V1"), _expr=cst(9.));
                      auto r = lhs.solve(_rhs=rhs, _solution=v, _rebuild=true);
                      if( !r.isConverged() )
                          Feel::cout << "solve for mu = " << mu.toString()
                                     << " failed to converge" << std::endl;
                      return r.isConverged();
                  };

    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-7; //SI unit : H.m-1 = m.kg.s-2.A-2
    auto mu1 = 58e6;
    auto p1 = vec(cst(0.),cst(0.),cst(0.));
    auto dist = inner( p1-psi, p1-psi,
                       mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
    auto ex1 = -mu0*coeff*mu1*cross(trans(gradv(v)*Jinv),p1-psi)/(dist*dist*dist);
    int comp = 2;
    // auto ex1 = idv(v);
    // int comp = 0;
    eq.addExpression(ex1, mu, v, lambda, comp );

    // auto ex2 = idv(v);
    // eq.addExpression(ex2, mu, v, lambda );

    int err = eq.offline();

    if( !err )
    {
        auto trainset = Dmu->sampling();
        trainset->randomize(int(doption("parameters.a")));
        trainset->addElement(eq.sample(0));
        for(auto const& m : *trainset )
        {
            // mu << -2.5, 1;
            // mu = eq.sample(0);
            mu = m;
            Feel::cout << "mu = " << mu.toString() << std::endl;
            tic();
            auto d1 = eq.evaluate(mu, 0);
            Feel::cout << std::setprecision(10) << "d1 = " << d1 << std::endl;
            toc("evaluate 1");
            tic();
            auto c1 = eq.evaluateOffline(mu, 0);
            Feel::cout << std::setprecision(10) << "c1 = " << c1 << std::endl;
            toc("evaluateOff 1");
            tic();
            auto a1 = integrate(_range=r, _expr=ex1).evaluate()(comp,0);
            Feel::cout << std::setprecision(10) << "a1 = " << a1 << std::endl;
            toc("integrate 1");
            // tic();
            // auto d2 = eq.evaluate(mu, 1);
            // Feel::cout << std::setprecision(10) << "d2 = " << d2 << std::endl;
            // toc("evaluate 2");
            // tic();
            // auto c2 = eq.evaluateOffline(mu, 1);
            // Feel::cout << std::setprecision(10) << "c2 = " << c2 << std::endl;
            // toc("evaluateOff 2");
            // tic();
            // auto a2 = integrate(_range=r, _expr=ex2).evaluate()(0,0);
            // Feel::cout << std::setprecision(10) << "a2 = " << a2 << std::endl;
            // toc("integrate 2");
        }
    }

    return 0;
}
