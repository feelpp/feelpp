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

std::string alphaRef(int nb)
{
    using boost::math::binomial_coefficient;

    double zb = 0.5;//-0.08651;
    double zt = 2.5;//0.08651;
    auto paramNames = std::vector<std::string>();
    for( int i = 1; i <= nb; ++i )
        paramNames.push_back("p"+std::to_string(i));

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

std::string alphaPrimeRef(int nb)
{
    using boost::math::binomial_coefficient;

    double zb = 0.5;//-0.08651;
    double zt = 2.5;//0.08651;
    auto paramNames = std::vector<std::string>();
    for( int i = 1; i <= nb; ++i )
        paramNames.push_back("p"+std::to_string(i));
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
    // auto r = markedelements(mesh,"Cu");
    auto r = elements(mesh);

    auto eq = EmpiricalQuadrature(r,Dmu);

    auto Xh = Pch<1>(mesh, r);
    auto u = Xh->element();
    auto v = Xh->element();

    int nbParams = doption("parameters.b");
    auto paramNames = std::vector<std::string>();
    for( int i = 1; i <= nbParams; ++i )
        paramNames.push_back("p"+std::to_string(i));
    std::vector<Expr<Feel::vf::Cst<boost::reference_wrapper<double> > >> paramRefs;
    for(int i = 0; i < paramNames.size(); ++i )
        paramRefs.push_back(cst_ref(mu.parameterNamed(paramNames[i])));
    auto alphaStr = alphaRef(nbParams);
    auto alphaExpr = expr(alphaStr, paramNames, paramRefs);
    auto alphaPrimeStr = alphaPrimeRef(nbParams);
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
                      int nbParams = doption("parameters.b");
                      auto paramNames = std::vector<std::string>();
                      for( int i = 1; i <= nbParams; ++i )
                          paramNames.push_back("p"+std::to_string(i));
                      std::vector<Expr<Feel::vf::Cst<boost::reference_wrapper<const double> > >> paramRefs;
                      for(int i = 0; i < paramNames.size(); ++i )
                          paramRefs.push_back(cst_ref(mu.parameterNamed(paramNames[i])));
                      auto alphaStr = alphaRef(nbParams);
                      auto alphaExpr = expr(alphaStr, paramNames, paramRefs);
                      auto alphaPrimeStr = alphaPrimeRef(nbParams);
                      auto alphaPrimeExpr = expr(alphaPrimeStr, paramNames, paramRefs);
                      auto Jinv = mat<3,3>( cos(alphaExpr), sin(alphaExpr), alphaPrimeExpr*(cos(alphaExpr)*Py()-sin(alphaExpr)*Px()),
                                            -sin(alphaExpr), cos(alphaExpr), -alphaPrimeExpr*(cos(alphaExpr)*Px()+sin(alphaExpr)*Py()),
                                            cst(0.), cst(0.), cst(1.) );
                      lhs = integrate(_range=r, _expr=inner(gradt(v)*Jinv,grad(v)*Jinv) );
                      lhs += on(_rhs=rhs, _element=v, _range=markedfaces(mesh,"base"), _expr=cst(0));
                      lhs += on(_rhs=rhs, _element=v, _range=markedfaces(mesh,"top"), _expr=cst(9.));
                      auto r = lhs.solve(_rhs=rhs, _solution=v, _rebuild=true);
                      if( !r.isConverged() )
                          Feel::cout << "solve for mu = " << mu.toString()
                                     << " failed to converge" << std::endl;
                      return r.isConverged();
                  };

    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-7; //SI unit : H.m-1 = m.kg.s-2.A-2
    auto mu1 = 58e6;
    auto p1 = vec(cst(0.),cst(0.),cst(1.5));
    auto dist = inner( p1-psi, p1-psi,
                       mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
    auto ex1 = -mu0*coeff*mu1*cross(trans(gradv(v)*Jinv),p1-psi)/(dist*dist*dist);
    int comp = 2;
    // auto ex1 = idv(v);
    // int comp = 0;
    eq.addExpression(ex1, mu, v, lambda, comp );

    // auto p2 = vec(cst(0.1),cst(0.1),cst(1.5));
    // auto dist2 = inner( p2-psi, p2-psi,
    //                     mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
    // auto ex2 = -mu0*coeff*mu1*cross(trans(gradv(v)*Jinv),p2-psi)/(dist2*dist2*dist2);
    // eq.addExpression(ex2, mu, v, lambda, comp );

    // auto p3 = vec(cst(-0.15),cst(0.03),cst(1.5));
    // auto dist3 = inner( p3-psi, p3-psi,
    //                     mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
    // auto ex3 = -mu0*coeff*mu1*cross(trans(gradv(v)*Jinv),p3-psi)/(dist3*dist3*dist3);
    // eq.addExpression(ex3, mu, v, lambda, comp );

    // auto p4 = vec(cst(-0.15),cst(-0.15),cst(1.5));
    // auto dist4 = inner( p4-psi, p4-psi,
    //                     mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
    // auto ex4 = -mu0*coeff*mu1*cross(trans(gradv(v)*Jinv),p4-psi)/(dist4*dist4*dist4);
    // eq.addExpression(ex4, mu, v, lambda, comp );

    // auto ex2 = idv(v);
    // eq.addExpression(ex2, mu, v, lambda );

    // auto ex1 = cst_ref(mu(0))*cst_ref(mu(1))*Pz();
    // std::string ex1str = soption("functions.e");
    // std::string ex1str = "(z>0.5)*(z<2.5)*(3*((z-0.5)/2)^1*(1-(z-0.5)/2)^2*p1+3*((z-0.5)/2)^2*(1-(z-0.5)/2)^1*p2):z:p1:p2";
    // Feel::cout << ex1str << std::endl;
    // Feel::cout << alphaStr << std::endl;
    // auto ex1 = expr(ex1str, paramNames, paramRefs);
    // auto ex1 = cst(2.)*(Pz()-0.5)/2*(1-(Pz()-0.5)/2)*cst_ref(mu.parameterNamed("p1"));
    // int comp = 0;
    // eq.addExpression(alphaExpr, mu, comp);

    // eq.addExpression(alphaPrimeExpr, mu, comp);

    std::tuple expressions = std::make_tuple(ex1);//,ex2,ex3,ex4);
    auto rangeExpr = hana::make_range(hana::int_c<0>, hana::int_c<1>);//4>);

    int err = eq.offline();

    if( !err )
    {
        int size = int(doption("parameters.a"));
        std::vector<std::vector<double> > errs(eq.nbExpressions());
        auto trainset = Dmu->sampling();
        trainset->randomize(size);
        // trainset->addElement(eq.sample(0));
        for(auto const& m : *trainset )
        {
            mu = m;
            hana::for_each(rangeExpr,
                           [&](auto i) {
            // for(int i = 0; i < eq.nbExpressions(); ++i )
            // {
                // Feel::cout << "mu = " << mu.toString() << std::endl;
                auto d1 = eq.evaluate(mu, i);
                // Feel::cout << std::setprecision(10) << "d1 = " << d1 << std::endl;
                // tic();
                // auto c1 = eq.evaluateOffline(mu, i);
                // Feel::cout << std::setprecision(10) << "c1 = " << c1 << std::endl;
                // toc("evaluateOff 1");
                tic();
                auto a1 = integrate(_range=r, _expr=std::get<i>(expressions)).evaluate()(comp,0);
                Feel::cout << std::setprecision(10) << "a1 = " << a1 << std::endl;
                toc("integrate", eq.verbosity() > 0);
                errs[i].push_back(abs(d1-a1));
            // }
                           });
        }
        std::vector<double> min(eq.nbExpressions()), max(eq.nbExpressions()), mean(eq.nbExpressions());
        for( int i = 0; i < eq.nbExpressions(); ++i )
        {
            min[i] = *std::min_element(errs[i].begin(), errs[i].end());
            max[i] = *std::max_element(errs[i].begin(), errs[i].end());
            double s = std::accumulate(errs[i].begin(), errs[i].end(), 0.0);
            mean[i] = s/size;
        }
        fs::ofstream cvgStat( "stat.dat" );
        if( cvgStat && Environment::isMasterRank() )
        {
            cvgStat << std::setw(5) << "M" << std::setw(25) << "min" << std::setw(25) << "max"
                    << std::setw(25) << "mean" << std::endl;
            for(int n = 0; n < eq.nbExpressions(); ++n)
                cvgStat << std::setw(5) << n+1 << std::setw(25) << min[n] << std::setw(25) << max[n]
                        << std::setw(25) << mean[n] << std::endl;
            cvgStat.close();
        }
        Feel::cout << std::setw(5) << "M" << std::setw(25) << "min" << std::setw(25) << "max"
                   << std::setw(25) << "mean" << std::endl;
        for(int n = 0; n < eq.nbExpressions(); ++n)
            Feel::cout << std::setw(5) << n+1 << std::setw(25) << min[n] << std::setw(25) << max[n]
                       << std::setw(25) << mean[n] << std::endl;
    }

    std::ofstream os ( "timers.md" );
    Environment::saveTimersMD(os);
    os.close();

    return 0;
}
