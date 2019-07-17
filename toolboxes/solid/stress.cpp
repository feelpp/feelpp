#include <feel/feelmodels/solid/solidmechanics.hpp>

namespace Feel
{

template <uint16_type OrderDisp>
void
runApplicationSolid()
{
    typedef FeelModels::SolidMechanics< Simplex<FEELPP_DIM,1>,
                                        Lagrange<OrderDisp, Vectorial,Continuous,PointSetFekete> > model_type;
    auto SM = model_type::New("solid");

    SM->init();
    SM->printAndSaveInfo();
    SM->solve();
    SM->updateStressCriterions();
    auto u = SM->fieldDisplacement();
    auto s = SM->fieldsPrincipalStresses();
    auto vm = SM->fieldVonMisesCriterions();
    auto tr = SM->fieldTrescaCriterions();

    auto mesh = SM->mesh();

#if defined(ELASTICITY_ANA)
    auto mp = SM->mechanicalProperties();

    // analytical solution
    double pc = mp->cstCoeffPoisson( "Omega" );
    double yc = mp->cstYoungModulus( "Omega" );
    double r_int = 1.;
    double r_ext = 2.;
    double bz1 = 10;
    double bz2 = 30;
    double jt = 5000;

    auto r = sqrt( Px()*Px() + Py()*Py() );
    auto cos_theta = Px()/r;
    auto sin_theta = Py()/r;

    auto bterm_1 = (bz1-bz2)/((r_ext/r_int) - 1);
    auto bterm_2 = bz1 + bterm_1;

    // compute constant C1 and C2
    auto factor_C1 = ( ((1+pc)*(1-2*pc))/(yc*(1-pc)) )*jt*r_int;
    auto coeff_C1 = pow(r_ext,2)/( pow(r_ext,2) - pow(r_int,2)  );
    auto term1_C1 = bterm_2*((2 - pc)/3)*( 1 - coeff_C1*( 1 - (r_ext/r_int) ) );
    auto term2_C1 = bterm_1*((2*pc - 3)/8)*( 1 - coeff_C1*( 1 - pow(r_ext/r_int,2) ) );
    auto C1 = factor_C1*(term1_C1 + term2_C1);

    auto factor_C2 = -(  (pow(r_int,3)*pow(r_ext,2)*jt*(1+pc))/((pow(r_ext,2)-pow(r_int,2))*yc*(1-pc))  );
    auto term1_C2 = bterm_2*((2 - pc)/3)*(1 - (r_ext/r_int));
    auto term2_C2 = bterm_1*((2*pc - 3)/8)*(1 - pow(r_ext/r_int,2));
    auto C2 = factor_C2*(term1_C2 + term2_C2);

    // compute particular solution (up)
    auto factor_up = factor_C1;
    auto term1_up = -(r_int/3)*bterm_2*pow(r/r_int, 2);
    auto term2_up = (r_int/8)*bterm_1*pow(r/r_int, 3);
    auto up = factor_up*(term1_up + term2_up);

    auto u_cyl = C1*r + C2/r + up;

    auto Xh_vec = SM->functionSpaceDisplacement();
    auto anal = vf::project( Xh_vec, elements(mesh), cos_theta*u_cyl*oneX() + sin_theta*u_cyl*oneY() + cst(0.)*oneZ() );
#endif

    auto e = exporter( _mesh=mesh, _name="solidMec" );
    e->add( "displacement", u );
#if defined(ELASTICITY_ANA)
    e->add( "analytic", anal );
#endif
    e->add( "vonMises", vm );
    e->add( "tresca", tr );
    for( auto i : range(FEELPP_DIM))
        e->add( (boost::format( "principal-stress-%1%") %i).str(), *(s[i]) );
    e->save();
}

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;
	po::options_description solidmecoptions( "application solid-mechanics options" );
    solidmecoptions.add( toolboxes_options("solid") );
    solidmecoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P1,P2 ")
        ("solve-quasi-static", Feel::po::value<bool>()->default_value(false), "solve-quasi-static ")
        ("solve-quasi-static.variable-initial", Feel::po::value<double>()->default_value(0.0), "solve-quasi-static.variable-initial")
        ("solve-quasi-static.variable-step", Feel::po::value<double>()->default_value(0.1), "solve-quasi-static.variable-step")
        ("solve-quasi-static.variable-symbol", Feel::po::value<std::string>()->default_value(""), "solve-quasi-static.variable-symbol")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=solidmecoptions,
                     _about=about(
#if defined(ELASTICITY_ANA)
                         _name="application_solenoid",
#else
                         _name="application_solid_stress",
#endif
                         _author="Feel++ Consortium",
                         _email="feelpp-devel@feelpp.org"));

    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P1" )
        runApplicationSolid<1>();
    else if ( feapprox == "P2" )
        runApplicationSolid<2>();
    else CHECK( false ) << "invalid feapprox " << feapprox;


    return 0;
}
