#include "../feel/feelmodels/hdg/mixedpoisson.hpp"

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "mixed-poisson-model" ,
                     "mixed-poisson-model" ,
                     "0.1",
                     "Mixed-Poisson-Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Romain Hild", "developer", "", "" );
    about.addAuthor( "Daniele Prada", "developer", "", "" );
    return about;
}

template<int nDim, int OrderT, int GOrder = 1>
int
runApplicationMixedPoisson( std::string  const& prefix )
{
    using namespace Feel;

    typedef FeelModels::MixedPoisson<nDim,OrderT,GOrder> mp_type;

    std::string p = "hdg.poisson";
    if( !prefix.empty() )
        p += "."+prefix;
    auto MP = mp_type::New(p);
    auto mesh = loadMesh( _mesh=new typename mp_type::mesh_type );
    decltype( IPtr( _domainSpace=Pdh<OrderT>(mesh), _imageSpace=Pdh<OrderT>(mesh) ) ) Idh ;
    decltype( IPtr( _domainSpace=Pdhv<OrderT>(mesh), _imageSpace=Pdhv<OrderT>(mesh) ) ) Idhv;
    if ( soption( "gmsh.submesh" ).empty() )
        MP -> init(mesh);
    else
    {
        Feel::cout << "Using submesh: " << soption("gmsh.submesh") << std::endl;
        auto cmesh = createSubmesh( _mesh=mesh, _range=markedelements(mesh,soption("gmsh.submesh")) );

        Idh = IPtr( _domainSpace=Pdh<OrderT>(cmesh), _imageSpace=Pdh<OrderT>(mesh) );
        Idhv = IPtr( _domainSpace=Pdhv<OrderT>(cmesh), _imageSpace=Pdhv<OrderT>(mesh) );
        MP -> init( cmesh, mesh );
    }

	// Feel::cout << "Stationary: " << MP -> isStationary() << std::endl;
	// Feel::cout << "boption steady: " << boption("ts.steady") << std::endl;

    if ( MP -> isStationary() )
    {
        if( boption("use-picard") )
        {
            int maxit = ioption("picard.maxit");
            double tol = doption("picard.tol");
            double error;
            int i = 0;
            do
            {
                MP->setMatricesAndVectorToZero();
                auto oldPotential = MP->potentialField();
                MP->assembleCstPart();
                MP->modelProperties().parameters().updateParameterValues();
                MP->updateConductivityTerm(true);
                MP->assembleRHS();
                MP->assembleRhsBoundaryCond();
                MP->solve();
                auto p = MP->potentialField();
                auto np = normL2(_range=elements(mesh),_expr=idv(p));
                error = normL2(_range=elements(mesh),_expr=idv(p)-idv(oldPotential))/np;
                Feel::cout << "error[" << i++ << "] = " << error << std::endl;
            }
            while( error > tol && i < maxit );
            MP->exportResults( mesh, Idh, Idhv );
        }
        else
        {
            MP->assembleAll();
            MP->solve();
            MP->exportResults( mesh, Idh, Idhv );
        }
    }
    else
    {
        //MP->assembleCstPart();
        for ( ; !MP->timeStepBase()->isFinished() ; MP->updateTimeStep() )
        {
            Feel::cout << "============================================================\n";
            Feel::cout << "time simulation: " << MP->time() << "s \n";
            Feel::cout << "============================================================\n";
            // MP->assembleNonCstPart();
            MP->assembleAll();
            MP->solve();
            MP->exportResults( mesh, Idh, Idhv );
        }
    }

    // MP->computeError();
    return !MP->checkResults();
}

int main(int argc, char *argv[])
{
    using namespace Feel;

    po::options_description mpoptions( "hdg.poisson options" );
    mpoptions.add( FeelModels::makeMixedPoissonOptions("","hdg.poisson") );
    mpoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ("use-picard", Feel::po::value<bool>()->default_value(false), "use picard to solve non linear problem" )
        ("picard.tol", Feel::po::value<double>()->default_value(1e-8), "picard tolerance" )
        ("picard.maxit", Feel::po::value<int>()->default_value(50), "picard maximum number of iteration")
        ;

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=mpoptions,
                           _desc_lib=FeelModels::makeMixedPoissonLibOptions("","hdg.poisson").add(feel_options())
                           );

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1>, hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2>, hana::int_c<1> ),
                                             hana::make_tuple("P3", hana::int_c<3>, hana::int_c<1> ),
                                             hana::make_tuple("P1G2", hana::int_c<1>, hana::int_c<2> ),
                                             hana::make_tuple("P2G2", hana::int_c<2>, hana::int_c<2> ),
                                             hana::make_tuple("P3G2", hana::int_c<3>, hana::int_c<2> ) );
#elif FEELPP_INSTANTIATION_ORDER_MAX >= 2
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1>, hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2>, hana::int_c<1> ),
                                             hana::make_tuple("P1G2", hana::int_c<1>, hana::int_c<2> ),
                                             hana::make_tuple("P2G2", hana::int_c<2>, hana::int_c<2> ) );
#else
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1>, hana::int_c<1> ),
                                             hana::make_tuple("P1G2", hana::int_c<1>, hana::int_c<2> ) );
#endif

    int status = 1;
    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)),
                    [&discretization,&dimension,&status]( auto const& d )
                        {
                            constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                            std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                            constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                            constexpr int _gorder = std::decay_t<decltype(hana::at_c<2>( hana::at_c<1>(d) ))>::value;
                            if ( dimension == _dim && discretization == _discretization )
                                status = runApplicationMixedPoisson<_dim,_torder,_gorder>( "" );
                        } );

    return status;
}
