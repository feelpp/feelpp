#include "../feel/feelmodels/hdg/mixedelasticity.hpp"

using namespace Feel;


inline
AboutData
makeAbout()
{
    AboutData about( "mixed-elasticity-model" ,
                     "mixed-elasticity-model" ,
                     "0.1",
                     "Mixed-Elasticity-Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Romain Hild", "developer", "", "" );
    about.addAuthor( "Lorenzo Sala", "developer", "", "" );
    return about;
}

template<int nDim, int OrderT>
void
runApplicationMixedElasticity()
{
    using namespace Feel;

    typedef FeelModels::MixedElasticity<nDim,OrderT> me_type;

    auto ME = me_type::New("mixedelasticity");
    auto mesh = loadMesh( _mesh=new typename me_type::mesh_type );
    
    decltype( IPtr( _domainSpace=Pdhv<OrderT>(mesh), _imageSpace=Pdhv<OrderT>(mesh) ) ) Idh ;
    decltype( IPtr( _domainSpace=Pdhms<OrderT>(mesh), _imageSpace=Pdhms<OrderT>(mesh) ) ) Idhv;

    std::list<std::string> listSubmesh;

    std::string unseparatedList = soption( "mixedelasticity.gmsh.submesh");

    char help;
    std::string nameSubmesh;
    for ( int i = 0; i < unseparatedList.size(); i++)
    {

        help = unseparatedList[i];
        if ( help == ',' || i == unseparatedList.size()-1 )
        {        
            if ( i ==  unseparatedList.size()-1)
                nameSubmesh.push_back(help);
            listSubmesh.push_back(nameSubmesh);
            nameSubmesh.erase();
        }
        else
        {
            nameSubmesh.push_back(help);
        }
    }

    if ( listSubmesh.empty() )
        ME -> init(mesh);
    else
    {
        Feel::cout << "Using submesh: " << listSubmesh << std::endl;
        auto cmesh = createSubmesh( mesh, markedelements( mesh, listSubmesh ), Environment::worldComm() );
        Idh = IPtr( _domainSpace=Pdhv<OrderT>(cmesh), _imageSpace=Pdhv<OrderT>(mesh) );
        Idhv = IPtr( _domainSpace=Pdhms<OrderT>(cmesh), _imageSpace=Pdhms<OrderT>(mesh) );
        ME -> init( cmesh, mesh );
    }



/*	
	if ( soption("gmsh.submesh").empty() )
	{
		ME -> init(mesh);
	} 
	else
	{
		Feel::cout << "Using submesh: " << soption("gmsh.submesh") << std::endl;
		std::list<std::string> listSubmeshes;
		listSubmeshes.push_back( soption("gmsh.submesh") );
		if ( !soption("gmsh.submesh2").empty() )
		{
			Feel::cout << "Using submesh 2: " << soption("gmsh.submesh2") << std::endl;
			listSubmeshes.push_back( soption("gmsh.submesh2") );
		}
		auto cmesh = createSubmesh( mesh, markedelements(mesh,listSubmeshes), Environment::worldComm() );
		Idh = IPtr( _domainSpace=Pdhv<OrderT>(cmesh), _imageSpace=Pdhv<OrderT>(mesh) );
    	Idhv = IPtr( _domainSpace=Pdhms<OrderT>(cmesh), _imageSpace=Pdhms<OrderT>(mesh) );
    	ME -> init( cmesh, mesh );
	}
*/	 
    
	if ( ME -> isStationary() )
    {
		ME->assembleCst();
		ME->assembleNonCst();
        ME->solve();
        ME->exportResults( mesh );
        ME->exportTimers(); 
    }
    else
    {
		//ME->assembleCst();
    	for ( ; !ME->timeStepBase()->isFinished() ; ME->updateTimeStep() )
        {
        	Feel::cout << "============================================================\n";
        	Feel::cout << "time simulation: " << ME->time() << "s \n";
        	Feel::cout << "============================================================\n";
			ME->assembleCst();
			ME->assembleNonCst();
        	ME->solve();
        	ME->exportResults( mesh , Idh, Idhv );
        }
     }

#if 0

    ME->geometricTest();

#endif
}

int main(int argc, char *argv[])
{
    using namespace Feel;

    po::options_description meoptions( "mixedelasticity options" );
    meoptions.add( FeelModels::makeMixedElasticityOptions("mixedelasticity") );
    meoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ;

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=meoptions,
                           _desc_lib=FeelModels::makeMixedElasticityLibOptions("mixedelasticity").add(feel_options())
                           );


    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ),
                                             hana::make_tuple("P3", hana::int_c<3> ) );
#elif FEELPP_INSTANTIATION_ORDER_MAX >= 2
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ) );
#else
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ) );
#endif

    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension]( auto const& d )
                    {
                        constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                        std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                        constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                        if ( dimension == _dim && discretization == _discretization )
                            runApplicationMixedElasticity<_dim,_torder>();
                    } );


    return 0;
}



