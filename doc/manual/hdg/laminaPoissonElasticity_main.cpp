#include <mixedpoissonelasticity.hpp>

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "mixed-poisson-elasticity-model" ,
                     "mixed-poisson-elasticyt-model" ,
                     "0.1",
                     "Mixed-Poisson-Elasticity-Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Romain Hild", "developer", "", "" );
    about.addAuthor( "Lorenzo Sala", "developer", "", "" );
    return about;
}




int main(int argc, char *argv[])
{
    using namespace Feel;
    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=FeelModels::makeMixedPoissonElasticityOptions(),
                           _desc_lib=FeelModels::makeMixedPoissonElasticityLibOptions("mixedpoissonelasticity").add(feel_options())
                           );


    typedef FeelModels::MixedPoissonElasticity<FEELPP_DIM,FEELPP_ORDER> mpe_type;
    


    auto mesh = loadMesh( _mesh=new mpe_type::mesh_type );

	mpe_type MPE( mesh );

	/*
    decltype( IPtr( _domainSpace=Pdh<FEELPP_ORDER>(mesh), _imageSpace=Pdh<1>(mesh) ) ) Idh_poi ;
    decltype( IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(mesh), _imageSpace=Pdhv<1>(mesh) ) ) Idhv_poi;

    decltype( IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(mesh), _imageSpace=Pdhv<1>(mesh) ) ) Idh_el ;
    decltype( IPtr( _domainSpace=Pdhms<FEELPP_ORDER>(mesh), _imageSpace=Pdhms<1>(mesh) ) ) Idhv_el;
	*/

	Feel::cout << __LINE__ << std::endl;

	MPE.run();

	Feel::cout << __LINE__ << std::endl;
	/*
	// solve the elasticity stationary
	Feel::cout << "STARTING SOLVE THE ELASTICITY PART . . ." << std::endl;
	ME->init( mesh );
	ME->solve();
	ME->exportResults( mesh );
	Feel::cout << "ELASTICITY PART SOLVED" << std::endl;

	// solve the poisson stationary
	Feel::cout << "STARTING SOLVE THE DARCY PART . . ." << std::endl;
	MP->init( mesh );
	MP->solve();
	MP->exportResults( mesh );
	Feel::cout << "DARCY PART SOLVED" << std::endl;
	

	
    // if ( soption( "mixedpoisson.gmsh.submesh" ).empty() )
        MP -> init(mesh);
    // else
    // {
    //     auto cmesh = createSubmesh( mesh, markedelements(mesh,soption("mixedpoisson.gmsh.submesh")), Environment::worldComm() );
    //     Idh = IPtr( _domainSpace=Pdh<FEELPP_ORDER>(cmesh), _imageSpace=Pdh<1>(mesh) );
    //     Idhv = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmesh), _imageSpace=Pdhv<1>(mesh) );
    //     MP -> init( cmesh, 0, 0, mesh );
    // }
    
    if ( MP -> isStationary() )
    {
        MP->solve();
        MP->exportResults( mesh );
    }
    // else
    // {
    //     for ( ; !MP->timeStepBase()->isFinished() ; MP->updateTimeStep() )
    //     {
    //         Feel::cout << "============================================================\n";
    //         Feel::cout << "time simulation: " << MP->time() << "s \n";
    //         Feel::cout << "============================================================\n";
    //         MP->solve();
    //         MP->exportResults( mesh, Idh, Idhv );
    //     }
    // }

 
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
		Idh = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmesh), _imageSpace=Pdhv<1>(mesh) );
    	Idhv = IPtr( _domainSpace=Pdhms<FEELPP_ORDER>(cmesh), _imageSpace=Pdhms<1>(mesh) );
    	ME -> init( cmesh, mesh );
	}
 
    if ( ME -> isStationary() )
    {
        ME->solve();
        ME->exportResults( mesh );
    }
    else
    {
    	for ( ; !ME->timeStepBase()->isFinished() ; ME->updateTimeStep() )
        {
        	Feel::cout << "============================================================\n";
        	Feel::cout << "time simulation: " << ME->time() << "s \n";
        	Feel::cout << "============================================================\n";
        	ME->solve();
        	ME->exportResults( mesh , Idh, Idhv );
        }
     }
	*/
   
	
	return 0;
}
