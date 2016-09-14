#include <mixedelasticity.hpp>

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


int main(int argc, char *argv[])
{

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=FeelModels::makeMixedElasticityOptions("mixedelasticity"),
                           _desc_lib=FeelModels::makeMixedElasticityLibOptions("mixedelasticity").add(feel_options())
                           );


    typedef FeelModels::MixedElasticity<FEELPP_DIM,FEELPP_ORDER> me_type;

    auto ME = me_type::New("mixedelasticity");
    auto mesh = loadMesh( _mesh=new me_type::mesh_type );
    
    decltype( IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(mesh), _imageSpace=Pdhv<1>(mesh) ) ) Idh ;
    decltype( IPtr( _domainSpace=Pdhms<FEELPP_ORDER>(mesh), _imageSpace=Pdhms<1>(mesh) ) ) Idhv;
    ME -> init (mesh);
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
		Idh = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmesh), _imageSpace=Pdhv<1>(mesh) );
    	Idhv = IPtr( _domainSpace=Pdhms<FEELPP_ORDER>(cmesh), _imageSpace=Pdhms<1>(mesh) );
    	ME -> init( cmesh, mesh );
	}
	*/ 
    
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


    return 0;
}



