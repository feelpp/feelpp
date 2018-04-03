#include "mixedpoissonelasticity.hpp"

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

	decltype( IPtr( _domainSpace=Pdh<FEELPP_ORDER>(mesh), _imageSpace=Pdh<FEELPP_ORDER>(mesh) ) ) Idh_poi ;
    decltype( IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(mesh), _imageSpace=Pdhv<FEELPP_ORDER>(mesh) ) ) Idhv_poi;

    decltype( IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(mesh), _imageSpace=Pdhv<FEELPP_ORDER>(mesh) ) ) Idh_el ;
    decltype( IPtr( _domainSpace=Pdhms<FEELPP_ORDER>(mesh), _imageSpace=Pdhms<FEELPP_ORDER>(mesh) ) ) Idhv_el;

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

		
	if ( soption( "gmsh.submesh" ).empty() && soption("mixedelasticity.gmsh.submesh").empty() )
	{
		mpe_type MPE( mesh, mesh, mesh, mesh);
		MPE.run( Idh_el, Idhv_el, Idh_poi, Idhv_poi );
	}
	else if ( soption( "gmsh.submesh" ).empty() )
    {
		Feel::cout << "Using submesh for Elasticity: " << soption("mixedelasticity.gmsh.submesh") << std::endl;
        auto cmeshElasticity = createSubmesh( mesh, markedelements(mesh,listSubmesh ), Environment::worldComm() );
    	
        Idh_el = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmeshElasticity), _imageSpace=Pdhv<FEELPP_ORDER>(mesh) );
		Idhv_el = IPtr( _domainSpace=Pdhms<FEELPP_ORDER>(cmeshElasticity), _imageSpace=Pdhms<FEELPP_ORDER>(mesh) );
		
		mpe_type MPE( mesh, cmeshElasticity, cmeshElasticity, mesh ); 
		MPE.run( Idh_el, Idhv_el, Idh_poi, Idhv_poi );

    }
	else if ( soption("mixedelasticity.gmsh.submesh").empty() )
    {
		Feel::cout << "Using submesh for Poisson: " << soption("gmsh.submesh") << std::endl;
        auto cmeshPoisson = createSubmesh( mesh, markedelements(mesh,soption("gmsh.submesh")), Environment::worldComm() );

		Idh_poi = IPtr( _domainSpace=Pdh<FEELPP_ORDER>(cmeshPoisson), _imageSpace=Pdh<FEELPP_ORDER>(mesh) );
    	Idhv_poi = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmeshPoisson), _imageSpace=Pdhv<FEELPP_ORDER>(mesh) );
		
        mpe_type MPE( cmeshPoisson, mesh, cmeshPoisson, mesh ); 
		MPE.run( Idh_el, Idhv_el, Idh_poi, Idhv_poi );
	} 
    else
    {
		Feel::cout << "Using submesh for Poisson: " << soption("gmsh.submesh") << std::endl;
		Feel::cout << "Using submesh for Elasticity: " << soption("mixedelasticity.gmsh.submesh") << std::endl;
	    
        auto cmeshPoisson = createSubmesh( mesh, markedelements(mesh,soption("gmsh.submesh")), Environment::worldComm() );
        auto cmeshElasticity = createSubmesh( mesh, markedelements(mesh,listSubmesh ), Environment::worldComm() );

		Idh_poi = IPtr( _domainSpace=Pdh<FEELPP_ORDER>(cmeshPoisson), _imageSpace=Pdh<FEELPP_ORDER>(mesh) );
    	Idhv_poi = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmeshPoisson), _imageSpace=Pdhv<FEELPP_ORDER>(mesh) );

    	Idh_el = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmeshElasticity), _imageSpace=Pdhv<FEELPP_ORDER>(mesh) );
		Idhv_el = IPtr( _domainSpace=Pdhms<FEELPP_ORDER>(cmeshElasticity), _imageSpace=Pdhms<FEELPP_ORDER>(mesh) );
		
        auto meshCommon = (soption("gmsh.submesh")<soption("mixedelasticity.gmsh.submesh")) ? cmeshPoisson : cmeshElasticity ;

		mpe_type MPE( cmeshPoisson, cmeshElasticity, meshCommon, mesh); 
		MPE.run( Idh_el, Idhv_el, Idh_poi, Idhv_poi );
    }

	
	return 0;
}
