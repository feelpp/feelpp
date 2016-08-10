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

	decltype( IPtr( _domainSpace=Pdh<FEELPP_ORDER>(mesh), _imageSpace=Pdh<1>(mesh) ) ) Idh_poi ;
    decltype( IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(mesh), _imageSpace=Pdhv<1>(mesh) ) ) Idhv_poi;

    decltype( IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(mesh), _imageSpace=Pdhv<1>(mesh) ) ) Idh_el ;
    decltype( IPtr( _domainSpace=Pdhms<FEELPP_ORDER>(mesh), _imageSpace=Pdhms<1>(mesh) ) ) Idhv_el;
		
	if ( soption( "gmsh.submesh" ).empty() )
	{
		mpe_type MPE( mesh );
		MPE.run( mesh, Idh_el, Idhv_el, Idh_poi, Idhv_poi );
	}
	else
	{
		Feel::cout << "Using submesh: " << soption("gmsh.submesh") << std::endl;
	    auto cmesh = createSubmesh( mesh, markedelements(mesh,soption("gmsh.submesh")), Environment::worldComm() );

		Idh_poi = IPtr( _domainSpace=Pdh<FEELPP_ORDER>(cmesh), _imageSpace=Pdh<1>(mesh) );
    	Idhv_poi = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmesh), _imageSpace=Pdhv<1>(mesh) );

    	Idh_el = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmesh), _imageSpace=Pdhv<1>(mesh) );
		Idhv_el = IPtr( _domainSpace=Pdhms<FEELPP_ORDER>(cmesh), _imageSpace=Pdhms<1>(mesh) );
		
		mpe_type MPE( cmesh, mesh); 
		MPE.run( mesh, Idh_el, Idhv_el, Idh_poi, Idhv_poi );
	}

	
	return 0;
}
