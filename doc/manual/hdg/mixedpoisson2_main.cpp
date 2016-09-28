#include <mixedpoisson2.hpp>

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

int main(int argc, char *argv[])
{
    using namespace Feel;
    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=FeelModels::makeMixedPoissonOptions("mixedpoisson"),
                           _desc_lib=FeelModels::makeMixedPoissonLibOptions("mixedpoisson").add(feel_options())
                           );


    typedef FeelModels::MixedPoisson<FEELPP_DIM,FEELPP_ORDER> mp_type;

    auto MP = mp_type::New("mixedpoisson");
    auto mesh = loadMesh( _mesh=new mp_type::mesh_type );
    decltype( IPtr( _domainSpace=Pdh<FEELPP_ORDER>(mesh), _imageSpace=Pdh<1>(mesh) ) ) Idh ;
    decltype( IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(mesh), _imageSpace=Pdhv<1>(mesh) ) ) Idhv;
    if ( soption( "mixedpoisson.gmsh.submesh" ).empty() )
        MP -> init(mesh);
    else
    {
        Feel::cout << "Using submesh: " << soption("mixedpoisson.gmsh.submesh") << std::endl;
		auto cmesh = createSubmesh( mesh, markedelements(mesh,soption("mixedpoisson.gmsh.submesh")), Environment::worldComm() );
        //Idh = IPtr( _domainSpace=Pdh<FEELPP_ORDER>(cmesh), _imageSpace=Pdh<1>(mesh) );
        //Idhv = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmesh), _imageSpace=Pdhv<1>(mesh) );
        MP -> init( cmesh, mesh );
    }

    if ( MP -> isStationary() )
    {
        MP->solve();
        MP->exportResults( mesh, Idh, Idhv );
    }
    else
    {
        for ( ; !MP->timeStepBase()->isFinished() ; MP->updateTimeStep() )
        {
            Feel::cout << "============================================================\n";
            Feel::cout << "time simulation: " << MP->time() << "s \n";
            Feel::cout << "============================================================\n";
            MP->solve();
            MP->exportResults( mesh, Idh, Idhv );
        }
    }

    // MP->computeError();

    return 0;
}
