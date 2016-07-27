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


    typedef FeelModels::MixedElasticity<2,1> me_type;

    auto ME = me_type::New("mixedelasticity");
    auto mesh = loadMesh( _mesh=new me_type::mesh_type );
    
    ME -> init(mesh);
 
    return 0;
}



