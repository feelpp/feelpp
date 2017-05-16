#include "thermoelectric.hpp"

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "thermoelectric-hdg-model" ,
                     "thermoelectric-hdg-model" ,
                     "0.1",
                     "ThermoElectricHDG-Model",
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
                           _desc=makeThermoElectricHDGOptions(),
                           _desc_lib=feel_options()
                           );

    using thermoelectric_type = ThermoElectricHDG<FEELPP_DIM, FEELPP_ORDER, FEELPP_ORDER>;
    using thermoelectric_ptrtype = boost::shared_ptr<thermoelectric_type>;

    auto te = boost::make_shared<thermoelectric_type>();
    te->run();

    return 0;
}
