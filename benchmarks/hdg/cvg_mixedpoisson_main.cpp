#include <cvg_mixedpoisson.hpp>

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( (boost::format("conv-mixed-poisson_%1%DP%2%G%3%E%4%") % FEELPP_DIM % FEELPP_ORDER % FEELPP_GEO_ORDER % FEELPP_EXP_ORDER).str(),
                     (boost::format("conv-mixed-poisson_%1%DP%2%G%3%E%4%") % FEELPP_DIM % FEELPP_ORDER % FEELPP_GEO_ORDER % FEELPP_EXP_ORDER).str(),
                     "0.1",
                     "Convergence Mixed Poisson Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Christophe Trophime", "developer", "christophe.trophime@lncmi.cnrs.fr", "" );
    about.addAuthor( "Daniele Prada", "developer", "", "" );
    return about;

}

int main(int argc, char *argv[])
{
    using namespace Feel;
    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=makeConvOptions(),
                           _desc_lib=makeConvLibOptions().add(feel_options()) );

    typedef ConvergenceTest<FEELPP_DIM,FEELPP_ORDER,FEELPP_GEO_ORDER,FEELPP_EXP_ORDER> cv_type;

    cv_type CV;
    int status = CV.run();
    return status;
}
