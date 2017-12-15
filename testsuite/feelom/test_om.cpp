#include "feel/feelfmi/fmu.hpp"


int main( int argc, char* argv[] )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv );

    FMU my_fmu;
    int ret = my_fmu.load();
    std::cout << "all ok for now ! \n";
    my_fmu.simulate();
}
