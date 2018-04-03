#include "feel/feelfmi/fmu.hpp"

int main( int argc, char* argv[] )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv );

    FMU my_fmu;
    my_fmu.load();
    my_fmu.printModelInfo();

    my_fmu.simulate();
}
