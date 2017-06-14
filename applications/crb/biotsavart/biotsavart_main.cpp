#include "biotsavart.hpp"
#include "thermoelectric.hpp"

using namespace Feel;

int main( int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=biotsavartOptions()
                     .add(crbOptions())
                     .add(crbSEROptions())
                     .add(eimOptions())
                     .add(podOptions())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("ThermoElectricCRB"))
                     .add(makeOptions()) );

    BiotSavartCRB<Thermoelectric> bs = BiotSavartCRB<Thermoelectric>();
    bs.runBS();

    return 0;
}
