#include "feel/feelfmi/fmu.hpp"


int main( int argc, char* argv[] )
{
    using namespace Feel;

    FMU my_fmu;
    int ret = my_fmu.load("test_om/test_om.fmu");
    std::cout << "all ok for now ! \n";
}
