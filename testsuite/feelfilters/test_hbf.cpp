#define BOOST_TEST_MODULE hbf

#include <feel/feelfilters/gmsh.hpp>

int test_Write_Read()
{
    holo3_image<int_8> x (2,2);
    std::string const * s;
    writeHBF(s, x);
    holo3_image<int_8> y = readHBF(s);

    int res = (x==y);

    return res;
}

int main()
{
    int res = test_Write_Read();
    return res;
}