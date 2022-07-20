#include "view_factor.hpp"
#include <math.h>
#define BOOST_TEST_MODULE view_factor_tests tests
using namespace Feel;
#include <feel/feelcore/testsuite.hpp>
inline
po::options_description
makeOptions()
{
    po::options_description code1doptions( "test options" );
    return code1doptions.add( Feel::feel_options() );
}

//_____________________________________________________________________________________________________//
//_____________________________________________________________________________________________________//

inline
AboutData
makeAbout()
{
    AboutData about( "Test_view_factor_circle" ,
                     "Test_view_factor_circle" ,
                     "0.1",
                     "Test : Test_view_factor_circle",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2022 Universite de Strasbourg" );

    about.addAuthor( "Luca Berti", "developer", "", "" );
    return about;

}

double PointtoCircle_exact(double H, double R)
{
    double h = H/R;
    return 1./(1.+h*h);
}


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( view_factor_testsuite )
BOOST_AUTO_TEST_CASE( view_factor_test_circle )
{
    using namespace Feel;

    double H=2.0;
    double R=1.0;                    
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>); 
    int N = 100;

    view_factor<Mesh<Simplex<3>>> vf(mesh,N);
    
    Eigen::Vector3d patch_origin,patch_normal;
    patch_origin << 0.,0.,0;
    patch_normal << 0.,0.,-1.;
 
    auto M_view_factor = vf.computeViewFactorPatchMeshControlledRays(patch_origin,patch_normal); 
    

    double vf_exact = PointtoCircle_exact(H,R);    

    std::cout << M_view_factor[2] << " "<< M_view_factor[1] << " "<< M_view_factor[0] << " "<< vf_exact << std::endl;

    BOOST_CHECK_SMALL((vf_exact-M_view_factor[2])/vf_exact,5e-2);    

    
}
BOOST_AUTO_TEST_SUITE_END()