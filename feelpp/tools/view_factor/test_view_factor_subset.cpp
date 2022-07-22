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
// Compute the view factor for a subset of markers of the given mesh
inline
AboutData
makeAbout()
{
    AboutData about( "Test_view_factor_subset" ,
                     "Test_view_factor_subset" ,
                     "0.1",
                     "Test : Test_view_factor_subset",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2022 Universite de Strasbourg" );

    about.addAuthor( "Luca Berti", "developer", "", "" );
    return about;

}

double view_factor_cylinder_bottom_top(double height, double radius)
{
    double view_factor_cylinder_bottom_top=0.;
    double b=radius/height;
    double rho=(sqrt(4*b*b+1)-1)/b;
    view_factor_cylinder_bottom_top = 1-rho/(2*b);
    
    return view_factor_cylinder_bottom_top;
}

double view_factor_cylinder_bottom_lateral(double height, double radius)
{

    double view_factor_cylinder_bottom_lateral=0.;
    double b=radius/height;
    double rho=(sqrt(4*b*b+1)-1)/b;
    view_factor_cylinder_bottom_lateral = rho/(2*b);
    
    return view_factor_cylinder_bottom_lateral;
}
    

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( view_factor_testsuite )
BOOST_AUTO_TEST_CASE( view_factor_test_subset )
{
    using namespace Feel;

    double height=2.0;
    double radius=1.0;

    auto vf_exact_top = view_factor_cylinder_bottom_top(height,radius);
    auto vf_exact_lateral = view_factor_cylinder_bottom_lateral(height,radius);
    int constexpr nDim=3;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<nDim>>); 
    
    view_factor<Mesh<Simplex<nDim>>> view_f_topBottom(mesh,{"TopDisk","BottomDisk"},10);
    view_factor<Mesh<Simplex<nDim>>> view_f_lateralBottom(mesh,{"LateralSurface","BottomDisk"},10);


    auto viewFactors_tb = view_f_topBottom.computeViewFactors();
    auto viewFactors_lb = view_f_lateralBottom.computeViewFactors();

    // std::cout << viewFactors_tb << " " << vf_exact_top << std::endl;
    // std::cout << viewFactors_lb << " " << vf_exact_lateral << std::endl;

    BOOST_CHECK_SMALL((vf_exact_top-viewFactors_tb(0,1))/vf_exact_top,5e-2);
    BOOST_CHECK_SMALL((vf_exact_lateral-viewFactors_lb(0,1))/vf_exact_lateral,5e-2);
    
}
BOOST_AUTO_TEST_SUITE_END()