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
    AboutData about( "Test_view_factor_cube" ,
                     "Test_view_factor_cube" ,
                     "0.1",
                     "Test : Test_view_factor_cube",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2022 Universite de Strasbourg" );

    about.addAuthor( "Luca Berti", "developer", "", "" );
    return about;

}

double view_factor_parallel_walls_exact(double length, double width,double separation)
{
    double view_factor_bottom_to_top_wall=0.;
    double a=length;
    double b=width;
    double c=separation;
    double X = a/c;
    double Y = b/c;
    view_factor_bottom_to_top_wall = log(sqrt((1+X*X)*(1+Y*Y)/(1+X*X+Y*Y)));
    view_factor_bottom_to_top_wall += X*sqrt(1+Y*Y)*atan(X/sqrt(1+Y*Y));
    view_factor_bottom_to_top_wall += Y*sqrt(1+X*X)*atan(Y/sqrt(1+X*X));
    view_factor_bottom_to_top_wall += -X*atan(X)-Y*atan(Y);
    view_factor_bottom_to_top_wall *= 2/(M_PI*X*Y);

    return view_factor_bottom_to_top_wall;
}

double view_factor_perp_walls_exact(double length, double width,double separation)
{
    double view_factor_bottom_to_side_wall=0.;
    double a=separation;
    double b=length;
    double c=width;
    double h = a/c;
    double w = b/c;
    view_factor_bottom_to_side_wall = h*atan(1/h)+w*atan(1/w);
    view_factor_bottom_to_side_wall += -sqrt(h*h+w*w)*atan(1/sqrt(h*h+w*w));
    double fact1 = (1+h*h)*(1+w*w)/(1+h*h+w*w);
    double fact2 = w*w*(1+h*h+w*w)/(((1+w*w))*(h*h+w*w));
    double fact3 = h*h*(1+h*h+w*w)/(((1+h*h))*(h*h+w*w));
    view_factor_bottom_to_side_wall +=0.25*log(fact1*pow(fact2,w*w)*pow(fact3,h*h));    
    view_factor_bottom_to_side_wall *= 1/(M_PI*w);

    return view_factor_bottom_to_side_wall;
}

double PointtoCircle_exact(double H, double R)
{
    double h = H/R;
    return 1./(1.+h*h);
}


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( view_factor_testsuite )
BOOST_AUTO_TEST_CASE( view_factor_test_cube )
{
    using namespace Feel;

    double length=1.0;
    double width=1.0;
    double separation=1.0; 

    auto vf_exact_parall = view_factor_parallel_walls_exact(length,width,separation);
    auto vf_exact_perp = view_factor_perp_walls_exact(length,width,separation);
    int constexpr nDim=3;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<nDim>>); 
    
    view_factor<Mesh<Simplex<nDim>>> view_f(mesh,10);

    auto viewFactors = view_f.computeViewFactors();

    for (std::string i: view_f.markerNames())
        std::cout << i << ' ';
    std::cout << "" << std::endl;    
    std::cout << viewFactors << std::endl;

    std::cout << "Theoretical values (theoretical values not known for the last line, from cube walls): " << std::endl;

    Eigen::MatrixXd theoretical_vf(view_f.markerNames().size(),view_f.markerNames().size());

    theoretical_vf.row(0) <<0. , vf_exact_parall , 4*vf_exact_perp;
    theoretical_vf.row(1) << vf_exact_parall ,0., 4*vf_exact_perp;
    theoretical_vf.row(2) << 0.,0.,0.;

    std::cout << theoretical_vf << std::endl;

    std::cout << "Factors sum to 1? " << std::endl;
    std::cout << viewFactors.rowwise().sum() << std::endl;

    BOOST_CHECK_SMALL((vf_exact_parall-viewFactors(0,1))/vf_exact_parall,5e-2);
    BOOST_CHECK_SMALL((4*vf_exact_perp-viewFactors(0,2))/(4*vf_exact_perp),5e-2);
    BOOST_CHECK_SMALL((vf_exact_parall-viewFactors(1,0))/vf_exact_parall,5e-2);
    BOOST_CHECK_SMALL((4*vf_exact_perp-viewFactors(1,2))/(4*vf_exact_perp),5e-2);

    BOOST_CHECK_EQUAL(viewFactors.rowwise().sum().sum(),3);
    
}
BOOST_AUTO_TEST_SUITE_END()