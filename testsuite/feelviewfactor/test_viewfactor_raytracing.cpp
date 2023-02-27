/**
 * @file test_viewfactor.cpp
 * @author Christophe Prud'homme (christophe.prudhomme@cemosis.fr)
 * @brief various tests for view factors computations
 * @version 0.1
 * @date 2022-07-29
 * 
 * @copyright Copyright (c) 2022 Feel++ Consortium
 * @copyright copyright (c) 2022 Université de Strasbourg
 * 
 */

#define BOOST_TEST_MODULE test_viewfactor_raytracing
#include <feel/feelcore/testsuite.hpp>

#include <fmt/ostream.h>
#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
//#include <feel/feelviewfactor/unobstructedplanarviewfactor.hpp>
#include <feel/feelviewfactor/raytracingviewfactor.hpp>

/** use Feel namespace */
using namespace Feel;
using Feel::project;

inline
AboutData
makeAbout()
{
    AboutData about( "test_vf" ,
                     "test_vf" ,
                     "0.2",
                     "nD(n=2,3) test vf",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2022 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@cemosis.fr", "" );
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

std::map<std::string,double> cylinder_vf_map(double radius, double height)
{
    double r = radius/height;
    double rho = (sqrt(4*math::pow(r,2)+1)-1)/r;
    std::map<std::string,double> vf{{"BasisToBasis",1-rho/(2*r)},\
                                    {"LateralToLateral",1-rho*0.5},\
                                    {"BasisToLateral",rho/(2*r)},
                                    {"LateralToBasis",rho*0.25}};
    return vf;
    
}

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description opts( "Test Environment options" );
    opts.add_options()
        ( "N", Feel::po::value<int>()->default_value( 2 ), "Number of problems to solve successively" )
    ;
    opts.add( Feel::feel_options(  ) );
    opts.add(feel_options( "cube-raytracing" ) );
    opts.add(feel_options( "cylinder-raytracing" ) );
    return opts;
}

template<typename MeshType>
void checkViewFactorRaytracing(std::string const& prefix)
{   
    auto mesh = loadMesh(_mesh=new MeshType, _filename =Environment::expand(soption(fmt::format("{}.gmsh.filename",prefix) ) )  );
    std::cout << Environment::expand(soption(fmt::format("{}.gmsh.filename",prefix) ) )<< std::endl;
    auto jsons = vsoption( fmt::format("{}.json.filename", prefix ) );
    BOOST_TEST_MESSAGE( fmt::format( "reading json: {}", Environment::expand(jsons[0]) ) );
    nl::json j;
    std::ifstream f(  Environment::expand(jsons[0]) );
    f >> j;
    BOOST_TEST_MESSAGE(fmt::format("j: {}",j.dump(1)));
 
    RayTracingViewFactor<MeshType> rtvf( mesh, j );
    
    rtvf.compute();
    BOOST_TEST_MESSAGE( fmt::format("Max dev reciprocity {}", rtvf.maxDevReciprocity()));
    
    BOOST_TEST_MESSAGE( fmt::format("{}", rtvf.viewFactors() ) );
    auto row_sum_vf = rtvf.viewFactors().rowwise().sum();
    auto exact_vf = eigen_vector_x_col_type<double>::Ones(rtvf.viewFactors().rows()) ;
    auto difference_infNorm = (exact_vf-row_sum_vf).template lpNorm<Eigen::Infinity>();
    BOOST_TEST_MESSAGE( fmt::format("View factors sum to one {}; infinity norm of 1 - rowwise sum {}; rowwise sum {} ", difference_infNorm<1e-3, difference_infNorm,row_sum_vf) );
        
    if(prefix=="cube-raytracing")
    {
        BOOST_CHECK_MESSAGE( difference_infNorm<5e-2, fmt::format("Infinity norm of (1 - rowwise) sum is larger than 4e-2" ) );
    }
    else if(prefix=="cylinder-raytracing")
    {
        BOOST_CHECK_MESSAGE( difference_infNorm<5e-2, fmt::format("Infinity norm of (1 - rowwise) sum is larger than 4e-2" ) );
    }

    BOOST_CHECK_MESSAGE(rtvf.maxDevReciprocity()<2e-2, fmt::format("Max dev reciprocity less than 2e-2, {}",rtvf.maxDevReciprocity()));

    if(prefix=="cube-raytracing")
    {
        auto vf_parallel_walls = view_factor_parallel_walls_exact(1.,1.,1.); // cube of side 1.
        auto vf_perp_walls = view_factor_perp_walls_exact(1.,1.,1.); // cube of side 1.
        // std::cout << "vf_parallel_walls" << vf_parallel_walls << std::endl;
        // std::cout << "vf_perp_walls" << vf_perp_walls << std::endl;
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(0,1)-vf_parallel_walls)/vf_parallel_walls <4e-2, fmt::format("Relative error view factors between parallel walls 0 1 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(2,3)-vf_parallel_walls)/vf_parallel_walls <4e-2, fmt::format("Relative error view factors between parallel walls 2 3 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(4,5)-vf_parallel_walls)/vf_parallel_walls <4e-2, fmt::format("Relative error view factors between parallel walls 4 5 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(0,2)-vf_perp_walls)/vf_perp_walls <4e-2, fmt::format("Relative error view factors between perp walls 0 2 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(0,3)-vf_perp_walls)/vf_perp_walls <4e-2, fmt::format("Relative error view factors between perp walls 0 3 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(0,4)-vf_perp_walls)/vf_perp_walls <4e-2, fmt::format("Relative error view factors between perp walls 0 4 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(0,5)-vf_perp_walls)/vf_perp_walls <4e-2, fmt::format("Relative error view factors between perp walls 0 5 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(1,3)-vf_perp_walls)/vf_perp_walls <4e-2, fmt::format("Relative error view factors between perp walls 1 3 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(1,4)-vf_perp_walls)/vf_perp_walls <4e-2, fmt::format("Relative error view factors between perp walls 1 4 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(1,5)-vf_perp_walls)/vf_perp_walls <4e-2, fmt::format("Relative error view factors between perp walls 1 5 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(2,4)-vf_perp_walls)/vf_perp_walls <4e-2, fmt::format("Relative error view factors between perp walls 2 4 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(2,5)-vf_perp_walls)/vf_perp_walls <4e-2, fmt::format("Relative error view factors between perp walls 2 5 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(3,5)-vf_perp_walls)/vf_perp_walls <4e-2, fmt::format("Relative error view factors between perp walls 3 5 is less than 4e-2 ") );
    }
    else if(prefix=="cylinder-raytracing")
    {
        auto vf_cylinder = cylinder_vf_map(1.,2.);  
        std::cout <<       rtvf.viewFactors() << std::endl;
        // std::cout <<       vf_cylinder["BasisToBasis"] << std::endl;
        // std::cout <<       vf_cylinder["BasisToLateral"]<< std::endl;
        // std::cout <<      vf_cylinder["LateralToBasis"] << std::endl;
        // std::cout <<       vf_cylinder["LateralToLateral"] << std::endl;
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(0,1)-vf_cylinder["BasisToBasis"])/vf_cylinder["BasisToBasis"] <4e-2, fmt::format("Relative error view factors between bases 0 1 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(0,2)-vf_cylinder["BasisToLateral"])/vf_cylinder["BasisToLateral"] <4e-2, fmt::format("Relative error view factors between bottom basis to lateral walls 0 2 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(1,0)-vf_cylinder["BasisToBasis"])/vf_cylinder["BasisToBasis"] <4e-2, fmt::format("Relative error view factors between bases 1 0 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(1,2)-vf_cylinder["BasisToLateral"])/vf_cylinder["BasisToLateral"] <4e-2, fmt::format("Relative error view factors between bottom basis to lateral walls 1 2 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(2,0)-vf_cylinder["LateralToBasis"])/vf_cylinder["LateralToBasis"]  <4e-2, fmt::format("Relative error view factors between lateral wall to top basis 2 0 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(2,1)-vf_cylinder["LateralToBasis"])/vf_cylinder["LateralToBasis"]  <4e-2, fmt::format("Relative error view factors between lateral wall to bottom basis 2 1 is less than 4e-2 ") );
        BOOST_CHECK_MESSAGE( abs(rtvf.viewFactors()(2,2)-vf_cylinder["LateralToLateral"])/vf_cylinder["LateralToLateral"]  <4e-2, fmt::format("Relative error view factors between lateral wall to itself 2 2 is less than 4e-2: {}",(rtvf.viewFactors()(2,2)-vf_cylinder["LateralToLateral"])/vf_cylinder["LateralToLateral"]) );
    }
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( viewfactor_raytracing )

BOOST_AUTO_TEST_CASE( test_cube_raytracing )
{
    checkViewFactorRaytracing<Mesh<Simplex<3>>>( "cube-raytracing" );
}
BOOST_AUTO_TEST_CASE( test_cylinder_raytracing )
{
    checkViewFactorRaytracing<Mesh<Simplex<3>>>( "cylinder-raytracing" );
}

// Tests for the raytracing part

// Function tested:
// get_random_direction
// Test if the algorithm computing the random directions on a sphere or circle 
// provides a uniform distribution over the manifold using a chi-squared test
#if 1
BOOST_AUTO_TEST_CASE( test_random_direction )
{
    std::vector<double> three_dimensional_direction(3);
    std::vector<double> two_dimensional_direction(2);
    std::vector<double> two_dimensional_bins(4);
    std::vector<double> three_dimensional_bins(8);
    Eigen::VectorXd normal3d(3),normal2d(2);
    std::mt19937 gen1;
    std::mt19937 gen2;
    normal3d << 0,0,1;
    normal2d << 0,1;

    for(int i=0; i<1e4; i++)
    {
        get_random_direction(three_dimensional_direction,gen1,gen2,normal3d);
        get_random_direction(two_dimensional_direction,gen1,gen2,normal2d);

        BOOST_CHECK_SMALL(math::pow(two_dimensional_direction[0],2)+math::pow(two_dimensional_direction[1],2)-1,1e-6);
        BOOST_CHECK_SMALL(math::pow(three_dimensional_direction[0],2)+math::pow(three_dimensional_direction[1],2)+math::pow(three_dimensional_direction[2],2)-1,1e-6);

        // Build bins for the 8 spherical sectors and 4 circular sectors to check if the distribution is uniform
        // 2d case - 4 sectors
        if(two_dimensional_direction[0]>=0 && two_dimensional_direction[1]>=0)
        {
            two_dimensional_bins[0]+=1;
        }
        else if (two_dimensional_direction[0]<0 && two_dimensional_direction[1]>=0)
        {
            two_dimensional_bins[1]+=1;
        }
        else if (two_dimensional_direction[0]>=0 && two_dimensional_direction[1]<0)
        {
            two_dimensional_bins[2]+=1;
        }
        else if (two_dimensional_direction[0]<0 && two_dimensional_direction[1]<0)
        {
            two_dimensional_bins[3]+=1;
        }
        else
        {
            throw std::logic_error( "Problem for the 2d random direction" );
        }

        // 3d case - 8 sectors
        if(three_dimensional_direction[0]>=0 && three_dimensional_direction[1]>=0 && three_dimensional_direction[2]>=0)
        {
            three_dimensional_bins[0]+=1;
        }
        else if (three_dimensional_direction[0]>=0 && three_dimensional_direction[1]>=0 && three_dimensional_direction[2]<0)
        {
            three_dimensional_bins[1]+=1;
        }
        else if (three_dimensional_direction[0]>=0 && three_dimensional_direction[1]<0 && three_dimensional_direction[2]<0)
        {
            three_dimensional_bins[2]+=1;
        }
        else if (three_dimensional_direction[0]>=0 && three_dimensional_direction[1]<0 && three_dimensional_direction[2]>=0)
        {
            three_dimensional_bins[3]+=1;
        }
        else if (three_dimensional_direction[0]<0 && three_dimensional_direction[1]>=0 && three_dimensional_direction[2]>=0)
        {
            three_dimensional_bins[4]+=1;
        }
        else if (three_dimensional_direction[0]<0 && three_dimensional_direction[1]>=0 && three_dimensional_direction[2]<0)
        {
            three_dimensional_bins[5]+=1;
        }
        else if (three_dimensional_direction[0]<0 && three_dimensional_direction[1]<0 && three_dimensional_direction[2]>=0)
        {
            three_dimensional_bins[6]+=1;
        }
        else if (three_dimensional_direction[0]<0 && three_dimensional_direction[1]<0 && three_dimensional_direction[2]<0)
        {
            three_dimensional_bins[7]+=1;
        }
        else
        {
            throw std::logic_error( "Problem for the 3d random direction" );
        }
    }
    double expected_bin_content_3d = 1./8.;
    double expected_bin_content_2d = 1./4.;
    double chi_squared_three_dim = 0.;
    double chi_squared_two_dim = 0.;

    for(int i=0;i<4;i++)
    {
        chi_squared_two_dim += math::pow(two_dimensional_bins[i]/1e4-expected_bin_content_2d,2);
    }
    for(int i=0;i<8;i++)
    {
        chi_squared_three_dim += math::pow(three_dimensional_bins[i]/1e4-expected_bin_content_3d,2);
    }

    chi_squared_two_dim = chi_squared_two_dim/expected_bin_content_2d;
    chi_squared_three_dim = chi_squared_three_dim/expected_bin_content_3d;    
    std::cout << "Chisquared" << chi_squared_two_dim << " " << chi_squared_three_dim << std::endl;
    BOOST_CHECK(chi_squared_two_dim < 9.48); // p-value larger than 0.05 for chi squared distrib(4)
    BOOST_CHECK(chi_squared_three_dim < 15.50); // p-value larger than 0.05 for chi squared distrib(8)

}
#endif 
// Function tested:
// isOnSurface, elementArea
// Test if the computed area is correct and if a point belongs on a triangle
BOOST_AUTO_TEST_CASE( test_element_area_and_point_belonging )
{
    Eigen::VectorXd point1_2d(2),point2_2d(2),barycenter_2d(2);
    Eigen::VectorXd point1_3d(3),point2_3d(3),point3_3d(3),barycenter_3d(3);
    double area2d,area3d;
    bool isOnSurf2d,isOnSurf3d;

    point1_2d << 0.,3.;
    point2_2d << 3.,3.;    

    point1_3d << 0.,3.,1.;
    point2_3d << 3.,3.,1.;
    point3_3d << 3.,0.,1.;

    barycenter_2d = (point1_2d+point2_2d)/2.;
    barycenter_3d = (point1_3d+point2_3d+point3_3d)/3.;

    area2d = element_area(barycenter_2d,point1_2d,point2_2d);
    area3d = element_area(barycenter_3d,point1_3d,point2_3d,point3_3d);

    BOOST_CHECK_SMALL(area2d-3.,1e-6);
    BOOST_CHECK_SMALL(area3d-4.5,1e-6);

    isOnSurf2d = isOnSurface(barycenter_2d,point1_2d,point2_2d);
    isOnSurf3d = isOnSurface(barycenter_3d,point1_3d,point2_3d,point3_3d);

    BOOST_CHECK(isOnSurf2d==true);
    BOOST_CHECK(isOnSurf3d==true);

    isOnSurf2d = isOnSurface(point1_2d,point1_2d,point2_2d);
    isOnSurf3d = isOnSurface(point1_3d,point1_3d,point2_3d,point3_3d);

    BOOST_CHECK(isOnSurf2d==true);
    BOOST_CHECK(isOnSurf3d==true);

    Eigen::VectorXd new_point_2d(2);
    Eigen::VectorXd new_point_3d(3);

    new_point_2d << 4.,4.;
    new_point_3d << 4.,4.,4.;

    isOnSurf2d = isOnSurface(new_point_2d,point1_2d,point2_2d);
    isOnSurf3d = isOnSurface(new_point_3d,point1_3d,point2_3d,point3_3d);

    BOOST_CHECK(isOnSurf2d==false);
    BOOST_CHECK(isOnSurf3d==false);


}

// Function tested:
// get_random_point
BOOST_AUTO_TEST_CASE( test_random_point )
{
    RayTracingViewFactor<Mesh<Simplex<2>>> rtv_class_2d;
    RayTracingViewFactor<Mesh<Simplex<3>>> rtv_class_3d;
    
    Eigen::VectorXd point1_2d(2),point2_2d(2),random_point_2d(2);
    Eigen::VectorXd point1_3d(3),point2_3d(3),point3_3d(3),random_point_3d(3);
    double area2d,area3d;
    bool isOnSurf2d,isOnSurf3d;

    ublas::matrix<double, ublas::column_major> matrix_2d(2,2),matrix_3d(3,3);

    point1_2d << 0.,3.;
    point2_2d << 3.,3.;

    point1_3d << 0.,3.,1.;
    point2_3d << 3.,3.,1.;
    point3_3d << 3.,0.,1.;

    for(int j =0;j<3;j++)
    {
        matrix_3d(j,0) = point1_3d(j);
        matrix_3d(j,1) = point2_3d(j);
        matrix_3d(j,2) = point3_3d(j);
    }
    for(int j=0;j<2;j++)
    {
        matrix_2d(j,0) = point1_2d(j);
        matrix_2d(j,1) = point2_2d(j);
    }

    random_point_2d = rtv_class_2d.get_random_point(matrix_2d);

    random_point_3d = rtv_class_3d.get_random_point(matrix_3d);

    isOnSurf2d = isOnSurface(random_point_2d,point1_2d,point2_2d);
    isOnSurf3d = isOnSurface(random_point_3d,point1_3d,point2_3d,point3_3d);

    BOOST_CHECK(isOnSurf2d==true);
    BOOST_CHECK(isOnSurf3d==true);
}
#if 0
// Functions tested:
// checkIntersectionWithSegment, checkIntersectionWithTriangle
BOOST_AUTO_TEST_CASE( test_ray_intersections )
{
    
    RayTracingViewFactor<Mesh<Simplex<2>>> rtv_class_2d;
    RayTracingViewFactor<Mesh<Simplex<3>>> rtv_class_3d;

    Eigen::VectorXd point1_2d(2),point2_2d(2),random_point_2d(2);
    Eigen::VectorXd point1_3d(3),point2_3d(3),point3_3d(3),random_point_3d(3);
    double area2d,area3d;
    bool isOnSurf2d,isOnSurf3d;

    ublas::matrix<double, ublas::column_major> matrix_2d(2,2),matrix_3d(3,3);

    point1_2d << 0.,3.;
    point2_2d << 3.,3.;

    point1_3d << 0.,3.,1.;
    point2_3d << 3.,3.,1.;
    point3_3d << 3.,6.,1.;

    for(int j =0;j<3;j++)
    {
        matrix_3d(j,0) = point1_3d(j);
        matrix_3d(j,1) = point2_3d(j);
        matrix_3d(j,2) = point3_3d(j);
    }
    for(int j=0;j<2;j++)
    {
        matrix_2d(j,0) = point1_2d(j);
        matrix_2d(j,1) = point2_2d(j);
    }

    Eigen::VectorXd ray2d_orig(2),ray2d_dir(2);
    Eigen::VectorXd ray3d_orig(3),ray3d_dir(3);

    ray2d_orig << 0.,1.;
    ray3d_orig << 0.,3.,0.5;
        
    for(int i=0;i<100;i++)
    {   
        // Check intersection between the segment (0,3)--(3,3) and the rays (0,1) + t(cos(alpha),sin(alpha))
        // There should be intersection when atan(2/3) <= alpha <= pi/2
        double alpha = i/100.*2*pi;
        ray2d_dir << math::cos(alpha),math::sin(alpha);
        RayTracingViewFactor<Mesh<Simplex<2>>>::Ray ray2d(ray2d_orig,ray2d_dir);
        auto [check2d,point2d] = rtv_class_2d.checkIntersectionWithSegment(matrix_2d,ray2d);
        if(alpha>=math::atan(2./3.) && alpha <= pi/2.)
            {
                BOOST_CHECK(check2d==true);
            }
        else
            {
                BOOST_CHECK(check2d==false);
            }
    }

    for(int i=0;i<100;i++)
    {   
        double alpha = i/100.*2*pi;
        for(int j=0;j<100;j++)
        {
            // Check intersection between the triangle (0,3,1)--(3,3,1)--(3,6,1) and the 
            // rays (0,3,0.5) + t(cos(alpha)*sin(beta),sin(alpha)*sin(beta),cos(beta))
            // There should be intersection when 0 <= alpha <= pi/4            
            // and (0.5/3*cos(alpha) <= beta <= pi/2)
            double beta = j/100.*pi;
            ray3d_dir << math::cos(alpha)*math::sin(beta),math::sin(alpha)*math::sin(beta),math::cos(beta);
            RayTracingViewFactor<Mesh<Simplex<3>>>::BVHRay ray3d(ray3d_orig,ray3d_dir);
            auto [check3d,point3d] = rtv_class_3d.checkIntersectionWithTriangle(matrix_3d,ray3d);
            if((beta<=math::atan(3/0.5/math::cos(alpha)) && beta > 0. && alpha>=0 && alpha <=pi/4. )||beta==0)
            {
                // BOOST_CHECK_MESSAGE(check3d==true);
                BOOST_CHECK_MESSAGE(check3d==true,"CHECK TRUE NOT WORKING" << alpha << " " << beta << " " << point3d);
            }
            else
            {
                BOOST_CHECK_MESSAGE(check3d==false,"CHECK FALSE NOT WORKING" << alpha << " " << beta << " " << point3d);
            }
        }
    }
}
#endif

BOOST_AUTO_TEST_SUITE_END()