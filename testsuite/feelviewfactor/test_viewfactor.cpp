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

#define BOOST_TEST_MODULE test_viewfactor
#include <feel/feelcore/testsuite.hpp>

#include <fmt/ostream.h>
#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelviewfactor/unobstructedplanarviewfactor.hpp>
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

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description opts( "Test Environment options" );
    opts.add_options()
        ( "N", Feel::po::value<int>()->default_value( 2 ), "Number of problems to solve successively" )
    ;
    return opts.add( Feel::feel_options(  ) ).add( Feel::feel_options( "square" ) ).add( feel_options( "cube" ) );
}

template<typename MeshType>
void checkViewFactorEnclosure(std::string const& prefix)
{
    //auto mesh = loadMesh( _mesh=new MeshType, _filename=Environment::expand(soption(fmt::format("{}.gmsh.filename",prefix) ) ) );
    // auto mesh = createGMSHMesh( _mesh=new MeshType,
    //                                     _desc=domain( _name=Environment::expand(soption(fmt::format("{}.gmsh.filename",prefix) ) )  ,
    //                                                   _dim=(prefix=="square"?2:3) ) );
    
    auto mesh = loadMesh(_mesh=new MeshType, _filename =Environment::expand(soption(fmt::format("{}.gmsh.filename",prefix) ) )  );
    std::cout << Environment::expand(soption(fmt::format("{}.gmsh.filename",prefix) ) )<< std::endl;
    auto jsons = vsoption( fmt::format("{}.json.filename", prefix ) );
    BOOST_TEST_MESSAGE( fmt::format( "reading json: {}", Environment::expand(jsons[0]) ) );
    nl::json j;
    std::ifstream f(  Environment::expand(jsons[0]) );
    f >> j;
    BOOST_TEST_MESSAGE(fmt::format("j: {}",j.dump(1)));
 
    UnobstructedPlanarViewFactor<MeshType> upvf( mesh, j );
    std::cout << "ok till here" << std::endl;
    upvf.compute();
    std::cout << "ok till here2" << std::endl;
    BOOST_TEST_MESSAGE( fmt::format("{}", upvf.viewFactors() ) );
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( viewfactor )


BOOST_AUTO_TEST_CASE( test_square )
{
    checkViewFactorEnclosure<Mesh<Simplex<2>>>("square");
}
BOOST_AUTO_TEST_CASE( test_cube )
{
    checkViewFactorEnclosure<Mesh<Simplex<3>>>( "cube" );
}

// Tests for the raytracing part

// Function tested:
// get_random_direction
// Test if the algorithm computing the random directions on a sphere or circle 
// provides a uniform distribution over the manifold using a chi-squared test
BOOST_AUTO_TEST_CASE( test_random_direction )
{
    std::vector<double> three_dimensional_direction(3);
    std::vector<double> two_dimensional_direction(2);
    std::vector<double> two_dimensional_bins(4);
    std::vector<double> three_dimensional_bins(8);
    for(int i=0; i<1e4; i++)
    {
        get_random_direction(three_dimensional_direction);
        get_random_direction(two_dimensional_direction);

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
        if(two_dimensional_direction[0]>=0 && two_dimensional_direction[1]>=0 && two_dimensional_direction[2]>=0)
        {
            three_dimensional_bins[0]+=1;
        }
        else if (two_dimensional_direction[0]>=0 && two_dimensional_direction[1]>=0 && two_dimensional_direction[2]<0)
        {
            three_dimensional_bins[1]+=1;
        }
        else if (two_dimensional_direction[0]>=0 && two_dimensional_direction[1]<0 && two_dimensional_direction[2]<0)
        {
            three_dimensional_bins[2]+=1;
        }
        else if (two_dimensional_direction[0]>=0 && two_dimensional_direction[1]<0 && two_dimensional_direction[2]>=0)
        {
            three_dimensional_bins[3]+=1;
        }
        else if (two_dimensional_direction[0]<0 && two_dimensional_direction[1]>=0 && two_dimensional_direction[2]>=0)
        {
            three_dimensional_bins[4]+=1;
        }
        else if (two_dimensional_direction[0]<0 && two_dimensional_direction[1]>=0 && two_dimensional_direction[2]<0)
        {
            three_dimensional_bins[5]+=1;
        }
        else if (two_dimensional_direction[0]<0 && two_dimensional_direction[1]<0 && two_dimensional_direction[2]>=0)
        {
            three_dimensional_bins[6]+=1;
        }
        else if (two_dimensional_direction[0]<0 && two_dimensional_direction[1]<0 && two_dimensional_direction[2]<0)
        {
            three_dimensional_bins[7]+=1;
        }
        else
        {
            throw std::logic_error( "Problem for the 2d random direction" );
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
    BOOST_CHECK(chi_squared_two_dim < 9.48); // p-value larger than 0.05 for chi squared distrib(4)
    BOOST_CHECK(chi_squared_three_dim < 15.50); // p-value larger than 0.05 for chi squared distrib(8)

}
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

// Functions tested:
// check_intersection_with_segment, check_intersection_with_triangle
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
        auto [check2d,point2d] = rtv_class_2d.check_intersection_with_segment(matrix_2d,ray2d);
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
            RayTracingViewFactor<Mesh<Simplex<3>>>::Ray ray3d(ray3d_orig,ray3d_dir);
            auto [check3d,point3d] = rtv_class_3d.check_intersection_with_triangle(matrix_3d,ray3d);
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
    
// BOOST_AUTO_TEST_CASE( initialize_raytracing )
// {    
    
//     auto mesh3d =  createGMSHMesh( _mesh=new Mesh<Simplex<3>>,
//                            _desc=geo( _filename="${cfgdir}/raytracing/view_factor_cube.geo",
//                                       _dim=3,
//                                       _order=1,
//                                       _h=0.2 ) );
//     //loadMesh(_mesh=new Mesh<Simplex<2>>,_filename=Environment::expand(soption("square.gmsh.filename"))); 
//     auto mesh2d = loadMesh(_mesh=new Mesh<Simplex<2>>,_filename=Environment::expand(soption("square.gmsh.filename"))); 

//     RayTracingViewFactor<Mesh<Simplex<2>>> rtv_class_2d(mesh2d,10);
//     RayTracingViewFactor<Mesh<Simplex<3>>> rtv_class_3d(mesh3d,10);

//     std::cout << "Markers 2d" << rtv_class_2d.markerNames() <<std::endl;
//     std::cout << "Markers 3d" << rtv_class_3d.markerNames() <<std::endl;

// }


BOOST_AUTO_TEST_SUITE_END()