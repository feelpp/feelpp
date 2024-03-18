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

#define BOOST_TEST_MODULE test_viewfactor_quadrature
#include <feel/feelcore/testsuite.hpp>

#include <fmt/ostream.h>
#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelviewfactor/unobstructedplanarviewfactor.hpp>
// #include <feel/feelviewfactor/raytracingviewfactor.hpp>

/** use Feel namespace */
using namespace Feel;
using Feel::project;

inline
AboutData
makeAbout()
{
    AboutData about( "test_vf_quadrature" ,
                     "test_vf_quadrature",
                     "0.2",
                     "nD(n=2,3) test vf",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2022-2024 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@cemosis.fr", "" );
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
    opts.add( Feel::feel_options( "square" ) );
    opts.add(feel_options( "cube" ) );
    opts.add(feel_options( "cylinder" ) );
    // opts.add(feel_options( "cube-raytracing" ) );
    // opts.add(feel_options( "cylinder-raytracing" ) );
    return opts;
}

template<typename MeshType>
void checkViewFactorEnclosure(std::string const& prefix)
{   
    auto mesh = loadMesh(_mesh=new MeshType, _filename =Environment::expand(soption(fmt::format("{}.gmsh.filename",prefix) ) )  );
    std::cout << Environment::expand(soption(fmt::format("{}.gmsh.filename",prefix) ) )<< std::endl;
    auto jsons = vsoption( fmt::format("{}.json.filename", prefix ) );
    BOOST_TEST_MESSAGE( fmt::format( "reading json: {}", Environment::expand(jsons[0]) ) );
    nl::json j;
    std::ifstream f(  Environment::expand(jsons[0]) );
    f >> j;
    BOOST_TEST_MESSAGE(fmt::format("j: {}",j.dump(1)));
 
    UnobstructedPlanarViewFactor<MeshType> upvf( mesh, j );
    
    upvf.compute();
    BOOST_TEST_MESSAGE( fmt::format("Max dev reciprocity {}", upvf.maxDevReciprocity()));
    
    BOOST_TEST_MESSAGE( fmt::format("{}", upvf.viewFactors() ) );
    auto row_sum_vf = upvf.viewFactors().rowwise().sum();
    auto exact_vf = eigen_vector_x_col_type<double>::Ones(upvf.viewFactors().rows()) ;
    auto difference_infNorm = (exact_vf-row_sum_vf).template lpNorm<Eigen::Infinity>();
    BOOST_TEST_MESSAGE( fmt::format("View factors sum to one {}; infinity norm of 1 - rowwise sum {}; rowwise sum {} ", difference_infNorm<1e-3, difference_infNorm,row_sum_vf) );
    
    if(prefix=="square")
        BOOST_CHECK_MESSAGE( difference_infNorm<1e-3, fmt::format("Infinity norm of (1 - rowwise) sum is larger than 1e-3" ) );
    else if(prefix=="cube")
    {
        BOOST_CHECK_MESSAGE( difference_infNorm<5e-2, fmt::format("Infinity norm of (1 - rowwise) sum is larger than 4e-2" ) );
    }
    else if(prefix=="cylinder")
    {
        BOOST_CHECK_MESSAGE( difference_infNorm<5e-2, fmt::format("Infinity norm of (1 - rowwise) sum is larger than 4e-2" ) );
    }

    BOOST_CHECK_MESSAGE(upvf.maxDevReciprocity()<1e-6, fmt::format("Max dev reciprocity less than 1e-6, {}",upvf.maxDevReciprocity()));
    double eps = 8e-2;
    if(prefix=="cube")
    {
        auto vf_parallel_walls = view_factor_parallel_walls_exact(1.,1.,1.); // cube of side 1.
        auto vf_perp_walls = view_factor_perp_walls_exact(1.,1.,1.); // cube of side 1.
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(0,1)-vf_parallel_walls)/vf_parallel_walls < eps, fmt::format("Relative error view factors between parallel walls 0 1 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(2,3)-vf_parallel_walls)/vf_parallel_walls < eps, fmt::format("Relative error view factors between parallel walls 2 3 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(4,5)-vf_parallel_walls)/vf_parallel_walls < eps, fmt::format("Relative error view factors between parallel walls 4 5 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(0,2)-vf_perp_walls)/vf_perp_walls < eps, fmt::format("Relative error view factors between perp walls 0 2 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(0,3)-vf_perp_walls)/vf_perp_walls < eps, fmt::format("Relative error view factors between perp walls 0 3 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(0,4)-vf_perp_walls)/vf_perp_walls < eps, fmt::format("Relative error view factors between perp walls 0 4 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(0,5)-vf_perp_walls)/vf_perp_walls < eps, fmt::format("Relative error view factors between perp walls 0 5 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(1,3)-vf_perp_walls)/vf_perp_walls < eps, fmt::format("Relative error view factors between perp walls 1 3 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(1,4)-vf_perp_walls)/vf_perp_walls < eps, fmt::format("Relative error view factors between perp walls 1 4 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(1,5)-vf_perp_walls)/vf_perp_walls < eps, fmt::format("Relative error view factors between perp walls 1 5 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(2,4)-vf_perp_walls)/vf_perp_walls < eps, fmt::format("Relative error view factors between perp walls 2 4 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(2,5)-vf_perp_walls)/vf_perp_walls < eps, fmt::format("Relative error view factors between perp walls 2 5 is less than {}", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(3,5)-vf_perp_walls)/vf_perp_walls < eps, fmt::format("Relative error view factors between perp walls 3 5 is less than {}", eps) );
    }
    else if(prefix=="cylinder")
    {
        auto vf_cylinder = cylinder_vf_map(1.,2.);  
        std::cout <<       upvf.viewFactors() << std::endl;
        // std::cout <<       vf_cylinder["BasisToBasis"] << std::endl;
        // std::cout <<       vf_cylinder["BasisToLateral"]<< std::endl;
        // std::cout <<      vf_cylinder["LateralToBasis"] << std::endl;
        // std::cout <<       vf_cylinder["LateralToLateral"] << std::endl;
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(0,1)-vf_cylinder["BasisToBasis"])/vf_cylinder["BasisToBasis"] < eps, fmt::format("Relative error view factors between bases 0 1 is less than {} ", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(0,2)-vf_cylinder["BasisToLateral"])/vf_cylinder["BasisToLateral"] < eps, fmt::format("Relative error view factors between bottom basis to lateral walls 0 2 is less than {} ", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(1,0)-vf_cylinder["BasisToBasis"])/vf_cylinder["BasisToBasis"] < eps, fmt::format("Relative error view factors between bases 1 0 is less than {} ", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(1,2)-vf_cylinder["BasisToLateral"])/vf_cylinder["BasisToLateral"] < eps, fmt::format("Relative error view factors between bottom basis to lateral walls 1 2 is less than {} ", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(2,0)-vf_cylinder["LateralToBasis"])/vf_cylinder["LateralToBasis"]  < eps, fmt::format("Relative error view factors between lateral wall to top basis 2 0 is less than {} ", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(2,1)-vf_cylinder["LateralToBasis"])/vf_cylinder["LateralToBasis"]  < eps, fmt::format("Relative error view factors between lateral wall to bottom basis 2 1 is less than {} ", eps) );
        BOOST_CHECK_MESSAGE( abs(upvf.viewFactors()(2,2)-vf_cylinder["LateralToLateral"])/vf_cylinder["LateralToLateral"]  < eps, fmt::format("Relative error view factors between lateral wall to itself 2 2 is less than {}: {}",eps, (upvf.viewFactors()(2,2)-vf_cylinder["LateralToLateral"])/vf_cylinder["LateralToLateral"]) );
    }
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( viewfactor_quadrature )


BOOST_AUTO_TEST_CASE( test_square )
{
    checkViewFactorEnclosure<Mesh<Simplex<2>>>("square");
}
BOOST_AUTO_TEST_CASE( test_cube )
{
    checkViewFactorEnclosure<Mesh<Simplex<3>>>( "cube" );
}
BOOST_AUTO_TEST_CASE( test_cylinder )
{
    checkViewFactorEnclosure<Mesh<Simplex<3>>>( "cylinder" );
}

BOOST_AUTO_TEST_SUITE_END()