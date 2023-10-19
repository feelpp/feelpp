/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-*/

#define BOOST_TEST_MODULE test_markerhdf5
#include <feel/feelcore/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace test_markerhdf5
{

using namespace Feel;
using namespace Feel::vf;

typedef Mesh<Simplex<CHECKH5_DIM>> mesh_type;
typedef std::shared_ptr<mesh_type> mesh_ptrtype;

void checkMarker(mesh_ptrtype mesh)
{
    std::map<std::string,double> M_map;
    M_map["elements"] =  integrate( _range = elements( mesh ),_expr = cst(1.) ).evaluate()(0,0);
    M_map["markedelements(Omega)"] =  integrate( _range = markedelements( mesh,"Omega" ),_expr = cst(1.) ).evaluate()(0,0);
    M_map["markedFaces(Border)"] =  integrate( _range = markedfaces( mesh, "Border" ),_expr = cst(1.) ).evaluate()(0,0);
    M_map["boundaryFaces"] =  integrate( _range = boundaryfaces( mesh ),_expr = cst(1.) ).evaluate()(0,0);

    BOOST_CHECK_EQUAL(M_map["elements"], M_map["markedelements(Omega)"]);
    BOOST_CHECK_EQUAL(M_map["markedFaces(Border)"], M_map["boundaryFaces"]);
}
void
run()
{
    fs::path p(
#if CHECKH5_DIM==2
        "markerhdf5_2D.geo"
#else
        "markerhdf5_3D.geo"
#endif
        );
    /// Load mesh given by gmsh.filename
    std::cout << "Check with " << p << std::endl;
    auto mesh1 = loadMesh(_mesh = new mesh_type, _filename=p.string(),_savehdf5=true);
    checkMarker(mesh1);
    /// Load the generated msh file
    p.replace_extension("msh");
    std::cout << "Check with " << p << std::endl;
    auto mesh2 = loadMesh(_mesh = new mesh_type, _filename=(boost::format( "%1%/%2%" )
                                                                          % fs::current_path().string()
                                                                          % p.string()
                                                                          ).str());
    checkMarker(mesh2);
    /// Load the generated json metadata file
    p.replace_extension("json");
    std::cout << "Check with " << p << std::endl;
    auto mesh3 = loadMesh(_mesh = new mesh_type, _filename=(boost::format( "%1%/%2%" )
                                                                          % fs::current_path().string()
                                                                          % p.string()
                                                                          ).str());
    checkMarker(mesh3);
} // run

} //namespace test_matching

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::about(Feel::_name="test_markerHDF5", //strcat(strcat((char*)"test_markerHDF5_",TOSTRING(CHECKH5_DIM)),(char*)"D"),
                                             Feel::_author="Feel++ Consortium",
                                             Feel::_email="feelpp-devel@feelpp.org"),
                                 Feel::feel_options().add( Feel::backend_options("test_markerhdf5") ) )
//FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( markerhdf5 )

typedef Feel::Application Application_type;
typedef std::shared_ptr<Application_type> Application_ptrtype;

BOOST_AUTO_TEST_CASE( markerhdf5 )
{
    //auto theApp = Application_ptrtype( new Application_type );

    test_markerhdf5::run( /*theApp*/ );

}
BOOST_AUTO_TEST_SUITE_END()
