/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-*/

#define BOOST_TEST_MODULE test_matching
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feeldiscr/pch.hpp>

namespace test_matching
{

using namespace Feel;
using namespace Feel::vf;

typedef Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

void run( Application_ptrtype & theApp )
{
    theApp->changeRepository( boost::format( "testsuite/feeldiscr/%1%/h_%2%/" )
                              % theApp->about().appName()
                              % option(_name="gmsh.hsize2").as<double>() );

    auto mesh1 = createGMSHMesh( _mesh=new Mesh<Hypercube<2,1,2> >,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                 _desc=domain( _name="mesh1", _shape="hypercube",
                                               _dim=2, _h=option(_name="gmsh.hsize2"). as<double>(),
                                               _convex="Hypercube",_structured=2,
                                               _xmin=0., _xmax=1., _ymin=0., _ymax=1. ) );

    auto mesh2 = createGMSHMesh( _mesh=new Mesh<Hypercube<2,1,2> >,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                 _desc=domain( _name="mesh2", _shape="hypercube",
                                               _dim=2, _h=option(_name="gmsh.hsize2"). as<double>(),
                                               _convex="Hypercube",_structured=2,
                                               _xmin=0., _xmax=1., _ymin=1., _ymax=2. ) );

    auto mesh_1 = createSubmesh(mesh1, markedfaces(mesh1,(boost::any)4),Environment::worldComm() );
    auto mesh_2 = createSubmesh(mesh2, markedfaces(mesh2,(boost::any)2),Environment::worldComm() );

    auto Xh = Pch<1>(mesh_1);
    auto Yh = Pch<1>(mesh_2);

    auto u = Xh->element();
    auto v = Yh->element();

    auto backend = backend_type::build( theApp->vm() );
    auto M = backend->newMatrix( _test=Xh, _trial=Yh );
    form2( _trial=Xh, _test=Yh, _matrix=M ) = integrate( _range=elements(mesh_2), _expr=idt(u)*id(v) );

    auto g = Px()+Py();
    auto gproj = vf::project( _space=Xh, _range=elements( mesh_1 ),_expr=g );
    auto F = backend->newVector( _test=Yh );
    form1( _test=Yh, _vector=F ) = integrate( _range=elements(mesh_2), _expr=g*id(v) );

    backend->solve(_matrix=M, _solution=u, _rhs=F,_pcfactormatsolverpackage="mumps");

    auto expo = exporter(_mesh=mesh_1, _name="Exporter");
    expo->add( "g", gproj );
    expo->add( "u", u );
    expo->save();

} // run

} //namespace test_matching

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( matching )

typedef Feel::Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

BOOST_AUTO_TEST_CASE( matching )
{
    auto theApp = Application_ptrtype( new Application_type );

    test_matching::run( theApp );

}
BOOST_AUTO_TEST_SUITE_END()
