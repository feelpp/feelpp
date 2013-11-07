/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-*/

#define TWOLAP_USE_BOOST_TEST 0

#if TWOLAP_USE_BOOST_TEST
#define BOOST_TEST_MODULE test_twolaplaciansdistributed
#include <testsuite/testsuite.hpp>
#endif

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>

namespace test_twolaplaciansdistributed
{

using namespace Feel;
using namespace Feel::vf;

typedef Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_twolaplaciansdistributed options" );
    desc_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.02 ), "mesh size" )
    ( "1d-hsize", po::value<double>()->default_value( 0.02 ), "mesh size 1d" )
    ;
    return desc_options.add( Feel::feel_options() ).add( backend_options( "backend1" ) );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "Test_Twolaplaciansdistributed" ,
                     "Test_Twolaplaciansdistributed" ,
                     "0.1",
                     "test bilinear and linear form for several trial/test mesh ",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;
}


/*_________________________________________________*
 * About
 *_________________________________________________*/

void run( Application_ptrtype & theApp )
{
    typedef Mesh< Simplex<2,1,2> > mesh_type;
    double meshSize = theApp->vm()["hsize"].as<double>();

    const int domainColor1=1;
    const int domainColor2=2;
    const int nTotalProc = Environment::worldComm().size();

    auto worldCommDomain1 = Environment::worldComm();
    auto worldCommDomain2 = Environment::worldComm();

    if ( nTotalProc>1)
        {
            std::vector<int> MapWorld(nTotalProc);
            if (nTotalProc==1) // a particular test case
                {
                    MapWorld[0]=domainColor2;
                    MapWorld[1]=domainColor1;
                    MapWorld[2]=domainColor2;
                    MapWorld[3]=domainColor2;
                    MapWorld[4]=domainColor1;
                }
            else
                {
                    for (int proc = 0 ; proc < nTotalProc; ++proc)
                        {
                            if (proc < nTotalProc/2 )
                                MapWorld[proc] = domainColor1;
                            else
                                MapWorld[proc] = domainColor2;
                        }
                }
            auto worldCommColored = WorldComm(MapWorld);
            //worldCommColored.showMe();
            worldCommDomain1 = worldCommColored.subWorldComm(domainColor1);
            worldCommDomain2 = worldCommColored.subWorldComm(domainColor2);
        }
    //worldCommDomain1.showMe();
    //worldCommDomain2.showMe();
    //Environment::worldComm().showMe();
    //---------------------------------------------------------------------------------------//
    // meshes
    GeoTool::Node x11( 0,0 );
    GeoTool::Node x12( 1,1 );
    GeoTool::Rectangle Omega1( meshSize,"Omega",x11,x12 );
    Omega1.setMarker( _type="line",_name="Paroi",_markerAll=true );
    Omega1.setMarker( _type="surface",_name="Omega",_markerAll=true );
    auto mesh1 = Omega1.createMesh(_mesh=new mesh_type,
                                   _name ="omega1_"+ mesh_type::shape_type::name(),
                                   _partitions=worldCommDomain1.localSize(),
                                   _worldcomm=worldCommDomain1 );

    GeoTool::Node x21( -1,0.3 );
    GeoTool::Node x22( -1+0.5,0.3 );
    GeoTool::Circle Omega2( meshSize,"Omega",x21,x22 );
    Omega2.setMarker( _type="line",_name="Paroi",_markerAll=true );
    Omega2.setMarker( _type="surface",_name="Omega",_markerAll=true );
    auto mesh2 = Omega2.createMesh(_mesh=new mesh_type,
                                   _name ="omega2_"+ mesh_type::shape_type::name(),
                                   _partitions=worldCommDomain2.localSize(),
                                   _worldcomm=worldCommDomain2 );
    //---------------------------------------------------------------------------------------//
    // functionspaces
    auto Xh1 = Pch<1>( mesh1 );
    auto Xh2 = Pch<1>( mesh2 );
    auto u1 = Xh1->element();
    auto u2 = Xh2->element();
    //---------------------------------------------------------------------------------------//
    // backends
    auto backend1 = backend(_rebuild=true,_worldcomm=Xh1->worldComm());
    auto backend2 = backend(_rebuild=true,_worldcomm=Xh2->worldComm());
    //---------------------------------------------------------------------------------------//
    // init matrix and vectors
    auto A1 = backend1->newMatrix(_test=Xh1,_trial=Xh1);
    auto A2 = backend2->newMatrix(_test=Xh2,_trial=Xh2);
    auto F1 = backend1->newVector(Xh1);
    auto F2 = backend2->newVector(Xh2);
    //---------------------------------------------------------------------------------------//
    // assembly
    form2( _test=Xh1, _trial=Xh1, _matrix=A1 ) +=
        integrate( _range=elements(mesh1), _expr=gradt(u1)*trans(grad(u1)) );
    form1( _test=Xh1, _vector=F1 ) +=
        integrate( _range=elements(mesh1), _expr=id(u1) );

    if (Xh1->worldComm().isActive())  //(todo vincent : simplifier)
        {
            //A1->close();F1->close();
            form2( _test=Xh1, _trial=Xh1, _matrix=A1 ) +=
                on( _range=boundaryfaces(mesh1),
                    _element=u1,_rhs=F1,
                    _expr=cst(0.) );
        }

    form2( _test=Xh2, _trial=Xh2, _matrix=A2 ) +=
        integrate( _range=elements(mesh2), _expr=gradt(u2)*trans(grad(u2)) );
    form1( _test=Xh2, _vector=F2 ) +=
        integrate( _range=elements(mesh2), _expr=id(u2) );

    if (Xh2->worldComm().isActive())  //(todo vincent : simplifier)
        {
            //A2->close();F2->close();
            form2( _test=Xh2, _trial=Xh2, _matrix=A2 ) +=
                on( _range=boundaryfaces(mesh2),
                    _element=u2,_rhs=F2,
                    _expr=cst(0.) );
        }

    //---------------------------------------------------------------------------------------//
    // solve (todo vincent : simplifier (+rajouter le worldcomm automatiquement dans prec defaut de la fonction solve et nlsolve  )
    if (Xh1->worldComm().isActive())
        {
#if 0
            auto prec1 = preconditioner(_pc=(PreconditionerType) backend1->pcEnumType() /*LU_PRECOND*/,
                                        _matrix=A1,
                                        _backend= backend1,
                                        _pcfactormatsolverpackage=(MatSolverPackageType) backend1->matSolverPackageEnumType(),
                                        _worldcomm=backend1->comm(),
                                        _prefix=backend1->prefix() );
#endif
            backend1->solve(_matrix=A1,_solution=u1,_rhs=F1/*,_prec=prec1*/);
        }
    if (Xh2->worldComm().isActive())
        {
#if 0
            auto prec2 = preconditioner(_pc=(PreconditionerType) backend2->pcEnumType() /*LU_PRECOND*/,
                                        _matrix=A2,
                                        _backend= backend2,
                                        _pcfactormatsolverpackage=(MatSolverPackageType) backend2->matSolverPackageEnumType(),
                                        _worldcomm=backend2->comm(),
                                        _prefix=backend2->prefix() );
#endif
            backend2->solve(_matrix=A2,_solution=u2,_rhs=F2/*,_prec=prec2*/);
        }

    //---------------------------------------------------------------------------------------//
    // exports
    auto myexporter1 = exporter( _mesh=mesh1, _name="MyExportDomain1" );
    myexporter1->step(0)->add( "u1", u1 );
    myexporter1->save();
    auto myexporter2 = exporter( _mesh=mesh2, _name="MyExportDomain2" );
    myexporter2->step(0)->add( "u2", u2 );
    myexporter2->save();

}

} //namespace test_twolaplaciansdistributed


/*_________________________________________________*
 * Main
 *_________________________________________________*/

#if TWOLAP_USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_twolaplaciansdistributed::makeAbout(),
                                 test_twolaplaciansdistributed::makeOptions() )

BOOST_AUTO_TEST_SUITE( twolaplaciansdistributed )

typedef Feel::Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;


BOOST_AUTO_TEST_CASE( twolaplaciansdistributed1 )
{
    auto theApp = Application_ptrtype( new Application_type );

    if ( theApp->vm().count( "help" ) )
    {
        std::cout << theApp->optionsDescription() << "\n";
        exit( 0 );
    }

    theApp->changeRepository( boost::format( "/testsuite/feeldiscr/%1%/" )
                              % theApp->about().appName() );

    test_twolaplaciansdistributed::run( theApp );

}
BOOST_AUTO_TEST_SUITE_END()


#else
int
main( int argc, char** argv )
{
    typedef Feel::Application Application_type;
    typedef boost::shared_ptr<Application_type> Application_ptrtype;
    Feel::Environment env( Feel::_argc=argc,Feel::_argv=argv,
                           Feel::_about=test_twolaplaciansdistributed::makeAbout(),
                           Feel::_desc=test_twolaplaciansdistributed::makeOptions() );
    auto theApp = Application_ptrtype( new Application_type );



    if ( theApp->vm().count( "help" ) )
    {
        std::cout << theApp->optionsDescription() << "\n";
        exit( 0 );
    }

    theApp->changeRepository( boost::format( "/testsuite/feeldiscr/%1%/" )
                              % theApp->about().appName() );

    test_twolaplaciansdistributed::run( theApp );

}
#endif

