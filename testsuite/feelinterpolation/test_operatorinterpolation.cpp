
#define BOOST_TEST_MODULE test_operatorinterpolation
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/dh.hpp>

#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>


using namespace Feel;

namespace test_operatorinterpolation
{

typedef Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_operatorinterpolation options" );
    desc_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ;
    return desc_options.add( Feel::feel_options() );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "Test_OperatorInterpolation" ,
                     "Test_OperatorInterpolation" ,
                     "0.1",
                     "test composant of a field vectorial",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;
}


template <uint32_type OrderGeo>
void
test2dTo1d()
{
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_2d_type;

    double meshSize = doption(_name="hsize");
    BOOST_TEST_MESSAGE( "meshSize=" << meshSize );

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 2,1 );
    GeoTool::Rectangle C( meshSize,"OMEGA",x1,x2 );
    C.setMarker( _type="line",_name="Sortie",_markerAll=true );
    C.setMarker( _type="surface",_name="OmegaFluide",_markerAll=true );
    auto mesh2d = C.createMesh( _mesh=new mesh_2d_type,
                                _name="test2dTo1d_domain"+mesh_2d_type::shape_type::name() );

#if 0
    GeoTool::Node x3( 0,0 );
    GeoTool::Node x4( 2,1 );
    GeoTool::Line L( meshSize, "Line",x3,x4 );
    L.setMarker( _type="point",_name="Sortie",_markerAll=true );
    L.setMarker( _type="line",_name="Omega1d",_markerAll=true );
    auto mesh1d = L.createMesh( _mesh=new mesh_1d_type,
                                _name="test2dTo1d_domain1d"+mesh_1d_type::shape_type::name(),_straighten=0 );
#else
    // create the 1d mesh from a submesh of 2d mesh because gmsh don't partion 1d mesh
    GeoTool::Node x3( 0,0 );
    GeoTool::Node x4( 2,1 );
    GeoTool::Node x5( 0,1 );
    GeoTool::Triangle T( meshSize,"OMEGA",x3,x4,x5);
    T.setMarker( _type="line",_name="Boundary1",_marker1=true );
    T.setMarker( _type="line",_name="Boundary2",_marker2=true,_marker3=true );
    T.setMarker( _type="surface",_name="Omega2dFor1d",_markerAll=true );
    auto mesh2dFor1d = T.createMesh( _mesh=new mesh_2d_type,
                                     _name="test2dTo1d_domain2dFor1d"+mesh_2d_type::shape_type::name(), _straighten=0 );
    auto mesh1d = createSubmesh(mesh2dFor1d,markedfaces(mesh2dFor1d,"Boundary1"));
#endif

    //-----------------------------------------------------------//

    auto Xh1d = Pchv<5,PointSetFekete>( mesh1d );
    auto Xh2d = Pchv<3,PointSetFekete>( mesh2d );

    auto exprProj = vec( cos( M_PI*Px() ),sin( M_PI*Py() ) );
    //auto exprProj = vec( 2*Px()*Py(), 0.25*(1.0-Py())*Px() );

    auto u1dinterp = Xh1d->element();
    auto u2dproj = vf::project( _space=Xh2d,
                                _range=elements( mesh2d ),
                                _expr=exprProj );
    auto u1dproj = vf::project( _space=Xh1d,
                                _range=elements( mesh1d ),
                                _expr=exprProj );

    //-----------------------------------------------------------//

    auto opI=opInterpolation( _domainSpace=Xh2d,
                              _imageSpace=Xh1d );
    opI->apply( u2dproj,u1dinterp );

    auto s = integrate( _range=elements( mesh1d ),
                        _expr=inner( idv( u1dinterp )-idv( u1dproj ), idv( u1dinterp )-idv( u1dproj ), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s,1e-8 );
    BOOST_TEST_MESSAGE( "s=" << s );

} // test2dTo1d

//---------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------//

template <uint32_type OrderGeo>
void
test2dTo2d()
{
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;

    //-------------------------------------------------------
    //case 1 : same mesh
    //-------------------------------------------------------
    WorldComm myWorldComm = Environment::worldComm();
    auto meshSize = doption(_name="hsize");
    LOG(INFO) << "meshSize=" << meshSize << "\n";
    BOOST_TEST_MESSAGE( "meshSize=" << meshSize );
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 0.6,0 );
    GeoTool::Circle C( meshSize,"OMEGA",x1,x2 );
    C.setMarker( _type="line",_name="Sortie",_markerAll=true );
    C.setMarker( _type="surface",_name="OmegaFluide",_markerAll=true );
    auto mesh = C.createMesh( _mesh = new mesh_type,
                              _name="test2dTo2d_domain"+mesh_type::shape_type::name(),
                              _partitions=myWorldComm.localSize() );

    auto Xh1 = Pchv<3,PointSetFekete>( mesh );
    auto Xh2 = Pchv<4,PointSetFekete>( mesh );
    auto u1 = Xh1->element();
    auto u2 = Xh2->element();
    auto u2a = Xh2->element();

    auto exprProj = vec( cos( M_PI*Px() ),sin( M_PI*Py() ) );
    u1 = vf::project( _space=Xh1,
                      _range=elements( mesh ),
                      _expr=exprProj );

    auto mybackend = backend(_rebuild=true);

    auto opI = opInterpolation( _domainSpace=Xh1,
                                _imageSpace=Xh2,
                                _backend=mybackend );
    opI->apply( u1,u2 );

    auto s1 = integrate( _range=elements( mesh ),
                         _expr=inner( idv( u1 ) , idv( u1 ),mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()( 0,0 );
    auto s2 = integrate( _range=elements( mesh ),
                         _expr=inner( idv( u2 ) ,idv( u2 ), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( s1, s2,1e-8 );
    BOOST_TEST_MESSAGE( "s1=" << s1 << " s2=" << s2 << " s1-s2 = " << s1-s2 );

    auto opIa=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2,
                               _range=boundaryfaces( mesh ),
                               _backend=mybackend );
    opIa->apply( u1,u2a );

    auto s2a = integrate( _range=boundaryfaces( mesh ),
                          _expr=inner( idv( u1 )-idv( u2a ) , idv( u1 )-idv( u2a ), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s2a,1e-8 );
    BOOST_TEST_MESSAGE( "s2a=" << s2a );

    if ( OrderGeo == 1 )
    {
    //-----------------------------------------------------
    // Lagrange <-> Nedelec
    auto XhNed = Ned1h<0>( mesh );
    auto opINed=opInterpolation( _domainSpace=Xh1,
                                 _imageSpace=XhNed );
    auto uNed = XhNed->element();
    opINed->apply( u1,uNed );
    auto sNed = integrate( _range=elements( mesh ),
                           _expr=inner( idv( uNed ) , idv( uNed ), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( std::abs(s1-sNed),1e-2 );
    BOOST_TEST_MESSAGE( "sNed=" << sNed << "(vs s1=" << s1 << ")" );
    auto opINed2=opInterpolation( _domainSpace=XhNed,
                                  _imageSpace=Xh1 );
    uNed.on(_range=elements(mesh),_expr=vec(-sin(Py()),-sin(Px()) ) );
    auto u1Ned = Xh1->element();
    opINed2->apply( uNed, u1Ned );
    double sNed2 = integrate( _range=elements( mesh ),
                              _expr=inner( idv( u1Ned ) , idv( u1Ned ), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()( 0,0 );
    //BOOST_CHECK_SMALL( std::abs(sNed-sNed2),1e-2 );
    BOOST_TEST_MESSAGE( "sNed2=" << sNed << "(vs sNed=" << sNed << ")" );

    //-----------------------------------------------------
    // Lagrange <-> Raviart-Thomas
    auto XhRT = Dh<0>( mesh );
    auto opIRT = opInterpolation( _domainSpace=Xh1,
                                  _imageSpace=XhRT );
    auto uRT = XhRT->element();
    opIRT->apply( u1,uRT );
    double sRT = integrate( _range=elements( mesh ),
                            _expr=inner( idv( uRT ) , idv( uRT ), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( std::abs(s1-sRT),1e-2 );
    BOOST_TEST_MESSAGE( "sRT=" << sRT << "(vs s1=" << s1 << ")" );
    auto opIRT2 = opInterpolation( _domainSpace=XhRT,
                                   _imageSpace=Xh1 );
    uRT.on(_range=elements(mesh),_expr=vec(-sin(Py()),-sin(Px()) ) );
    auto u1RT = Xh1->element();
    opIRT2->apply( uRT, u1RT );
    double sRT2 = integrate( _range=elements( mesh ),
                              _expr=inner( idv( u1RT ) , idv( u1RT ), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()( 0,0 );
    //BOOST_CHECK_SMALL( std::abs(sRT-sRT2),1e-3 );
    BOOST_TEST_MESSAGE( "sRT2=" << sRT << "(vs sRT=" << sRT << ")" );
    }

    //-------------------------------------------------------//
    //case 2 : with interpolation tool
    //-------------------------------------------------------//
    GeoTool::Node x3(0.4,0);
    GeoTool::Circle C2( meshSize/2.,"OMEGA",x1,x3);
    C2.setMarker(_type="line",_name="Boundary",_markerAll=true);
    C2.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh2 = C2.createMesh(_mesh=new mesh_type,
                               _name="test2dTo2d_domain2"+mesh_type::shape_type::name(),
                               _partitions=myWorldComm.localSize() );
    auto Xh2bis = Pchv<4,PointSetFekete>(mesh2);
    auto u2bis = Xh2bis->element();
    auto u2bisbis = Xh2bis->element();
    auto u2bisproj = vf::project( _space=Xh2bis,
                                  _range=elements( mesh2 ),
                                  _expr=exprProj );
    // opInterp on elements
    auto opI2=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=elements( mesh2 ),
                               _backend=mybackend );
    opI2->apply( u1,u2bis );
    auto s3 = integrate( _range=elements( mesh2 ),
                         _expr=inner( idv( u2bis )-idv( u2bisproj ), idv( u2bis )-idv( u2bisproj ), mpl::int_<InnerProperties::IS_SAME>() )
                         /*_geomap=GeomapStrategyType::GEOMAP_HO*/ ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s3,1e-8 );
    BOOST_TEST_MESSAGE( "s3=" << s3 );

    // opInterp on faces
    auto opI3=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=boundaryfaces( mesh2 ),
                               _backend=mybackend );
    opI3->apply( u1,u2bisbis );
    auto s4 = integrate( _range=boundaryfaces( mesh2 ),
                         _expr=inner( idv( u2bisbis )-idv( u2bisproj ), idv( u2bisbis )-idv( u2bisproj ), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s4,1e-8 );
    BOOST_TEST_MESSAGE( "s4=" << s4 );

} // test2dTo2d

template <uint16_type OrderGeo>
boost::shared_ptr<Mesh<Simplex<2,OrderGeo> > >
buildMeshSMD( mpl::int_<2> /**/)
{
    typedef Mesh<Simplex<2,OrderGeo> > mesh_type;
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 5,1 );
    double meshSize = doption(_name="hsize");
    GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );
    R.setMarker(_type="line",_name="BoundaryInterp",_marker1=true);
    R.setMarker(_type="line",_name="BoundaryOther",_marker2=true,_marker3=true,_marker4=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh(_mesh=new mesh_type,
                             _name=(boost::format("rectangleSMD-%1%")%OrderGeo).str() );
    return mesh;
}
template <uint16_type OrderGeo>
boost::shared_ptr<Mesh<Simplex<3,OrderGeo> > >
buildMeshSMD( mpl::int_<3> /**/)
{
    typedef Mesh<Simplex<3,OrderGeo> > mesh_type;
#if 1
    GeoTool::Node Center(0,0,0);
    GeoTool::Node Radius( 0.5);
    GeoTool::Node Dir(1,0,0);
    GeoTool::Node Lg(2,0,0);
    double meshSize = doption(_name="hsize")*3;
    GeoTool::Cylindre C( meshSize,"Cyl",Center,Dir,Radius,Lg);
    C.setMarker(_type="surface",_name="Inlet",_marker1=true);
    C.setMarker(_type="surface",_name="Outlet",_marker2=true);
    C.setMarker(_type="surface",_name="BoundaryInterp",_marker3=true);
    C.setMarker(_type="volume",_name="Omega",_markerAll=true);
    auto mesh = C.createMesh(_mesh=new mesh_type,
                             _name=(boost::format("cylinderSMD-%1%")%OrderGeo).str(),
                             _hmax=meshSize);
#else
    GeoTool::Node x1( -1,-1,-1 );
    GeoTool::Node x2(  1,-1,-1 );
    GeoTool::Node x3(  0, 1,-1 );
    GeoTool::Node x4(  0, 0, 1 );
    double meshSize = doption(_name="hsize");//1;//2;//
    GeoTool::Tetrahedron R( meshSize,"OMEGA",x1,x2,x3,x4 );
    R.setMarker(_type="surface",_name="BoundaryInterp",_marker1=true);
    R.setMarker(_type="surface",_name="Boundary2",_marker2=true);
    R.setMarker(_type="surface",_name="Boundary3",_marker3=true);
    R.setMarker(_type="surface",_name="Boundary4",_marker4=true);
    R.setMarker(_type="volume",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh(_mesh=new mesh_type,
                             _name="domainTetra",
                             _hmax=meshSize );

#endif

    return mesh;
}

template <uint16_type Dim,uint16_type OrderGeo>
void
testSMD()
{
    BOOST_TEST_MESSAGE( "start test SMD "<< Dim << "d Geo" << OrderGeo );
    auto mesh = buildMeshSMD<OrderGeo>( mpl::int_<Dim>() );
    BOOST_TEST_MESSAGE( "mesh done" );
    auto submesh = createSubmesh( mesh, markedfaces(mesh,"BoundaryInterp") );
    BOOST_TEST_MESSAGE( "submesh done" );
    auto Xh1 = Pch<3/*,PointSetFekete*/>(mesh);
    BOOST_TEST_MESSAGE( "spaces Xh1 done" );
    auto Xh2 = Pch<5/*,PointSetFekete*/>(submesh);
    BOOST_TEST_MESSAGE( "spaces Xh2 done" );
    auto u1 = Xh1->element();
    auto u2 = Xh2->element();

    auto thedispexpr = Px()*(5.0-Px())*0.1;//*oneY();
    //auto thedispexpr = Px();//*(5.0-Px())*0.1;//*oneY();
    auto u1Base = vf::project(_space=Xh1,_expr=thedispexpr);
    auto u2Base = vf::project(_space=Xh2,_expr=thedispexpr);
#if 1
    //-------------------------------------------------------------------------------------//
    auto opIa = opInterpolation(_domainSpace=Xh2,
                               _imageSpace=Xh1,
                               _range=markedfaces(mesh,"BoundaryInterp") );
    opIa->apply(u2Base,u1);
    BOOST_TEST_MESSAGE( "op Xh2->Xh1 done" );

    auto opIb = opInterpolation(_domainSpace=Xh1,
                               _imageSpace=Xh2,
                               _range=elements(submesh) );
    opIb->apply(u1Base,u2);
    BOOST_TEST_MESSAGE( "op Xh1->Xh2 done" );

    //-------------------------------------------------------------------------------------//
    /*double s1 = integrate(_range=elements(submesh),
                          _expr=norm2( idv(u1)-idv(u1Base) ) ).evaluate()(0,0);
    double s2 = integrate(_range=elements(submesh),
                          _expr=norm2( idv(u2)-idv(u2Base) ) ).evaluate()(0,0);
    BOOST_TEST_MESSAGE( "s1=" << s1 );
    BOOST_TEST_MESSAGE( "s2=" << s2 );*/
    double s1a = integrate(_range=elements(submesh),
                           _expr=inner( idv(u1)-idv(u2Base),idv(u1)-idv(u2Base), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()(0,0);
    double s2a = integrate(_range=elements(submesh),
                           _expr=inner( idv(u2)-idv(u1Base),idv(u2)-idv(u1Base), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()(0,0);
    BOOST_TEST_MESSAGE( "s1a=" << s1a );
    BOOST_TEST_MESSAGE( "s2a=" << s2a );
    BOOST_CHECK_SMALL( s1a,1e-16 );
    BOOST_CHECK_SMALL( s2a,1e-16 );
    double s1b = integrate(_range=markedfaces(mesh,"BoundaryInterp"),
                           _expr=inner( idv(u1)-idv(u2Base),idv(u1)-idv(u2Base), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()(0,0);
    double s2b = integrate(_range=markedfaces(mesh,"BoundaryInterp"),
                           _expr=inner( idv(u2)-idv(u1Base),idv(u2)-idv(u1Base), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()(0,0);
    BOOST_TEST_MESSAGE( "s1b=" << s1b );
    BOOST_TEST_MESSAGE( "s2b=" << s2b );
    BOOST_CHECK_SMALL( s1b,1e-16 );
    BOOST_CHECK_SMALL( s2b,1e-16 );
    double s1c = integrate(_range=markedfaces(mesh,"BoundaryInterp"),
                           _expr=inner( idv(u1)-idv(u1Base),idv(u1)-idv(u1Base), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()(0,0);
    double s2c = integrate(_range=elements(submesh),
                           _expr=inner( idv(u2)-idv(u2Base),idv(u2)-idv(u2Base), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()(0,0);
    BOOST_TEST_MESSAGE( "s1c=" << s1c );
    BOOST_TEST_MESSAGE( "s2c=" << s2c );
    BOOST_CHECK_SMALL( s1c,1e-16 );
    BOOST_CHECK_SMALL( s2c,1e-16 );
    double s1d = integrate(_range=markedfaces(mesh,"BoundaryInterp"),
                           _expr=inner( idv(u1Base)-idv(u2Base),idv(u1Base)-idv(u2Base), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()(0,0);
    double s2d = integrate(_range=elements(submesh),
                           _expr=inner( idv(u1Base)-idv(u2Base),idv(u1Base)-idv(u2Base), mpl::int_<InnerProperties::IS_SAME>() ) ).evaluate()(0,0);
    BOOST_TEST_MESSAGE( "s1d=" << s1d );
    BOOST_TEST_MESSAGE( "s2d=" << s2d );
    BOOST_CHECK_SMALL( s1d,1e-16 );
    BOOST_CHECK_SMALL( s2d,1e-16 );
    /*double s1e = integrate(_range=markedfaces(mesh,"BoundaryInterp"),
                           _expr=idv(u1Base) ).evaluate()(0,0);
    double s2e = integrate(_range=elements(submesh),
                           _expr=idv(u2Base) ).evaluate()(0,0);
    BOOST_TEST_MESSAGE( "s1e=" << s1e );
    BOOST_TEST_MESSAGE( "s2e=" << s2e );*/
#endif

#if 0
    for ( auto const& face : markedfaces(mesh,"BoundaryInterp") )
    {
      // subMeshToMesh
      auto elt = submesh->element( submesh->meshToSubMesh(face.id()) );
      for ( uint16_type iloc = 0; iloc < Xh1->dof()->nLocalDofOnFace(true) ; ++iloc )
            for ( uint16_type comp = 0; comp < 1; ++comp )
            {
              size_type i1 =  boost::get<0>( Xh1->dof()->localToGlobal( face, iloc, comp ) );
              auto const imageGlobDofPt = Xh1->dof()->dofPoint( i1 ).template get<0>();
              size_type thelocDofToFind=iloc;
              bool find=false;
              for (uint16_type iloc2 = 0; iloc2 < Xh1->dof()->nLocalDofOnFace() && !find ; ++iloc2)
                {
                  size_type i2b =  boost::get<0>( Xh2->dof()->localToGlobal( elt, iloc2, comp ) );
                  auto const domainGlobDofPt = Xh2->dof()->dofPoint( i2b ).template get<0>();
                  bool find2=true;
                  for (uint16_type d=0;d< Dim;++d)
                    {
                      find2 = find2 && (std::abs( imageGlobDofPt[d]-domainGlobDofPt[d] )<1e-9);
                    }
                  if (find2) { thelocDofToFind=iloc2;find=true; }
                }
              CHECK( find ) << "not find a compatible dof\n ";

              size_type i2 =  boost::get<0>( Xh2->dof()->localToGlobal( elt, thelocDofToFind/*iloc*/, comp ) );
              std::cout << "iloc " << iloc << " u1-u2 " << u1Base(i1) -u2Base(i2) << " u1 " << u1Base(i1) << " u2 " << u2Base(i2)
                        << " " << imageGlobDofPt << " " << Xh2->dof()->dofPoint( i2 ).template get<0>()
                        << std::endl;
            }
    }
#endif

#if 0
    for ( auto const& elt : elements(submesh) )
    {
      for ( uint16_type iloc = 0; iloc < Xh2->dof()->nLocalDof(true) ; ++iloc )
        for ( uint16_type comp = 0; comp < 1; ++comp )
          {
              size_type i1 =  boost::get<0>( Xh2->dof()->localToGlobal( elt, iloc, comp ) );
              std::cout << "iloc " << iloc << " u2 " << u2Base(i1) << " "
                        << Xh2->dof()->dofPoint( i1 ).template get<0>()
                        << std::endl;
          }
    }
#endif
#if 0
    for ( auto const& elt : elements(mesh) )
    {
      for ( uint16_type iloc = 0; iloc < Xh1->dof()->nLocalDof(true) ; ++iloc )
        for ( uint16_type comp = 0; comp < 1; ++comp )
          {
              size_type i1 =  boost::get<0>( Xh1->dof()->localToGlobal( elt, iloc, comp ) );
              std::cout << "iloc " << iloc << " u1 " << u1Base(i1) << " "
                        << Xh1->dof()->dofPoint( i1 ).template get<0>()
                        << std::endl;
          }
    }
#endif
    //-------------------------------------------------------------------------------------//
#if 0
    auto e1 = exporter( _mesh=mesh,_name="Export1" );
    e1->add( "u1", u1 );
    e1->save();
    auto e2 = exporter( _mesh=submesh,_name="Export2" );
    e2->add( "u2", u2 );
    e2->save();
#endif
}

#if defined(FEELPP_HAS_VTK)
template <uint32_type OrderGeo>
void
test2dOpLagrangeP1()
{
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;
    double meshSize = doption(_name="hsize");
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 0.6,0 );
    GeoTool::Circle C( meshSize,"OMEGA",x1,x2 );
    C.setMarker( _type="line",_name="Sortie",_markerAll=true );
    C.setMarker( _type="surface",_name="OmegaFluide",_markerAll=true );
    auto mesh = C.createMesh( _mesh = new mesh_type,
                              _name="test2dOpLagrangeP1_domain"+mesh_type::shape_type::name(),
                              _partitions=Environment::worldComm().localSize() );

    auto Xh = Pchv<3,PointSetFekete>( mesh );
    auto exprProj = vec( cos( M_PI*Px() ),sin( M_PI*Py() ) );
    auto u = vf::project( _space=Xh,
                          _range=elements( mesh ),
                          _expr=exprProj );

    auto mybackend = backend(_rebuild=true);

    auto opLagP1 = lagrangeP1(_space=Xh);
    auto meshLagP1 = opLagP1->mesh();

    auto XhLagP1 = Pchv<1>( meshLagP1 );
    auto uLagP1interp = XhLagP1->element();

    auto uLagP1proj = vf::project( _space=XhLagP1,
                          _range=elements( meshLagP1 ),
                          _expr=exprProj );

    auto opI=opInterpolation( _domainSpace=Xh,
                              _imageSpace=XhLagP1,
                              _range=elements( meshLagP1 ) );
    opI->apply( u,uLagP1interp );

    auto s1 = integrate(_range=elements(meshLagP1),
                        _expr=inner( idv(uLagP1interp)-idv(uLagP1proj) , idv(uLagP1interp)-idv(uLagP1proj) ) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s1,1e-6);

#if 0
    auto myexporter = exporter( _mesh=meshLagP1, _name="test2dOpLagrangeP1_MyExport" );
    myexporter->add( "test2dOpLagrangeP1_uLagP1", uLagP1 );
    myexporter->save();
#endif
}


#endif // FEELPP_HAS_VTK

} // namespace test_operatorinterpolation

/**
 * main code
 */
FEELPP_ENVIRONMENT_WITH_OPTIONS( test_operatorinterpolation::makeAbout(),
                                 test_operatorinterpolation::makeOptions() )

BOOST_AUTO_TEST_SUITE( interp_operatorinterpolation )

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_2d_2d_geo1 )
{
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_2d_geo1" );
    using namespace test_operatorinterpolation;

    Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/2d2dgeo1" )
                                   % Environment::about().appName() );

    test_operatorinterpolation::test2dTo2d<1>();
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_2d_geo1 done" );
}

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_2d_2d_geo2 )
{
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_2d_geo2" );
    using namespace test_operatorinterpolation;

    Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/2d2dgeo2" )
                                   % Environment::about().appName() );

    test_operatorinterpolation::test2dTo2d<2>();
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_2d_geo2 done" );
}

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_2d_1d_geo1 )
{
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_1d_geo1" );
    using namespace test_operatorinterpolation;

    Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/2d1dgeo1" )
                                   % Environment::about().appName() );

    //test_operatorinterpolation::test2dTo2d<3>(test_app);
    test_operatorinterpolation::test2dTo1d<1>();
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_1d_geo1 done" );
}

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_2d_1d_geo2 )
{
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_1d_geo2" );
    using namespace test_operatorinterpolation;

    Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/2d1dgeo2" )
                                   % Environment::about().appName() );

    test_operatorinterpolation::test2dTo1d<2>();
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_1d_geo2 done" );
}

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_smd )
{
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_smd" );
    using namespace test_operatorinterpolation;

    Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/3d3dgeo1smd" )
                                   % Environment::about().appName() );

    //test_operatorinterpolation::testSMD<2,1>();
    test_operatorinterpolation::testSMD<3,1>();
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_smd done" );
}

#if defined(FEELPP_HAS_VTK)
BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_oplagp1_geo1 )
{
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_oplagp1_geo1" );
    using namespace test_operatorinterpolation;

    Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/oplagp1geo1" )
                                   % Environment::about().appName() );


    test_operatorinterpolation::test2dOpLagrangeP1<1>();
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_oplagp1_geo1 done" );
}
#endif // FEELPP_HAS_VTK

BOOST_AUTO_TEST_SUITE_END()

