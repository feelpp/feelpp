#define BOOST_TEST_MODULE test_exporter_disc
#include <feel/feelcore/testsuite.hpp>


#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feelmesh/meshmover.hpp>

#if __has_include(<vtkEnSightGoldBinaryReader.h>)
# define have_vtk_ensightgoldbinaryreader_h
#include <vtkEnSightGoldBinaryReader.h>
#endif

/** use Feel namespace */
using namespace Feel;


FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( exporter_disc )


BOOST_AUTO_TEST_CASE( test_1 )
{
    auto mesh = unitSquare();

    auto e = exporter(_mesh=mesh,
                      _path=(fs::path(Environment::exportsRepository())/fs::path("test_1")/fs::path(soption("exporter.format"))).string() );
    auto it_ets = e->beginTimeSet();
    auto en_ets = e->endTimeSet();
    BOOST_CHECK( std::distance(it_ets,en_ets) == 1 );
    //auto ets = *it_ets;

    auto XhDisp =Pchv<1>( mesh );
    auto disp = XhDisp->elementPtr();

    std::vector<std::vector<double>> timeByTimeSet(2);
    double t=3;
    for (;t<3.9;t+=0.1 )
        timeByTimeSet[0].push_back( t );
    for (;t<4.6;t+=0.1 )
        timeByTimeSet[1].push_back( t );

    for ( int tsIndex=0, k=0;k<1/*2*/;++k )
    {
        if ( k > 0 )
        {
            tsIndex = e->addTimeSet();
            e->timeSet(tsIndex)->setMesh( mesh );
        }
        for ( double t : timeByTimeSet[k] )
        {
            disp->scale(-1);
            meshMove( mesh, *disp );
            disp->on(_range=elements(mesh),_expr=(t-timeByTimeSet[0][0])*P() );
            meshMove( mesh, *disp );
            //std::cout << "k=" << k << " tsIndex="<<tsIndex << " t=" << t << std::endl;
            auto currentStep = e->step(t,tsIndex);
            BOOST_CHECK( !currentStep->isOnDisk() );
            for ( auto it_step = e->timeSet(tsIndex)->beginStep(); it_step != e->timeSet(tsIndex)->endStep(); ++ it_step )
            {
                auto s = *it_step;
                CHECK( s->hasData() );
                if ( s->index() < currentStep->index() )
                {
                    BOOST_CHECK( !s->isInMemory() );
                    BOOST_CHECK( s->isOnDisk() );
                }
                else
                    BOOST_CHECK( !s->isOnDisk() );
            }

            e->step(t)->add( "disp", disp, std::set<std::string>{"nodal","element"} );
            e->step(t)->add( "v", t*Px(), elements(mesh), std::set<std::string>{"nodal","element"} );

            BOOST_CHECK( currentStep->isInMemory() );
            BOOST_CHECK( !currentStep->isOnDisk() );

            e->save();
            BOOST_CHECK( !currentStep->isInMemory() );
            BOOST_CHECK( currentStep->isOnDisk() );
        }
    }

}

typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( test_2, T, dim_types )
{
    static const uint16_type nDim = T::value;
    typedef Mesh<Simplex<nDim> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    //auto mesh = loadMesh( _mesh=new mesh_type);

    mesh_ptrtype mesh;
    double meshSize = doption(_name="gmsh.hsize");
    bool keepInterface = true;
    if constexpr ( nDim == 2 )
    {
        GeoTool::Node x1a( 0,0 );
        GeoTool::Node x2a( 0.5,1 );
        GeoTool::Rectangle Ra( meshSize,"Omega1",x1a,x2a );
        Ra.setMarker(_type="line",_name="Boundary1",_marker1=true,_marker3=true,_marker4=true);
        Ra.setMarker(_type="line",_name="InternalInterface",_marker2=true);
        Ra.setMarker(_type="surface",_name="Omega1",_markerAll=true);

        GeoTool::Node x1b( 0.5,0 );
        GeoTool::Node x2b( 1.5,1 );
        GeoTool::Rectangle Rb( meshSize,"Omega2",x1b,x2b );
        Rb.setMarker(_type="line",_name="Boundary1",_marker1=true,_marker3=true);
        Rb.setMarker(_type="line",_name="InternalInterface",_marker2=true,_marker4=true);
        Rb.setMarker(_type="surface",_name="Omega2",_markerAll=true);

        GeoTool::Node x1e( 1.5,0 );
        GeoTool::Node x2e( 2,1 );
        GeoTool::Rectangle Rc( meshSize,"Omega3",x1e,x2e );
        Rc.setMarker(_type="line",_name="Boundary1",_marker1=true,_marker3=true,_marker2=true);
        Rc.setMarker(_type="line",_name="InternalInterface",_marker4=true);
        Rc.setMarker(_type="surface",_name="Omega3",_markerAll=true);

        mesh = (Ra+Rb+Rc).
            fusion(Ra,2,Rb,4,keepInterface).
            fusion(Rb,2,Rc,4,keepInterface).
            createMesh(_mesh=new mesh_type,
                       _name="test_2_domain2d" );
    }
    else
    {
        GeoTool::Node x1a(0,0,0);
        GeoTool::Node x2a(0.5,1,0.5);
        GeoTool::Cube Ca( meshSize,"Cube1",x1a,x2a);
        Ca.setMarker(_type="surface",_name="Boundary1",_marker1=true,_marker2=true,_marker3=true,_marker5=true,_marker6=true );
        Ca.setMarker(_type="surface",_name="InternalInterface_1_2",_marker4=true);
        Ca.setMarker(_type="volume",_name="Omega1",_markerAll=true);

        GeoTool::Node x1b(0.5,0,0);
        GeoTool::Node x2b(1.5,1,0.5);
        GeoTool::Cube Cb( meshSize,"Cube2",x1b,x2b);
        Cb.setMarker(_type="surface",_name="Boundary1",_marker1=true,_marker2=true,_marker3=true,_marker5=true );
        Cb.setMarker(_type="surface",_name="InternalInterface_2_3",_marker4=true);
        Cb.setMarker(_type="surface",_name="InternalInterface_1_2",_marker6=true);
        Cb.setMarker(_type="volume",_name="Omega2",_markerAll=true);

        GeoTool::Node x1c(1.5,0,0);
        GeoTool::Node x2c(2,1,0.5);
        GeoTool::Cube Cc( meshSize,"Cube3",x1c,x2c);
        Cc.setMarker(_type="surface",_name="Boundary1",_marker1=true,_marker2=true,_marker3=true,_marker5=true,_marker4=true );
        Cc.setMarker(_type="surface",_name="InternalInterface_2_3",_marker6=true);
        Cc.setMarker(_type="volume",_name="Omega3",_markerAll=true);

        mesh = (Ca+Cb+Cc).
            fusion(Ca,4,Cb,6,keepInterface).
            fusion(Cb,4,Cc,6,keepInterface).
            createMesh(_mesh=new mesh_type,
                       _name="test_2_domain3d" );
    }

    auto VhScalar = Pch<2>( mesh );
    auto uScalar = VhScalar->element();
    auto VhVectorial = Pchv<2>( mesh );
    auto uVectorial = VhVectorial->element();
    auto VhTensor2 = Pchm<2>( mesh );
    auto uTensor2 = VhTensor2->element();
    auto VhTensor2Symm = Pchms<2>( mesh );
    auto uTensor2Symm = VhTensor2Symm->element();
    auto WhScalar = Pdh<2>( mesh );
    auto wScalar = WhScalar->element();
    auto WhVectorial = Pdhv<2>( mesh );
    auto wVectorial = WhVectorial->element();
    auto WhTensor2 = Pdhm<2>( mesh );
    auto wTensor2 = WhTensor2->element();
    auto WhTensor2Symm = Pdhms<2>( mesh );
    auto wTensor2Symm = WhTensor2Symm->element();

    uScalar.on(_range=elements(mesh),_expr=Px()*Py());
    wScalar.on(_range=elements(mesh),_expr=Px()*Py());

    if constexpr ( nDim == 2 )
    {
        uVectorial.on(_range=elements(mesh),_expr=vec(Px(),Py()));
        uTensor2.on(_range=elements(mesh),_expr=mat<2,2>(Px(),Px()*Py(),-Px()*Py(),Py()) );
        uTensor2Symm.on(_range=elements(mesh),_expr=mat<2,2>(Px(),Px()*Py(),Px()*Py(),Py()) );
        wVectorial.on(_range=elements(mesh),_expr=vec(Px(),Py()));
        wTensor2.on(_range=elements(mesh),_expr=mat<2,2>(Px(),Px()*Py(),-Px()*Py(),Py()) );
        wTensor2Symm.on(_range=elements(mesh),_expr=mat<2,2>(Px(),Px()*Py(),Px()*Py(),Py()) );
    }
    else
    {
        uVectorial.on(_range=elements(mesh),_expr=vec(Px(),Py(),Pz()));
        uTensor2.on(_range=elements(mesh),_expr=mat<3,3>(Px(),Px()*Py(),Px()*Pz(),
                                                   -Px()*Py(),Py(),Py()*Pz(),
                                                   -Px()*Pz(),-Py()*Pz(),Pz()
                                                   ) );
        uTensor2Symm.on(_range=elements(mesh),_expr=mat<3,3>(Px(),Px()*Py(),Px()*Pz(),
                                                             Px()*Py(),Py(),Py()*Pz(),
                                                             Px()*Pz(),Py()*Pz(),Pz()
                                                             ) );
        wVectorial.on(_range=elements(mesh),_expr=vec(Px(),Py(),Pz()));
        wTensor2.on(_range=elements(mesh),_expr=mat<3,3>(Px(),Px()*Py(),Px()*Pz(),
                                                   -Px()*Py(),Py(),Py()*Pz(),
                                                   -Px()*Pz(),-Py()*Pz(),Pz()
                                                   ) );
        wTensor2Symm.on(_range=elements(mesh),_expr=mat<3,3>(Px(),Px()*Py(),Px()*Pz(),
                                                             Px()*Py(),Py(),Py()*Pz(),
                                                             Px()*Pz(),Py()*Pz(),Pz()
                                                             ) );
    }

    auto e = exporter(_mesh=mesh,
                      _path=(fs::path(Environment::exportsRepository())/fs::path("test_2")/fs::path(soption("exporter.format"))/std::to_string(nDim) ).string() );
    //    auto e = exporter( _mesh = mesh );
    e->addRegions();
    e->add( "uScalar", uScalar );
    e->add( "uVectorial", uVectorial );
    //e->add( "uTensor2", uTensor2 );
    e->add( "uTensor2Symm", uTensor2Symm );

    e->add( "wScalar", wScalar );
    e->add( "wVectorial", wVectorial );
    //e->add( "wTensor2", wTensor2 );
    e->add( "wTensor2Symm", wTensor2Symm );

    e->add( "expr1Scalar", inner(P()), markedelements(mesh,"Omega1") );
    e->add( "expr1Scalar", inner(P()), markedelements(mesh,"Omega2") );
    e->add( "expr1Scalar", inner(P()), markedelements(mesh,"Omega3") );

    e->add( "expr2Scalar", cst(1.), markedelements(mesh,"Omega1"), "element" );
    e->add( "expr2Scalar", cst(2.), markedelements(mesh,"Omega2"), "element" );
    e->add( "expr2Scalar", cst(3.), markedelements(mesh,"Omega3"), "element" );

    e->add( "expr3Vectorial", P(), std::set<std::string>{"nodal","element"} );
    e->save();


#if defined(FEELPP_HAS_VTK) && defined(have_vtk_ensightgoldbinaryreader_h)
    Environment::worldComm().barrier();
    if ( Environment::isMasterRank() )
    {
        vtkSmartPointer<vtkEnSightGoldBinaryReader> reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
        reader->SetCaseFileName( (fs::path(e->path())/ (e->prefix() + ".case")).string().c_str() );
        reader->Update();
        BOOST_CHECK_EQUAL( reader->GetNumberOfVariables() , 11 );

        vtkSmartPointer<vtkMultiBlockDataSet> mbds;
        mbds = reader->GetOutput();
        BOOST_CHECK_EQUAL(  mbds->GetNumberOfBlocks() , 3 );

        std::set<std::string> eltMarkerNames = { "Omega1", "Omega2", "Omega3" };
        double tol = 1e-4;
        for (int k=0;k< mbds->GetNumberOfBlocks();++k)
        {
            std::string blockName = mbds->GetMetaData(k)->Get(vtkCompositeDataSet::NAME() );
            BOOST_CHECK( eltMarkerNames.find( blockName ) != eltMarkerNames.end() );
            vtkUnstructuredGrid* grid = vtkUnstructuredGrid::SafeDownCast( mbds->GetBlock( k ) );

            vtkPoints * points = grid->GetPoints();
            vtkPointData * pointData = grid->GetPointData();
            vtkCellData * cellData = grid->GetCellData();

            int nPoint = points->GetNumberOfPoints();
            int nCell =  grid->GetNumberOfCells();

            BOOST_CHECK( pointData->HasArray("uScalar") );
            vtkDataArray * pointDataScalar = pointData->GetScalars("uScalar");
            BOOST_CHECK_EQUAL( pointDataScalar->GetNumberOfComponents() , 1 );
            BOOST_CHECK_EQUAL( pointDataScalar->GetDataSize() , nPoint );
            BOOST_CHECK( pointData->HasArray("uVectorial") );
            vtkDataArray * pointDataVector = pointData->GetVectors("uVectorial");
            BOOST_CHECK_EQUAL( pointDataVector->GetNumberOfComponents() , 3 );
            BOOST_CHECK_EQUAL( pointDataVector->GetDataSize() , 3*nPoint );
            BOOST_CHECK( pointData->HasArray("uTensor2Symm") );
            vtkDataArray * pointDataTensor2Symm = pointData->GetTensors( "uTensor2Symm" );
            BOOST_CHECK_EQUAL( pointDataTensor2Symm->GetNumberOfComponents() , 6 );
            BOOST_CHECK_EQUAL( pointDataTensor2Symm->GetDataSize() , 6*nPoint );

            BOOST_CHECK( pointData->HasArray("expr1Scalar") );
            vtkDataArray * pointDataExpr1Scalar = pointData->GetScalars("expr1Scalar");
            BOOST_CHECK_EQUAL( pointDataExpr1Scalar->GetNumberOfComponents() , 1 );
            BOOST_CHECK_EQUAL( pointDataExpr1Scalar->GetDataSize() , nPoint );
            BOOST_CHECK( pointData->HasArray("expr3Vectorial_n") );
            vtkDataArray * pointDataExpr3Vector = pointData->GetVectors("expr3Vectorial_n");
            BOOST_CHECK_EQUAL( pointDataExpr3Vector->GetNumberOfComponents() , 3 );
            BOOST_CHECK_EQUAL( pointDataExpr3Vector->GetDataSize() , 3*nPoint );


            double coord[3];
            for ( vtkIdType p=0;p<nPoint;++p )
            {
                points->GetPoint ( p, coord);

                double scalarValue = pointDataScalar->GetTuple1( p );
                BOOST_CHECK_SMALL( scalarValue - coord[0]*coord[1], tol );

                double * vectorialValue = pointDataVector->GetTuple3( p );
                for ( int c=0;c<nDim;++c )
                    BOOST_CHECK_SMALL( vectorialValue[c] - coord[c], tol );
                for ( int c=nDim;c<3;++c )
                    BOOST_CHECK_SMALL( vectorialValue[c], tol );

                double * tensor2SymmValue = pointDataTensor2Symm->GetTuple6( p );
                for ( int c=0;c<nDim;++c )
                    BOOST_CHECK_SMALL( tensor2SymmValue[c] - coord[c], tol );
                for ( int c=nDim;c<3;++c )
                    BOOST_CHECK_SMALL( tensor2SymmValue[c], tol );
                BOOST_CHECK_SMALL( tensor2SymmValue[3] - coord[0]*coord[1], tol );
                if ( nDim == 3 )
                {
                    BOOST_CHECK_SMALL( tensor2SymmValue[4] - coord[1]*coord[2], tol );
                    BOOST_CHECK_SMALL( tensor2SymmValue[5] - coord[0]*coord[2], tol );
                }

                 double scalarValueExpr1 = pointDataExpr1Scalar->GetTuple1( p );
                 if (nDim == 2 )
                     BOOST_CHECK_SMALL( scalarValueExpr1 - ( std::pow(coord[0],2) + std::pow(coord[1],2) ), tol );
                 else
                     BOOST_CHECK_SMALL( scalarValueExpr1 - ( std::pow(coord[0],2) + std::pow(coord[1],2) + std::pow(coord[2],2) ), tol );
            }

            BOOST_CHECK( cellData->HasArray("wScalar") );
            vtkDataArray * cellDataScalar = cellData->GetScalars("wScalar");
            BOOST_CHECK_EQUAL( cellDataScalar->GetNumberOfComponents() , 1 );
            BOOST_CHECK_EQUAL( cellDataScalar->GetDataSize() , nCell );
            BOOST_CHECK( cellData->HasArray("wVectorial") );
            vtkDataArray * cellDataVector = cellData->GetVectors("wVectorial");
            BOOST_CHECK_EQUAL( cellDataVector->GetNumberOfComponents() , 3 );
            BOOST_CHECK_EQUAL( cellDataVector->GetDataSize() , 3*nCell );
            BOOST_CHECK( cellData->HasArray("wTensor2Symm") );
            vtkDataArray * cellDataTensor2Symm = cellData->GetTensors( "wTensor2Symm" );
            BOOST_CHECK_EQUAL( cellDataTensor2Symm->GetNumberOfComponents() , 6 );
            BOOST_CHECK_EQUAL( cellDataTensor2Symm->GetDataSize() , 6*nCell );

            BOOST_CHECK( cellData->HasArray("expr2Scalar") );
            vtkDataArray * cellDataExpr2Scalar = cellData->GetScalars("expr2Scalar");
            BOOST_CHECK_EQUAL( cellDataExpr2Scalar->GetNumberOfComponents() , 1 );
            BOOST_CHECK_EQUAL( cellDataExpr2Scalar->GetDataSize() , nCell );
            BOOST_CHECK( cellData->HasArray("expr3Vectorial_e") );
            vtkDataArray * cellDataExpr3Vector = cellData->GetVectors("expr3Vectorial_e");
            BOOST_CHECK_EQUAL( cellDataExpr3Vector->GetNumberOfComponents() , 3 );
            BOOST_CHECK_EQUAL( cellDataExpr3Vector->GetDataSize() , 3*nCell );

            double valueExpr2 = 1;
            if ( blockName == "Omega2" )
                valueExpr2 = 2;
            else if ( blockName == "Omega3" )
                valueExpr2 = 3;
            for ( vtkIdType p=0;p<nCell;++p )
            {
                double scalarValue = cellDataExpr2Scalar->GetTuple1( p );
                BOOST_CHECK_SMALL( scalarValue - valueExpr2, tol );
            }
        }

    }
#endif

}

BOOST_AUTO_TEST_CASE( test_3 )
{
    auto mesh = unitSquare();
    auto Xh = FunctionSpace<Mesh<Simplex<2> >, bases<Lagrange<1>,Lagrange<1> > >::New(mesh);
    auto V = Xh->element();
    auto e = exporter(_mesh=mesh);
    e->add("V",V);
    e->save();
}

BOOST_AUTO_TEST_SUITE_END()
