/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feells/disttoentityrange.hpp>
#include <feel/feelmesh/remesh.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/options.hpp>
namespace Feel
{
struct Body
{
    std::vector<std::string> v_markers;
    std::vector<std::string> f_markers;
};

template <typename DistFieldT, typename DistToRangeExpr, typename = std::enable_if_t<is_functionspace_element_v<DistFieldT>>>
auto flatThenIncreaseAroundEntityRange( DistFieldT const& d, DistToRangeExpr const& e, int coefexp = 3 )
{
    return ( idv( d ) < e ) + ( idv( d ) > e - cst( 1e-8 ) ) * exp( coefexp * idv( d ) / h() );
}

template <typename type_element>
auto computeRatioRadii( type_element elt )
{
    auto a = elt.faceMeasure( 0 );
    auto b = elt.faceMeasure( 1 );
    auto c = elt.faceMeasure( 2 );
    auto rho = a * b * c * ( a + b + c ) / 2.0 / ( 4.0 * elt.measure() * elt.measure() );
    return rho;
}

template <typename type_mesh, typename type_tolerance, typename quality_field_type>
bool mesh_quality( type_mesh const& mesh, type_tolerance tolerance, quality_field_type quality_field )
{
    bool toRemesh = false;
    //Cycle on the elements of mesh
    for ( auto const& eltWrap : elements( mesh ) )
    {
        auto const& elt = unwrap_ref( eltWrap );
        auto RatioRadii = computeRatioRadii( elt );

        // If the half ratio of radii is larger than the tolerance (>>1), remesh
        if ( RatioRadii * 0.5 >= tolerance )
        {
            std::cout << elt.barycenter() << " " << RatioRadii << std::endl;
            toRemesh = true;
            return toRemesh;
        }
    }
    // See if the
    return toRemesh;
}
template <typename typemesh, typename typealemesh>
typemesh const rem( typemesh const& mesh, typealemesh const& alemesh, Body const& body, typemesh parent = {} )
{
    bool first_material = true;
    //auto mats = rangeMeshElementsByMaterial();
    //for ( auto const& matPair : mats )

    //auto meshMarkers = mat.meshMarkers();
    auto submesh = mesh; //createSubmesh(mesh,markedelements(mesh,meshMarkers));
    auto Xh = Pch<1>( submesh );
    auto metric = Xh->element();
    auto phi = distToEntityRange( Xh, markedfaces( Xh->mesh(), body.f_markers ) );
    auto nlayers = ioption( _name = "remesh.metric.layers" );
    auto [havg, hmin, hmax] = hMeasures( submesh );
    auto expr_magnitude_disp = ( sqrt( inner( idv( alemesh->displacement() ), idv( alemesh->displacement() ) ) ) + cst( 1 ) );
    if(soption("remesh.strategy")=="progressive")
    {
        metric.on( _range = elements( mesh ),
                   _expr = cst( havg ) );
    }
    else if(soption("remesh.strategy")=="constant")
    {
        metric.on(_range=elements(mesh),_expr=cst(doption("remesh.strategy.constant.value")));
    }
    else
    {
        std::cout << "problem on the mesh metric, not defined" << std::endl;
    }
    ///_expr=flatThenIncreaseAroundEntityRange(phi,nlayers*h())*havg);//(hmin+0.5*havg)/(expr_magnitude_disp));
    //metric.on(_range=elements(mesh),_expr=flatThenIncreaseAroundEntityRange(phi,nlayers*h(),nlayers)*(hmin+havg)*0.5);//cst(havg));
    //auto r = remesher( submesh, mesh->markerName( "Swimmer" ), alemesh->aleFactory()->flagSet( "moving" ) );
    auto r = remesher( submesh, body.v_markers, body.f_markers, parent );

    r.setMetric( metric );
    auto new_mesh_remeshed = r.execute();
    new_mesh_remeshed->updateForUse();
    new_mesh_remeshed->setMarkerNames( submesh->markerNames() );
    /*
        auto r2 = remesher(new_mesh_remeshed2,mesh->markerName("Fluid"), alemesh->aleFactory()->flagSet("moving") );
        
        r2.setMetric(metric);
        auto new_mesh_remeshed = r2.execute();
        new_mesh_remeshed->updateForUse();
        new_mesh_remeshed->setMarkerNames(new_mesh_remeshed2->markerNames());
    */

    return new_mesh_remeshed; //mesh.reset(new typemesh(new_mesh_remeshed));

    std::cout << "We need to remesh!!" << std::endl;
}
#if 0
void updateMeasureFile( std::string fileName, std::string rootRepository )
{
    std::string filename = rootRepository + "/" + fileName;
    std::ofstream of_c( filename, std::ios_base::binary | std::ios_base::app );
    std::string measureName = rootRepository + "/fluid.measures.csv";
    std::ifstream if_b( measureName, std::ios_base::binary );
    of_c.seekp( 0, std::ios_base::end );
    of_c << if_b.rdbuf();
}
#endif
template <int nDim, uint16_type OrderVelocity, uint16_type OrderPressure, uint16_type OrderGeo = 1>
int runApplicationFluid( Body const& body )
{
    using namespace Feel;

    typedef FeelModels::FluidMechanics<Simplex<nDim, OrderGeo>,
                                       Lagrange<OrderVelocity, Vectorial, Continuous, PointSetFekete>,
                                       Lagrange<OrderPressure, Scalar, Continuous, PointSetFekete>>
        model_type;
    using trace_space_component_type = typename model_type::space_trace_velocity_type::component_functionspace_type;
    std::shared_ptr<model_type> FM_ref( new model_type( "fluid" ) );
    //auto FM_ref = model_type::New("fluid");
    FM_ref->init();
    double tolerance = 8;
    auto global_mesh = FM_ref->mesh();
    auto Xh_ref_swimmer = Pchv<1>( FM_ref->mesh() );
    auto bdf_global = bdf( _space = FM_ref->functionSpaceVelocity(), _name = "mybdf", _prefix = "mybdf" );
    FeelModels::ModelMeasuresIO measures( "fluid-global.measures.csv", Environment::worldCommPtr() );    
    auto e = exporter( _mesh = global_mesh, _name = "move", _geo = "change" );
    auto eSave = [&e]( double t, auto alemesh, auto vel, auto press, auto w ) {
        e->step( t )->setMesh( alemesh->movingMesh() );
        e->step( t)->add("disp",alemesh->displacement());
        e->step( t )->add( "vel", vel );
        e->step( t )->add( "press",press);
        e->step( t )->add( "test_interp", w );
        e->save(); };
    auto eRef = exporter( _mesh = global_mesh, _name = "reference", _geo = "change" );
    auto eRefSave = [&eRef]( double t, auto alemesh ) {
        eRef->step( t )->setMesh( alemesh->referenceMesh() );
        eRef->step( t )->add( "disp", alemesh->displacement() );
        eRef->save();
    };
    bool remesh_boolean = true;
    std::shared_ptr<model_type> FM( new model_type( "fluid" ) );
    //typedef typename model_type::space_mesh_disp_type Disp_Space_type;
    //Disp_Space_type Disp_Space;
    // Quality field
    typedef Lagrange<0, Scalar, Discontinuous> basis_quality_type;
    typedef FunctionSpace<Mesh<Simplex<nDim, OrderGeo>>, bases<basis_quality_type>> space_quality_type;
    typedef std::shared_ptr<space_quality_type> space_quality_ptrtype;
    space_quality_ptrtype Xh_Quality;
    typedef typename space_quality_type::element_type element_quality_type;
    typedef std::shared_ptr<element_quality_type> element_quality_ptrtype;
    element_quality_type quality_field;
    double delta_CM, delta_CM_init;
    int number_of_remesh = 0;

    /*typedef Lagrange<1, Vectorial,Discontinuous> basis_disp_type;
    typedef FunctionSpace<Mesh<Simplex<nDim,OrderGeo>>, bases<basis_disp_type>> space_disp_type;
    typedef std::shared_ptr<space_disp_type> space_disp_ptrtype;
    space_disp_ptrtype Xh_disp;
    typedef typename space_disp_type::element_type element_disp_type;
    typedef std::shared_ptr<element_disp_type> element_disp_ptrtype;*/
    typedef std::shared_ptr<trace_space_component_type> space_disp_ptrtype;
    space_disp_ptrtype Xh_disp;
    typedef typename trace_space_component_type::element_type element_disp_type;
    typedef std::shared_ptr<element_disp_type> element_disp_ptrtype;
    element_disp_type PosInit;
    element_disp_type Disp;

    for ( bdf_global->start(); !bdf_global->isFinished(); bdf_global->next() )
    {
        if ( remesh_boolean == true )
        {
            // PREPARATION OF THE INITIAL CONFIGURATION

            //auto alemesh_init = FeelModels::meshale( _mesh=FM_ref->mesh(),_prefix="alemesh_init" );
            //alemesh_init->init();
            auto disp = FM_ref->meshALE()->functionSpace()->element();
            if constexpr ( nDim == 2 )
            {
                auto dispExpr = expr<nDim, 1>( "{0,0}:x" );
                disp.on( _range = markedfaces( FM_ref->meshALE()->movingMesh(), body.f_markers ), _expr = dispExpr );
            }
            else if constexpr ( nDim == 3 )
            {
                auto dispExpr = expr<nDim, 1>( "{0,0,0}:x" );
                disp.on( _range = markedfaces( FM_ref->meshALE()->movingMesh(), body.f_markers ), _expr = dispExpr );
            }
            FM_ref->meshALE()->update( disp );
            auto new_mesh = rem( FM_ref->meshALE()->movingMesh(), FM_ref->meshALE(), body );
            CHECK( FM_ref->meshALE()->movingMesh()->isParentMeshOf( new_mesh) );
            new_mesh->saveHDF5( "init_mesh.json" );

            auto remeshing_time = FM_ref->time();
            saveGMSHMesh( _mesh = new_mesh, _filename = "init_mesh.msh" );
            FM = model_type::New( "fluid" );
            FM->setMesh( new_mesh );
            FM->setTimeInitial( remeshing_time );
            FM->addParameterInModelProperties( "T0", 0 );
            FM->init( true );

            Xh_Quality = space_quality_type::New( _mesh = FM->mesh() );
            quality_field = Xh_Quality->element();
            FM->timeStepBase()->setTimeInitial( bdf_global->time() - bdf_global->timeStep() );
            FM->printAndSaveInfo();
            FM->startTimeStep();
            auto center_of_mass_init = integrate( _range = markedelements( FM->mesh(), body.v_markers ), _expr = P() ).evaluate();
            auto mass = integrate( _range = markedelements( FM->mesh(), body.v_markers ), _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
            center_of_mass_init /= mass;
            delta_CM_init = center_of_mass_init( 0, 0 );
            remesh_boolean = false;
        }

        if ( FM_ref->worldComm().isMasterRank() )
        {
            std::cout << "============================================================\n";
            std::cout << "bdf_global: " << bdf_global->time() << "s \n";
            std::cout << "FM_time" << FM->time() << "s \n";
            std::cout << "============================================================\n";
        }

        //FM->solve();

        if ( FM->worldComm().isMasterRank() )
        {
            std::cout << "============================================================\n";
            std::cout << "time simulation: " << FM->time() << "s \n";
            std::cout << "============================================================\n";
            std::cout << "CHECK IF THE MESH QUALITY IS GOOD\n";

            bool toRemesh = mesh_quality( FM->meshALE()->movingMesh(), tolerance, quality_field );
            if ( toRemesh )
            {
                //Disp_Space.elementsPtr(number_of_remesh).push_back(FM->meshALE()->displacement());
                //FM->meshALE()->functionSpace()->elementsPtr(number_of_remesh).push_back(FM->meshALE()->displacement());
                /*auto submesh_swimmer = createSubmesh(FM->mesh(),markedelements(FM->mesh(),"Swimmer"));
                Xh_disp = space_disp_type::New(_mesh=submesh_swimmer);
                auto PosCurrent = Xh_disp->element();
                PosCurrent.on(_range=elements(support(Xh_disp)),_expr=P());*/
                number_of_remesh += 1;

                delta_CM = delta_CM_init;
                std::cout << "Need to remesh now!  \n";
                //bool toRemesh = mesh_quality(FM->meshALE()->movingMesh(),tolerance,quality_field);
                auto new_mesh = rem( FM->meshALE()->movingMesh(), FM->meshALE(), body, FM_ref->meshALE()->movingMesh() );
                new_mesh->saveHDF5( "new_mesh.json" );
                auto u = FM->fieldVelocity();
                u.saveHDF5( "new_initialCondition.json" );
                auto p = FM->fieldPressure();
                p.saveHDF5( "new_Pressure.json" );
                auto remeshing_time = FM->time();
                saveGMSHMesh( _mesh = new_mesh, _filename = "new_mesh.msh" );

                std::shared_ptr<model_type> FM_after_remeshing( new model_type( "fluid" ) );
                FM_after_remeshing->setMesh( new_mesh );
                FM_after_remeshing->setTimeInitial( remeshing_time - FM->timeStep() );
                FM_after_remeshing->addParameterInModelProperties( "T0", remeshing_time - FM->timeStep() );
                auto center_of_mass_cur = integrate( _range = markedelements( FM->mesh(), body.v_markers ), _expr = P() ).evaluate();
                auto mass = integrate( _range = markedelements( FM->mesh(), body.v_markers ), _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
                center_of_mass_cur /= mass;
                delta_CM -= center_of_mass_cur( 0, 0 );
                std::cout << "delta_CM " << delta_CM << std::endl;

                FM_after_remeshing->init( true );
                FM_after_remeshing->updateParameterValues();
                auto Yh = Pch<1>( new_mesh );
                auto w = Yh->element();
                w.on( _range = elements( support( Yh ) ), _expr = Px() );

                std::cout << "OK after init" << std::endl;
                FM_after_remeshing->printAndSaveInfo();
                //eSave(bdf_global->time()+FM->timeStep(),FM_after_remeshing->meshALE(),FM->fieldVelocity(),FM->fieldPressure(),w);
                std::cout << "OK after print" << std::endl;
                //std::cout << FM_after_remeshing->functionSpaceVelocity() << std::endl;
                //std::cout << FM_after_remeshing->fieldVelocity() << std::endl;

                auto u_new = FM_after_remeshing->fieldVelocity();
                std::cout << "OK after vel" << std::endl;
                auto p_new = FM_after_remeshing->fieldPressure();
                std::cout << "OK before VelOpinterp" << std::endl;

                toRemesh = false;

                std::vector<std::string> markersInterpolate;
                markersInterpolate = { "Fluid" };
                FM_after_remeshing->init( FM, markersInterpolate );
                FM_after_remeshing->startTimeStep();
                FM = FM_after_remeshing;

                Feel::cout << "FM " << FM << "FM_after_remeshing " << FM_after_remeshing << std::endl;             
            }
            // Picard sub iterations 
            for(int i_picard = 0; i_picard < ioption("picard.iterations"); ++ i_picard )
            {
                FM->updateParameterValues();
                Feel::cout << "ok fino export results" << std::endl;
                Feel::cout << " is MoveDomain " << FM->isMoveDomain() << std::endl;
                FM->setApplyMovingMeshBeforeSolve( false );

                auto expr_swimming_f = [&FM](double t){
                    if constexpr ( nDim == 2 )
                    {
                        auto ud = expr(expr<nDim,1>("{udxt_0,udxt_1}:udxt_0:udxt_1"),FM->symbolsExpr());
                        ud.setParameterValues({{"t",t}});
                        return ud;
                    }
                    else
                    {
                        auto ud = expr(expr<nDim,1>("{udxt_0,udxt_1,udxt_2}:udxt_0:udxt_1:udxt_2"),FM->symbolsExpr());
                        ud.setParameterValues({{"t",t}});
                        return ud;
                    }

                };
                auto expr_swimming_tnp1 =  expr_swimming_f(FM->time());
                auto expr_swimming_tnp05 =  expr_swimming_f(FM->time()-0.5*FM->timeStep());
                auto expr_swimming_tn =  expr_swimming_f(FM->time()-FM->timeStep());
                expr_swimming_tnp1.setParameterValues( FM->modelProperties().parameters().toParameterValues() );
                expr_swimming_tnp05.setParameterValues( FM->modelProperties().parameters().toParameterValues() );
                expr_swimming_tn.setParameterValues( FM->modelProperties().parameters().toParameterValues() );
                auto ud = Xh_ref_swimmer->element( expr_swimming_tnp1 );
                auto integ_ud = Xh_ref_swimmer->element( 1.0/6.0*(expr_swimming_tn+4.0*expr_swimming_tnp05+expr_swimming_tnp1) );
                auto integral_1 = integrate(_range=elements(Xh_ref_swimmer->mesh()),_expr=expr_swimming_tnp1).evaluate();
                auto integral_2 = integrate(_range=elements(Xh_ref_swimmer->mesh()),_expr=idv(ud)).evaluate();
                auto integral_3 = integrate(_range=elements(Xh_ref_swimmer->mesh()),_expr=idv(integ_ud)).evaluate();
                std::cout << "Integral1 "  << integral_1 << " Integral2 " << integral_2 <<" Integral3 " << integral_3 << std::endl;            
                auto se = Feel::vf::symbolsExpr( FM->symbolsExpr(), 
                                                 symbolExpr( "ud", idv( ud ),SymbolExprComponentSuffix( nDim,1 ) ),
                                                 symbolExpr( "integ_ud", idv( integ_ud ),SymbolExprComponentSuffix( nDim,1 ) )  );

                FM->updateALEmesh( se );
                FM->solve();

            } // end of Picard
            FM->exportResults();
            measures.setMeasures( FM->postProcessMeasuresIO().currentMeasures() );
            measures.exportMeasures();
            // updateMeasureFile( "fluid_glob.measures.csv", FM->rootRepository() );
            {
                auto Yh = Pch<1>( FM->mesh() );
                auto w = Yh->element();
                w.on( _range = elements( support( Yh ) ), _expr = Px() );
                eSave( bdf_global->time(), FM->meshALE(), FM->fieldVelocity(), FM->fieldPressure(), w );
                eRefSave( bdf_global->time(), FM->meshALE() );
            }
            Feel::cout << "ok post export results" << std::endl;
            FM->updateTimeStep();
        }
    }
    return 0;
}

} // namespace Feel

int main( int argc, char** argv )
{
    using namespace Feel;
    po::options_description fluidmecoptions( "application fluid-mechanics options" );
    fluidmecoptions.add( toolboxes_options( "fluid" ) );
    fluidmecoptions.add( bdf_options( "mybdf" ) ).add( ts_options( "mybdf" ) );
    fluidmecoptions.add( backend_options( "Ip" ) );
    fluidmecoptions.add( backend_options( "Iv" ) );
    fluidmecoptions.add( backend_options( "Idisp" ) );
    fluidmecoptions.add( alemesh_options( "alemesh_init" ) );
    // clang-format off
    fluidmecoptions.add_options()
      ( "case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension" )
      ( "case.discretization", Feel::po::value<std::string>()->default_value( "P2P1G1" ), "discretization : P2P1G1,P2P1G2" )
      ( "export.matlab", po::value<bool>()->default_value( true ), "export matrix and vector to matlab" )
      ( "remesh.metric.layers", po::value<int>()->default_value( 2 ), "number of remeshing layers" )
      ( "body.markers.volume", po::value<std::vector<std::string> >()->multitoken(), "list of volume markers for the moving body" )
      ( "body.markers.facet", po::value<std::vector<std::string> >()->multitoken(), "list of facet markers for the moving body" )
      ( "picard.iterations", Feel::po::value<int>()->default_value( 1 ), "number of picard iterations" )
      ( "remesh.strategy" ,Feel::po::value<std::string>()->default_value("progressive"),"remesh strategy")
      ( "remesh.strategy.constant.value" ,Feel::po::value<double>()->default_value(1),"remesh constant value")
      ( "remesh.tolerance" ,Feel::po::value<double>()->default_value(4),"remesh constant value")
      ;
    // clang-format on

    Environment env( _argc = argc, _argv = argv,
                     _desc = fluidmecoptions,
                     _about = about( _name = "application_fluid",
                                     _author = "Feel++ Consortium",
                                     _email = "feelpp-devel@feelpp.org" ) );

    int dimension = ioption( _name = "case.dimension" );
    std::string discretization = soption( _name = "case.discretization" );
    if ( discretization == "P2P1" )
        discretization = "P2P1G1";

    auto dimt = hana::make_tuple( hana::int_c<2>, hana::int_c<3> );

    auto discretizationt = hana::make_tuple( hana::make_tuple( "P2P1G1", hana::make_tuple( hana::int_c<2>, hana::int_c<1>, hana::int_c<1> ) ) );

    
    Body body;
    body.v_markers = vsoption("body.markers.volume");
    body.f_markers = vsoption("body.markers.facet");
    std::cout << fmt::format("body volume makers: {}", body.v_markers) << std::endl;
    std::cout << fmt::format("body facet makers: {}", body.f_markers) << std::endl;
    int status = 0;
    hana::for_each( hana::cartesian_product( hana::make_tuple( dimt, discretizationt ) ), [&discretization, &dimension, &status, &body]( auto const& d ) {
        constexpr int _dim = std::decay_t<decltype( hana::at_c<0>( d ) )>::value;
        std::string const& _discretization = hana::at_c<0>( hana::at_c<1>( d ) );
        constexpr int _uorder = std::decay_t<decltype( hana::at_c<0>( hana::at_c<1>( hana::at_c<1>( d ) ) ) )>::value;
        constexpr int _porder = std::decay_t<decltype( hana::at_c<1>( hana::at_c<1>( hana::at_c<1>( d ) ) ) )>::value;
        constexpr int _gorder = std::decay_t<decltype( hana::at_c<2>( hana::at_c<1>( hana::at_c<1>( d ) ) ) )>::value;
        if ( dimension == _dim && discretization == _discretization )
            status = runApplicationFluid<_dim, _uorder, _porder, _gorder>( body );
    } );
    return status;
}
