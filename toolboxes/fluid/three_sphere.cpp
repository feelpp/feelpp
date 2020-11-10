/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <fstream>

namespace Feel
{

    template<class FluidT> 
    void testPmatrix(std::shared_ptr<FluidT>& FM)
    {
        auto P = FM->matP();
        auto u = FM->fieldVelocity();
        
        for(auto const& [name,body] : FM->bodySet())
        {
            if(body.name()=="SphereCenter")
            {
                for(int i=0;i<=2;i++)
                {
                    (*body.fieldTranslationalVelocityPtr())(i) = 1;
                }
                
            }   
        }
        
        FM->buildBlockVector();
        auto U = FM->blockVectorSolution().vectorMonolithic()->clone();
        * U = *FM->blockVectorSolution().vectorMonolithic();
        P->multVector(U,FM->blockVectorSolution().vectorMonolithic());
    #if 0
        auto explictPartOfSolution = FM->algebraicFactory()->explictPartOfSolution();
        FM->blockVectorSolution().vectorMonolithic()->add( 1.0, explictPartOfSolution );
    #endif
        FM->blockVectorSolution().localize();

        for(auto const& [name,body] : FM->bodySet())
        {
            std::cout << body.name() << " ";
            for(int i=0;i<=2;i++)
            {
                std::cout << (*body.fieldTranslationalVelocityPtr())(i);
            }
        
            std::cout << " " << std::endl;
        }

        FM->exportResults();
        
    }

    template <int nDim,uint16_type OrderVelocity,uint16_type OrderPressure, uint16_type OrderGeo = 1>
    int
    runApplicationFluid()
    {
        using namespace Feel;

        typedef FeelModels::FluidMechanics< Simplex<nDim,OrderGeo>,
                                            Lagrange<OrderVelocity, Vectorial,Continuous,PointSetFekete>,
                                            Lagrange<OrderPressure, Scalar,Continuous,PointSetFekete> > model_type;
        auto FM = model_type::New("fluid");

        // Initialize the files where the spheres' speeds and positions are stored
        std::ofstream outfile;
        std::ofstream outfile2;
        outfile.open("/ssd/berti/phd/threeSphereCM.csv", std::ios_base::out); // here the path is absolute -> should be changed to relative
        outfile2.open("/ssd/berti/phd/threeSphereSPEED.csv", std::ios_base::out);
        if(outfile.is_open() )
        {
            std::cout << "the file outfile is open" << std::endl;
        }
        else
        {
            std::cout << "outfile not open" <<std::endl;
        }
        if(outfile2.is_open() )
        {
            std::cout << "the file outfile2 is open" << std::endl;
        }
        else
        {
            std::cout << "outfile not open" <<std::endl;
        }
        // Solve the fluid problem and store the results
        FM->init();
        FM->printAndSaveInfo();
        if ( FM->isStationary() )
        {
            FM->solve();
            FM->exportResults();
        }
        else
        {
            if ( !FM->doRestart() )
                FM->exportResults(FM->timeInitial());
            
            for ( FM->startTimeStep(); !FM->timeStepBase()->isFinished(); FM->updateTimeStep() )
            {
                if (FM->worldComm().isMasterRank())
                {
                    std::cout << "============================================================\n";
                    std::cout << "time simulation: " << FM->time() << "s \n";
                    std::cout << "============================================================\n"; 
                }
                FM->solve();
                FM->exportResults();

                //testPmatrix(FM);

                auto u = FM->fieldVelocity();
                auto p = FM->fieldPressure();
                auto Id = eye<nDim,nDim>();
                auto defv = sym(gradv(u));
                auto sigmav = -idv(p)*Id + 2*defv;
                
                auto fluid_mesh = createSubmesh(FM->mesh(),markedelements(FM->mesh(),"Fluid"));

                auto force_center = integrate(_range=markedfaces(fluid_mesh,"SphereCenter"),_expr=sigmav*N()).evaluate();
                auto force_left =integrate(_range=markedfaces(fluid_mesh,"SphereLeft"),_expr=sigmav*N()).evaluate();
                auto force_right =integrate(_range=markedfaces(fluid_mesh,"SphereRight"),_expr=sigmav*N()).evaluate();

                auto speed_center = mean(_range=markedfaces(fluid_mesh,"SphereCenter"),_expr=idv(u));
                auto speed_left =mean(_range=markedfaces(fluid_mesh,"SphereLeft"),_expr=idv(u));
                auto speed_right =mean(_range=markedfaces(fluid_mesh,"SphereRight"),_expr=idv(u));               
                
                Feel::cout << force_left << " " << force_center << " " << force_right << std::endl;
                auto cm_central_sphere = integrate(_range=markedfaces(FM->mesh(),"SphereCenter"),_expr=P()).evaluate();
                auto cm_left_sphere = integrate(_range=markedfaces(FM->mesh(),"SphereLeft"),_expr=P()).evaluate();
                auto cm_right_sphere = integrate(_range=markedfaces(FM->mesh(),"SphereRight"),_expr=P()).evaluate();
                auto meas_c = measure(_range=markedfaces(FM->mesh(),"SphereCenter"));
                auto meas_l = measure(_range=markedfaces(FM->mesh(),"SphereLeft"));
                auto meas_r = measure(_range=markedfaces(FM->mesh(),"SphereRight")); 
                cm_central_sphere /= meas_c;
                cm_left_sphere /= meas_l;
                cm_right_sphere /= meas_r;
                if (FM->worldComm().isMasterRank())
                {
                    // store the CM data (calculated above)
                    outfile << FM->time()<< ", "<< cm_central_sphere(0,0)<<", " << cm_central_sphere(1,0) <<", "<< cm_central_sphere(2,0) <<", "
                    << cm_left_sphere(0,0)<<", " << cm_left_sphere(1,0) <<", "<< cm_left_sphere(2,0) <<", "
                    << cm_right_sphere(0,0)<<", " << cm_right_sphere(1,0) <<", "<< cm_right_sphere(2,0) <<", "
                    << cm_left_sphere(0,0)-cm_central_sphere(0,0)<<", " << cm_left_sphere(1,0)-cm_central_sphere(1,0) <<", "<< cm_left_sphere(2,0)-cm_central_sphere(2,0) <<
                    ", "<< cm_central_sphere(0,0)-cm_right_sphere(0,0)<<", " << cm_central_sphere(1,0)-cm_right_sphere(1,0) <<", "<< cm_central_sphere(2,0)-cm_right_sphere(2,0) <<std::endl;
                    // Store the translational velocities (computed above)
                    outfile2 << FM->time() << "," << speed_center(0,0) << ", " << speed_center(1,0) << ", " << speed_center(2,0) << ", " 
                                                << speed_left(0,0) << ", " << speed_left(1,0) << ", " << speed_left(2,0) << ", " 
                                                << speed_right(0,0) << ", " << speed_right(1,0) << ", " << speed_right(2,0) ;
                    for(auto const& [name,body] : FM->bodySet())
                    {
                        // Store the translational velocities for comparison
                        for(int i=0;i<=2;i++)
                        {
                            outfile2 << ", " << (*body.fieldTranslationalVelocityPtr())(i);
                        }
                    }
                    outfile2 << std::endl;   
                }            
            }
            outfile.close();
            outfile2.close();
        }

        return !FM->checkResults();

    }

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;
	po::options_description fluidmecoptions( "application fluid-mechanics options" );
    fluidmecoptions.add( toolboxes_options("fluid") );
    fluidmecoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P2P1G1" ), "discretization : P2P1G1,P2P1G2")
        ( "export.matlab", po::value<bool>()->default_value( true ), "export matrix and vector to matlab" )
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=fluidmecoptions,
                   _about=about(_name="application_fluid",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");
    if ( discretization == "P2P1" )
        discretization = "P2P1G1";

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);

    auto discretizationt = hana::make_tuple( hana::make_tuple("P2P1G1", hana::make_tuple( hana::int_c<2>,hana::int_c<1>,hana::int_c<1>) ));
                                             //hana::make_tuple("P2P1G2", hana::make_tuple( hana::int_c<2>,hana::int_c<1>,hana::int_c<2>) ),
                                             //hana::make_tuple("P1P1G1", hana::make_tuple( hana::int_c<1>,hana::int_c<1>,hana::int_c<1>) ) );

    int status = 0;
    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension,&status]( auto const& d )
                    {
                        constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                        std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                        constexpr int _uorder = std::decay_t<decltype(hana::at_c<0>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        constexpr int _porder = std::decay_t<decltype(hana::at_c<1>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        constexpr int _gorder = std::decay_t<decltype(hana::at_c<2>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        if ( dimension == _dim && discretization == _discretization )
                            status = runApplicationFluid<_dim,_uorder,_porder,_gorder>();
                    } );
    return status;
}


