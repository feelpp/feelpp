/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <fstream>

namespace Feel
{

template <int nDim,uint16_type OrderVelocity,uint16_type OrderPressure, uint16_type OrderGeo = 1>
int
runApplicationFluid()
{
    using namespace Feel;

    typedef FeelModels::FluidMechanics< Simplex<nDim,OrderGeo>,
                                        Lagrange<OrderVelocity, Vectorial,Continuous,PointSetFekete>,
                                        Lagrange<OrderPressure, Scalar,Continuous,PointSetFekete> > model_type;
    auto FM = model_type::New("fluid");

    std::ofstream outfile;
    outfile.open("/ssd/berti/phd/threeSphereCM.csv", std::ios_base::out);
    if(outfile.is_open() )
    {
        std::cout << "the file outfile is open" << std::endl;
    }
    else
    {
        std::cout << "outfile not open" <<std::endl;
    }
    double sphere_radius = 1;
    FM->init();
    FM->printAndSaveInfo();
    double single_step = 100*FM->timeStep();
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
            auto u = FM->fieldVelocity();
            auto p = FM->fieldPressure();
            auto Id = eye<nDim,nDim>();
            auto defv = sym(gradv(u));
            auto sigmav = -idv(p)*Id + 2*1e-3*defv;
            
            auto force_center = integrate(_range=markedfaces(FM->mesh(),"SphereCenter"),_expr=sigmav*N()).evaluate();
            auto force_left =integrate(_range=markedfaces(FM->mesh(),"SphereLeft"),_expr=sigmav*N()).evaluate();
            auto force_right =integrate(_range=markedfaces(FM->mesh(),"SphereRight"),_expr=sigmav*N()).evaluate();

            /*auto force_center_newmethod = integrate(_range=markedfaces(FM->mesh(),"SphereCenter"),_expr=leftfacev(sigmav*N()*(emarker()==FM->mesh()->markerId("Fluid")))+
                                                    rightfacev(sigmav*N()*(emarker()==FM->mesh()->markerId("Fluid")))).evaluate();
            auto force_left_newmethod =integrate(_range=markedfaces(FM->mesh(),"SphereLeft"),_expr=leftfacev(sigmav*N()*(emarker()==FM->mesh()->markerId("Fluid")))+
                                                 rightfacev(sigmav*N()*(emarker()==FM->mesh()->markerId("Fluid")))).evaluate();
            auto force_right_newmethod =integrate(_range=markedfaces(FM->mesh(),"SphereRight"),_expr=leftfacev(sigmav*N()*(emarker()==FM->mesh()->markerId("Fluid")))+
             rightfacev(sigmav*N()*(emarker()==FM->mesh()->markerId("Fluid")))).evaluate();*/

            auto l = form1(_test=FM->fieldVelocity().functionSpace());
            l = integrate(_range=markedfaces(FM->mesh(),{"SphereCenter","SphereLeft","SphereRight"}),_expr=inner(sigmav*N(),id(u)));
            auto vitesse = FM->fieldVelocity();
            vitesse.setZero();
            vitesse[ComponentType::X].setOnes();
            auto force_linform_X = l(vitesse);
            vitesse.setZero();
            vitesse[ComponentType::Y].setOnes();
            auto force_linform_Y = l(vitesse);
            vitesse.setZero();
            vitesse[ComponentType::Z].setOnes();
            auto force_linform_Z = l(vitesse);
            vitesse.setZero();
            
            
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
                outfile << FM->time()<< "Central sphere: ( "<< cm_central_sphere(0,0)<<", " << cm_central_sphere(1,0) <<", "<< cm_central_sphere(2,0) <<")\n";
                outfile << FM->time()<< "Left sphere: ( "<< cm_left_sphere(0,0)<<", " << cm_left_sphere(1,0) <<", "<< cm_left_sphere(2,0) <<")\n";
                outfile << FM->time()<< "Right sphere: ( "<< cm_right_sphere(0,0)<<", " << cm_right_sphere(1,0) <<", "<< cm_right_sphere(2,0) <<")\n";
                outfile << FM->time()<< "Left-center spheres distance: ( "<< cm_left_sphere(0,0)-cm_central_sphere(0,0)<<", " << cm_left_sphere(1,0)-cm_central_sphere(1,0) <<", "<< cm_left_sphere(2,0)-cm_central_sphere(2,0) <<")\n";
                outfile << FM->time()<< "Right-center spheres distance: ( "<< cm_central_sphere(0,0)-cm_right_sphere(0,0)<<", " << cm_central_sphere(1,0)-cm_right_sphere(1,0) <<", "<< cm_central_sphere(2,0)-cm_right_sphere(2,0) <<")\n";
                outfile << FM->time()<< "Central sphere force: ( "<< force_center(0,0)<<", " << force_center(1,0) <<", "<< force_center(2,0) <<")\n";
                outfile << FM->time()<< "Left sphere force: ( "<< force_left(0,0)<<", " << force_left(1,0) <<", "<< force_left(2,0) <<")\n";
                outfile << FM->time()<< "Right sphere force: ( "<< force_right(0,0)<<", " << force_right(1,0) <<", "<< force_right(2,0) <<")\n";
                outfile << FM->time()<< "Overall force: ( "<< force_right(0,0)+force_left(0,0)+force_center(0,0)<<", " << force_right(1,0)+force_left(1,0)+force_center(1,0) <<", "<< force_right(2,0)+force_left(2,0)+force_center(2,0) <<")\n";
                /* outfile << FM->time()<< "Central sphere force - new method: ( "<< force_center_newmethod(0,0)<<", " << force_center_newmethod(1,0) <<", "<< force_center_newmethod(2,0) <<")\n";
                outfile << FM->time()<< "Left sphere force - new method: ( "<< force_left_newmethod(0,0)<<", " << force_left_newmethod(1,0) <<", "<< force_left_newmethod(2,0) <<")\n";
                outfile << FM->time()<< "Right sphere force - new method: ( "<< force_right_newmethod(0,0)<<", " << force_right_newmethod(1,0) <<", "<< force_right_newmethod(2,0) <<")\n";
                 outfile << FM->time()<< "Overall force - new method: ( "<< force_right_newmethod(0,0)+force_left_newmethod(0,0)+force_center_newmethod(0,0)<<", " << force_right_newmethod(1,0)+force_left_newmethod(1,0)+force_center_newmethod(1,0) <<", "<< force_right_newmethod(2,0)+force_left_newmethod(2,0)+force_center_newmethod(2,0) <<")\n";*/
                outfile << FM->time()<< "Overall force - linspace: (" << force_linform_X <<", " << force_linform_Y << ", " <<force_linform_Z << ")\n";
                    
            }
            /*            if(FM->time()<=single_step)
            {
                 auto u = FM->fieldVelocity();
                u.on(_range=markedfaces(FM->mesh(),"SphereLeft"),_expr=oneX()*cst(4.0*sphere_radius/single_step));
            }
                //Impose Dirichlet first part
            else if (FM->time()>single_step && FM->time()<=2*single_step )
            {
                auto u = FM->fieldVelocity();
                u.on(_range=markedfaces(FM->mesh(),"SphereRight"),_expr=-oneX()*cst(4.0*sphere_radius/single_step));
            }//Impose Dirichlet second part
            else if (FM->time()>2*single_step && FM->time()<=3*single_step )
            {
                auto u = FM->fieldVelocity();
                u.on(_range=markedfaces(FM->mesh(),"SphereLeft"),_expr=-oneX()*cst(4.0*sphere_radius/single_step));
            }  //Impose Dirichlet third part
            else if (FM->time()>3*single_step && FM->time()<=4*single_step )
            {
                auto u = FM->fieldVelocity();
                 u.on(_range=markedfaces(FM->mesh(),"SphereRight"),_expr=oneX()*cst(4.0*sphere_radius/single_step));
            }    //Impose Dirichlet fourth part
             MPI_Barrier(FM->worldComm());*/
            FM->solve();
            FM->exportResults();
        }
        outfile.close();
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
