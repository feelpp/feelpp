/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/heat/heat.hpp>

template <typename ToolboxType>
int
runToolboxSimulation( std::shared_ptr<ToolboxType> toolbox )
{
    toolbox->init();
    toolbox->printAndSaveInfo();

    if ( toolbox->isStationary() )
    {
        toolbox->solve();
        toolbox->exportResults();
    }
    else
    {
        if ( !toolbox->doRestart() )
            toolbox->exportResults(toolbox->timeInitial());

        for ( toolbox->startTimeStep() ; !toolbox->timeStepBase()->isFinished(); toolbox->updateTimeStep() )
        {
            if (toolbox->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << toolbox->time() << "s \n";
                std::cout << "============================================================\n";
            }

            toolbox->solve();
            toolbox->exportResults();
        }
    }
    return !toolbox->checkResults();
}


template <int nDim,int OrderT>
int
runHeatSimulation()
{
    using namespace Feel;
    typedef Feel::FeelModels::Heat< Simplex<nDim,1>,
                                    Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_type;
    auto heat = model_type::New(_prefix="heat");
    return runToolboxSimulation( heat );
}
template <int nDim,int OrderT>
int
runHeatHConvergence()
{
    using namespace Feel;
    int status = 0;

    typedef FeelModels::Heat< Simplex<nDim,1>,
                                      Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_type;

    FeelModels::ModelMeasuresIO measures( "h-convergence.measures.csv", Environment::worldCommPtr() );
    CHECK( Environment::vm().count( "case.mode.h-convergence.hsize" ) ) << "require option : case.mode.h-convergence.hsize";
    std::vector<double> meshSizes = vdoption( _name="case.mode.h-convergence.hsize" );

    for ( int k=0;k<meshSizes.size();++k )
    {
        double meshSize = meshSizes[k];
        std::string runSubdir = (boost::format("run_h%1%")%k).str();
        Feel::Environment::setOptionValue( "heat.gmsh.hsize", meshSize );

        auto heat = model_type::New( _prefix="heat",
                                     _repository=Feel::FeelModels::ModelBaseRepository( (fs::path(Environment::appRepository())/runSubdir).string(), false, Environment::exprRepository() ) );
        status = runToolboxSimulation( heat );
        if ( status != 0 )
            return status;

        measures.setMeasure( "run_id", k );
        measures.setMeasure( "mesh_size", meshSize );
        measures.setMeasures( heat->postProcessMeasuresIO().currentMeasures() );
        measures.exportMeasures();
    }
    return status;
}

int
main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description heatoptions( "heat options" );
    heatoptions.add( toolboxes_options("heat") );
    heatoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ("case.mode", Feel::po::value<std::string>()->default_value( "simulation" ), "mode : simulation, h-convergence")
        ("case.mode.h-convergence.hsize", po::value<std::vector<double> >()->multitoken(), "mesh hsize used in h-convergence" )
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=heatoptions,
                     _about=about(_name="toolboxes_heat",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    Feel::FeelModels::printToolboxApplication( "heat" );

    std::string mode = soption(_name="case.mode");
    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ),
                                             hana::make_tuple("P3", hana::int_c<3> ) );
#elif FEELPP_INSTANTIATION_ORDER_MAX >= 2
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ) );
#else
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ) );
#endif
    int status = 0;

    if ( mode == "simulation" || mode == "h-convergence" )
    {
        hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension,&status,&mode]( auto const& d )
                        {
                            constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                            std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                            constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                            if ( dimension == _dim && discretization == _discretization )
                            {
                                if ( mode == "simulation" )
                                    status = runHeatSimulation<_dim,_torder>();
                                else
                                    status = runHeatHConvergence<_dim,_torder>();
                            }
                        } );
    }
    else if ( mode == "p-convergence" )
    {
        // TODO
    }
    return status;
}
