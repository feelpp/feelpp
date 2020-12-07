/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4
*/ 

#include <feel/feelmodels/solid/solidmechanics.hpp>

namespace Feel
{

template <int nDim,uint16_type OrderDisp>
int
runApplicationSolid()
{
    typedef FeelModels::SolidMechanics< Simplex<nDim,1>,
                                        Lagrange<OrderDisp, Vectorial,Continuous,PointSetFekete> > model_type;
    auto SM = model_type::New("solid");

    if ( SM->isStationary() )
    {
        bool algoQuasiStatic = boption(_name="solve-quasi-static");
        if ( !algoQuasiStatic )
        {
            SM->init();
            SM->printAndSaveInfo();
            SM->solve();
            SM->exportResults();
        }
        else
        {
            SM->init();
            SM->printAndSaveInfo();

            std::string variableSymbol = soption(_name="solve-quasi-static.variable-symbol");
            CHECK( SM->modelProperties().parameters().find(variableSymbol) != SM->modelProperties().parameters().end() ) << "not find symbol " << variableSymbol;
            double valueFinal = SM->modelProperties().parameters().find(variableSymbol)->second.value();
            double valueInitial = doption(_name="solve-quasi-static.variable-initial");
            double valueStep = doption(_name="solve-quasi-static.variable-step");
            int cptNeed = std::abs(valueFinal-valueInitial)/valueStep;
            double currentParam = valueInitial;

            std::string namefileVariableParameters = (fs::path(SM->rootRepository()) / fs::path("paramters-done.data")).string();

            std::string saveType = soption(_name="solve-quasi-static.save.file-format");

            int cptCurrent=1;
            if ( SM->doRestart() )
            {
                std::ifstream fileParameter(namefileVariableParameters.c_str());
                while ( ! fileParameter.eof() ) { fileParameter >> cptCurrent >> currentParam; }
                fileParameter.close();

                SM->fieldDisplacement().load(_path=SM->rootRepository()+"/"+(boost::format("uSol.field-%1%") %cptCurrent ).str(),
                                             _type=saveType );
                if ( SM->hasDisplacementPressureFormulation() )
                    SM->fieldPressure().load(_path=SM->rootRepository()+"/"+(boost::format("pSol.field-%1%") %cptCurrent ).str(),
                                             _type=saveType );

                SM->restartExporters(cptCurrent);
                ++cptCurrent;
            }
            else
            {
                std::ofstream fileParameter(namefileVariableParameters.c_str(), std::ios::out);
                fileParameter.close();
            }
            currentParam += valueStep;

            for ( ; cptCurrent <= cptNeed ; ++cptCurrent, currentParam += valueStep )
            {
                if ( Environment::isMasterRank() )
                {
                    std::cout << "============================================================\n";
                    std::cout << " value of parameter: " << currentParam << "\n";
                    std::cout << "============================================================\n";
                }

                SM->timerTool("Solve").setAdditionalParameter(variableSymbol,currentParam);
                SM->timerTool("PostProcessing").setAdditionalParameter(variableSymbol,currentParam);
                SM->timerTool("TimeStepping").setAdditionalParameter(variableSymbol,currentParam);

                SM->addParameterInModelProperties(variableSymbol,currentParam);
                SM->updateParameterValues();
                SM->solve();
                SM->exportResults(cptCurrent);

                SM->fieldDisplacement().save(_path=SM->rootRepository()+"/"+(boost::format("uSol.field-%1%") % cptCurrent ).str(),
                                             _type=saveType );
                if ( SM->hasDisplacementPressureFormulation() )
                    SM->fieldPressure().save(_path=SM->rootRepository()+"/"+(boost::format("pSol.field-%1%") % cptCurrent ).str(),
                                             _type=saveType );

                if (Environment::isMasterRank())
                {
                    std::ofstream fileParameter(namefileVariableParameters.c_str(), std::ios::out | std::ios::app);
                    fileParameter << cptCurrent << " " << currentParam << std::endl;
                    fileParameter.close();
                }
            }
        }
        bool saveSolution = boption(_name="save-solution");
        if ( saveSolution )
        {
            std::string saveType = soption(_name="save-solution.file-format");
            SM->fieldDisplacement().save(_path=(fs::path( SM->rootRepository() )/fs::path( "solution.displacement" )).string(),
                                         _type=saveType );
            if ( SM->hasDisplacementPressureFormulation() )
                SM->fieldPressure().save(_path=(fs::path( SM->rootRepository() )/fs::path( "solution.pressure" )).string(),
                                         _type=saveType );
        }
    }
    else
    {
        SM->init();
        SM->printAndSaveInfo();
        if ( !SM->doRestart() )
            SM->exportResults(SM->timeInitial());

        for ( SM->startTimeStep() ; !SM->timeStepBase()->isFinished(); SM->updateTimeStep() )
        {
            if (SM->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << SM->time() << "s \n";
                std::cout << "============================================================\n";
            }

            SM->solve();
            SM->exportResults();
        }
    }

     return !SM->checkResults();
}

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;
	po::options_description solidmecoptions( "application solid-mechanics options" );
    solidmecoptions.add( toolboxes_options("solid") );
    solidmecoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2")
        ("save-solution", Feel::po::value<bool>()->default_value(true), "save-solution")
#ifdef FEELPP_HAS_HDF5
        ("save-solution.file-format", Feel::po::value<std::string>()->default_value("hdf5"), "save-solution.file-format")
#else
        ("save-solution.file-format", Feel::po::value<std::string>()->default_value("binary"), "save-solution.file-format")
#endif
        ("solve-quasi-static", Feel::po::value<bool>()->default_value(false), "solve-quasi-static ")
        ("solve-quasi-static.variable-initial", Feel::po::value<double>()->default_value(0.0), "solve-quasi-static.variable-initial")
        ("solve-quasi-static.variable-step", Feel::po::value<double>()->default_value(0.1), "solve-quasi-static.variable-step")
        ("solve-quasi-static.variable-symbol", Feel::po::value<std::string>()->default_value(""), "solve-quasi-static.variable-symbol")
#ifdef FEELPP_HAS_HDF5
        ("solve-quasi-static.save.file-format", Feel::po::value<std::string>()->default_value("hdf5"), "solve-quasi-static.save.file-format")
#else
        ("solve-quasi-static.save.file-format", Feel::po::value<std::string>()->default_value("binary"), "solve-quasi-static.save.file-format")
#endif
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=solidmecoptions,
                     _about=about(_name="feelpp_toolbox_solid",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
#if FEELPP_INSTANTIATION_ORDER_MAX >= 2
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ) );
#else
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ) );
#endif
    int status = 0;
    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension,&status]( auto const& d )
                    {
                        constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                        std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                        constexpr int _dorder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                        if ( dimension == _dim && discretization == _discretization )
                            status = runApplicationSolid<_dim,_dorder>();
                    } );

    return status;
}
