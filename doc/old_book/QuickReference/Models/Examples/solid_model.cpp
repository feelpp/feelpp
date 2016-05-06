/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/solid/solidmechanics.hpp>

namespace Feel
{

template <uint16_type OrderDisp>
void
runApplicationSolid()
{
    typedef FeelModels::SolidMechanics< Simplex<FEELPP_DIM,1>,
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
                if ( SM->useDisplacementPressureFormulation() )
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
                SM->solve();
                SM->exportResults(cptCurrent);

                SM->fieldDisplacement().save(_path=SM->rootRepository()+"/"+(boost::format("uSol.field-%1%") % cptCurrent ).str(),
                                             _type=saveType );
                if ( SM->useDisplacementPressureFormulation() )
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
    }
    else
    {
        SM->init();
        SM->printAndSaveInfo();

        for ( ; !SM->timeStepBase()->isFinished(); SM->updateTimeStep() )
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
}

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;
	po::options_description solidmecoptions( "application solid-mechanics options" );
    solidmecoptions.add( feelmodels_options("solid") );
    solidmecoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P1,P2 ")
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
                     _about=about(_name="application_solid",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P1" )
        runApplicationSolid<1>();
    else if ( feapprox == "P2" )
        runApplicationSolid<2>();
    else CHECK( false ) << "invalid feapprox " << feapprox;


    return 0;
}
