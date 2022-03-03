#ifndef FEELPP_TOOLBOXES_CONVERGENCE_HPP
#define FEELPP_TOOLBOXES_CONVERGENCE_HPP 1

#include <feel/feelmath/polyfit.hpp>
#include <feel/feelmath/vector.hpp>


template <typename ToolboxType>
int
runHConvergence(std::string const& prefix,
               std::function<int(std::shared_ptr<ToolboxType>)> runToolboxSimulation)
{
    using namespace Feel;
    int status = 0;

    auto repbase = fs::path(Environment::appRepository());
    FeelModels::ModelMeasuresStorage measures( repbase.string(), Environment::worldCommPtr() );
    CHECK( Environment::vm().count( "case.mode.h-convergence.hsize" ) ) << "require option : case.mode.h-convergence.hsize";
    std::vector<double> meshSizes = vdoption( _name="case.mode.h-convergence.hsize" );

    std::vector<std::string> fitMeasuresNames;
    std::vector<double> fitMeasuresRefSlopes;
    if ( Environment::vm().count( "case.mode.h-convergence.measures" ) )
        fitMeasuresNames = vsoption( _name="case.mode.h-convergence.measures" );
    if ( Environment::vm().count( "case.mode.h-convergence.slopes" ) )
        fitMeasuresRefSlopes = vdoption( _name="case.mode.h-convergence.slopes" );
    CHECK( fitMeasuresNames.size() == fitMeasuresRefSlopes.size() ) << "invalid input";
    std::map<std::string,std::vector<double>> fitMeasuresValues;
    for ( int k=0;k<meshSizes.size();++k )
    {
        double meshSize = meshSizes[k];
        std::string runSubdir = (boost::format("run_h%1%")%k).str();
        Feel::Environment::setOptionValue( prefix+".gmsh.hsize", meshSize );

        auto tb = ToolboxType::New( _prefix=prefix,
                                    _repository=Feel::FeelModels::ModelBaseRepository( (fs::path(Environment::appRepository())/runSubdir).string(), false, Environment::exprRepository() ) );
        status = runToolboxSimulation( tb );
        if ( status != 0 )
            return status;

        for ( std::string const& fitMeasuresName : fitMeasuresNames )
            fitMeasuresValues[fitMeasuresName].push_back( tb->postProcessMeasures().value(fitMeasuresName) );

        measures.setValue( "h-convergence", "run_id", k );
        measures.setValue( "h-convergence", "mesh_size", meshSize );
        measures.setKeyValue( "h-convergence", tb->postProcessMeasures().values() );
        measures.save(k);
    }

    Feel::cout << "H convergence results are in " << (fs::path(Environment::appRepository())/"values.h-convergence.csv").string() << std::endl;

    if ( !fitMeasuresNames.empty() )
    {
        Feel::Table tableFitChecker;
        tableFitChecker.add_row( {"check","measure name","ref slope","measured slope","tolerance"} );
        tableFitChecker.format().setFirstRowIsHeader( true );

        for ( int k=0;k<fitMeasuresNames.size();++k )
        {
            std::string const& fitMeasuresName = fitMeasuresNames[k];
            double refSlope = fitMeasuresRefSlopes[k];
            double tol = 0.3;
            auto c = polyfit( log(meshSizes), log(fitMeasuresValues[fitMeasuresName]), 1 );
            double measuredSlope = c[1];
            double diffSlope = std::abs(measuredSlope-refSlope);
            bool checkIsOk = diffSlope < tol;
            if ( status == 0 && !checkIsOk )
                status = 1;

            std::string checkStr = checkIsOk? "[success]" : "[failure]";
            tableFitChecker.add_row( { checkStr, fitMeasuresName, refSlope, measuredSlope, tol } );
            tableFitChecker( tableFitChecker.nRow()-1,0).format().setFontColor( checkIsOk ? Font::Color::green :Font::Color::red );
        }

        auto tabInfoFitChecker = TabulateInformationsSections::New();
        tabInfoFitChecker->add( "Fit Checkers", TabulateInformations::New( tableFitChecker ) );
        if ( tableFitChecker.nRow() > 1 )
            Feel::cout << *tabInfoFitChecker << std::endl;
    }

    return status;
}

#endif
