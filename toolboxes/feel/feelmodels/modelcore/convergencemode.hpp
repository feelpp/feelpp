#ifndef FEELPP_TOOLBOXES_CONVERGENCE_HPP
#define FEELPP_TOOLBOXES_CONVERGENCE_HPP 1

template <typename ToolboxType>
int
runHConvergence(std::string const& prefix,
               std::function<int(std::shared_ptr<ToolboxType>)> runToolboxSimulation)
{
    using namespace Feel;
    int status = 0;
    CHECK( false ) << "TODO dir";
    FeelModels::ModelMeasuresStorage measures( "toto", Environment::worldCommPtr() );
    //FeelModels::ModelMeasuresIO measures( "h-convergence.measures.csv", Environment::worldCommPtr() );
    CHECK( Environment::vm().count( "case.mode.h-convergence.hsize" ) ) << "require option : case.mode.h-convergence.hsize";
    std::vector<double> meshSizes = vdoption( _name="case.mode.h-convergence.hsize" );

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

        measures.setValue( "run_id", k );
        measures.setValue( "mesh_size", meshSize );
        measures.setKeyValue( tb->postProcessMeasures().values() );
        measures.save(k);
    }
    Feel::cout << "H convergence results are in " << (fs::path(Environment::appRepository())/"h-convergence.measures.csv").string() << std::endl;
    return status;
}

#endif
