
#include <feel/feelmodels2/modelcore/applibase.hpp>

namespace Feel {

namespace FeelModels {


ModelBase::ModelBase( std::string _theprefix,
                      WorldComm const& _worldComm,
                      std::string subPrefix,
                      std::string appliShortRepository )
    :
    M_worldComm(_worldComm),
    M_worldsComm(std::vector<WorldComm>(1,_worldComm)),
    M_localNonCompositeWorldsComm(std::vector<WorldComm>(1,_worldComm)),
    M_prefix(_theprefix),
    M_subPrefix(subPrefix),
    M_appliShortRepository(appliShortRepository),
    M_verbose( boption(_name="verbose",_prefix=this->prefix()) ),
    M_verboseAllProc( boption(_name="verbose_allproc",_prefix=this->prefix()) ),
    M_filenameSaveInfo( prefixvm(this->prefix(),prefixvm(this->subPrefix(),"appli.info")) ),
    M_timersActivated( boption(_name="timers.activated",_prefix=this->prefix()) ),
    M_timersSaveFileMasterRank( boption(_name="timers.save-master-rank",_prefix=this->prefix()) ),
    M_timersSaveFileMax( boption(_name="timers.save-max",_prefix=this->prefix()) ),
    M_timersSaveFileMin( boption(_name="timers.save-min",_prefix=this->prefix()) ),
    M_timersSaveFileMean( boption(_name="timers.save-mean",_prefix=this->prefix()) ),
    M_timersSaveFileAll( boption(_name="timers.save-all",_prefix=this->prefix()) ),
    M_scalabilitySave( boption(_name="scalability-save",_prefix=this->prefix()) ),
    M_scalabilityReinitSaveFile( boption(_name="scalability-reinit-savefile",_prefix=this->prefix()) )
{
    if (Environment::vm().count(prefixvm(this->prefix(),"scalability-path")))
        M_scalabilityPath = Environment::vm()[prefixvm(this->prefix(),"scalability-path")].as< std::string >();
    else
        M_scalabilityPath = this->appliRepositoryWithoutNumProc();

    if (Environment::vm().count(prefixvm(this->prefix(),"scalability-filename")))
        M_scalabilityFilename = Environment::vm()[prefixvm(this->prefix(),"scalability-filename")].as< std::string >();
    else
        M_scalabilityFilename = this->prefix()+".scalibility";
}
ModelBase::~ModelBase()
{}

WorldComm const&
ModelBase::worldComm() const { return M_worldComm; }
std::vector<WorldComm> const&
ModelBase::worldsComm() const { return M_worldsComm; }
void
ModelBase::setWorldsComm(std::vector<WorldComm> const& _worldsComm) { M_worldsComm=_worldsComm; }
std::vector<WorldComm> const&
ModelBase::localNonCompositeWorldsComm() const { return M_localNonCompositeWorldsComm; }
void
ModelBase::setLocalNonCompositeWorldsComm(std::vector<WorldComm> const& _worldsComm) { M_localNonCompositeWorldsComm=_worldsComm; }
void
ModelBase::createWorldsComm() {}//warning

std::string
ModelBase::prefix() const { return M_prefix; }
std::string
ModelBase::subPrefix() const { return M_subPrefix; }

po::variables_map const&
ModelBase::vm() const { return Environment::vm(); }

std::string
ModelBase::appliRepositoryWithoutNumProc() const
{
    return Environment::rootRepository()+"/"+
        this->appliShortRepository();
}

std::string
ModelBase::appliRepository() const
{
    return Environment::rootRepository() + "/" + this->appliShortRepositoryWithNumProc();
}

std::string
ModelBase::appliShortRepository() const
{
    //return option(_name="exporter.directory").as<std::string>()+"/";
    return M_appliShortRepository+"/";
}

std::string
ModelBase::appliShortRepositoryWithNumProc() const
{
    return this->appliShortRepository() + (boost::format( "np_%1%" ) % Environment::numberOfProcessors() ).str() + "/";
}


// verbose
bool
ModelBase::verbose() const { return M_verbose; }
bool
ModelBase::verboseAllProc() const { return M_verboseAllProc; }

// info
std::string
ModelBase::filenameSaveInfo() const
{
    return M_filenameSaveInfo;
}
void
ModelBase::setFilenameSaveInfo(std::string s)
{
    M_filenameSaveInfo = s;
}
boost::shared_ptr<std::ostringstream>
ModelBase::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    return _ostr;
}
void
ModelBase::printInfo() const
{
    if ( this->verboseAllProc() )
        std::cout << this->getInfo()->str();
    else if (this->worldComm().isMasterRank() )
        std::cout << this->getInfo()->str();
}
void
ModelBase::saveInfo() const
{
    if (this->worldComm().isMasterRank() )
    {
        std::ofstream file(this->filenameSaveInfo().c_str(), std::ios::out);
        file << this->getInfo()->str();
        file.close();
    }
}
void
ModelBase::printAndSaveInfo() const
{
    this->printInfo();
    this->saveInfo();
}

// timer
TimerToolBase &
ModelBase::timerTool(std::string key)
{
    auto itFind = M_mapTimerTool.find( key );
    if ( itFind == M_mapTimerTool.end() )
    {
        this->addTimerTool(key,"applibase-defaultname.data");
        CHECK( M_mapTimerTool.find( key ) != M_mapTimerTool.end() ) << "key not exist";
        return *M_mapTimerTool[key];
    }
    else
        return *itFind->second;
}
void
ModelBase::addTimerTool(std::string key,std::string fileName)
{
    CHECK( M_mapTimerTool.find( key ) == M_mapTimerTool.end() ) << "key already exist";
    if ( M_timersActivated )
    {
        auto myTimerTool = std::make_shared<TimerTool>(fileName,this->worldComm());
        myTimerTool->setReinitSaveFile( this->scalabilityReinitSaveFile() );
        myTimerTool->setSaveFileMasterRank( M_timersSaveFileMasterRank || M_timersSaveFileAll );
        myTimerTool->setSaveFileMax( M_timersSaveFileMax || M_timersSaveFileAll );
        myTimerTool->setSaveFileMin( M_timersSaveFileMin || M_timersSaveFileAll );
        myTimerTool->setSaveFileMean( M_timersSaveFileMean || M_timersSaveFileAll );
        M_mapTimerTool.emplace(std::make_pair( key,myTimerTool ) );
    }
    else
    {
        auto myTimerTool = std::make_shared<TimerToolNull>();
        M_mapTimerTool.emplace(std::make_pair( key,myTimerTool ) );
    }
}

// save assembly/solver scalability
bool
ModelBase::scalabilitySave() const { return M_scalabilitySave; }
bool
ModelBase::scalabilityReinitSaveFile() const { return M_scalabilityReinitSaveFile; }
void
ModelBase::setScalabilitySave( bool b )  { M_scalabilitySave=b; }
std::string
ModelBase::scalabilityPath() const { return M_scalabilityPath; }
void
ModelBase::setScalabilityPath(std::string s) { M_scalabilityPath=s; }
std::string
ModelBase::scalabilityFilename() const { return M_scalabilityFilename; }
void
ModelBase::setScalabilityFilename(std::string s)  { M_scalabilityFilename=s; }



//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
// ModelAlgebraic
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//


ModelAlgebraic::ModelAlgebraic( std::string _theprefix,
                                          WorldComm const& _worldComm,
                                          std::string subPrefix,
                                          std::string appliShortRepository )
    :
    super_type( _theprefix,_worldComm,subPrefix,appliShortRepository),
    M_verboseSolverTimer( boption(_name="verbose_solvertimer",_prefix=this->prefix()) ),
    M_verboseSolverTimerAllProc( boption(_name="verbose_solvertimer_allproc",_prefix=this->prefix()) ),
    M_rebuildCstPartInLinearSystem( boption(_name="linearsystem-cst-update",_prefix=this->prefix()) ),
    M_useLinearJacobianInResidual( boption(_name="residual-uselinearjac",_prefix=this->prefix()) ),
    M_rebuildLinearPartInJacobian( boption(_name="jacobian-linear-update",_prefix=this->prefix()) ),
    M_rebuildCstPartInResidual(true), // not an option (just opitmisation with semi-implicit)
    M_useCstMatrix( boption(_name="use-cst-matrix",_prefix=this->prefix()) ),
    M_useCstVector( boption(_name="use-cst-vector",_prefix=this->prefix()) ),
    M_needToRebuildCstPart( false ),
    M_errorIfSolverNotConverged( boption(_name="error-if-solver-not-converged",_prefix=this->prefix()) ),
    M_printGraph( boption(_name="graph-print-python",_prefix=this->prefix()) )
{

    //-----------------------------------------------------------------------//
    //-----------------------------------------------------------------------//
    if (Environment::vm().count(prefixvm(this->prefix(),"graph-print-python-filename")))
        M_printGraphFileName = Environment::vm()[prefixvm(this->prefix(),"graph-print-python-filename")].as< std::string >();
    else
        M_printGraphFileName = this->prefix()+".graphPython.py";

}
ModelAlgebraic::~ModelAlgebraic()
{}



// verbose
bool
ModelAlgebraic::verboseSolverTimer() const { return M_verboseSolverTimer; }
bool
ModelAlgebraic::verboseSolverTimerAllProc() const { return M_verboseSolverTimerAllProc; }
// do rebuild cst part in linear/jacobian or use jac for residual
bool
ModelAlgebraic::rebuildCstPartInLinearSystem() const { return M_rebuildCstPartInLinearSystem; }
void
ModelAlgebraic::setRebuildCstPartInLinearSystem(bool b) { M_rebuildCstPartInLinearSystem=b; }
bool
ModelAlgebraic::useLinearJacobianInResidual() const { return M_useLinearJacobianInResidual; }
void
ModelAlgebraic::setUseLinearJacobianInResidual(bool b) { M_useLinearJacobianInResidual=b; }
bool
ModelAlgebraic::rebuildLinearPartInJacobian() const { return M_rebuildLinearPartInJacobian; }
void
ModelAlgebraic::setRebuildLinearPartInJacobian(bool b) { M_rebuildLinearPartInJacobian=b; }
// a utiliser avec precaution!!!
bool
ModelAlgebraic::rebuildCstPartInResidual() const { return M_rebuildCstPartInResidual; }
void
ModelAlgebraic::setRebuildCstPartInResidual(bool b) { M_rebuildCstPartInResidual=b; }
// define an other matrix/vector to store the cst part
bool
ModelAlgebraic::useCstMatrix() const { return M_useCstMatrix; }
void
ModelAlgebraic::setUseCstMatrix(bool b) { M_useCstMatrix = b; }
bool
ModelAlgebraic::useCstVector() const { return M_useCstVector; }
void
ModelAlgebraic::setUseCstVector(bool b) { M_useCstVector = b; }
// allow to rebuild cst part (once at next solve) if some parameters (model,time mode,..) change
bool
ModelAlgebraic::needToRebuildCstPart() const { return M_needToRebuildCstPart; }
void
ModelAlgebraic::setNeedToRebuildCstPart(bool b) { M_needToRebuildCstPart = b; }
// an option
bool
ModelAlgebraic::errorIfSolverNotConverged() const { return M_errorIfSolverNotConverged; }
void
ModelAlgebraic::setErrorIfSolverNotConverged( bool b ) { M_errorIfSolverNotConverged= b; }
// save a python script to view graph
bool
ModelAlgebraic::printGraph() const { return M_printGraph; }
void
ModelAlgebraic::setPrintGraph(bool b) { M_printGraph=b; }
std::string
ModelAlgebraic::printGraphFileName() const { return M_printGraphFileName; }
void
ModelAlgebraic::setPrintGraphFileName(std::string s) { M_printGraphFileName=s; }

/**
 * return false
 */
bool
ModelAlgebraic::hasExtendedPattern() const { return false; }

/**
 * return an empty blockPattern if not overhead
 */
ModelAlgebraic::block_pattern_type
ModelAlgebraic::blockPattern() const
{
    return block_pattern_type(0,0);
}

bool
ModelAlgebraic::buildMatrixPrecond() const
{
    return !( Environment::vm()[prefixvm(this->prefix(),"preconditioner.contribution")].as<std::string>() == "same_matrix" );
}

void
ModelAlgebraic::updatePreconditioner(const vector_ptrtype& X,
                                          sparse_matrix_ptrtype& A,
                                          sparse_matrix_ptrtype& A_extended,
                                          sparse_matrix_ptrtype& Prec) const
{
    std::string precType = option(_prefix=this->prefix(),_name="preconditioner.contribution").as<std::string>();

    if( precType =="same_matrix")
    {
        // only copy shrared_ptr (normally already done in constructor)
        Prec=A;
    }
    else if( precType =="standart")
    {
        // copy standart pattern
        Prec->zero();
        Prec->addMatrix(1.,A);
    }
    else if( precType =="extended" )
    {
        // copy standart and extended pattern
        Prec->zero();
        Prec->addMatrix(1.,A);
        if (hasExtendedPattern())
            Prec->addMatrix(1.,A_extended);
    }
}

ModelAlgebraic::graph_ptrtype
ModelAlgebraic::buildMatrixGraph() const
{
    return graph_ptrtype();
}

void
ModelAlgebraic::updateCLDirichlet(vector_ptrtype& U) const {} // const = 0;
void
ModelAlgebraic::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J , vector_ptrtype& R,
                                     bool BuildCstPart,
                                     sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                                     bool _doClose, bool _doBCStrongDirichlet) const {}// = 0;
void
ModelAlgebraic::updateResidual( const vector_ptrtype& X, vector_ptrtype& R,
                                     bool BuildCstPart, bool UseJacobianLinearTerms,
                                     bool _doClose, bool _doBCStrongDirichlet ) const {}// = 0;

void
ModelAlgebraic::updateLinearPDE(const vector_ptrtype& X,sparse_matrix_ptrtype& A , vector_ptrtype& F,
                                     bool _buildCstPart,
                                     sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                                     bool _doClose, bool _doBCStrongDirichlet ) const {}// = 0;




} // namespace FeelModels

} // namespace Feel
