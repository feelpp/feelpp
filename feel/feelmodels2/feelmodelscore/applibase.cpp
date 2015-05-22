
#include <feel/feelmodels2/feelmodelscore/applibase.hpp>

namespace Feel {

namespace FeelModels {


AppliBase::AppliBase( std::string _theprefix,
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
AppliBase::~AppliBase()
{}

WorldComm const&
AppliBase::worldComm() const { return M_worldComm; }
std::vector<WorldComm> const&
AppliBase::worldsComm() const { return M_worldsComm; }
void
AppliBase::setWorldsComm(std::vector<WorldComm> const& _worldsComm) { M_worldsComm=_worldsComm; }
std::vector<WorldComm> const&
AppliBase::localNonCompositeWorldsComm() const { return M_localNonCompositeWorldsComm; }
void
AppliBase::setLocalNonCompositeWorldsComm(std::vector<WorldComm> const& _worldsComm) { M_localNonCompositeWorldsComm=_worldsComm; }
void
AppliBase::createWorldsComm() {}//warning

std::string
AppliBase::prefix() const { return M_prefix; }
std::string
AppliBase::subPrefix() const { return M_subPrefix; }

po::variables_map const&
AppliBase::vm() const { return Environment::vm(); }

std::string
AppliBase::appliRepositoryWithoutNumProc() const
{
    return Environment::rootRepository()+"/"+
        this->appliShortRepository();
}

std::string
AppliBase::appliRepository() const
{
    return Environment::rootRepository() + "/" + this->appliShortRepositoryWithNumProc();
}

std::string
AppliBase::appliShortRepository() const
{
    //return option(_name="exporter.directory").as<std::string>()+"/";
    return M_appliShortRepository+"/";
}

std::string
AppliBase::appliShortRepositoryWithNumProc() const
{
    return this->appliShortRepository() + (boost::format( "np_%1%" ) % Environment::numberOfProcessors() ).str() + "/";
}


// verbose
bool
AppliBase::verbose() const { return M_verbose; }
bool
AppliBase::verboseAllProc() const { return M_verboseAllProc; }

// info
std::string
AppliBase::filenameSaveInfo() const
{
    return M_filenameSaveInfo;
}
void
AppliBase::setFilenameSaveInfo(std::string s)
{
    M_filenameSaveInfo = s;
}
boost::shared_ptr<std::ostringstream>
AppliBase::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    return _ostr;
}
void
AppliBase::printInfo() const
{
    if ( this->verboseAllProc() )
        std::cout << this->getInfo()->str();
    else if (this->worldComm().isMasterRank() )
        std::cout << this->getInfo()->str();
}
void
AppliBase::saveInfo() const
{
    if (this->worldComm().isMasterRank() )
    {
        std::ofstream file(this->filenameSaveInfo().c_str(), std::ios::out);
        file << this->getInfo()->str();
        file.close();
    }
}
void
AppliBase::printAndSaveInfo() const
{
    this->printInfo();
    this->saveInfo();
}

// timer
TimerToolBase &
AppliBase::timerTool(std::string key)
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
AppliBase::addTimerTool(std::string key,std::string fileName)
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
AppliBase::scalabilitySave() const { return M_scalabilitySave; }
bool
AppliBase::scalabilityReinitSaveFile() const { return M_scalabilityReinitSaveFile; }
void
AppliBase::setScalabilitySave( bool b )  { M_scalabilitySave=b; }
std::string
AppliBase::scalabilityPath() const { return M_scalabilityPath; }
void
AppliBase::setScalabilityPath(std::string s) { M_scalabilityPath=s; }
std::string
AppliBase::scalabilityFilename() const { return M_scalabilityFilename; }
void
AppliBase::setScalabilityFilename(std::string s)  { M_scalabilityFilename=s; }



//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
// AppliBaseMethodsNum
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//


AppliBaseMethodsNum::AppliBaseMethodsNum( std::string _theprefix,
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
AppliBaseMethodsNum::~AppliBaseMethodsNum()
{}



// verbose
bool
AppliBaseMethodsNum::verboseSolverTimer() const { return M_verboseSolverTimer; }
bool
AppliBaseMethodsNum::verboseSolverTimerAllProc() const { return M_verboseSolverTimerAllProc; }
// do rebuild cst part in linear/jacobian or use jac for residual
bool
AppliBaseMethodsNum::rebuildCstPartInLinearSystem() const { return M_rebuildCstPartInLinearSystem; }
void
AppliBaseMethodsNum::setRebuildCstPartInLinearSystem(bool b) { M_rebuildCstPartInLinearSystem=b; }
bool
AppliBaseMethodsNum::useLinearJacobianInResidual() const { return M_useLinearJacobianInResidual; }
void
AppliBaseMethodsNum::setUseLinearJacobianInResidual(bool b) { M_useLinearJacobianInResidual=b; }
bool
AppliBaseMethodsNum::rebuildLinearPartInJacobian() const { return M_rebuildLinearPartInJacobian; }
void
AppliBaseMethodsNum::setRebuildLinearPartInJacobian(bool b) { M_rebuildLinearPartInJacobian=b; }
// a utiliser avec precaution!!!
bool
AppliBaseMethodsNum::rebuildCstPartInResidual() const { return M_rebuildCstPartInResidual; }
void
AppliBaseMethodsNum::setRebuildCstPartInResidual(bool b) { M_rebuildCstPartInResidual=b; }
// define an other matrix/vector to store the cst part
bool
AppliBaseMethodsNum::useCstMatrix() const { return M_useCstMatrix; }
void
AppliBaseMethodsNum::setUseCstMatrix(bool b) { M_useCstMatrix = b; }
bool
AppliBaseMethodsNum::useCstVector() const { return M_useCstVector; }
void
AppliBaseMethodsNum::setUseCstVector(bool b) { M_useCstVector = b; }
// allow to rebuild cst part (once at next solve) if some parameters (model,time mode,..) change
bool
AppliBaseMethodsNum::needToRebuildCstPart() const { return M_needToRebuildCstPart; }
void
AppliBaseMethodsNum::setNeedToRebuildCstPart(bool b) { M_needToRebuildCstPart = b; }
// an option
bool
AppliBaseMethodsNum::errorIfSolverNotConverged() const { return M_errorIfSolverNotConverged; }
void
AppliBaseMethodsNum::setErrorIfSolverNotConverged( bool b ) { M_errorIfSolverNotConverged= b; }
// save a python script to view graph
bool
AppliBaseMethodsNum::printGraph() const { return M_printGraph; }
void
AppliBaseMethodsNum::setPrintGraph(bool b) { M_printGraph=b; }
std::string
AppliBaseMethodsNum::printGraphFileName() const { return M_printGraphFileName; }
void
AppliBaseMethodsNum::setPrintGraphFileName(std::string s) { M_printGraphFileName=s; }

/**
 * return false
 */
bool
AppliBaseMethodsNum::hasExtendedPattern() const { return false; }

/**
 * return an empty blockPattern if not overhead
 */
AppliBaseMethodsNum::block_pattern_type
AppliBaseMethodsNum::blockPattern() const
{
    return block_pattern_type(0,0);
}

bool
AppliBaseMethodsNum::buildMatrixPrecond() const
{
    return !( Environment::vm()[prefixvm(this->prefix(),"preconditioner.contribution")].as<std::string>() == "same_matrix" );
}

void
AppliBaseMethodsNum::updatePreconditioner(const vector_ptrtype& X,
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

AppliBaseMethodsNum::graph_ptrtype
AppliBaseMethodsNum::buildMatrixGraph() const
{
    return graph_ptrtype();
}

void
AppliBaseMethodsNum::updateCLDirichlet(vector_ptrtype& U) const {} // const = 0;
void
AppliBaseMethodsNum::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J , vector_ptrtype& R,
                                     bool BuildCstPart,
                                     sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                                     bool _doClose, bool _doBCStrongDirichlet) const {}// = 0;
void
AppliBaseMethodsNum::updateResidual( const vector_ptrtype& X, vector_ptrtype& R,
                                     bool BuildCstPart, bool UseJacobianLinearTerms,
                                     bool _doClose, bool _doBCStrongDirichlet ) const {}// = 0;

void
AppliBaseMethodsNum::updateLinearPDE(const vector_ptrtype& X,sparse_matrix_ptrtype& A , vector_ptrtype& F,
                                     bool _buildCstPart,
                                     sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                                     bool _doClose, bool _doBCStrongDirichlet ) const {}// = 0;




} // namespace FeelModels

} // namespace Feel
