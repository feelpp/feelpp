#include <feel/feelmodels2/modelcore/applibasenumericalsimulationtransitory.hpp>

namespace Feel
{
namespace FeelModels
{

    ModelNumerical::ModelNumerical(bool _isStationary, std::string _theprefix, WorldComm const& _worldComm, std::string subPrefix,
                                                                                   std::string appliShortRepository )
        :
        super_type( _theprefix, _worldComm, subPrefix, appliShortRepository ),
        M_rebuildMeshPartitions( option(_name="rebuild_mesh_partitions",_prefix=this->prefix()).as<bool>() ),
        M_isStationary(_isStationary),
        M_doRestart( option(_name="bdf.restart").as<bool>() ),
        M_restartPath( option(_name="bdf.restart.path").as<std::string>() ),
        M_restartAtLastSave( option(_name="bdf.restart.at-last-save").as<bool>() ),
        M_timeInitial( option(_name="bdf.time-initial").as<double>() ),
        M_timeFinal( option(_name="bdf.time-final").as<double>() ),
        M_timeStep( option(_name="bdf.time-step").as<double>() ),
        M_bdfSaveInFile( option(_name="bdf.save").as<bool>() ),
        M_bdfSaveFreq( option(_name="bdf.save.freq").as<int>() ),
        M_timeCurrent(M_timeInitial),
        M_modelProps( Environment::expand( soption( _name=prefixvm(this->prefix(),"filename")) ) ),
        M_parameters(std::vector<double>(FEELMODELS_FSIBASE_NUMBER_OF_PARAMETERS,0)),
        M_geoParameters(std::vector<double>(FEELMODELS_FSIBASE_NUMBER_OF_GEOPARAMETERS,0)),
        M_ginacExpr(std::vector<std::pair<std::string,std::string> >(FEELMODELS_FSIBASE_NUMBER_OF_GINACEXPR)),
        M_ginacExprCompilationDirectory( "" ),
        M_geotoolMeshIndex( option(_name="geotool-mesh-index",_prefix=this->prefix()).as<int>() ),
        M_geotoolSaveDirectory( option(_name="geotool-save-directory",_prefix=this->prefix()).as<std::string>() ),
        M_geotoolSaveName( option(_name="geotool-save-name",_prefix=this->prefix()).as<std::string>() ),
        M_row_startInMatrix(0),
        M_col_startInMatrix(0),
        M_row_startInVector(0),
        M_mshFileStr("FEELMODELS_WARNING_NODEFINE"),
        M_geoFileStr("FEELMODELS_WARNING_NODEFINE"),
        M_exporterPath( this->appliRepository()+"/"+prefixvm(this->prefix(), prefixvm(this->subPrefix(),"exports")) )
        //M_PsLogger( new PsLogger(prefixvm(this->prefix(),"PsLogger"),this->worldComm() ) )
    {
        //-----------------------------------------------------------------------//
        // init user cst parameters
        for ( uint16_type k=1;k<=M_parameters.size();++k )
        {
            double val = doption(_prefix=this->prefix(),_name=(boost::format("parameter%1%") %k ).str());
            this->setUserCstParameter(k,val);
        }
        // init user cst geo parameters
        for ( uint16_type k=1;k<=M_geoParameters.size();++k )
        {
            double val = doption(_prefix=this->prefix(),_name=(boost::format("geo-parameter%1%") %k ).str());
            this->setUserCstGeoParameter(k,val);
        }


        // init user ginac expr
        for ( uint16_type k=1;k<M_ginacExpr.size();++k )
            {
                std::string gexpr= option(_prefix=this->prefix(),_name=(boost::format("ginac-expr%1%") %k ).str()).as<std::string>();
                std::string gname= option(_prefix=this->prefix(),_name=(boost::format("ginac-name%1%") %k ).str()).as<std::string>();
                this->setUserGinacExpr(k,gexpr,gname);
            }
        if ( Environment::vm().count(prefixvm(this->prefix(),"ginac-expr-directory").c_str()) )
        {
            M_ginacExprCompilationDirectory=Environment::rootRepository()+"/"+soption(_name="ginac-expr-directory",_prefix=this->prefix());
        }
        //-----------------------------------------------------------------------//
        // mesh file : .msh
        if (Environment::vm().count(prefixvm(this->prefix(),"mshfile").c_str()))
            M_mshFileStr = Environment::vm()[prefixvm(this->prefix(),"mshfile")].as< std::string >();
        // mesh file : .geo
        if (Environment::vm().count(prefixvm(this->prefix(),"geofile").c_str()))
            M_geoFileStr = Environment::vm()[prefixvm(this->prefix(),"geofile")].as< std::string >();
        // mesh file : geotool with .mesh file
        if (M_geotoolSaveDirectory.empty()) M_geotoolSaveDirectory = this->appliShortRepository();//this->appliRepository();
        if (M_geotoolSaveName.empty()) M_geotoolSaveName = this->prefix();
        //-----------------------------------------------------------------------//
        if (Environment::vm()[prefixvm(this->prefix(),"geomap")].as<std::string>()=="opt")
            M_geomap=GeomapStrategyType::GEOMAP_OPT;
        else
            M_geomap=GeomapStrategyType::GEOMAP_HO;
        //-----------------------------------------------------------------------//
    }


    void
    ModelNumerical::setStationary(bool b)
    {
        if ( M_isStationary != b)
        {
            M_isStationary=b;
            this->setNeedToRebuildCstPart(true);
        }
    }

    // ginac expr
    Expr< GinacEx<2> >
    ModelNumerical::userGinacExpr(uint16_type i, std::map<std::string,double> const& mp ) const
    {
        CHECK( i >=1 && i <=M_ginacExpr.size() ) << "invalid index\n";
        std::string pathGinacExpr = (this->ginacExprCompilationDirectory().empty()) ?
            this->userGinacExprName(i) :
            this->ginacExprCompilationDirectory() + "/" + this->userGinacExprName(i) ;
        return expr( this->userGinacExprStr(i), mp , pathGinacExpr );
    }
    Expr< GinacEx<2> >
    ModelNumerical::userGinacExpr(uint16_type i, std::pair<std::string,double> const& mp ) const
    {
        return this->userGinacExpr( i,{ { mp.first, mp.second } } );
    }
    std::string
    ModelNumerical::userGinacExprStr(uint16_type i) const
    {
        CHECK( i >=1 && i <=M_ginacExpr.size() ) << "invalid index\n";
        return M_ginacExpr[i-1].first;
    }
    std::string
    ModelNumerical::userGinacExprName(uint16_type i) const
    {
        CHECK( i >=1 && i <=M_ginacExpr.size() ) << "invalid index\n";
        return M_ginacExpr[i-1].second;
    }
    void
    ModelNumerical::setUserGinacExpr(uint16_type i,std::string expr)
    {
        CHECK( i >=1 && i <=M_ginacExpr.size() ) << "invalid index\n";
        M_ginacExpr[i-1]=std::make_pair(expr,(boost::format("defaultNameGinacExpr%1%")%i).str());
    }
    void
    ModelNumerical::setUserGinacExpr(uint16_type i,std::string expr,std::string name)
    {
        CHECK( i >=1 && i <=M_ginacExpr.size() ) << "invalid index\n";
        M_ginacExpr[i-1]=std::make_pair(expr,name);
    }
    std::string
    ModelNumerical::ginacExprCompilationDirectory() const { return M_ginacExprCompilationDirectory; }

    void
    ModelNumerical::saveMSHfilePath(std::string namePath) const
    {
        if ( this->worldComm().isMasterRank() )
        {
            //std::string nameFile = prefixvm(this->prefix(),namePath);
            std::ofstream file(namePath.c_str(), std::ios::out);
            //M_mshFileStr = this->application()->vm()[prefixvm(this->prefix(),"mshfile")].as< std::string >() ;
            file << M_mshFileStr;
            file.close();
        }
    }


} // namespace FeelModels

} // namespace Feel

