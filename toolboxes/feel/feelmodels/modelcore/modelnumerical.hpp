/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>,
            Thibaut Metivet <thibaut.metivet@inria.fr>
 Date: 2012-01-19

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
 Copyright (C) 2020 Inria

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3.0 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
/**
 \file modelnumerical.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \author Thibaut Metivet <thibaut.metivet@inria.fr>
 \date 2012-01-19
 */

#ifndef FEELPP_MODELNUMERICAL_HPP
#define FEELPP_MODELNUMERICAL_HPP 1

#include <feel/feelmodels/modelcore/modelalgebraic.hpp>

#include <feel/feelpoly/geomap.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <feel/feelvf/ginac.hpp>

#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelmodels/modelcore/modelmeasures.hpp>
#include <feel/feelfit/fit.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>

#include <feel/feelcore/tuple_utils.hpp>

#include <feel/feelmodels/modelcore/modelcontext.hpp>

#include <feel/feelmodels/modelcore/modelmeasuresnormevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresstatisticsevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasurespointsevaluation.hpp>

namespace Feel
{
namespace FeelModels
{
/**
 * Handles some numerical model aspects: timestepping, mesh and properties
 */
class ModelNumerical : public ModelAlgebraic
    {
    public:
        typedef ModelAlgebraic super_type;

        static const bool is_class_null = false;

        typedef double value_type;
        using index_type = typename super_type::index_type;
        using size_type = typename super_type::size_type;
        typedef super_type::backend_type backend_type;
        typedef super_type::backend_ptrtype backend_ptrtype;

        typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
        typedef backend_type::vector_ptrtype vector_ptrtype;

        typedef backend_type::indexsplit_type indexsplit_type;
        typedef backend_type::indexsplit_ptrtype indexsplit_ptrtype;

        typedef vf::BlocksBase<size_type> block_pattern_type;


        ModelNumerical( std::string const& _theprefix, std::string const& keyword,
                        worldcomm_ptr_t const& _worldComm=Environment::worldCommPtr(),
                        std::string const& subPrefix="",
                        ModelBaseRepository const& modelRep = ModelBaseRepository(),
                        ModelBaseCommandLineOptions const& modelCmdLineOpt = ModelBaseCommandLineOptions() );
        ModelNumerical( std::string const& _theprefix, worldcomm_ptr_t const& _worldComm=Environment::worldCommPtr(),
                        std::string const& subPrefix="",
                        ModelBaseRepository const& modelRep = ModelBaseRepository(),
                        ModelBaseCommandLineOptions const& modelCmdLineOpt = ModelBaseCommandLineOptions() )
            :
            ModelNumerical( _theprefix, _theprefix, _worldComm, subPrefix, modelRep, modelCmdLineOpt )
            {}

        ModelNumerical( ModelNumerical const& app ) = default;

        virtual ~ModelNumerical() {};

        //std::shared_ptr<PsLogger> psLogger()  { return M_PsLogger; }
        //std::shared_ptr<PsLogger> const& psLogger() const { return M_PsLogger; }

        bool isStationary() const { return M_isStationary; }
        void setStationary(bool b);

        bool doRestart() const { return M_doRestart; }
        void setRestart(bool b) { M_doRestart=b; }
        std::string restartPath() const { return M_restartPath; }
        void setRestartPath(std::string const& s) { M_restartPath=s; }
        bool restartAtLastSave() const { return M_restartAtLastSave; }
        void setRestartAtLastSave( bool b) { M_restartAtLastSave=b; }

        double time() const { return this->currentTime(); }
        double currentTime() const { return M_timeCurrent; }
        void updateTime(double t);

        double timeInitial() const { return M_timeInitial; }
        double timeFinal() const { return M_timeFinal; }
        double timeStep() const { return M_timeStep; }
        int timeOrder() const { return M_timeOrder; }
        bool tsSaveInFile() const { return M_tsSaveInFile; }
        int tsSaveFreq() const { return M_tsSaveFreq; }
        void setTimeInitial(double v)  { M_timeInitial=v; }
        void setTimeFinal(double v)  { M_timeFinal=v; }
        void setTimeStep(double v)  { M_timeStep=v; }
        void setTimeOrder(int o) { M_timeOrder=o; }


        template<typename SpaceType>
        auto createBdf( std::shared_ptr<SpaceType> space, std::string const& name, int bdfOrder, int nConsecutiveSave, std::string const& myFileFormat ) const
            {
                std::string suffixName = "";
                if ( myFileFormat == "binary" )
                    suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
                fs::path saveTsDir = fs::path(this->rootRepository())/fs::path( prefixvm(this->prefix(),prefixvm(this->subPrefix(),"ts")) );

                double ti = this->timeInitial();
                double tf = this->timeFinal();
                double dt = this->timeStep();

                auto thebdf = bdf( _space=space,
                                   _vm=this->clovm(),
                                   _name=name+suffixName,
                                   _prefix=this->prefix(),
                                   _order=bdfOrder,
                                   // don't use the fluid.bdf {initial,final,step}time but the general bdf info, the order will be from fluid.bdf
                                   _initial_time=ti, _final_time=tf, _time_step=dt,
                                   _restart=this->doRestart(),
                                   _restart_path=this->restartPath(),
                                   _restart_at_last_save=this->restartAtLastSave(),
                                   _save=this->tsSaveInFile(), _format=myFileFormat, _freq=this->tsSaveFreq(),
                                   _n_consecutive_save=nConsecutiveSave );
                thebdf->setfileFormat( myFileFormat );
                thebdf->setPathSave( ( saveTsDir/name ).string() );
                return thebdf;
            }

        bool hasModelProperties() const { return (M_modelProps)? true : false; }
        std::shared_ptr<ModelProperties> modelPropertiesPtr() const { return M_modelProps; }
        ModelProperties const& modelProperties() const { return *M_modelProps; }
        ModelProperties & modelProperties() { return *M_modelProps; }
        void setModelProperties( std::shared_ptr<ModelProperties> modelProps ) { M_modelProps = modelProps; }
        void addParameterInModelProperties( std::string const& symbolName,double value );

        bool manageParameterValues() const { return M_manageParameterValues; }
        void setManageParameterValues( bool b ) { M_manageParameterValues = b; }
        bool manageParameterValuesOfModelProperties() const { return M_manageParameterValuesOfModelProperties; }
        void setManageParameterValuesOfModelProperties( bool b ) { M_manageParameterValuesOfModelProperties = b; }


        GeomapStrategyType geomap() const { return M_geomap; }

        //----------------------------------------------------------------------------------//


        std::string meshFile() const { return M_meshFile; }
        void setMeshFile(std::string const& file)  { M_meshFile=file; }
        std::string geoFile() const { return M_geoFile; }
        void setGeoFile(std::string const& file)  { M_geoFile=file; }
        bool hasMeshFile() const { return !M_meshFile.empty(); }
        bool hasGeoFile() const { return !M_geoFile.empty(); }

        void saveMeshFile( std::string const& fileSavePath, std::string const& meshPath = "" ) const;

        void setExporterPath(std::string const& s)  { M_exporterPath=s; }
        std::string exporterPath() const { return M_exporterPath; }

        bool hasPostProcessExportsField( std::string const& fieldName ) const { return hasPostProcessExportsField( "",fieldName ); }
        bool hasPostProcessExportsField(  std::string const& exportTag, std::string const& fieldName ) const
            {
                auto itFindTag = M_postProcessExportsFields.find( exportTag );
                if ( itFindTag == M_postProcessExportsFields.end() )
                    return false;
                return std::get<0>( itFindTag->second ).find( fieldName ) != std::get<0>( itFindTag->second ).end();
            }
        std::set<std::string> const& postProcessExportsFields( std::string const& exportTag = "" ) const { return std::get<0>( M_postProcessExportsFields.find( exportTag )->second ); }
        std::set<std::string> const& postProcessExportsAllFieldsAvailable( std::string const& exportTag = "" ) const { return std::get<1>( M_postProcessExportsFields.find( exportTag )->second ); }
        std::string const& postProcessExportsPidName( std::string const& exportTag = "" ) const { return std::get<2>( M_postProcessExportsFields.find( exportTag )->second ); }
        std::set<std::string> const& postProcessSaveFields() const { return M_postProcessSaveFields; }
        fs::path const& postProcessSaveRepository() const { return M_postProcessSaveRepository; }
        std::set<std::string> postProcessExportsFields( std::string const& tag, std::set<std::string> const& ifields, std::string const& prefix = "" ) const;
        std::set<std::string> postProcessSaveFields( std::set<std::string> const& ifields, std::string const& prefix = "" ) const;

        template <typename ExporterType,typename ModelFieldsType, typename SymbolsExprType = symbols_expression_empty_t, typename TupleExprOnRangeType = hana::tuple<> >
        void executePostProcessExports( std::shared_ptr<ExporterType> exporter, std::string const& tag, double time, ModelFieldsType const& tupleFields,
                                        SymbolsExprType const& symbolsExpr = symbols_expression_empty_t{}, TupleExprOnRangeType const& tupleExprOnRange = hana::make_tuple() );
        template <typename ExporterType,typename ModelFieldsType, typename SymbolsExprType = symbols_expression_empty_t, typename TupleExprOnRangeType = hana::tuple<> >
        void executePostProcessExports( std::shared_ptr<ExporterType> exporter, double time, ModelFieldsType const& tupleFields,
                                        SymbolsExprType const& symbolsExpr = symbols_expression_empty_t{}, TupleExprOnRangeType const& tupleExprOnRange = hana::make_tuple() )
            {
                this->executePostProcessExports(exporter,"",time,tupleFields,symbolsExpr,tupleExprOnRange);
            }
        template <typename ExporterType,typename ModelFieldsType, typename SymbolsExprType, typename TupleExprOnRangeType>
        bool updatePostProcessExports( std::shared_ptr<ExporterType> exporter, std::set<std::string> const& fields, double time, ModelFieldsType const& tupleFields,
                                       SymbolsExprType const& symbolsExpr, TupleExprOnRangeType const& tupleExprOnRange );

        template <typename MeshType, typename RangeType, typename MeasurePointEvalType, typename SymbolsExpr, typename ModelFieldsType, typename TupleQuantitiesType>
        void executePostProcessMeasures( double time, std::shared_ptr<MeshType> mesh, RangeType const& rangeMeshElements, std::shared_ptr<MeasurePointEvalType> measurePointsEvaluation, SymbolsExpr const& symbolsExpr, ModelFieldsType const& tupleFields, TupleQuantitiesType const& tupleQuantities );
        template<typename TupleQuantitiesType>
        bool updatePostProcessMeasuresQuantities( TupleQuantitiesType const& tupleQuantities );
        template <typename MeshType, typename RangeType, typename SymbolsExpr, typename ModelFieldsType>
        bool updatePostProcessMeasuresNorm( std::shared_ptr<MeshType> mesh, RangeType const& rangeMeshElements, SymbolsExpr const& symbolsExpr, ModelFieldsType const& tupleFields );
        template <typename MeshType, typename RangeType, typename SymbolsExpr, typename ModelFieldsType>
        bool updatePostProcessMeasuresStatistics( std::shared_ptr<MeshType> mesh, RangeType const& rangeMeshElements, SymbolsExpr const& symbolsExpr, ModelFieldsType const& tupleFields );
        template <typename MeasurePointEvalType, typename ModelFieldsType>
        bool updatePostProcessMeasuresPoint( std::shared_ptr<MeasurePointEvalType> measurePointsEvaluation, ModelFieldsType const& tupleFields );
        template <typename MeshType, typename RangeType, typename MeasurePointEvalType, typename SymbolsExpr, typename ModelFieldsType, typename TupleQuantitiesType>
        bool updatePostProcessMeasures( std::shared_ptr<MeshType> mesh, RangeType const& rangeMeshElements, std::shared_ptr<MeasurePointEvalType> measurePointsEvaluation, SymbolsExpr const& symbolsExpr, ModelFieldsType const& tupleFields, TupleQuantitiesType const& tupleQuantities );

        template <typename ModelFieldsType>
        void executePostProcessSave( uint32_type index, ModelFieldsType const& fields )
            {
                this->executePostProcessSave( this->postProcessSaveFields(), M_postProcessSaveFieldsFormat, index, fields );
            }
        template <typename ModelFieldsType>
        void executePostProcessSave( std::set<std::string> const& fieldsNamesToSave, std::string const& format, uint32_type index, ModelFieldsType const& fields );

        ModelMeasuresIO const& postProcessMeasuresIO() const { return M_postProcessMeasuresIO; }
        ModelMeasuresIO & postProcessMeasuresIO() { return M_postProcessMeasuresIO; }
        ModelMeasuresEvaluatorContext const& postProcessMeasuresEvaluatorContext() const { return M_postProcessMeasuresEvaluatorContext; }
        ModelMeasuresEvaluatorContext & postProcessMeasuresEvaluatorContext() { return M_postProcessMeasuresEvaluatorContext; }

        virtual bool checkResults() const;

    protected :

        void setPostProcessExportsAllFieldsAvailable( std::set<std::string> const& ifields ) { this->setPostProcessExportsAllFieldsAvailable( "", ifields ); }
        void setPostProcessExportsAllFieldsAvailable( std::string const& exportTag, std::set<std::string> const& ifields ) { std::get<1>( M_postProcessExportsFields[exportTag] ) = ifields; }
        void addPostProcessExportsAllFieldsAvailable( std::set<std::string> const& ifields ) { this->addPostProcessExportsAllFieldsAvailable( "", ifields ); }
        void addPostProcessExportsAllFieldsAvailable( std::string const& exportTag, std::set<std::string> const& ifields ) { std::get<1>( M_postProcessExportsFields[exportTag] ).insert( ifields.begin(), ifields.end() ); }

        void setPostProcessExportsPidName( std::string const& pidName ) { this->setPostProcessExportsPidName( "", pidName ); }
        void setPostProcessExportsPidName( std::string const& exportTag,std::string const& pidName ) { std::get<2>( M_postProcessExportsFields[exportTag] ) = pidName; }
        void setPostProcessSaveAllFieldsAvailable( std::set<std::string> const& ifields ) { M_postProcessSaveAllFieldsAvailable = ifields; }
        virtual void initPostProcess();

        auto symbolsExprParameter() const
            {
                if ( this->hasModelProperties() )
                    return this->modelProperties().parameters().symbolsExpr();
                else
                    return std::decay_t<decltype(this->modelProperties().parameters().symbolsExpr())>{};
            }

        template <typename ElementType, typename RangeType, typename SymbolsExpr>
        void
        updateInitialConditions( std::string const& fieldName,  RangeType const& defaultRange, SymbolsExpr const& symbolsExpr,
                                 std::vector< std::shared_ptr<ElementType> > & dataToUpdate )
            {
                ModelNumerical::updateInitialConditions( this->modelProperties().initialConditions().get( fieldName ),
                                                         defaultRange, symbolsExpr, this->geomap(), dataToUpdate );
            }
    private :
        template <typename ElementType, typename RangeType, typename SymbolsExpr>
        static
        void
        updateInitialConditions( ModelInitialConditionTimeSet const& icts, RangeType const& defaultRange, SymbolsExpr const& symbolsExpr,
                                 GeomapStrategyType geomapStrategy, std::vector< std::shared_ptr<ElementType> > & dataToUpdate );

    private :

        bool M_isStationary;
        bool M_doRestart;
        std::string M_restartPath;
        bool M_restartAtLastSave;

        double M_timeInitial,M_timeFinal,M_timeStep;
        int M_timeOrder;
        bool M_tsSaveInFile;
        int M_tsSaveFreq;
        double M_timeCurrent;

        std::shared_ptr<ModelProperties> M_modelProps;
        bool M_manageParameterValues, M_manageParameterValuesOfModelProperties;


        std::string M_meshFile, M_geoFile;

        std::string M_exporterPath;
        std::map<std::string,std::tuple< std::set<std::string>, std::set<std::string>, std::string > > M_postProcessExportsFields; // (fields, allFieldsAvailable,pidName)
        std::set<std::string> M_postProcessSaveFields, M_postProcessSaveAllFieldsAvailable;
        std::string M_postProcessSaveFieldsFormat;
        fs::path M_postProcessSaveRepository;
        ModelMeasuresIO M_postProcessMeasuresIO;
        ModelMeasuresEvaluatorContext M_postProcessMeasuresEvaluatorContext;

        //std::shared_ptr<PsLogger> M_PsLogger;

        GeomapStrategyType M_geomap;

        bool M_useChecker;
    };


template <typename ElementType, typename RangeType, typename SymbolsExpr>
void
ModelNumerical::updateInitialConditions( ModelInitialConditionTimeSet const& icts, RangeType const& defaultRange, SymbolsExpr const& symbolsExpr,
                                         GeomapStrategyType geomapStrategy, std::vector< std::shared_ptr<ElementType> > & dataToUpdate )
{
    if ( icts.empty() )
        return;

    CHECK( icts.size() == 1 ) << "TODO";
    CHECK( !dataToUpdate.empty() ) << "require a non empty vector";
    for( auto const& [time,icByType] : icts )
    {
        auto & u = *dataToUpdate[0];
        auto itFindIcFile = icByType.find( "File" );
        if ( itFindIcFile != icByType.end() )
        {
            for ( auto const& ic : itFindIcFile->second )
            {
                CHECK( ic.isFile() ) << "must be an ic file";
                CHECK( !ic.fileName().empty() || !ic.fileDirectory().empty() ) << "entry filename or directory is required";
                std::string fileName = ic.fileName();
                std::string fileType = ic.fileType();
                if ( fileType.empty() )
                    fileType = "default";
                std::string fileDirectory = ic.fileDirectory();
                if ( fileDirectory.empty() )
                    u.load(_path=fileName,_type=fileType);
                else
                {
                    if ( fileName.empty() )
                        fileName = u.name();
                    u.load(_path=fileDirectory,_type=fileType,_name=fileName);
                }
            }
        }

        auto itFindIcExpr = icByType.find( "Expression" );
        if ( itFindIcExpr != icByType.end() )
        {
            for ( auto const& ic : itFindIcExpr->second )
            {
                CHECK( ic.isExpression() ) << "must be an ic expression";
                if ( !ic.expression().template hasExpr<ElementType::nComponents1,ElementType::nComponents2>() )
                    CHECK( false ) << "must be a scalar expression";
                auto theExpr = expr( ic.expression().template expr<ElementType::nComponents1,ElementType::nComponents2>(), symbolsExpr );

                if ( ic.markers().empty() )
                    u.on(_range=defaultRange,_expr=theExpr,_geomap=geomapStrategy);
                else
                {
                    auto markersByEntity = Feel::FeelModels::detail::distributeMarkerListOnSubEntity(u.mesh(),ic.markers());
                    auto const& listMarkerFaces = std::get<0>( markersByEntity );
                    auto const& listMarkerEdges = std::get<1>( markersByEntity );
                    auto const& listMarkerPoints = std::get<2>( markersByEntity );
                    auto const& listMarkerElements = std::get<3>( markersByEntity );
                    if ( !listMarkerElements.empty() )
                        u.on(_range=markedelements(u.mesh(),listMarkerElements),_expr=theExpr,_geomap=geomapStrategy);
                    if ( !listMarkerFaces.empty() )
                        u.on(_range=markedfaces(u.mesh(),listMarkerFaces),_expr=theExpr,_geomap=geomapStrategy);
                    if ( !listMarkerEdges.empty() )
                        u.on(_range=markededges(u.mesh(),listMarkerEdges),_expr=theExpr,_geomap=geomapStrategy);
                    if ( !listMarkerPoints.empty() )
                        u.on(_range=markedpoints(u.mesh(),listMarkerPoints),_expr=theExpr,_geomap=geomapStrategy);
                }
            }
        }
    }

    for (int k=1;k<dataToUpdate.size();++k)
        *dataToUpdate[k] = *dataToUpdate[0];
}

template <typename ExporterType,typename ModelFieldsType,typename SymbolsExpr,typename TupleExprOnRangeType>
void
ModelNumerical::executePostProcessExports( std::shared_ptr<ExporterType> exporter, std::string const& tag, double time, ModelFieldsType const& tupleFields, SymbolsExpr const& symbolsExpr, TupleExprOnRangeType const& tupleExprOnRange )
{
    if ( !exporter ) return;
    if ( !exporter->doExport() ) return;

    std::set<std::string> /*const&*/ fieldsNamesToExport = this->postProcessExportsFields( tag );
#if 1
    std::map<std::string,std::vector<std::tuple<ModelExpression, elements_reference_wrapper_t<typename ExporterType::mesh_type>, std::set<std::string> > > > mapExportsExpr;
    if ( this->hasModelProperties() )
    {
        auto mesh = exporter->defaultTimeSet()->mesh();
        auto defaultRange = elements(mesh);
        for ( auto const& [_name,_expr,_markers,_rep,_tag] : this->modelProperties().postProcess().exports( this->keyword() ).expressions() )
        {
            if ( !_tag.empty() && _tag.find( tag ) == _tag.end() )
                continue;
            if ( _tag.empty() && tag != "" )
                continue;
            std::set<std::string> themarker = _markers;
            auto therange =  (themarker.empty())?defaultRange: markedelements(mesh,themarker);
            std::string nameUsed = prefixvm( "expr", _name );
            mapExportsExpr[nameUsed].push_back( std::make_tuple( _expr, therange, _rep ) );
            fieldsNamesToExport.insert( nameUsed );
        }
    }
#endif

    bool hasFieldToExport = this->updatePostProcessExports( exporter, fieldsNamesToExport, time, tupleFields, symbolsExpr,
                                                            hana::concat( hana::make_tuple(mapExportsExpr),tupleExprOnRange) );

    std::string const& pidName = this->postProcessExportsPidName( tag );
    if ( !pidName.empty() && fieldsNamesToExport.find( pidName ) != fieldsNamesToExport.end() )
    {
        exporter->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
        hasFieldToExport = true;
    }

    if ( hasFieldToExport )
    {
        exporter->save();
        this->upload( exporter->path() );
    }
}

template <typename ExporterType,typename ModelFieldsType,typename SymbolsExpr,typename TupleExprOnRangeType>
bool
ModelNumerical::updatePostProcessExports( std::shared_ptr<ExporterType> exporter, std::set<std::string> const& fieldsNamesToExport, double time, ModelFieldsType const& tupleFields, SymbolsExpr const& symbolsExpr, TupleExprOnRangeType const& tupleExprOnRange )
{
    if ( !exporter ) return false;
    if ( !exporter->doExport() ) return false;

    bool hasFieldToExport = false;
    hana::for_each( tupleFields.tuple(), [this,&exporter,&fieldsNamesToExport,&time,&hasFieldToExport]( auto const& e )
                    {
                        if constexpr ( is_iterable_v<decltype(e)> )
                            {
                                for ( auto const& mfield : e )
                                {
                                    std::string fieldName = mfield.nameWithPrefix();
                                    if ( fieldsNamesToExport.find( fieldName ) == fieldsNamesToExport.end() )
                                        continue;

                                    mfield.applyUpdateFunction();
                                    auto const& thefield = unwrap_ptr( mfield.field() );
                                    if constexpr ( decay_type<ExporterType>::mesh_type::nDim == decay_type<decltype(thefield)>::mesh_type::nDim )
                                                 {
                                                         exporter->step( time )->add( prefixvm(this->prefix(),fieldName),
                                                                                      prefixvm(this->prefix(),prefixvm(this->subPrefix(),fieldName)),
                                                                                      thefield );
                                                         hasFieldToExport = true;
                                                 }
                                }
                            }
                        else
                        {
#if 0
                            std::string const& fieldName = e.first;
                            auto const& fieldPtr = e.second;
                            if (!fieldPtr )
                                return;
                            if constexpr ( decay_type<ExporterType>::mesh_type::nDim == decay_type<decltype(fieldPtr)>::mesh_type::nDim )
                                         {
                                             if ( fieldPtr && fieldsNamesToExport.find( fieldName ) != fieldsNamesToExport.end() )
                                             {
                                                 exporter->step( time )->add( prefixvm(this->prefix(),fieldName),
                                                                              prefixvm(this->prefix(),prefixvm(this->subPrefix(),fieldName)),
                                                                              *fieldPtr );
                                                 hasFieldToExport = true;
                                             }
                                         }
#endif
                        }
                    });

    hana::for_each( tupleExprOnRange, [this,&exporter,&fieldsNamesToExport,&time,&hasFieldToExport,&symbolsExpr]( auto const& e )
                    {
                        for ( auto const& [fieldName,exprDatas] : e )
                        {
                            if ( fieldsNamesToExport.find( fieldName ) == fieldsNamesToExport.end() )
                                continue;
                            for ( auto const&[theexpr,range,reprs] : exprDatas )
                            {
                                if constexpr ( std::is_base_of_v<ExprBase,decay_type<decltype(theexpr)>> )
                                {
                                    using _expr_shape = typename std::decay_t<decltype(theexpr)>::template evaluator_t<typename  decay_type<ExporterType>::mesh_type::element_type>::shape;
                                    if constexpr ( _expr_shape::is_tensor2 && _expr_shape::M > 1 && _expr_shape::N > 1 ) // tensor2 asym is not supported with ParaView -> export each components in wating
                                        {
                                            for ( int i=0;i<_expr_shape::M;++i )
                                                for ( int j=0;j<_expr_shape::N;++j )
                                                {
                                                    std::string fieldName2 = (boost::format("%1%_%2%%3%") %fieldName %i %j).str();
                                                    exporter->step( time )->add( prefixvm(this->prefix(),fieldName2),
                                                                                 prefixvm(this->prefix(),prefixvm(this->subPrefix(),fieldName2)),
                                                                                 theexpr(i,j),range,reprs );
                                                }
                                        }
                                    else
                                    {
                                        exporter->step( time )->add( prefixvm(this->prefix(),fieldName),
                                                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),fieldName)),
                                                                     theexpr,range,reprs );
                                    }
                                    hasFieldToExport = true;
                                }
#if 1
                                else if constexpr ( std::is_same_v<decay_type<decltype(theexpr)>,ModelExpression> )
                                {
                                    auto const& fieldNameBIS = fieldName;
                                    auto const& theexprBIS = theexpr;
                                    auto const& rangeBIS = range;
                                    auto const& reprsBIS = reprs;
                                    hana::for_each( ModelExpression::expr_shapes, [this,&exporter,&time,&hasFieldToExport,&symbolsExpr,&fieldNameBIS,&theexprBIS,&rangeBIS,&reprsBIS]( auto const& e_ij )
                                                    {
                                                        constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                                                        constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                                                        if ( theexprBIS.template hasExpr<ni,nj>() )
                                                        {
                                                            if constexpr ( ni == nj && ni > 1 )  // tensor2 asym is not supported with ParaView -> export each components in wating 
                                                            {
                                                                for ( int i=0;i<ni;++i )
                                                                    for ( int j=0;j<nj;++j )
                                                                    {
                                                                        std::string fieldNameBIS2 = (boost::format("%1%_%2%%3%") %fieldNameBIS %i %j).str();
                                                                        exporter->step( time )->add( prefixvm(this->prefix(),fieldNameBIS2),
                                                                                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),fieldNameBIS2)),
                                                                                                     expr(theexprBIS.template expr<ni,nj>(),symbolsExpr)(i,j),rangeBIS,reprsBIS );
                                                                    }
                                                            }
                                                            else
                                                            {
                                                                exporter->step( time )->add( prefixvm(this->prefix(),fieldNameBIS),
                                                                                             prefixvm(this->prefix(),prefixvm(this->subPrefix(),fieldNameBIS)),
                                                                                             expr(theexprBIS.template expr<ni,nj>(),symbolsExpr),rangeBIS,reprsBIS );
                                                            }
                                                            hasFieldToExport = true;
                                                        }
                                                    });
                                } // is ModelExpression
#endif
                            }
                        }
                    });

    return hasFieldToExport;
}

template <typename MeshType, typename RangeType, typename MeasurePointEvalType, typename SymbolsExpr, typename ModelFieldsType, typename TupleQuantitiesType>
void 
ModelNumerical::executePostProcessMeasures( double time, std::shared_ptr<MeshType> mesh, RangeType const& range, std::shared_ptr<MeasurePointEvalType> measurePointsEvaluation, SymbolsExpr const& symbolsExpr, ModelFieldsType const& tupleFields, TupleQuantitiesType const& tupleQuantities )
{
    bool hasMeasure = this->updatePostProcessMeasures( mesh, range, measurePointsEvaluation, symbolsExpr, tupleFields, tupleQuantities );

    if ( hasMeasure )
    {
        if ( !this->isStationary() )
            this->postProcessMeasuresIO().setMeasure( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
        this->upload( this->postProcessMeasuresIO().pathFile() );
    }
}

template<typename TupleQuantitiesType>
bool 
ModelNumerical::updatePostProcessMeasuresQuantities( TupleQuantitiesType const& tupleQuantities )
{
    bool hasMeasure = false;
    std::set<std::string> const& quantitiesToMeasure = this->modelProperties().postProcess().measuresQuantities( this->keyword() ).quantities();
    Feel::for_each( tupleQuantities, [this,&hasMeasure,&quantitiesToMeasure]( auto const& mquantity )
            {
                if constexpr( is_iterable_v<decltype(mquantity)> )
                {
                    for( auto const& quantity : mquantity )
                    {
                        std::string quantityName = quantity.nameWithPrefix();
                        if( quantitiesToMeasure.find( quantityName ) != quantitiesToMeasure.end() )
                        {
                            if constexpr( is_iterable_v<decltype(quantity.value())> )
                            {
                                this->postProcessMeasuresIO().setMeasureComp( quantityName, quantity.value() );
                            }
                            else
                            {
                                this->postProcessMeasuresIO().setMeasure( quantityName, quantity.value() );
                            }
                            hasMeasure = true;
                        }
                    }
                }
                else
                {
                    std::string quantityName = mquantity.nameWithPrefix();
                    if( quantitiesToMeasure.find( quantityName ) != quantitiesToMeasure.end() )
                    {
                        if constexpr( is_iterable_v<decltype(mquantity.value())> )
                        {
                            this->postProcessMeasuresIO().setMeasureComp( quantityName, mquantity.value() );
                        }
                        else
                        {
                            this->postProcessMeasuresIO().setMeasure( quantityName, mquantity.value() );
                        }
                        hasMeasure = true;
                    }
                }
            });
    return hasMeasure;
}

template <typename MeshType, typename RangeType, typename SymbolsExpr, typename ModelFieldsType>
bool
ModelNumerical::updatePostProcessMeasuresNorm( std::shared_ptr<MeshType> mesh, RangeType const& rangeMeshElements, SymbolsExpr const& symbolsExpr, ModelFieldsType const& tupleFields )
{
    bool hasMeasure = false;
    for ( auto const& ppNorm : this->modelProperties().postProcess().measuresNorm( this->keyword() ) )
    {
        std::map<std::string,double> resPpNorms;
        measureNormEvaluation( mesh, rangeMeshElements, ppNorm, resPpNorms, symbolsExpr, tupleFields );
        for ( auto const& resPpNorm : resPpNorms )
        {
            this->postProcessMeasuresIO().setMeasure( resPpNorm.first, resPpNorm.second );
            hasMeasure = true;
        }
    }
    return hasMeasure;
}

template <typename MeshType, typename RangeType, typename SymbolsExpr, typename ModelFieldsType>
bool
ModelNumerical::updatePostProcessMeasuresStatistics( std::shared_ptr<MeshType> mesh, RangeType const& rangeMeshElements, SymbolsExpr const& symbolsExpr, ModelFieldsType const& tupleFields )
{
    bool hasMeasure = false;
    for ( auto const& ppStat : this->modelProperties().postProcess().measuresStatistics( this->keyword() ) )
    {
        std::map<std::string,double> resPpStats;
        measureStatisticsEvaluation( mesh, rangeMeshElements, ppStat, resPpStats, symbolsExpr, tupleFields );
        for ( auto const& resPpStat : resPpStats )
        {
            this->postProcessMeasuresIO().setMeasure( resPpStat.first, resPpStat.second );
            hasMeasure = true;
        }
    }
    return hasMeasure;
}

template <typename MeasurePointEvalType, typename ModelFieldsType>
bool
ModelNumerical::updatePostProcessMeasuresPoint( std::shared_ptr<MeasurePointEvalType> measurePointsEvaluation, ModelFieldsType const& mfields )
{
    if ( !measurePointsEvaluation )
        return false;
    bool hasMeasure = false;
    std::map<std::string,double> resPpPoints;
    measurePointsEvaluation->eval( this->modelProperties().postProcess().measuresPoint( this->keyword() ), resPpPoints, mfields );
    for ( auto const& resPpPoint : resPpPoints )
    {
        this->postProcessMeasuresIO().setMeasure( resPpPoint.first, resPpPoint.second );
        hasMeasure = true;
    }
    return hasMeasure;
}

template <typename MeshType, typename RangeType, typename MeasurePointEvalType, typename SymbolsExpr, typename ModelFieldsType, typename TupleQuantitiesType>
bool 
ModelNumerical::updatePostProcessMeasures( std::shared_ptr<MeshType> mesh, RangeType const& rangeMeshElements, std::shared_ptr<MeasurePointEvalType> measurePointsEvaluation,
                                           SymbolsExpr const& symbolsExpr, ModelFieldsType const& mfields, TupleQuantitiesType const& tupleQuantities )
{
    bool hasMeasureQuantities = this->updatePostProcessMeasuresQuantities( tupleQuantities );
    bool hasMeasureNorm = this->updatePostProcessMeasuresNorm( mesh, rangeMeshElements, symbolsExpr, mfields );
    bool hasMeasureStatistics = this->updatePostProcessMeasuresStatistics( mesh, rangeMeshElements, symbolsExpr, mfields );
    bool hasMeasurePoint = false;
    if( measurePointsEvaluation )
        hasMeasurePoint = this->updatePostProcessMeasuresPoint( measurePointsEvaluation, mfields );

    return ( hasMeasureQuantities || hasMeasureNorm || hasMeasureStatistics || hasMeasurePoint );
}

template <typename ModelFieldsType>
void
ModelNumerical::executePostProcessSave( std::set<std::string> const& fieldsNamesToSave, std::string const& format, uint32_type index, ModelFieldsType const& fieldTuple )
{
    std::string formatUsed = (format.empty())? "default" : format;
    hana::for_each( fieldTuple.tuple(), [this,&fieldsNamesToSave,&formatUsed,&index]( auto const& e )
                    {
                        if constexpr ( is_iterable_v<decltype(e)> )
                            {
                                for ( auto const& mfield : e )
                                {
                                    std::string fieldName = mfield.nameWithPrefix();
                                    if ( fieldsNamesToSave.find( fieldName ) == fieldsNamesToSave.end() )
                                        continue;

                                    mfield.applyUpdateFunction();

                                    std::string fieldNameSaved = fieldName;
                                    if ( index != invalid_uint32_type_value )
                                        fieldNameSaved = (boost::format("%1%_%2%")%fieldNameSaved%index).str();
                                    auto const& thefield = unwrap_ptr( mfield.field() );
                                    thefield.save(_path=this->postProcessSaveRepository().string(),_name=fieldNameSaved,_type=formatUsed );
                                }
                            }
                        else
                        {
#if 0
                            std::string const& fieldName = e.first;
                            auto const& fieldPtr = e.second;
                            if ( fieldPtr && fieldsNamesToSave.find( fieldName ) != fieldsNamesToSave.end() )
                            {
                                std::string fieldNameSaved = fieldName;
                                if ( index != invalid_uint32_type_value )
                                    fieldNameSaved = (boost::format("%1%_%2%")%fieldNameSaved%index).str();
                                fieldPtr->save(_path=this->postProcessSaveRepository().string(),_name=fieldNameSaved,_type=formatUsed );
                            }
#endif
                        }
                    });
}



} // namespace FeelModels
} // namespace feel

#endif // FEELPP_MODELNUMERICAL_HPP
