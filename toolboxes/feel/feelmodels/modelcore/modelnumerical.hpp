/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2012-01-19

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
                        ModelBaseRepository const& modelRep = ModelBaseRepository() );
        ModelNumerical( std::string const& _theprefix, worldcomm_ptr_t const& _worldComm=Environment::worldCommPtr(),
                        std::string const& subPrefix="",
                        ModelBaseRepository const& modelRep = ModelBaseRepository() )
            :
            ModelNumerical( _theprefix, _theprefix, _worldComm, subPrefix, modelRep )
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

        bool hasModelProperties() const { return (M_modelProps)? true : false; }
        std::shared_ptr<ModelProperties> modelPropertiesPtr() const { return M_modelProps; }
        ModelProperties const& modelProperties() const { return *M_modelProps; }
        ModelProperties & modelProperties() { return *M_modelProps; }
        void setModelProperties( std::shared_ptr<ModelProperties> modelProps ) { M_modelProps = modelProps; }
        void addParameterInModelProperties( std::string const& symbolName,double value );


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

        ModelMeasuresIO const& postProcessMeasuresIO() const { return M_postProcessMeasuresIO; }
        ModelMeasuresIO & postProcessMeasuresIO() { return M_postProcessMeasuresIO; }
        ModelMeasuresEvaluatorContext const& postProcessMeasuresEvaluatorContext() const { return M_postProcessMeasuresEvaluatorContext; }
        ModelMeasuresEvaluatorContext & postProcessMeasuresEvaluatorContext() { return M_postProcessMeasuresEvaluatorContext; }

        virtual bool checkResults() const;

    protected :
        template<typename SymbExprField>
        auto symbolsExprFit( SymbExprField const& sef ) const
            {
                typedef Expr< Fit<decltype(expr(scalar_field_expression<2>{},sef)),0> > fit_expr_type;
                std::vector<std::pair<std::string,fit_expr_type>> fitSymbs;
                if ( this->hasModelProperties() )
                {
                    for ( auto const& param : this->modelProperties().parameters() )
                    {
                        if ( param.second.type() != "fit" )
                            continue;
                        auto exprInFit = expr( param.second.expression(), sef );
                        fitSymbs.push_back( std::make_pair( param.first, fit( exprInFit, param.second.fitInterpolator() ) ) );
                    }
                }
                return Feel::vf::symbolsExpr( symbolExpr( fitSymbs ) );
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


        std::string M_meshFile, M_geoFile;

        std::string M_exporterPath;
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
                    fileType = "hdf5";
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

} // namespace FeelModels
} // namespace feel

#endif // FEELPP_MODELNUMERICAL_HPP
