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


        ModelNumerical( std::string const& _theprefix, WorldComm const& _worldComm=WorldComm(), std::string const& subPrefix="",
                        std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );

        ModelNumerical( ModelNumerical const& app ) = default;

        virtual ~ModelNumerical() {};


        boost::shared_ptr<PsLogger> psLogger()  { return M_PsLogger; }
        boost::shared_ptr<PsLogger> const& psLogger() const { return M_PsLogger; }


        bool isStationary() const { return M_isStationary; }
        void setStationary(bool b);

        bool doRestart() const { return M_doRestart; }
        void setRestart(bool b) { M_doRestart=b; }
        std::string restartPath() const { return M_restartPath; }
        void setRestartPath(std::string s) { M_restartPath=s; }
        bool restartAtLastSave() const { return M_restartAtLastSave; }
        void setRestartAtLastSave( bool b) { M_restartAtLastSave=b; }

        double time() const { return this->currentTime(); }
        double currentTime() const { return M_timeCurrent; }
        void updateTime(double t);

        double timeInitial() const { return M_timeInitial; }
        double timeFinal() const { return M_timeFinal; }
        double timeStep() const { return M_timeStep; }
        bool tsSaveInFile() const { return M_tsSaveInFile; }
        int tsSaveFreq() const { return M_tsSaveFreq; }
        void setTimeInitial(double v)  { M_timeInitial=v; }
        void setTimeFinal(double v)  { M_timeFinal=v; }
        void setTimeStep(double v)  { M_timeStep=v; }

        ModelProperties const& modelProperties() const { return *M_modelProps; }
        ModelProperties & modelProperties() { return *M_modelProps; }
        void addParameterInModelProperties( std::string const& symbolName,double value );

        std::string const& directoryLibSymbExpr() const { return M_directoryLibSymbExpr; }

        size_type rowStartInMatrix() const { return M_row_startInMatrix; }
        void setRowStartInMatrix(size_type r) { M_row_startInMatrix=r; }
        size_type colStartInMatrix() const { return M_col_startInMatrix; }
        void setColStartInMatrix(size_type c) { M_col_startInMatrix=c; }
        size_type rowStartInVector() const { return M_row_startInVector; }
        void setRowStartInVector(size_type r) { M_row_startInVector=r; }

        bool rebuildMeshPartitions() const { return (M_rebuildMeshPartitions && (this->worldComm().localSize()>1)); }
        void setRebuildMeshPartitions(bool b) { M_rebuildMeshPartitions=b; }

        GeomapStrategyType geomap() const { return M_geomap; }

        //----------------------------------------------------------------------------------//


        std::string mshfileStr() const { return M_mshFileStr; }
        void setMshfileStr(std::string file)  { M_mshFileStr=file; }

        std::string geofileStr() const { return M_geoFileStr; }
        void setGeofileStr(std::string file)  { M_geoFileStr=file; }

        bool hasMshfileStr() const { return M_mshFileStr!="FEELMODELS_WARNING_NODEFINE"; }
        bool hasGeofileStr() const { return M_geoFileStr!="FEELMODELS_WARNING_NODEFINE"; }


        void saveMSHfilePath( std::string const& fileSavePath, std::string const& meshPath = "" ) const;

        void setExporterPath(std::string s)  { M_exporterPath=s; }
        std::string exporterPath() const { return M_exporterPath; }

        ModelMeasuresIO const& postProcessMeasuresIO() const { return M_postProcessMeasuresIO; }
        ModelMeasuresIO & postProcessMeasuresIO() { return M_postProcessMeasuresIO; }
        ModelMeasuresEvaluatorContext const& postProcessMeasuresEvaluatorContext() const { return M_postProcessMeasuresEvaluatorContext; }
        ModelMeasuresEvaluatorContext & postProcessMeasuresEvaluatorContext() { return M_postProcessMeasuresEvaluatorContext; }

    private :


        bool M_rebuildMeshPartitions;

        bool M_isStationary;
        bool M_doRestart;
        std::string M_restartPath;
        bool M_restartAtLastSave;

        double M_timeInitial,M_timeFinal,M_timeStep;
        bool M_tsSaveInFile;
        int M_tsSaveFreq;
        double M_timeCurrent;

        std::string M_directoryLibSymbExpr;
        std::shared_ptr<ModelProperties> M_modelProps;

        size_type M_row_startInMatrix, M_col_startInMatrix, M_row_startInVector;

        std::string M_mshFileStr;
        std::string M_geoFileStr;

        std::string M_exporterPath;
        ModelMeasuresIO M_postProcessMeasuresIO;
        ModelMeasuresEvaluatorContext M_postProcessMeasuresEvaluatorContext;

        boost::shared_ptr<PsLogger> M_PsLogger;

        GeomapStrategyType M_geomap;

    };


} // namespace FeelModels
} // namespace feel

#endif // FEELPP_MODELNUMERICAL_HPP
