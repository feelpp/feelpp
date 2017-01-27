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
 \file modelnumerical.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2012-01-19
 */

#include <feel/feelmodels/modelcore/modelnumerical.hpp>

namespace Feel
{
namespace FeelModels
{


ModelNumerical::ModelNumerical( std::string const& _theprefix, WorldComm const& _worldComm, std::string const& subPrefix,
                                std::string const& rootRepository )
        :
        super_type( _theprefix, _worldComm, subPrefix, rootRepository ),
        M_rebuildMeshPartitions( boption(_name="rebuild_mesh_partitions",_prefix=this->prefix()) ),
        M_isStationary( /*false*/ boption(_name="ts.steady") ),
        M_doRestart( boption(_name="ts.restart") ),
        M_restartPath( soption(_name="ts.restart.path") ),
        M_restartAtLastSave( boption(_name="ts.restart.at-last-save") ),
        M_timeInitial( doption(_name="ts.time-initial") ),
        M_timeFinal( doption(_name="ts.time-final") ),
        M_timeStep( doption(_name="ts.time-step") ),
        M_tsSaveInFile( boption(_name="ts.save") ),
        M_tsSaveFreq( ioption(_name="ts.save.freq") ),
        M_timeCurrent(M_timeInitial),
        M_directoryLibSymbExpr( "" ),
        M_row_startInMatrix(0),
        M_col_startInMatrix(0),
        M_row_startInVector(0),
        M_mshFileStr("FEELMODELS_WARNING_NODEFINE"),
        M_geoFileStr("FEELMODELS_WARNING_NODEFINE"),
        M_exporterPath( this->rootRepository()+"/"+prefixvm(this->prefix(), prefixvm(this->subPrefix(),"exports")) ),
        M_postProcessMeasuresIO( this->rootRepository()+"/"+prefixvm(this->prefix(), prefixvm(this->subPrefix(),"measures.csv")),this->worldComm() )
        //M_PsLogger( new PsLogger(prefixvm(this->prefix(),"PsLogger"),this->worldComm() ) )
    {
        //-----------------------------------------------------------------------//
        // move in stationary mode if we have this relation
        if ( M_timeInitial + M_timeStep == M_timeFinal)
            M_isStationary=true;
        //-----------------------------------------------------------------------//
        if ( Environment::vm().count(prefixvm(this->prefix(),"symbolic-expr.directory").c_str()) )
        {
            M_directoryLibSymbExpr = soption(_name="symbolic-expr.directory",_prefix=this->prefix());
            if ( fs::path( M_directoryLibSymbExpr ).is_relative() )
                M_directoryLibSymbExpr = (fs::path(this->rootRepositoryWithoutNumProc())/fs::path(M_directoryLibSymbExpr)).string();
        }
        else
            M_directoryLibSymbExpr = (fs::path(this->rootRepositoryWithoutNumProc() )/fs::path("symbolic_expr")).string();
        //-----------------------------------------------------------------------//
        // mesh file : .msh
        if (Environment::vm().count(prefixvm(this->prefix(),"mshfile").c_str()))
            M_mshFileStr = Environment::expand( soption(_prefix=this->prefix(),_name="mshfile") );
        // mesh file : .geo
        if (Environment::vm().count(prefixvm(this->prefix(),"geofile").c_str()))
            M_geoFileStr = Environment::expand( soption(_prefix=this->prefix(),_name="geofile") );
        //-----------------------------------------------------------------------//
        if (soption(_prefix=this->prefix(),_name="geomap")=="opt")
            M_geomap=GeomapStrategyType::GEOMAP_OPT;
        else
            M_geomap=GeomapStrategyType::GEOMAP_HO;
        //-----------------------------------------------------------------------//
        M_modelProps = std::make_shared<ModelProperties>( Environment::expand( soption( _name=prefixvm(this->prefix(),"filename")) ),
                                                          M_directoryLibSymbExpr, this->worldComm() );
    }

   void
   ModelNumerical::addParameterInModelProperties( std::string const& symbolName,double value)
   {
        M_modelProps->parameters()[symbolName] = ModelParameter(symbolName,value);
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

    void
    ModelNumerical::updateTime(double t)
    {
        M_timeCurrent=t;
        M_modelProps->parameters()["t"] = ModelParameter("current_time",M_timeCurrent);
    }



    void
    ModelNumerical::saveMSHfilePath( std::string const& fileSavePath, std::string const& meshPath ) const
    {
        std::string meshPathUsed = (meshPath.empty())? this->mshfileStr() : meshPath;
        if (this->verbose()) FeelModels::Log(this->prefix()+".ModelNumerical","saveMSHfilePath",
                                             "fileSavePath :"+ fileSavePath + "\nwrite :\n" + meshPathUsed,
                                             this->worldComm(),this->verboseAllProc());

        if ( this->worldComm().isMasterRank() )
        {
            fs::path thedir = fs::path( fileSavePath ).parent_path();
            if ( !fs::exists(thedir))
                fs::create_directories(thedir);

            std::ofstream file(fileSavePath.c_str(), std::ios::out);
            file << meshPathUsed;
            file.close();
        }
    }


} // namespace FeelModels

} // namespace Feel

