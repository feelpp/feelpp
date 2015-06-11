/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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

#include <feel/feelmodels2/modelcore/modelnumerical.hpp>

namespace Feel
{
namespace FeelModels
{

ModelNumerical::ModelNumerical(/*bool _isStationary,*/ std::string _theprefix, WorldComm const& _worldComm, std::string subPrefix,
                                                                                   std::string appliShortRepository )
        :
        super_type( _theprefix, _worldComm, subPrefix, appliShortRepository ),
        M_rebuildMeshPartitions( boption(_name="rebuild_mesh_partitions",_prefix=this->prefix()) ),
        M_isStationary( false /*_isStationary*/),
        M_doRestart( boption(_name="ts.restart") ),
        M_restartPath( soption(_name="ts.restart.path") ),
        M_restartAtLastSave( boption(_name="ts.restart.at-last-save") ),
        M_timeInitial( doption(_name="ts.time-initial") ),
        M_timeFinal( doption(_name="ts.time-final") ),
        M_timeStep( doption(_name="ts.time-step") ),
        M_bdfSaveInFile( boption(_name="ts.save") ),
        M_bdfSaveFreq( ioption(_name="ts.save.freq") ),
        M_timeCurrent(M_timeInitial),
        M_modelProps( Environment::expand( soption( _name=prefixvm(this->prefix(),"filename")) ) ),
        M_parameters(std::vector<double>(FEELMODELS_FSIBASE_NUMBER_OF_PARAMETERS,0)),
        M_geoParameters(std::vector<double>(FEELMODELS_FSIBASE_NUMBER_OF_GEOPARAMETERS,0)),
        M_ginacExpr(std::vector<std::pair<std::string,std::string> >(FEELMODELS_FSIBASE_NUMBER_OF_GINACEXPR)),
        M_ginacExprCompilationDirectory( "" ),
        M_geotoolMeshIndex( ioption(_name="geotool-mesh-index",_prefix=this->prefix()) ),
        M_geotoolSaveDirectory( soption(_name="geotool-save-directory",_prefix=this->prefix()) ),
        M_geotoolSaveName( soption(_name="geotool-save-name",_prefix=this->prefix()) ),
        M_row_startInMatrix(0),
        M_col_startInMatrix(0),
        M_row_startInVector(0),
        M_mshFileStr("FEELMODELS_WARNING_NODEFINE"),
        M_geoFileStr("FEELMODELS_WARNING_NODEFINE"),
        M_exporterPath( this->appliRepository()+"/"+prefixvm(this->prefix(), prefixvm(this->subPrefix(),"exports")) )
        //M_PsLogger( new PsLogger(prefixvm(this->prefix(),"PsLogger"),this->worldComm() ) )
    {

        // move in stationary mode if we have this relation
        if ( M_timeInitial + M_timeStep == M_timeFinal)
            M_isStationary=true;

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
                std::string gexpr= soption(_prefix=this->prefix(),_name=(boost::format("ginac-expr%1%") %k ).str());
                std::string gname= soption(_prefix=this->prefix(),_name=(boost::format("ginac-name%1%") %k ).str());
                this->setUserGinacExpr(k,gexpr,gname);
            }
        if ( Environment::vm().count(prefixvm(this->prefix(),"ginac-expr-directory").c_str()) )
        {
            M_ginacExprCompilationDirectory=Environment::rootRepository()+"/"+soption(_name="ginac-expr-directory",_prefix=this->prefix());
        }
        else
        {
            M_ginacExprCompilationDirectory = (fs::path(this->appliRepositoryWithoutNumProc() )/fs::path("symbolic_expr")).string();
        }
        //-----------------------------------------------------------------------//
        // mesh file : .msh
        if (Environment::vm().count(prefixvm(this->prefix(),"mshfile").c_str()))
            M_mshFileStr = Environment::expand( soption(_prefix=this->prefix(),_name="mshfile") );
        // mesh file : .geo
        if (Environment::vm().count(prefixvm(this->prefix(),"geofile").c_str()))
            M_geoFileStr = Environment::expand( soption(_prefix=this->prefix(),_name="geofile") );
        // mesh file : geotool with .mesh file
        if (M_geotoolSaveDirectory.empty()) M_geotoolSaveDirectory = this->appliShortRepository();//this->appliRepository();
        if (M_geotoolSaveName.empty()) M_geotoolSaveName = this->prefix();
        //-----------------------------------------------------------------------//
        if (soption(_prefix=this->prefix(),_name="geomap")=="opt")
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
        if (this->verbose()) FeelModels::Log(this->prefix()+".ModelNumerical","saveMSHfilePath",
                                             "namePath "+ namePath + "\nwrite :\n" + M_mshFileStr,
                                             this->worldComm(),this->verboseAllProc());

        if ( this->worldComm().isMasterRank() )
        {
            fs::path thedir = fs::path( namePath ).parent_path();
            if ( !fs::exists(thedir))
                fs::create_directories(thedir);

            //std::string nameFile = prefixvm(this->prefix(),namePath);
            std::ofstream file(namePath.c_str(), std::ios::out);
            //M_mshFileStr = this->application()->vm()[prefixvm(this->prefix(),"mshfile")].as< std::string >() ;
            file << M_mshFileStr;
            file.close();
        }
    }


} // namespace FeelModels

} // namespace Feel

