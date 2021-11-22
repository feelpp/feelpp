/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2015-12-02

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
 \file modelpostprocessextra.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2015-12-02
 */

#ifndef FEELPP_MODELPOSTPROCESS_EXTRA_H
#define FEELPP_MODELPOSTPROCESS_EXTRA_H 1

#include <feel/feelmodels/modelpostprocess.hpp>

namespace Feel
{
namespace FeelModels
{

class ModelMeasuresIO
{
public :
    ModelMeasuresIO( std::string const& pathFile, worldcomm_ptr_t const& worldComm /*= Environment::worldComm()*/ );
    ModelMeasuresIO( ModelMeasuresIO const& ) = default;
    void clear();
    FEELPP_DEPRECATED void start() {}
    void restart( std::string const& paramKey, double val );
    void exportMeasures();
    FEELPP_DEPRECATED void setParameter(std::string const& key,double val);
    void setMeasure(std::string const& key,double val);

    template <typename T>
    void setMeasure( std::string const& key, std::vector<T> const& values ) { this->setMeasure( key, values.data(), values.size() ); }

    template<typename Derived>
    void setMeasure(std::string const& key, Eigen::MatrixBase<Derived> const& mat )
        {
            if ( mat.rows() == 1 && mat.cols() == 1 )
                this->setMeasure( key, mat(0,0) );
            else if ( mat.rows() == 1 || mat.cols() == 1 )
                for (int d=0;d<mat.rows()*mat.cols();++d)
                    this->setMeasure( (boost::format("%1%_%2%")%key %d).str(), mat(d) );
            else
                for (int i=0;i<mat.rows();++i)
                    for (int j=0;j<mat.cols();++j)
                        this->setMeasure( (boost::format("%1%_%2%%3%")%key %i %j).str(), mat(i,j) );
        }
    template <typename T>
    void setMeasure(std::string const& key, const T * data, int dim)
        {
            if ( dim == 1 )
                this->setMeasure( key, data[0] );
            else
                for (int d=0;d<dim;++d)
                    this->setMeasure( (boost::format("%1%_%2%")%key %d).str(), data[d] );
        }
    FEELPP_DEPRECATED void setMeasureComp( std::string const& key,std::vector<double> const& values );
    void setMeasures( std::map<std::string,double> const& m );
    FEELPP_DEPRECATED bool hasParameter( std::string const& key ) const { return this->hasMeasure( key ); }
    bool hasMeasure( std::string const& key ) const { return M_dataNameToIndex.find( key ) != M_dataNameToIndex.end(); }
    //! return measure from a key
    double measure( std::string const& key ) const;
    //! return the current measures in memory
    std::map<std::string,double> currentMeasures() const;
    //! return path of file where measures are stored
    std::string const& pathFile() const { return M_pathFile; }
    //! set path of file where measures are stored
    void setPathFile( std::string const& s ) { M_pathFile = s; }

    //! update measures values into the mapping of values \mp
    void updateParameterValues( std::map<std::string,double> & mp, std::string const& prefix_symbol ) const;
private :
    void writeHeader();
private :
    std::shared_ptr<WorldComm> M_worldComm;
    std::string M_pathFile;
    std::map<std::string,uint16_type> M_dataNameToIndex;
    std::vector<std::string> M_dataIndexToName;
    std::vector<double> M_data;
    bool M_addNewDataIsLocked;
};

class ModelMeasuresEvaluatorContext
{
public :
    ModelMeasuresEvaluatorContext() = default;
    ModelMeasuresEvaluatorContext( ModelMeasuresEvaluatorContext const& ) = default;

    std::map<std::string, std::map<int,std::string> > const& mapFieldToMapCtxIdToName() const { return M_mapFieldToMapCtxIdToName; }
    void add( std::string const& field, int ctxId, std::string const& name );
    bool has( std::string const& field ) const;
    bool has( std::string const& field, int ctxId ) const;
    std::string const& name( std::string const& field, int ctxId ) const;
    int ctxId( std::string const& field, std::string const& name ) const;
private :
    // for each field, store data names evaluted : field -> ( (ctxId1->dataName1), (ctxId2->dataName2),...)
    std::map<std::string, std::map<int,std::string> > M_mapFieldToMapCtxIdToName;
    std::string M_emptyString;
};

class ModelMeasuresForces
{
public :
    ModelMeasuresForces() = default;
    ModelMeasuresForces( ModelMeasuresForces const& ) = default;
    ModelMeasuresForces( ModelMeasuresForces&& ) = default;
    ModelMeasuresForces& operator=( ModelMeasuresForces const& ) = default;
    ModelMeasuresForces& operator=( ModelMeasuresForces && ) = default;

    std::string const& name() const { return M_name; }
    std::list<std::string> const& meshMarkers() const { return M_meshMarkers; }

    void setName( std::string const& s ) { M_name = s; }
    void addMarker( std::string const& mark )
        {
            if ( std::find( M_meshMarkers.begin(),M_meshMarkers.end(), mark ) == M_meshMarkers.end() )
                M_meshMarkers.push_back( mark );
        }

private :
    std::string M_name;
    std::list<std::string> M_meshMarkers;
    std::string M_direction;
};


class ModelMeasuresFlowRate
{
public :
    ModelMeasuresFlowRate();
    ModelMeasuresFlowRate( ModelMeasuresFlowRate const& ) = default;
    ModelMeasuresFlowRate( ModelMeasuresFlowRate&& ) = default;
    ModelMeasuresFlowRate& operator=( ModelMeasuresFlowRate const& ) = default;
    ModelMeasuresFlowRate& operator=( ModelMeasuresFlowRate && ) = default;

    std::string const& name() const { return M_name; }
    std::list<std::string> const& meshMarkers() const { return M_meshMarkers; }
    std::string const& direction() const { return M_direction; }
    bool useExteriorNormal() const { return M_direction == "exterior_normal"; }

    void setName( std::string const& s ) { M_name = s; }
    void addMarker( std::string const& mark );
    void setDirection( std::string const& dir );

    void setup( pt::ptree const& _pt, std::string const& name );

private :
    std::string M_name;
    std::list<std::string> M_meshMarkers;
    std::string M_direction;
};


class ModelMeasuresNormalFluxGeneric
{
public :
    ModelMeasuresNormalFluxGeneric() = default;
    ModelMeasuresNormalFluxGeneric( ModelMeasuresNormalFluxGeneric const& ) = default;
    ModelMeasuresNormalFluxGeneric( ModelMeasuresNormalFluxGeneric&& ) = default;
    ModelMeasuresNormalFluxGeneric& operator=( ModelMeasuresNormalFluxGeneric const& ) = default;
    ModelMeasuresNormalFluxGeneric& operator=( ModelMeasuresNormalFluxGeneric && ) = default;

    std::string const& name() const { return M_name; }
    std::set<std::string> const& markers() const { return M_markers; }
    std::string const& direction() const { return M_direction; }
    bool isOutward() const { return M_direction == "outward"; }

    void setup( pt::ptree const& _pt, std::string const& name, ModelIndexes const& indexes );

private :
    std::string M_name;
    ModelMarkers M_markers;
    std::string M_direction;
};



} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELPOSTPROCESS_EXTRA_H
