/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2014-03-21

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
 \file viscositymodeldescription.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2014-03-21
 */

#ifndef FEELPP_FLUIDMECHANICS_DENSITYVISCOSITYMODEL_H
#define FEELPP_FLUIDMECHANICS_DENSITYVISCOSITYMODEL_H 1

#include <feel/feelmodels/fluid/dynamicviscositymodel.hpp>

namespace Feel
{
namespace FeelModels
{

template <class SpaceType>
class DensityViscosityModel : public DynamicViscosityModel<SpaceType>
{
    typedef DensityViscosityModel<SpaceType> self_type;

  public:
    typedef SpaceType space_type;
    typedef boost::shared_ptr<SpaceType> space_ptrtype;
    typedef typename SpaceType::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef DynamicViscosityModel<space_type> super_type;

    DensityViscosityModel( std::string prefix )
        : super_type( prefix )
    {
        M_cstDensity[self_type::defaultMaterialName()] = doption( _name = "rho", _prefix = prefix );
        M_cstCinematicViscosity[self_type::defaultMaterialName()] = this->cstMu() / this->cstDensity();
    }
    DensityViscosityModel( DensityViscosityModel const& app ) = default;

    void initFromMesh( mesh_ptrtype const& mesh, bool useExtendedDofTable )
    {
        super_type::initFromMesh( mesh, useExtendedDofTable );
        M_fieldDensity = this->dynamicViscositySpace()->elementPtr( cst( this->cstDensity() ) );
        M_fieldCinematicViscosity = this->dynamicViscositySpace()->elementPtr( cst( this->cstCinematicViscosity() ) );
    }

    void initFromSpace( space_ptrtype const& space )
    {
        super_type::initFromSpace( space );
        M_fieldDensity = this->dynamicViscositySpace()->elementPtr( cst( this->cstDensity() ) );
        M_fieldCinematicViscosity = this->dynamicViscositySpace()->elementPtr( cst( this->cstCinematicViscosity() ) );
    }

    double cstRho( std::string const& marker = "" ) const { return this->cstDensity( marker ); }
    double cstDensity( std::string const& marker = "" ) const
    {
        std::string markerUsed = ( marker.empty() ) ? self_type::defaultMaterialName() : marker;
        auto itFindMarker = M_cstDensity.find( markerUsed );
        CHECK( itFindMarker != M_cstDensity.end() ) << "invalid marker not registered " << markerUsed;
        return itFindMarker->second;
    }
    void setCstDensity( double d, std::string const& marker = "" )
    {
        std::string markerUsed = ( marker.empty() ) ? self_type::defaultMaterialName() : marker;
        M_cstDensity[markerUsed] = d;
        this->updateDensity( cst( d ), marker );
        this->updateCinematicViscosity( marker );
    }
    element_type const& fieldRho() const { return this->fieldDensity(); }
    element_type const& fieldDensity() const { return *M_fieldDensity; }
    element_ptrtype const& fieldDensityPtr() const { return M_fieldDensity; }

    double cstNu( std::string const& marker = "" ) const { return this->cstCinematicViscosity( marker ); }
    double cstCinematicViscosity( std::string const& marker = "" ) const
    {
        std::string markerUsed = ( marker.empty() ) ? self_type::defaultMaterialName() : marker;
        auto itFindMarker = M_cstCinematicViscosity.find( markerUsed );
        CHECK( itFindMarker != M_cstCinematicViscosity.end() ) << "invalid marker not registered " << markerUsed;
        return itFindMarker->second;
    }
    void setCstCinematicViscosity( double d, std::string const& marker = "" )
    {
        std::string markerUsed = ( marker.empty() ) ? self_type::defaultMaterialName() : marker;
        M_cstCinematicViscosity[markerUsed] = d;
        this->updateCinematicViscosity( cst( d ), marker );
    }
    element_type const& fieldNu() const { return this->fieldCinematicViscosity(); }
    element_type const& fieldCinematicViscosity() const { return *M_fieldCinematicViscosity; }
    element_ptrtype const& fieldCinematicViscosityPtr() const { return M_fieldCinematicViscosity; }

    void setCstDynamicViscosity( double d, std::string const& marker = "" )
    {
        super_type::setCstDynamicViscosity( d, marker );
        this->updateCinematicViscosity( marker );
    }

    template <typename ExprT>
    void updateDensity( vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
    {
        if ( !M_fieldDensity ) return;
        if ( marker.empty() )
            M_fieldDensity->on( _range = elements( M_fieldDensity->mesh() ), _expr = __expr );
        else
            M_fieldDensity->on( _range = markedelements( M_fieldDensity->mesh(), marker ), _expr = __expr );
    }
    template <typename ExprT>
    void updateCinematicViscosity( vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
    {
        if ( !M_fieldCinematicViscosity ) return;
        if ( marker.empty() )
            M_fieldCinematicViscosity->on( _range = elements( M_fieldCinematicViscosity->mesh() ), _expr = __expr );
        else
            M_fieldCinematicViscosity->on( _range = markedelements( M_fieldCinematicViscosity->mesh(), marker ), _expr = __expr );
    }
    template <typename ExprT>
    void updateDynamicViscosity( vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
    {
        super_type::updateDynamicViscosity( __expr, marker );
        this->updateCinematicViscosity( marker );
    }

    void updateCinematicViscosity( std::string const& marker = "" )
    {
        std::string markerUsed = ( marker.empty() ) ? self_type::defaultMaterialName() : marker;
        M_cstCinematicViscosity[markerUsed] = this->cstMu( markerUsed ) / this->cstRho( markerUsed );
        if ( M_fieldCinematicViscosity )
        {
            if ( marker.empty() )
                M_fieldCinematicViscosity->on( _range = elements( M_fieldCinematicViscosity->mesh() ),
                                               _expr = idv( this->fieldMu() ) / idv( this->fieldRho() ) );
            else
                M_fieldCinematicViscosity->on( _range = markedelements( M_fieldCinematicViscosity->mesh(), marker ),
                                               _expr = idv( this->fieldMu() ) / idv( this->fieldRho() ) );
        }
    }

    void updateFromModelMaterials( ModelMaterials const& mat )
    {
        if ( mat.empty() ) return;
        super_type::updateFromModelMaterials( mat );
        for ( auto const& m : mat )
        {
            auto const& mat = m.second;
            auto const& matmarker = m.first;
            this->setCstDensity( mat.rho(), matmarker );
        }
    }

  private:
    std::map<std::string, double> M_cstDensity;            // rho
    element_ptrtype M_fieldDensity;                        // rho
    std::map<std::string, double> M_cstCinematicViscosity; // nu
    element_ptrtype M_fieldCinematicViscosity;             // nu
};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_FLUIDMECHANICS_DENSITYVISCOSITYMODEL_H
