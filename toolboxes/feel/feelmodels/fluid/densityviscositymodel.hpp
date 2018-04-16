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

template<class SpaceType>
class DensityViscosityModel : public DynamicViscosityModel<SpaceType>
{
    typedef DensityViscosityModel<SpaceType> self_type;
public :
    typedef SpaceType space_type;
    typedef boost::shared_ptr<SpaceType> space_ptrtype;
    typedef typename SpaceType::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef DynamicViscosityModel<space_type> super_type;

    static std::string defaultMaterialName() { return std::string("FEELPP_DEFAULT_MATERIAL_NAME"); }

    DensityViscosityModel( std::string prefix )
        :
        super_type( prefix )
        {
            M_cstDensity[self_type::defaultMaterialName()] = doption(_name="rho",_prefix=prefix);
            M_cstCinematicViscosity[self_type::defaultMaterialName()] = this->cstMu()/this->cstDensity();
        }
    DensityViscosityModel( DensityViscosityModel const& app  ) = default;

    void updateForUse( mesh_ptrtype const& mesh , ModelMaterials const& mats, bool useExtendedDofTable )
    {
        super_type::updateForUse( mesh,mats,useExtendedDofTable );

        M_fieldDensity = this->dynamicViscositySpace()->elementPtr( cst( this->cstDensity( self_type::defaultMaterialName() ) ) );
        M_fieldCinematicViscosity = this->dynamicViscositySpace()->elementPtr( cst( this->cstCinematicViscosity( self_type::defaultMaterialName() ) ) );

        for( auto const& m : mats )
        {
            auto const& mat = m.second;
            for ( std::string const& matmarker : mat.meshMarkers() )
            {
                if ( this->markers().find( matmarker ) == this->markers().end() )
                    continue;

                if ( mat.hasPropertyExprScalar("rho") )
                    this->setDensity( mat.propertyExprScalar("rho"),matmarker );
                else
                    this->setCstDensity( mat.propertyConstant("rho"),matmarker );
            }
        }
    }
    void updateForUse( space_ptrtype const& space, ModelMaterials const& mats )
    {
        super_type::updateForUse( space, mats );
        M_fieldDensity = this->dynamicViscositySpace()->elementPtr( cst( this->cstDensity( self_type::defaultMaterialName() ) ) );
        M_fieldCinematicViscosity = this->dynamicViscositySpace()->elementPtr( cst( this->cstCinematicViscosity( self_type::defaultMaterialName() ) ) );

        for( auto const& m : mats )
        {
            auto const& mat = m.second;
            for ( std::string const& matmarker : mat.meshMarkers() )
            {
                if ( this->markers().find( matmarker ) == this->markers().end() )
                    continue;

                if ( mat.hasPropertyExprScalar("rho") )
                    this->setDensity( mat.propertyExprScalar("rho"),matmarker );
                else
                    this->setCstDensity( mat.propertyConstant("rho"),matmarker );
            }
        }
    }

    double cstRho( std::string const& marker = "" ) const { return this->cstDensity(marker); }
    double cstDensity( std::string const& marker = "" ) const
    {
        std::string markerUsed = ( marker.empty() )?
            ( ( this->markers().empty() )? self_type::defaultMaterialName() : *this->markers().begin() ) :
            marker;
        // std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
        auto itFindMarker = M_cstDensity.find( markerUsed );
        CHECK( itFindMarker != M_cstDensity.end() ) << "invalid marker not registered " << markerUsed;
        return itFindMarker->second;
    }
    void setCstDensity(double d, std::string const& marker = "", bool update = true )
    {
        std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
        M_cstDensity[markerUsed]=d;
        if ( update )
        {
            this->updateDensity( cst(d),marker);
            this->updateCinematicViscosity( marker );
        }
    }
    element_type const& fieldRho() const { return this->fieldDensity(); }
    element_type const& fieldDensity() const { return *M_fieldDensity; }
    element_ptrtype const& fieldDensityPtr() const { return M_fieldDensity; }

    double cstNu( std::string const& marker = "" ) const { return this->cstCinematicViscosity(marker); }
    double cstCinematicViscosity( std::string const& marker = "" ) const
    {
        std::string markerUsed = ( marker.empty() )?
            ( ( this->markers().empty() )? self_type::defaultMaterialName() : *this->markers().begin() ) :
            marker;
        // std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
        auto itFindMarker = M_cstCinematicViscosity.find( markerUsed );
        CHECK( itFindMarker != M_cstCinematicViscosity.end() ) << "invalid marker not registered " << markerUsed;
        return itFindMarker->second;
    }
    void setCstCinematicViscosity( double d, std::string const& marker = "" )
    {
        std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
        M_cstCinematicViscosity[markerUsed]=d;
        this->updateCinematicViscosity(cst(d),marker);
    }
    element_type const& fieldNu() const { return this->fieldCinematicViscosity(); }
    element_type const& fieldCinematicViscosity() const { return *M_fieldCinematicViscosity; }
    element_ptrtype const& fieldCinematicViscosityPtr() const { return M_fieldCinematicViscosity; }


    void setCstDynamicViscosity( double d, std::string const& marker = "" )
    {
        super_type::setCstDynamicViscosity(d,marker);
        this->updateCinematicViscosity(marker);
    }

    template < typename ExprT >
    void setDensity( vf::Expr<ExprT> const& vfexpr, std::string const& marker = "" )
    {
        this->updateDensity( vfexpr,marker);
        if ( M_fieldDensity )
            this->setCstDensity( M_fieldDensity->min(), marker, false );
    }
    template < typename ExprT >
    void updateDensity( vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
    {
        if ( !M_fieldDensity ) return;
        auto rangeEltUsed = ( marker.empty() )? elements(M_fieldDensity->mesh()) : markedelements(M_fieldDensity->mesh(),marker);
        M_fieldDensity->on(_range=rangeEltUsed,_expr=__expr );
    }
    template < typename ExprT >
    void updateCinematicViscosity( vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
    {
        if ( !M_fieldCinematicViscosity ) return;
        auto rangeEltUsed = ( marker.empty() )? elements(M_fieldCinematicViscosity->mesh()) : markedelements(M_fieldCinematicViscosity->mesh(),marker);
        M_fieldCinematicViscosity->on(_range=rangeEltUsed,_expr=__expr );
    }
    template < typename ExprT >
    void updateDynamicViscosity( vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
    {
        super_type::updateDynamicViscosity(__expr,marker);
        this->updateCinematicViscosity(marker);
    }

    void updateCinematicViscosity( std::string const& marker = "" )
    {
        std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
        M_cstCinematicViscosity[markerUsed] = this->cstMu(markerUsed)/this->cstRho(markerUsed);
        this->updateCinematicViscosity( idv(this->fieldMu())/idv(this->fieldRho()), marker );
    }

    boost::shared_ptr<std::ostringstream>
    getInfoMaterialParameters() const
        {
            boost::shared_ptr<std::ostringstream> ostr( new std::ostringstream() );
            *ostr << "\n   Materials parameters";
            std::set<std::string> matmarkerset = this->markers();
            if (  matmarkerset.empty() ) matmarkerset.insert(std::string(""));
            *ostr << "\n     -- number of materials : " << matmarkerset.size();
            for ( std::string const& matmarker : matmarkerset)
            {
                std::string matmarkertag = matmarker.empty()? std::string("") : (boost::format("[%1%] ")%matmarker).str();
                *ostr << "\n     -- " << matmarkertag << "rho : " << this->cstRho(matmarker);
                *ostr << this->getInfo( matmarker )->str();
            }
            return ostr;
        }


private :

    std::map<std::string,double> M_cstDensity;// rho
    element_ptrtype M_fieldDensity;// rho
    std::map<std::string,double> M_cstCinematicViscosity;// nu
    element_ptrtype M_fieldCinematicViscosity;// nu

};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_FLUIDMECHANICS_DENSITYVISCOSITYMODEL_H
