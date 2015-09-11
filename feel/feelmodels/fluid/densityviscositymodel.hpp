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
public :
    typedef SpaceType space_type;
    typedef boost::shared_ptr<SpaceType> space_ptrtype;
    typedef typename SpaceType::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef DynamicViscosityModel<space_type> super_type;

    DensityViscosityModel( std::string prefix )
        :
        super_type( prefix ),
        M_cstDensity( doption(_name="rho",_prefix=prefix) ),
        M_cstCinematicViscosity( this->cstMu()/M_cstDensity )
        {}
    DensityViscosityModel( DensityViscosityModel const& app  ) = default;

    void initFromMesh( mesh_ptrtype const& mesh, bool useExtendedDofTable )
    {
        super_type::initFromMesh(mesh,useExtendedDofTable);
        M_fieldDensity = this->dynamicViscositySpace()->elementPtr( cst( this->cstDensity() ) );
        M_fieldCinematicViscosity = this->dynamicViscositySpace()->elementPtr( cst( this->cstCinematicViscosity() ) );
    }

    void initFromSpace( space_ptrtype const& space )
    {
        super_type::initFromSpace(space);
        M_fieldDensity = this->dynamicViscositySpace()->elementPtr( cst( this->cstDensity() ) );
        M_fieldCinematicViscosity = this->dynamicViscositySpace()->elementPtr( cst( this->cstCinematicViscosity() ) );
    }


    double cstRho() const { return this->cstDensity(); }
    double cstDensity() const { return M_cstDensity; }
    void setCstDensity(double d) { M_cstDensity=d;this->updateDensity(cst(d)); this->updateCinematicViscosity(); }
    element_type const& fieldRho() const { return this->fieldDensity(); }
    element_type const& fieldDensity() const { return *M_fieldDensity; }
    element_ptrtype const& fieldDensityPtr() const { return M_fieldDensity; }

    double cstNu() const { return this->cstCinematicViscosity(); }
    double cstCinematicViscosity() const { return M_cstCinematicViscosity; }
    void setCstCinematicViscosity(double d) { M_cstCinematicViscosity=d;this->updateCinematicViscosity(cst(d)); }
    element_type const& fieldNu() const { return this->fieldCinematicViscosity(); }
    element_type const& fieldCinematicViscosity() const { return *M_fieldCinematicViscosity; }
    element_ptrtype const& fieldCinematicViscosityPtr() const { return M_fieldCinematicViscosity; }


    void setCstDynamicViscosity(double d)
    {
        super_type::setCstDynamicViscosity(d);
        this->updateCinematicViscosity();
    }

    template < typename ExprT >
    void updateDensity(vf::Expr<ExprT> const& __expr)
    {
        if ( M_fieldDensity )
            M_fieldDensity->on(_range=elements(M_fieldDensity->mesh()),_expr=__expr );
    }
    template < typename ExprT >
    void updateCinematicViscosity(vf::Expr<ExprT> const& __expr)
    {
        if ( M_fieldCinematicViscosity )
            M_fieldCinematicViscosity->on(_range=elements(M_fieldCinematicViscosity->mesh()),_expr=__expr );
    }
    template < typename ExprT >
    void updateDynamicViscosity(vf::Expr<ExprT> const& __expr)
    {
        super_type::updateDynamicViscosity(__expr);
        this->updateCinematicViscosity();
    }

    void updateCinematicViscosity()
    {
        M_cstCinematicViscosity = this->cstMu()/this->cstRho();
        if ( M_fieldCinematicViscosity )
            M_fieldCinematicViscosity->on(_range=elements(M_fieldCinematicViscosity->mesh()),_expr=idv(this->fieldMu())/idv(this->fieldRho()) );
    }

    void updateFromModelMaterials( ModelMaterials const& mat )
    {
        if ( mat.empty() ) return;
        CHECK( mat.size() == 1 ) << "TODO multi-mat";
        super_type::updateFromModelMaterials( mat );
        for( auto const& m : mat )
        {
            auto const& mat = m.second;
            auto const& matmarker = m.first;
            this->setCstDensity( mat.rho() );
        }
    }


private :

    double M_cstDensity;// rho
    element_ptrtype M_fieldDensity;// rho
    double M_cstCinematicViscosity;// nu
    element_ptrtype M_fieldCinematicViscosity;// nu

};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_FLUIDMECHANICS_DENSITYVISCOSITYMODEL_H
