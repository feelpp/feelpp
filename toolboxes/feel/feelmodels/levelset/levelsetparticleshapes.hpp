/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <metivet@math.unistra.fr>
 Date: 2019-03-04

 Copyright (C) 2019 Université de Strasbourg

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
 \file levelsetparticleshapes.hpp
 \author Thibaut Metivet <metivet@math.unistra.fr>
 \date 2019-03-04
 */
#ifndef _LEVELSETPARTICLESHAPES_HPP
#define _LEVELSETPARTICLESHAPES_HPP 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelmodels/levelset/parameter_map.hpp>

namespace Feel {
namespace FeelModels {

enum class LevelSetShapeType {
    SPHERE, ELLIPSE
};
static const std::map<std::string, LevelSetShapeType> LevelSetShapeTypeIdMap = {
    { "sphere", LevelSetShapeType::SPHERE },
    { "ellipse", LevelSetShapeType::ELLIPSE }
};

template< typename FunctionSpaceType >
class LevelSetParticleShapes
{
public:
    typedef LevelSetParticleShapes<FunctionSpaceType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef FunctionSpaceType functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::value_type value_type;
    //--------------------------------------------------------------------//
    // Mesh
    static const uint16_type nDim = functionspace_type::nDim;
    typedef typename functionspace_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // Levelset element
    typedef typename functionspace_type::element_type element_type;
    typedef typename functionspace_type::element_ptrtype element_ptrtype;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//

    //--------------------------------------------------------------------//
    // Constructor
    LevelSetParticleShapes( functionspace_ptrtype const& space );

    //--------------------------------------------------------------------//
    // Accessors
    functionspace_ptrtype const& functionSpace() const { return M_functionSpace; }
    //--------------------------------------------------------------------//
    // Parameters
    bool redistShape() const { return M_redistShape; }
    void setRedistShape( bool b ) { M_redistShape = b; }

    //--------------------------------------------------------------------//
    // Create
    element_type create( LevelSetShapeType shape, parameter_map const& params, bool redist ) const;
    element_type create( LevelSetShapeType shape, parameter_map const& params ) const;
    element_type createSphere( parameter_map const& params ) const;
    element_type createEllipse( parameter_map const& params ) const;

    //--------------------------------------------------------------------//
    // Redistantiation
    element_type redistanciate( element_type const& phi ) const;

    //--------------------------------------------------------------------//
    // Read parameters
    parameter_map readShapeParams( LevelSetShapeType shape, boost::property_tree::ptree const& pt ) const;
    parameter_map readShapeParams( std::string const& shape, boost::property_tree::ptree const& pt ) const;
    parameter_map readSphereParams( boost::property_tree::ptree const& pt ) const;
    parameter_map readEllipseParams( boost::property_tree::ptree const& pt ) const;

private:
    //--------------------------------------------------------------------//
    // Levelset
    functionspace_ptrtype M_functionSpace;
    //--------------------------------------------------------------------//
    // Parameters
    bool M_redistShape;

};

#define LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS \
    template< typename FunctionSpaceType > \
    /**/
#define LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE \
    LevelSetParticleShapes< FunctionSpaceType > \
    /**/

LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS
LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::LevelSetParticleShapes( functionspace_ptrtype const& space )
    : M_functionSpace( space ),
      M_redistShape( false )
{
}

LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::element_type
LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::create( LevelSetShapeType shape, parameter_map const& params, bool redist ) const
{
    element_type phi;
    switch( shape )
    {
        case LevelSetShapeType::SPHERE:
        {
            phi = this->createSphere( params );
        }
        break;
        case LevelSetShapeType::ELLIPSE:
        {
            phi = this->createEllipse( params );
        }
        break;
    }
    if( redist )
    {
        return this->redistanciate( phi );
    }
    else
    {
        return phi;
    }
}

LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::element_type
LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::create( LevelSetShapeType shape, parameter_map const& params ) const
{
    return this->create( shape, params, this->redistShape() );
}

LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::element_type
LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::createSphere( parameter_map const& params ) const
{
    auto X = Px() - params.dget("xc");
    auto Y = Py() - params.dget("yc");
    auto Z = Pz() - params.dget("zc"); 
    auto R = params.dget("radius");
    return vf::project(
            _space=this->functionSpace(),
            _range=elements(this->functionSpace()->mesh()),
            _expr=sqrt(X*X+Y*Y+Z*Z)-R
            );
}

LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::element_type
LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::createEllipse( parameter_map const& params ) const
{
    auto X = Px() - params.dget("xc");
    auto Y = Py() - params.dget("yc");
    auto Z = Pz() - params.dget("zc");
    double A = params.dget("a");
    double B = params.dget("b");
    double C = params.dget("c");
    double psi = params.dget("psi");
    double theta = params.dget("theta");
    // Apply inverse ZYX TaitBryan rotation
    double cosPsi = std::cos(psi); double sinPsi = std::sin(psi);
    double cosTheta = std::cos(theta); double sinTheta = std::sin(theta);
    auto Xp = cosTheta*(cosPsi*X+sinPsi*Y) + sinTheta*Z;
    auto Yp = -sinPsi*X + cosPsi*Y;
    auto Zp = -sinTheta*(cosPsi*X+sinPsi*Y) + cosTheta*Z;
    // Project
    return vf::project(
            _space=this->functionSpace(),
            _range=elements(this->functionSpace()->mesh()),
            _expr=sqrt(Xp*Xp+Yp*Yp*(A*A)/(B*B)+Zp*Zp*(A*A)/(C*C))-A
            );
}

LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::element_type
LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::redistanciate( element_type const& phi ) const
{
    auto phiILP = vf::project(
            _space=this->functionSpace(),
            _range=this->functionSpace()->template rangeElements<0>(),
            _expr=idv(phi)/sqrt( inner( gradv(phi), gradv(phi) ) )
            );
    //return this->reinitializerFM()->run( phiILP );
    return phiILP;
}

LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS
parameter_map
LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::readShapeParams( LevelSetShapeType shape, boost::property_tree::ptree const& pt ) const
{
    switch( shape )
    {
        case LevelSetShapeType::SPHERE:
        {
            return readSphereParams( pt );
        }
        break;
        case LevelSetShapeType::ELLIPSE:
        {
            return readEllipseParams( pt );
        }
        break;
    }
}

LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS
parameter_map
LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::readShapeParams( std::string const& shape, boost::property_tree::ptree const& pt ) const
{
    CHECK( LevelSetShapeTypeIdMap.count( shape ) ) << shape << " is not in the list of supported particle shapes\n";
    return this->readShapeParams( LevelSetShapeTypeIdMap.at( shape ), pt );
}

LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS
parameter_map
LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::readSphereParams( boost::property_tree::ptree const& pt ) const
{
    parameter_map shapeParams;
    
    double xc = pt.get( "xc", 0. );
    double yc = pt.get( "yc", 0. );
    double zc = pt.get( "zc", 0. );
    double radius = pt.get( "radius", 1. );

    shapeParams["xc"] = xc;
    shapeParams["yc"] = yc;
    shapeParams["zc"] = zc;
    shapeParams["radius"] = radius;
    return shapeParams;
}

LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS
parameter_map
LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE::readEllipseParams( boost::property_tree::ptree const& pt ) const
{
    parameter_map shapeParams;
    
    double xc = pt.get( "xc", 0. );
    double yc = pt.get( "yc", 0. );
    double zc = pt.get( "zc", 0. );
    double a = pt.get( "a", 1. );
    double b = pt.get( "b", 1. );
    double c = pt.get( "c", 1. );
    double psi = pt.get( "psi", 0. );
    double theta = pt.get( "theta", 0. );

    shapeParams["xc"] = xc;
    shapeParams["yc"] = yc;
    shapeParams["zc"] = zc;
    shapeParams["a"] = a;
    shapeParams["b"] = b;
    shapeParams["c"] = c;
    shapeParams["psi"] = psi;
    shapeParams["theta"] = theta;
    return shapeParams;
}

#undef LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_TYPE
#undef LEVELSETPARTICLESHAPES_CLASS_TEMPLATE_DECLARATIONS

} // namespace FeelModels
} // namespace Feel

#endif // _LEVELSETPARTICLESHAPES_HPP


