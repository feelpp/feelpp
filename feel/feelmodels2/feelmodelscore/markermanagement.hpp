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
 \file markermanagement.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2012-01-19
 */

#ifndef MARKERMANAGEMENT_HPP
#define MARKERMANAGEMENT_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {
namespace FeelModels {

class MarkerManagementDirichletBC
{
public :

    MarkerManagementDirichletBC();
    MarkerManagementDirichletBC( MarkerManagementDirichletBC const& op ) = default;

    void clearMarkerDirichletBC();

    void setMarkerDirichletBCByNameId( std::string type,std::string markerNameId,std::list<std::string> const& markers, ComponentType ct = ComponentType::NO_COMPONENT );
    void addMarkerDirichletBC(std::string type,std::string markerNameId, ComponentType ct = ComponentType::NO_COMPONENT);

    bool hasMarkerDirichletBC( std::string type, ComponentType ct = ComponentType::NO_COMPONENT) const;
    bool hasMarkerDirichletBCelimination( ComponentType ct ) const;
    bool hasMarkerDirichletBCnitsche( ComponentType ct ) const;
    bool hasMarkerDirichletBClm( ComponentType ct ) const;
    bool hasMarkerDirichletBCelimination() const;
    bool hasMarkerDirichletBCnitsche() const;
    bool hasMarkerDirichletBClm() const;

    std::map<std::string,std::list<std::string> > const& markerDirichletBCByType( ComponentType ct = ComponentType::NO_COMPONENT ) const;
    std::list<std::string> const& markerDirichletBCByNameId(std::string type,std::string markerNameId, ComponentType ct = ComponentType::NO_COMPONENT ) const;
    std::list<std::string> const& markerDirichletBCelimination( ComponentType ct = ComponentType::NO_COMPONENT ) const;
    std::list<std::string> const& markerDirichletBCnitsche( ComponentType ct = ComponentType::NO_COMPONENT ) const;
    std::list<std::string> const& markerDirichletBClm( ComponentType ct = ComponentType::NO_COMPONENT ) const;

    std::string getInfoDirichletBC() const;

private :
    void updateForUseMarkerDirichletBC();

    std::map<ComponentType,std::map<std::string,std::map<std::string,std::list<std::string> > > > M_dirichletBCType;
    std::map<ComponentType,std::map<std::string,std::list<std::string> > > M_dirichletBCMarkersListByType;
    std::list<std::string> M_listMarkerEmpty;
};


class MarkerManagementNeumannBC
{
public :
    enum NeumannBCShape { SCALAR = 0, VECTORIAL = 1 };

    MarkerManagementNeumannBC();
    MarkerManagementNeumannBC( MarkerManagementNeumannBC const& op ) = default;

    void clearMarkerNeumannBC();

    void setMarkerNeumannBC( NeumannBCShape shape, std::string markerNameId,std::list<std::string> const& markers );
    void addMarkerNeumannBC( NeumannBCShape shape, std::string markerNameId);

    std::map<std::string,std::list<std::string> > const& markerNeumannBC( NeumannBCShape shape ) const;
    std::list<std::string> const& markerNeumannBC( NeumannBCShape shape, std::string markerNameId ) const;

    std::string getInfoNeumannBC() const;

private :
    std::map<NeumannBCShape,std::map<std::string,std::list<std::string> > > M_containerMarkers;
};

class MarkerManagementALEMeshBC
{
public :

    MarkerManagementALEMeshBC();
    MarkerManagementALEMeshBC( MarkerManagementALEMeshBC const& op ) = default;

    void clearMarkerALEMeshBC();

    void setMarkerALEMeshBC( std::string type, std::list<std::string> const& markers );
    void addMarkerALEMeshBC(std::string type, std::string markerName);

    std::map<std::string,std::list<std::string> > const& markerALEMeshBC() const;
    std::list<std::string> const& markerALEMeshBC( std::string type ) const;

    std::string getInfoALEMeshBC() const;

private :
    std::map<std::string,std::list<std::string> > M_containerMarkers;
};

} // namespace FeelModels
} // namespace Feel


#endif // MARKERMANAGEMENT_HPP
