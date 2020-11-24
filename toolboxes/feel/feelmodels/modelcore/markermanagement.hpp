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
 \file markermanagement.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2012-01-19
 */

#ifndef FEELPP_TOOLBOXES_CORE_MARKERMANAGEMENT_HPP
#define FEELPP_TOOLBOXES_CORE_MARKERMANAGEMENT_HPP 1

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

    void setMarkerDirichletBCByNameId( std::string const& type,std::string const& name,std::set<std::string> const& markers, ComponentType ct = ComponentType::NO_COMPONENT );
    void addMarkerDirichletBC(std::string const& type, std::string const& name, std::string const& marker, ComponentType ct = ComponentType::NO_COMPONENT);
    void addMarkerDirichletBC(std::string const& type, std::string const& name, std::set<std::string> const& markers, ComponentType ct = ComponentType::NO_COMPONENT);

    bool hasMarkerDirichletBC( std::string const& type, ComponentType ct = ComponentType::NO_COMPONENT) const;
    bool hasMarkerDirichletBCelimination( ComponentType ct ) const;
    bool hasMarkerDirichletBCnitsche( ComponentType ct ) const;
    bool hasMarkerDirichletBClm( ComponentType ct ) const;
    bool hasMarkerDirichletBCelimination() const;
    bool hasMarkerDirichletBCnitsche() const;
    bool hasMarkerDirichletBClm() const;

    std::map<std::string,std::set<std::string> > const& markerDirichletBCByType( ComponentType ct = ComponentType::NO_COMPONENT ) const;
    std::set<std::string> const& markerDirichletBCByNameId(std::string const& type,std::string const& markerNameId, ComponentType ct = ComponentType::NO_COMPONENT ) const;
    std::set<std::string> const& markerDirichletBCelimination( ComponentType ct = ComponentType::NO_COMPONENT ) const;
    std::set<std::string> const& markerDirichletBCnitsche( ComponentType ct = ComponentType::NO_COMPONENT ) const;
    std::set<std::string> const& markerDirichletBClm( ComponentType ct = ComponentType::NO_COMPONENT ) const;

    std::string getInfoDirichletBC() const;
    void updateInformationObjectDirichletBC( pt::ptree & p ) const;

private :
    void updateForUseMarkerDirichletBC();

    std::map<ComponentType,std::map<std::string,std::map<std::string,std::set<std::string> > > > M_dirichletBCType;
    std::map<ComponentType,std::map<std::string,std::set<std::string> > > M_dirichletBCMarkersListByType;
    std::set<std::string> M_listMarkerEmpty;
};


class MarkerManagementNeumannBC
{
public :
    enum NeumannBCShape { SCALAR = 0, VECTORIAL = 1, TENSOR2 = 2 };

    MarkerManagementNeumannBC();
    MarkerManagementNeumannBC( MarkerManagementNeumannBC const& op ) = default;

    void clearMarkerNeumannBC();

    void setMarkerNeumannBC( NeumannBCShape shape, std::string const& name,std::set<std::string> const& markers );
    void addMarkerNeumannBC( NeumannBCShape shape, std::string const& name, std::string const& marker);
    void addMarkerNeumannBC( NeumannBCShape shape, std::string const& name, std::set<std::string> const& markers);

    std::map<std::string,std::set<std::string> > const& markerNeumannBC( NeumannBCShape shape ) const;
    std::set<std::string> const& markerNeumannBC( NeumannBCShape shape, std::string const& markerNameId ) const;

    std::string getInfoNeumannBC() const;
    void updateInformationObjectNeumannBC( pt::ptree & p ) const;

private :
    std::map<NeumannBCShape,std::map<std::string,std::set<std::string> > > M_containerMarkers;
    std::set<std::string> M_listMarkerEmpty;
};

class MarkerManagementNeumannEulerianFrameBC
{
public :
    enum NeumannEulerianFrameBCShape { SCALAR = 0, VECTORIAL = 1, TENSOR2 = 2 };

    MarkerManagementNeumannEulerianFrameBC();
    MarkerManagementNeumannEulerianFrameBC( MarkerManagementNeumannEulerianFrameBC const& op ) = default;

    void clearMarkerNeumannEulerianFrameBC();

    void setMarkerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape, std::string const& name,std::set<std::string> const& markers );
    void addMarkerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape, std::string const& name, std::string const& marker);
    void addMarkerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape, std::string const& name, std::set<std::string> const& markers);

    std::map<std::string,std::set<std::string> > const& markerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape ) const;
    std::set<std::string> const& markerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape, std::string const& markerNameId ) const;

    std::string getInfoNeumannEulerianFrameBC() const;
    void updateInformationObjectNeumannEulerianFrameBC( pt::ptree & p ) const;
private :
    std::map<NeumannEulerianFrameBCShape,std::map<std::string,std::set<std::string> > > M_containerMarkers;
    std::set<std::string> M_listMarkerEmpty;
};

class MarkerManagementALEMeshBC
{
public :

    MarkerManagementALEMeshBC();
    MarkerManagementALEMeshBC( MarkerManagementALEMeshBC const& op ) = default;

    void clearMarkerALEMeshBC();

    void setMarkerALEMeshBC( std::string const& type, std::set<std::string> const& markers );
    void addMarkerALEMeshBC(std::string const& type, std::string const& markerName);
    void addMarkerALEMeshBC( std::string const& type, std::set<std::string> const& markers );

    std::map<std::string,std::set<std::string> > const& markerALEMeshBC() const;
    std::set<std::string> const& markerALEMeshBC( std::string const& type ) const;

    std::string getInfoALEMeshBC() const;
    void updateInformationObjectALEMeshBC( pt::ptree & p ) const;

private :
    std::map<std::string,std::set<std::string> > M_containerMarkers;
    std::set<std::string> M_listMarkerEmpty;
};

class MarkerManagementSlipBC
{
public :

    MarkerManagementSlipBC();
    MarkerManagementSlipBC( MarkerManagementSlipBC const& op ) = default;
    void clearMarkerSlipBC();
    void setMarkerSlipBC( std::set<std::string> const& markers );
    void addMarkerSlipBC( std::string const& markerName);
    void addMarkerSlipBC( std::set<std::string> const& markers );
    std::set<std::string> const& markerSlipBC() const;
    std::string getInfoSlipBC() const;
    void updateInformationObjectSlipBC( pt::ptree & p ) const;
private :
    std::set<std::string> M_containerMarkers;
    std::set<std::string> M_listMarkerEmpty;
};

class MarkerManagementPressureBC
{
public :

    MarkerManagementPressureBC();
    MarkerManagementPressureBC( MarkerManagementPressureBC const& op ) = default;
    void clearMarkerPressureBC();
    void setMarkerPressureBC( std::string const& name, std::set<std::string> const& markers );
    void addMarkerPressureBC( std::string const& name, std::string const& marker );
    void addMarkerPressureBC( std::string const& name, std::set<std::string> const& markers );
    std::set<std::string> const& markerPressureBC() const;
    std::set<std::string> const& markerPressureBC( std::string const& markerNameId ) const;
    bool hasMarkerPressureBC() const;
    std::string getInfoPressureBC() const;
    void updateInformationObjectPressureBC( pt::ptree & p ) const;
private :
    std::map<std::string,std::set<std::string> > M_containerMarkers;
    std::set<std::string> M_listMarkers;
    std::set<std::string> M_listMarkerEmpty;
};

class MarkerManagementRobinBC
{
public :

    MarkerManagementRobinBC();
    MarkerManagementRobinBC( MarkerManagementRobinBC const& op ) = default;
    void clearMarkerRobinBC();
    void setMarkerRobinBC( std::string const& name, std::set<std::string> const& markers );
    void addMarkerRobinBC( std::string const& name, std::string const& marker );
    void addMarkerRobinBC( std::string const& name, std::set<std::string> const& markers );
    std::map<std::string,std::set<std::string> > const& markerRobinBC() const;
    std::set<std::string> const& markerRobinBC( std::string const& markerNameId ) const;
    std::string getInfoRobinBC() const;
    void updateInformationObjectRobinBC( pt::ptree & p ) const;
private :
    std::map<std::string,std::set<std::string> > M_containerMarkers;
    std::set<std::string> M_listMarkerEmpty;
};

class MarkerManagementFluidStructureInterfaceBC
{
public :

    MarkerManagementFluidStructureInterfaceBC();
    MarkerManagementFluidStructureInterfaceBC( MarkerManagementFluidStructureInterfaceBC const& op ) = default;
    void clearMarkerFluidStructureInterfaceBC();
    void setMarkerFluidStructureInterfaceBC( std::set<std::string> const& markers );
    void addMarkerFluidStructureInterfaceBC( std::string const& markerName );
    void addMarkerFluidStructureInterfaceBC( std::set<std::string> const& markers );
    std::set<std::string> const& markerFluidStructureInterfaceBC() const;
    std::string getInfoFluidStructureInterfaceBC() const;
    void updateInformationObjectFluidStructureInterfaceBC( pt::ptree & p ) const;
private :
    std::set<std::string> M_containerMarkers;
    std::set<std::string> M_listMarkerEmpty;
};


namespace detail
{

template <typename MeshType>
std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> >
distributeMarkerListOnSubEntity( std::shared_ptr<MeshType> const& mesh, std::set<std::string> const& listMarker )
{
    std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > res;
    for ( std::string const& marker : listMarker )
    {
        if ( !mesh->hasMarker( marker ) ) continue;

        if ( mesh->hasFaceMarker( marker ) )
        {
            std::get<0>( res ).insert( marker );
            //std::cout << "has face marker " << marker << "\n";
        }
        else if ( mesh->hasEdgeMarker( marker ) )
        {
            std::get<1>( res ).insert( marker );
            //std::cout << "has edge marker " << marker << "\n";
        }
        else if ( mesh->hasPointMarker( marker ) )
        {
            std::get<2>( res ).insert( marker );
            //std::cout << "has point marker " << marker << "\n";
        }
        else
        {
            std::get<3>( res ).insert( marker );
            //std::cout << "has element marker " << marker << "\n";
        }
    }
    return res;
}

// template <typename MeshType>
// std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> >
// distributeMarkerListOnSubEntity( std::shared_ptr<MeshType> const& mesh, std::set<std::string> const& listMarker )
// {
//     std::set<std::string> setOfMarkers(listMarker.begin(),listMarker.end());
//     return distributeMarkerListOnSubEntity( mesh,setOfMarkers );
// }

template <typename MeshType>
std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> >
distributeMarkerListOnSubEntity( std::shared_ptr<MeshType> const& mesh, std::initializer_list< std::set<std::string> > const& listOflistMarker )
{
    std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > res;
    std::set<std::string> setOfMarkers;
    for ( auto const& listMarker : listOflistMarker )
    {
        for (std::string const& mark : listMarker )
            setOfMarkers.insert( mark );
    }
    return distributeMarkerListOnSubEntity( mesh,setOfMarkers );
}


} // namespace detail


} // namespace FeelModels
} // namespace Feel


#endif // MARKERMANAGEMENT_HPP
