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

    void setMarkerDirichletBCByNameId( std::string type,std::string name,std::list<std::string> const& markers, ComponentType ct = ComponentType::NO_COMPONENT );
    void addMarkerDirichletBC(std::string type, std::string name, std::string marker, ComponentType ct = ComponentType::NO_COMPONENT);
    void addMarkerDirichletBC(std::string type, std::string name, std::list<std::string> const& markers, ComponentType ct = ComponentType::NO_COMPONENT);

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
    enum NeumannBCShape { SCALAR = 0, VECTORIAL = 1, TENSOR2 = 2 };

    MarkerManagementNeumannBC();
    MarkerManagementNeumannBC( MarkerManagementNeumannBC const& op ) = default;

    void clearMarkerNeumannBC();

    void setMarkerNeumannBC( NeumannBCShape shape, std::string name,std::list<std::string> const& markers );
    void addMarkerNeumannBC( NeumannBCShape shape, std::string name, std::string marker);
    void addMarkerNeumannBC( NeumannBCShape shape, std::string name, std::list<std::string> markers);

    std::map<std::string,std::list<std::string> > const& markerNeumannBC( NeumannBCShape shape ) const;
    std::list<std::string> const& markerNeumannBC( NeumannBCShape shape, std::string markerNameId ) const;

    std::string getInfoNeumannBC() const;

private :
    std::map<NeumannBCShape,std::map<std::string,std::list<std::string> > > M_containerMarkers;
    std::list<std::string> M_listMarkerEmpty;
};

class MarkerManagementNeumannEulerianFrameBC
{
public :
    enum NeumannEulerianFrameBCShape { SCALAR = 0, VECTORIAL = 1, TENSOR2 = 2 };

    MarkerManagementNeumannEulerianFrameBC();
    MarkerManagementNeumannEulerianFrameBC( MarkerManagementNeumannEulerianFrameBC const& op ) = default;

    void clearMarkerNeumannEulerianFrameBC();

    void setMarkerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape, std::string name,std::list<std::string> const& markers );
    void addMarkerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape, std::string name, std::string marker);
    void addMarkerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape, std::string name, std::list<std::string> markers);

    std::map<std::string,std::list<std::string> > const& markerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape ) const;
    std::list<std::string> const& markerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape, std::string markerNameId ) const;

    std::string getInfoNeumannEulerianFrameBC() const;
private :
    std::map<NeumannEulerianFrameBCShape,std::map<std::string,std::list<std::string> > > M_containerMarkers;
    std::list<std::string> M_listMarkerEmpty;
};

class MarkerManagementALEMeshBC
{
public :

    MarkerManagementALEMeshBC();
    MarkerManagementALEMeshBC( MarkerManagementALEMeshBC const& op ) = default;

    void clearMarkerALEMeshBC();

    void setMarkerALEMeshBC( std::string type, std::list<std::string> const& markers );
    void addMarkerALEMeshBC(std::string type, std::string markerName);
    void addMarkerALEMeshBC( std::string type, std::list<std::string> const& markers );

    std::map<std::string,std::list<std::string> > const& markerALEMeshBC() const;
    std::list<std::string> const& markerALEMeshBC( std::string type ) const;

    std::string getInfoALEMeshBC() const;

private :
    std::map<std::string,std::list<std::string> > M_containerMarkers;
    std::list<std::string> M_listMarkerEmpty;
};

class MarkerManagementSlipBC
{
public :

    MarkerManagementSlipBC();
    MarkerManagementSlipBC( MarkerManagementSlipBC const& op ) = default;
    void clearMarkerSlipBC();
    void setMarkerSlipBC( std::list<std::string> const& markers );
    void addMarkerSlipBC( std::string markerName);
    void addMarkerSlipBC( std::list<std::string> const& markers );
    std::list<std::string> const& markerSlipBC() const;
    std::string getInfoSlipBC() const;
private :
    std::list<std::string> M_containerMarkers;
    std::list<std::string> M_listMarkerEmpty;
};

class MarkerManagementPressureBC
{
public :

    MarkerManagementPressureBC();
    MarkerManagementPressureBC( MarkerManagementPressureBC const& op ) = default;
    void clearMarkerPressureBC();
    void setMarkerPressureBC( std::string const& name, std::list<std::string> const& markers );
    void addMarkerPressureBC( std::string const& name, std::string const& marker );
    void addMarkerPressureBC( std::string const& name, std::list<std::string> const& markers );
    std::list<std::string> const& markerPressureBC() const;
    std::list<std::string> const& markerPressureBC( std::string const& markerNameId ) const;
    bool hasMarkerPressureBC() const;
    std::string getInfoPressureBC() const;
private :
    std::map<std::string,std::list<std::string> > M_containerMarkers;
    std::list<std::string> M_listMarkers;
    std::list<std::string> M_listMarkerEmpty;
};

class MarkerManagementRobinBC
{
public :

    MarkerManagementRobinBC();
    MarkerManagementRobinBC( MarkerManagementRobinBC const& op ) = default;
    void clearMarkerRobinBC();
    void setMarkerRobinBC( std::string const& name, std::list<std::string> const& markers );
    void addMarkerRobinBC( std::string const& name, std::string const& marker );
    void addMarkerRobinBC( std::string const& name, std::list<std::string> const& markers );
    std::map<std::string,std::list<std::string> > const& markerRobinBC() const;
    std::list<std::string> const& markerRobinBC( std::string const& markerNameId ) const;
    std::string getInfoRobinBC() const;
private :
    std::map<std::string,std::list<std::string> > M_containerMarkers;
    std::list<std::string> M_listMarkerEmpty;
};

class MarkerManagementFluidStructureInterfaceBC
{
public :

    MarkerManagementFluidStructureInterfaceBC();
    MarkerManagementFluidStructureInterfaceBC( MarkerManagementFluidStructureInterfaceBC const& op ) = default;
    void clearMarkerFluidStructureInterfaceBC();
    void setMarkerFluidStructureInterfaceBC( std::list<std::string> const& markers );
    void addMarkerFluidStructureInterfaceBC( std::string markerName );
    void addMarkerFluidStructureInterfaceBC( std::list<std::string> const& markers );
    std::list<std::string> const& markerFluidStructureInterfaceBC() const;
    std::string getInfoFluidStructureInterfaceBC() const;
private :
    std::list<std::string> M_containerMarkers;
    std::list<std::string> M_listMarkerEmpty;
};


namespace detail
{

template <typename MeshType>
std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string>,std::list<std::string> >
distributeMarkerListOnSubEntity( std::shared_ptr<MeshType> const& mesh, std::set<std::string> const& listMarker )
{
    std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string>,std::list<std::string> > res;
    for ( std::string const& marker : listMarker )
    {
        if ( !mesh->hasMarker( marker ) ) continue;

        if ( mesh->hasFaceMarker( marker ) )
        {
            std::get<0>( res ).push_back( marker );
            //std::cout << "has face marker " << marker << "\n";
        }
        else if ( mesh->hasEdgeMarker( marker ) )
        {
            std::get<1>( res ).push_back( marker );
            //std::cout << "has edge marker " << marker << "\n";
        }
        else if ( mesh->hasPointMarker( marker ) )
        {
            std::get<2>( res ).push_back( marker );
            //std::cout << "has point marker " << marker << "\n";
        }
        else
        {
            std::get<3>( res ).push_back( marker );
            //std::cout << "has element marker " << marker << "\n";
        }
    }
    return res;
}

template <typename MeshType>
std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string>,std::list<std::string> >
distributeMarkerListOnSubEntity( std::shared_ptr<MeshType> const& mesh, std::list<std::string> const& listMarker )
{
    std::set<std::string> setOfMarkers(listMarker.begin(),listMarker.end());
    return distributeMarkerListOnSubEntity( mesh,setOfMarkers );
}

template <typename MeshType>
std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string>,std::list<std::string> >
distributeMarkerListOnSubEntity( std::shared_ptr<MeshType> const& mesh, std::initializer_list< std::list<std::string> > const& listOflistMarker )
{
    std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string>,std::list<std::string> > res;
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
