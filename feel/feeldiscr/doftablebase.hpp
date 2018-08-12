//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Vincent Chabannes <vincent.chabannes@feelpp.org>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_DOFTABLEBASE_HPP
#define FEELPP_DOFTABLEBASE_HPP 1

#include <feel/feelalg/datamap.hpp>
#include <feel/feeldiscr/dof.hpp>
#include <feel/feelmesh/meshsupportbase.hpp>

namespace Feel
{

struct DofTableInfos
{
    DofTableInfos() = default;
    DofTableInfos( DofTableInfos const& ) = default;
    DofTableInfos( DofTableInfos && ) = default;
    DofTableInfos& operator=( DofTableInfos const& ) = default;
    DofTableInfos& operator=( DofTableInfos && ) = default;

    uint16_type nOrder;
    uint16_type nDim;
    uint16_type nRealDim;
    uint16_type Shape;
    uint16_type nComponents;
    uint16_type nComponents1;
    uint16_type nComponents2;


    bool is_continuous;
    bool is_discontinuous_locally;
    bool is_discontinuous_totally;

    bool is_scalar;
    bool is_vectorial;
    bool is_tensor2;
    bool is_tensor2symm;
    bool is_modal;
    bool is_product;
    uint16_type nRealComponents;

    bool is_p0_continuous;

    bool is_hdiv_conforming;
    bool is_hcurl_conforming;

    uint16_type nDofPerEdge;
    uint16_type nDofPerElement;

    bool is_periodic;

    uint16_type nDofComponents;

    std::string feFamilyName;
};

class DofTableBase : public DataMap
{
    typedef DataMap super;
public:
    typedef Dof global_dof_type;
    typedef FaceDof global_dof_fromface_type;
    typedef std::shared_ptr<MeshSupportBase> mesh_support_base_ptrtype;

    DofTableBase( WorldComm const& _worldComm )
        :
        super( _worldComm.clone() )
        {}
    DofTableBase( worldcomm_ptr_t& _worldComm )
        :
        super( _worldComm )
        {}
    DofTableBase( DofTableBase const& ) = default;
    DofTableBase( DofTableBase && ) = default;
    DofTableBase& operator=( DofTableBase const& ) = default;
    DofTableBase& operator=( DofTableBase && ) = default;

    virtual DofTableInfos infos() const = 0;

    virtual global_dof_type const& localToGlobal( const size_type ElId,
                                                  const uint16_type localNode,
                                                  const uint16_type c = 0 ) const = 0;

    virtual global_dof_fromface_type const& faceLocalToGlobal( const size_type ElId,
                                                               const uint16_type localNode,
                                                               const uint16_type c = 0 ) const = 0;

    virtual mesh_support_base_ptrtype meshSupportBase() const = 0;

};

}
#endif
