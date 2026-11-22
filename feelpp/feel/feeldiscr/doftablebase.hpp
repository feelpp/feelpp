/*
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2017-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

//! @file
//! @author Vincent Chabannes <vincent.chabannes@feelpp.org>
//! @date 15 Jun 2017

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

    uint16_type nOrder = 0;
    uint16_type nDim = 0;
    uint16_type nRealDim = 0;
    uint16_type Shape = 0;
    uint16_type nComponents = 0;
    uint16_type nComponents1 = 0;
    uint16_type nComponents2 = 0;


    bool is_continuous = false;
    bool is_discontinuous_locally = false;
    bool is_discontinuous_totally = false;

    bool is_scalar = false;
    bool is_vectorial = false;
    bool is_tensor2 = false;
    bool is_tensor2symm = false;
    bool is_modal = false;
    bool is_product = false;
    uint16_type nRealComponents = 0;

    bool is_p0_continuous = false;

    bool is_hdiv_conforming = false;
    bool is_hcurl_conforming = false;

    uint16_type nDofPerEdge = 0;
    uint16_type nDofPerElement = 0;

    bool is_periodic = false;

    uint16_type nDofComponents = 0;

    std::string feFamilyName;
};

template<typename SizeT=uint32_type>
class DofTableBase : public DataMap<SizeT>
{
    typedef DataMap<SizeT> super;
public:
    using size_type = SizeT;
    typedef Dof<size_type> global_dof_type;
    
    typedef FaceDof<size_type> global_dof_fromface_type;
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

    global_dof_fromface_type const& facetLocalToGlobal( const size_type ElId,
                                                        const uint16_type localNode,
                                                        const uint16_type c = 0 ) const
        {
            return this->faceLocalToGlobal( ElId, localNode, c );
        }

    virtual mesh_support_base_ptrtype meshSupportBase() const = 0;

};

}
#endif
