/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 17 May 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_PARTITIONERMETIS_HPP
#define FEELPP_PARTITIONERMETIS_HPP 1

#include <feel/feelpartition/partitioner.hpp>

namespace Feel
{

/**
 * @brief Metis partitioner
 */
template <typename MeshType>
class PartitionerMetis : public Partitioner<MeshType>
{
  public:
    using super = Partitioner<MeshType>;
    using mesh_type = typename super::mesh_type;
    using mesh_ptrtype = typename super::mesh_ptrtype;
    using clone_ptrtype = typename super::clone_ptrtype;
    using element_type = typename super::element_type;
    /**
     * Constructor.
     */
    PartitionerMetis() {}

    /**
     * Creates a new partitioner of this type 
     */
    virtual clone_ptrtype clone() const override
    {
        return std::make_unique<PartitionerMetis<mesh_type>>();
    }

    void attachWeights( std::vector<double> const& weights ) override
    {
        this->M_weights = weights;
    }

  protected:
    /**
     * Partition the \p MeshBase into \p n subdomains.
     */
    void partitionImpl( mesh_ptrtype mesh, rank_type n ) override;
};

} // Feel

#include <feel/feelpartition/partitionermetis_impl.hpp>

#endif
