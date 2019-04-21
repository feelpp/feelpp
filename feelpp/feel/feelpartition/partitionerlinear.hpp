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
#ifndef FEELPP_PARTITIONERLINEAR_HPP
#define FEELPP_PARTITIONERLINEAR_HPP 1

namespace Feel {

/**
 * @brief split mesh into equal size parts associated to each processor
 * @note use with care
 */
template<typename MeshType>
class PartitionerLinear : public Partitioner<MeshType>
{
public:
    using super = Partitioner<MeshType>;
    using clone_ptrtype = typename super::clone_ptrtype;
    using mesh_ptrtype = typename super::mesh_ptrtype;
    using element_type = typename super::element_type;
    /**
     * Constructor.
     */
    PartitionerLinear () {}

    /**
     * Creates a new partitioner of this type 
     */
    virtual clone_ptrtype clone () const override
        {
            return std::make_unique<PartitionerLinear>();
        }

protected:
    /**
     * Partition the \p MeshBase<> into \p n subdomains.
     */
    virtual void partitionImpl ( mesh_ptrtype& mesh, rank_type n );

private:

};



} // Feel

#include <feel/feelpartition/partitionerlinear_impl.hpp>


#endif
