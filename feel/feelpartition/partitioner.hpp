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
#ifndef FEELPP_PARTITIONER_HPP
#define FEELPP_PARTITIONER_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feeldiscr/mesh.hpp>

namespace Feel
{

/**
 * The \p Partitioner class provides a uniform interface for
 * partitioning algorithms.  It takes a reference to a \p MeshBase
 * object as input, which it will partition into a number of
 * subdomains.
 */
template <typename MeshType>
class Partitioner
{
  public:
    using mesh_type = decay_type<MeshType>;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;
    using partitioner_type = Partitioner<mesh_type>;
    using clone_ptrtype = std::unique_ptr<partitioner_type>;

    using element_type = typename mesh_type::element_type;
    using face_type = typename mesh_type::face_type;
    using point_type = typename mesh_type::point_type;

    /**
     * Constructor.
     */
    Partitioner()
        : M_weights() {}

    /**
     * Destructor. Virtual so that we can derive from this class.
     */
    virtual ~Partitioner() = default;

    /**
     * Creates a new partitioner of this type and returns it in
     * an \p std::unique_ptr.
     * This is used when copying meshes, and must be overloaded in the
     * derived classes.
     */
    virtual clone_ptrtype clone() const = 0;

    /**
     * Partition the \p Mesh into \p n parts.
     * The partitioner currently does not modify the marker1 but rather marker3
     * of each element.  \p marker1 is reserved for things like
     * material properties, etc.
     */
    void partition( mesh_ptrtype mesh, rank_type n );

    /**
     * Partition the \p Mesh into \p Environment::numberOfProcessors() parts.
     * The partitioner currently does not modify the subdomain_id
     * of each element.  This number is reserved for things like
     * material properties, etc.
     */
    void partition( mesh_ptrtype mesh );

    /**
     * Repartitions the \p Mesh into \p n parts.  This
     * is required since some partitoning algorithms can repartition
     * more efficiently than computing a new partitioning from scratch.
     * The default behavior is to simply call this->partition(mesh,n)
     */
    void repartition( mesh_ptrtype mesh, rank_type n );

    /**
     * Repartitions the \p Mesh into \p Environment::numberOfProcessors() parts.
     * This is required since some partitoning algorithms can repartition more
     * efficiently than computing a new partitioning from scratch.
     */
    void repartition( mesh_ptrtype mesh );

    /**
     * This function is called after partitioning to set the processor IDs
     * for the nodes.  By definition, a Node's processor ID is the minimum
     * processor ID for all of the elements which share the node.
     */
    static void setNodeProcessorIds( mesh_ptrtype mesh );

    /**
     * Attach weights that can be used for partitioning.  This ErrorVector should be
     * _exactly_ the same on every processor and should have mesh->max_elem_id()
     * entries.
     */
    virtual void attachWeights( std::vector<double> const& weights )
    {
        CHECK( 0 ) << "attachWeights() not implemented in base class Partitioner";
    }

  protected:
    /**
     * Trivially "partitions" the mesh for one processor.
     * Simply loops through the elements and assigns all of them
     * to processor 0.  Is is provided as a separate function
     * so that derived classes may use it without reimplementing it.
     */
    void singlePartition( mesh_ptrtype mesh );

    /**
     * This is the actual partitioning method which must be overloaded
     * in derived classes.  It is called via the public partition()
     * method above by the user.
     */
    virtual void partitionImpl( mesh_ptrtype mesh, rank_type n ) = 0;

    /**
     * This is the actual re-partitioning method which can be overloaded
     * in derived classes.  Note that the default behavior is to simply
     * call the partition function.
     */
    virtual void repartitionImpl( mesh_ptrtype mesh, rank_type n )
    {
        this->partitionImpl( mesh, n );
    }

    /**
     * The blocksize to use when doing blocked parallel communication.  This limits the
     * maximum vector size which can be used in a single communication step.
     */
    static const dof_id_type communication_blocksize;

    /**
     * The weights that might be used for partitioning.
     */
    std::vector<double> M_weights;
};

} // namespace libMesh

#include <feel/feelpartition/partitioner_impl.hpp>

#endif
