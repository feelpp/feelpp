/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-11-09

  Copyright (C) 2005,2006 EPFL

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
// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
/**
   \file partitioner.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-11-09
 */
#ifndef __partitioner_h__
#define __partitioner_h__

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

namespace Feel
{
#if defined( FEELPP_INSTANTIATION_MODE )
void test_partitioner();
#endif // FEELPP_INSTANTIATION_MODE
/**
 * The \p Partitioner class provides a uniform interface for
 * partitioning algorithms.  It takes a reference to a \p mesh_type
 * object as input, which it will partition into a number of
 * subdomains.
 */
template<typename Mesh>
class Partitioner
{
public:

    typedef Mesh mesh_type;
    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::element_iterator element_iterator;

    struct Factory
    {
        typedef Feel::Singleton< Feel::Factory< Partitioner<Mesh>, std::string > > type;
    };
    /**
     * Constructor.
     */
    Partitioner ();

    /**
     * Destructor. Virtual so that we can derive from this class.
     */
    virtual ~Partitioner() {}

    /**
     * \return the communicator
     */
    mpi::communicator const& comm() const
    {
        return M_comm;
    }

    /**
     * create object \p partitioner from factory
     */
    static Partitioner<Mesh>* New( std::string const& partitioner )
    {
#if defined(FEELPP_INSTANTIATION_MODE)
        test_partitioner();
#endif
        return Factory::type::instance().createObject( partitioner );
    }

    /**
     * Partition the \p mesh_type into \p n parts corresponding to the number of
     * processors.  If the user does not specify a number of pieces into which
     * the mesh should be partitioned, then the default behavior of the
     * partitioner is to partition according to the number of processors defined
     * in libMesh::n_processors().  The partitioner currently does not modify
     * the subdomain_id of each element.  This number is reserved for things
     * like material properties, etc.
     */
    void partition ( Mesh& mesh )
    {
        partition( mesh, M_comm.size() );
    }

    /**
     * Partition the \p mesh_type into \p n parts.  If the user does not specify
     * a number of pieces into which the mesh should be partitioned, then the
     * default behavior of the partitioner is to partition according to the
     * number of processors defined in libMesh::n_processors().  The partitioner
     * currently does not modify the subdomain_id of each element.  This number
     * is reserved for things like material properties, etc.
     */
    void partition ( Mesh& mesh, size_type n_parts );

    /**
     * Repartitions the \p mesh_type into \p n parts corresponding to the number
     * of processors.  This is required since some partitoning algorithms can
     * repartition more efficiently than computing a new partitioning from
     * scratch.  The default behavior is to simply call this->partition(n)
     */
    void repartition ( mesh_type& mesh )
    {
        repartition( mesh, M_comm.size() );
    }


    /**
     * Repartitions the \p mesh_type into \p n parts.  This
     * is required since some partitoning algorithms can repartition
     * more efficiently than computing a new partitioning from scratch.
     * The default behavior is to simply call this->partition(n)
     */
    void repartition ( mesh_type& mesh, size_type n_parts );


protected:

    /**
     * Trivially "partitions" the mesh for one processor.
     * Simply loops through the elements and assigns all of them
     * to processor 0.  Is is provided as a separate function
     * so that derived classes may use it without reimplementing it.
     */
    void singlePartition ( mesh_type& mesh );

    /**
     * This is the actual partitioning method which must be overloaded
     * in derived classes.  It is called via the public partition()
     * method above by the user.
     */
    virtual void doPartition( mesh_type& mesh,
                              const size_type n ) = 0;

    /**
     * This is the actual re-partitioning method which can be overloaded
     * in derived classes.  Note that the default behavior is to simply
     * call the partition function.
     */
    virtual void doRepartition ( mesh_type& mesh,
                                 const size_type n )
    {
        this->doPartition ( mesh, n );
    }

private:

    mpi::communicator M_comm;
};

template<typename Mesh>
void
Partitioner<Mesh>::partition ( mesh_type& mesh,
                               size_type n )
{
    // Set the number of partitions in the mesh
    mesh.setNumberOfPartitions( n );

    // Call the partitioning function
    this->doPartition( mesh, n );
}



template<typename Mesh>
void
Partitioner<Mesh>::repartition ( mesh_type& mesh,
                                 size_type n )
{
    // Set the number of partitions in the mesh
    mesh.setNumberOfPartitions( n );

    // Call the partitioning function
    this->doRepartition( mesh, n );
}



template<typename Mesh>
void
Partitioner<Mesh>::singlePartition ( mesh_type& mesh )
{
    element_iterator       elem_it  = mesh.beginElement();
    const element_iterator elem_end = mesh.endElement();

    for ( ; elem_it != elem_end; ++elem_it )
    {
        element_type elt = *elem_it;
        elt.setProcessId( 0 );
        mesh.elements().replace( elem_it, elt );
    }
}

} // Feel

#if !defined( FEELPP_INSTANTIATION_MODE )
# include <feel/feeldiscr/partitioner.cpp>
#endif
#endif
