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
/**
   \file meshbase.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-11-09
 */
#ifndef __MeshBase_H
#define __MeshBase_H 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/context.hpp>
#include <feel/feelcore/environment.hpp>

namespace Feel
{
/**
 * Components of a mesh that can be enabled or disabled when calling
 * \c updateForUse()
 */
enum MeshComponents
{
    MESH_UPDATE_EDGES     = ( 1 << 0 ),
    MESH_UPDATE_FACES     = ( 1 << 1 ),
    MESH_CHECK            = ( 1 << 2 ),
    MESH_PARTITION        = ( 1 << 3 ),
    MESH_RENUMBER         = ( 1 << 4 )
};
const uint16_type MESH_ALL_COMPONENTS = MESH_UPDATE_EDGES | MESH_UPDATE_FACES | MESH_CHECK | MESH_PARTITION | MESH_RENUMBER;
const uint16_type MESH_COMPONENTS_DEFAULTS = MESH_RENUMBER | MESH_CHECK;

/**
 * \class MeshBase
 * \brief base mesh class
 *
 *  @author Christophe Prud'homme
 *  @see
 */
class MeshBase
{
public:


    /** @name Typedefs
     */
    //@{

    /**
     * Tuple that contains
     *
     * -# the index of the face
     *
     * -# the processor id the face belongs to
     */
    typedef boost::tuple<size_type, size_type> face_processor_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * Default constructor
     */
    MeshBase( WorldComm const& worldComm = Environment::worldComm() );

    /**
     * copy constructor
     */
    MeshBase( MeshBase const& );

    /**
     * destructor. make it virtual for derived classes
     */
    virtual ~MeshBase();

    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     * copy operator
     */
    MeshBase& operator=( MeshBase const& m );

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return \p true if the mesh is ready for use, \p false
     * otherwise
     */
    bool isUpdatedForUse() const
    {
        return M_is_updated;
    }

    /**
     * \return the number of elements
     */
    virtual size_type numElements() const = 0;

    /**
     * \return the number of faces
     */
    virtual size_type numFaces() const = 0;

    /**
     * \return the number of Points
     */
    virtual size_type numPoints() const = 0;

    /**
     * \return the number of vertices
     */
    size_type numVertices() const
    {
        return M_n_vertices;
    }

    /**
     * Returns the number of partitions.
     */
    uint16_type numberOfPartitions() const
    {
        return M_n_parts;
    }

    /**
     * \return \c true if mesh is partitioned, \c false otherwise
     */
    bool isPartitioned() const;

    /**
     * \return an integer(stored  in a \p Context) that encodes the components to be updated by
     * the mesh data structure.
     * \sa Context, MeshComponents
     */
    Context const& components() const
    {
        return M_components;
    }

    /**
     * \return an integer(stored  in a \p Context) that encodes the components to be updated by
     * the mesh data structure.
     * \sa Context, MeshComponents
     */
    Context&       components()
    {
        return M_components;
    }

    /**
     * \return the measure of the mesh
     */
    virtual double measure() const = 0;

    /**
     * \return true if the mesh has parametric nodes
     */
    bool isParametric() const
    {
        return M_is_parametric;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the number of partitions
     */
    void setNumberOfPartitions( uint16_type n )
    {
        M_n_parts = n;
    }

    /**
     * set the number of vertices
     */
    void setNumVertices( size_type n )
    {
        M_n_vertices = n ;
    }

    /**
     * set the components to be updated by \c updateForUse()
     * \sa updateForUse
     */
    void setComponents( size_type components = MESH_ALL_COMPONENTS )
    {
        M_components = components;
    }

    /**
     * set if the mesh is parametric ( e.g. has parametric nodes )
     */
    void setParametric( bool x )
    {
        M_is_parametric = x;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Empty all containers of the mesh. Must be redefined by derived
     * classes.
     */
    virtual void clear();


    /**
     * update all info for this mesh.
     */
    virtual void updateForUse() = 0;

    /**
     * update all info for this mesh according the \p components
     * \sa components(), setComponents()
     */
    virtual void updateForUse( size_type components );


    /**
     * Call the default partitioner (currently \p metis_partition()).
     */
    virtual void partition ( const uint16_type n_parts ) = 0;

    /**
     * \return the world comm
     */
    WorldComm const& worldComm() const
    {
        return M_worldComm;
    }

    virtual void setWorldComm( WorldComm const& _worldComm ) = 0;

    void setWorldCommMeshBase( WorldComm const& _worldComm )
    {
        M_worldComm = _worldComm;
    }

    mpi::communicator const& comm() const
    {
        return M_worldComm.localComm();
    }

    //@}



protected:

    /**
     * set to the flag whether the mesh is updated for proper use
     */
    void setUpdatedForUse( bool u )
    {
        M_is_updated = u;
    }

    /**
     * After loading/defining a mesh, we want to have as much locality
     * as possible (elements/faces/nodes to be contiguous). In order
     * to do that the mesh elements/faces/nodes are renumbered. That
     * will be then most helpful when generating the \p Dof table.
     * This procedure should work also with
     * \p comm().size() == 1
     *
     */
    virtual void renumber() = 0;

    /**
     * update the entities of co-dimension 1
     */
    virtual void updateEntitiesCoDimensionOne() = 0;

    /**
     * update the entities of co-dimension 2
     */
    virtual void updateEntitiesCoDimensionTwo() = 0;

    /**
     * check mesh connectivity
     */
    virtual void check() const = 0;

    /**
     * check elements orientation and fix it if needed
     */
    virtual void checkAndFixPermutation() = 0;

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & M_components;
            ar & M_is_updated;
            ar & M_is_parametric;
            ar & M_n_vertices;
            ar & M_n_parts;
        }
private:

    /**
     * encodes components that should be updated
     */
    Context M_components;

    /**
     * \p true if mesh ready to be used, \p false otherwise
     */
    bool M_is_updated;

    /**
     * \p true if the mesh is parametric (e.g. has parametric nodes), \p false otherwise
     */
    bool M_is_parametric;

    /**
     * number of vertices
     */
    size_type M_n_vertices;

    /**
     * The number of partitions the mesh has.  This is set by
     * the partitioners, and may not be changed directly by
     * the user.
     * \note The number of partitions *need not* equal
     * M_comm.size(), consider for example the case
     * where you simply want to partition a mesh on one
     * processor and view the result in GMV.
     */
    uint16_type M_n_parts;

    WorldComm M_worldComm;
};
}
#endif /* __MeshBase_H */
