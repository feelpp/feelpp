/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-28

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file dofmap.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-02-28
 */
#ifndef __DataMap_H
#define __DataMap_H 1

#include <vector>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/worldcomm.hpp>

namespace Feel
{
/**
 * \class DataMap
 *  \brief data layout in a multi-processor environnement
 *
 *  @author Christophe Prud'homme
 *  @see
 */
class DataMap
{

public:


    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    DataMap(WorldComm const& _worldComm = WorldComm() );

    /**
     * \param n total size of the vector
     * \param n_local local size of the vector on the curent processor
     */
    DataMap( size_type n, size_type n_local, WorldComm const& _worldComm = WorldComm() );

    /**
     * \param n total size of the vector
     * \param firstdof array of size n_processors containing the first index on each processor
     * \param lastdof array of size n_processors containing the last index on each processor
     */
    DataMap( size_type n, std::vector<int> const& firstdof, std::vector<int> const& lastdof );

    DataMap( DataMap const & dm );

    virtual ~DataMap();

    //@}

    /** @name Operator overloads
     */
    //@{

    DataMap& operator=( DataMap const& dm );

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @return the total number of degrees of freedom in the problem.
     */
    size_type nDof() const { return _M_n_dofs; }

    /**
     * @return the number of degrees of freedom on this processor(with ghosts).
     */
    size_type nLocalDof () const
    { return this->nLocalDofWithGhost()/*nDofOnProcessor (M_comm.rank())*/; }

    /**
     * @return the number of degrees of freedom on this processor without ghosts.
     */
    size_type nLocalDofWithoutGhost () const
    { return _M_n_localWithoutGhost_df[this->worldComm().rank()]; }

    /**
     * @return the number of degrees of freedom on this processor with ghosts.
     */
    size_type nLocalDofWithGhost () const
    { return _M_n_localWithGhost_df[this->worldComm().rank()]; }

    /**
     * @return the number of degrees of freedom on this processor.
     */
    size_type nMyDof () const
    { return this->nDofOnProcessor (M_worldComm.rank()); }

    /**
     * @return the number of degrees of freedom on subdomain \p proc.
     */
    size_type nDofOnProcessor(const size_type proc) const
    {
        FEELPP_ASSERT(proc < _M_first_df.size())( proc )( _M_first_df.size() ).error( "invalid proc id or dof table" );
        return _M_n_localWithoutGhost_df[proc];
        //return ( _M_last_df[proc] - _M_first_df[proc]+1);
    }

    size_type nProcessors() const
    {
        return ( M_worldComm.size() );
    }

    /**
     * @return the first dof index that is in  local subdomain
     */
    size_type firstDof() const
    {
        size_type proc = M_worldComm.rank();
        FEELPP_ASSERT(proc < _M_first_df.size())( proc )( _M_first_df.size() ).error( "invalid proc id or dof table" );
        return _M_first_df[proc];
    }
    /**
     * @return the first dof index that is local to subdomain \p proc.
     */
    size_type firstDof(const size_type proc) const
    {
        FEELPP_ASSERT(proc < _M_first_df.size())( proc )( _M_first_df.size() ).error( "invalid proc id or dof table" );
        return _M_first_df[proc];
    }

    size_type firstDofGlobalCluster() const
    {
        size_type proc = M_worldComm.rank();
        FEELPP_ASSERT(proc < _M_first_df_globalcluster.size())( proc )( _M_first_df_globalcluster.size() ).error( "invalid proc id or dof table" );
        return _M_first_df_globalcluster[proc];
    }

    size_type firstDofGlobalCluster(uint16_type proc) const
    {
        FEELPP_ASSERT(proc < _M_first_df_globalcluster.size())( proc )( _M_first_df_globalcluster.size() ).error( "invalid proc id or dof table" );
        return _M_first_df_globalcluster[proc];
    }

    /**
     * Returns the last dof index that is in local  subdomain
     */
    size_type lastDof() const
    {
        size_type proc = M_worldComm.rank();
        FEELPP_ASSERT(proc < _M_last_df.size())( proc )( _M_last_df.size() ).error( "invalid proc id or dof table" );
        return _M_last_df[proc];
    }
    /**
     * Returns the last dof index that is local to subdomain \p proc.
     */
    size_type lastDof(const unsigned int proc) const
    {
        FEELPP_ASSERT(proc < _M_last_df.size())( proc )( _M_last_df.size() ).error( "invalid proc id or dof table" );
        return _M_last_df[proc];
    }

    /**
     * Returns the last dof index that is in local  subdomain
     */
    size_type lastDofGlobalCluster() const
    {
        size_type proc = M_worldComm.rank();
        FEELPP_ASSERT(proc < _M_last_df_globalcluster.size())( proc )( _M_last_df_globalcluster.size() ).error( "invalid proc id or dof table" );
        return _M_last_df_globalcluster[proc];
    }

    size_type lastDofGlobalCluster(uint16_type proc) const
    {
        FEELPP_ASSERT(proc < _M_last_df_globalcluster.size())( proc )( _M_last_df_globalcluster.size() ).error( "invalid proc id or dof table" );
        return _M_last_df_globalcluster[proc];
    }

    uint16_type procOnGlobalCluster( size_type globDof ) const
    {
        uint16_type proc=0,res=0;
        bool find=false;
        while ( (!find) && (proc<this->nProcessors()) )
            {
                if ( (globDof <= _M_last_df_globalcluster[proc] ) && (globDof >= _M_first_df_globalcluster[proc] ) )
                    {
                        res = proc;
                        find=true;
                    }
                else
                    ++proc;
            }
        return res;
    }

    bool dofGlobalClusterIsOnProc( size_type globDof ) const
    {
        return dofGlobalClusterIsOnProc(globDof, this->worldComm().globalRank());
    }

    bool dofGlobalClusterIsOnProc( size_type globDof, int proc ) const
    {
        return ((globDof <= _M_last_df_globalcluster[proc] ) && (globDof >= _M_first_df_globalcluster[proc] ) );
    }

    //! Returns local ID of global ID, return invalid_size_type_value if not found on this processor.
    size_type  lid(size_type GID) const
    {
        uint16_type pid = M_worldComm.rank();
        if ( GID >= firstDof( pid ) &&
             GID <= lastDof( pid ) )
            return GID - firstDof( pid );
        return invalid_size_type_value;
    }

    //! Returns global ID of local ID, return -1 if not found on this processor.
    size_type gid( size_type LID) const
    {
        uint16_type pid = M_worldComm.rank();
        if ( LID < ( lastDof( pid )-firstDof( pid ) + 1 ) )
            return firstDof( pid ) + LID;
        return invalid_size_type_value;
    }

    //! Returns true if the GID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  myGID(size_type GID) const {return(lid(GID)!=invalid_size_type_value);}

    //! Returns true if the LID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  myLID(size_type LID) const {return(gid(LID)!=invalid_size_type_value);}

    //!Returns the minimum global ID across the entire map.
    size_type  minAllGID() const {return(firstDof( 0 ));}

    //! Returns the maximum global ID across the entire map.
    size_type  maxAllGID() const {return(lastDof( M_worldComm.size()-1 ) );}

    //! Returns the maximum global ID owned by this processor.
    size_type  minMyGID() const {return firstDof( M_worldComm.rank() );}

    //! Returns the maximum global ID owned by this processor.
    size_type  maxMyGID() const {return lastDof( M_worldComm.rank() );};

    //!  The minimum local index value on the calling processor.
    size_type  minLID() const {return 0;};

    //! The maximum local index value on the calling processor.
    size_type  maxLID() const {return lastDof( M_worldComm.rank() )-firstDof( M_worldComm.rank() );};

    //! number of elements across all processors.
    size_type nGlobalElements() const {return _M_n_dofs;};

    //! number of elements on the calling processor.
    //size_type nMyElements() const {return nLocalDofWithGhost();};
    size_type nMyElements() const {return nLocalDof();};

    //! Puts list of global elements on this processor size_typeo the user-provided array.
    std::vector<size_type> const& myGlobalElements() const;

    std::vector<size_type> const& mapGlobalProcessToGlobalCluster() const { return M_mapGlobalProcessToGlobalCluster; }
    std::vector<size_type> const& mapGlobalClusterToGlobalProcess() const { return M_mapGlobalClusterToGlobalProcess; }
    size_type mapGlobalProcessToGlobalCluster(size_type i) const { return M_mapGlobalProcessToGlobalCluster[i]; }
    size_type mapGlobalClusterToGlobalProcess(size_type i) const { return M_mapGlobalClusterToGlobalProcess[i]; }

    void setNDof(size_type ndof);

    void setNLocalDofWithoutGhost(const size_type proc, const size_type n, bool inWorld=true);
    void setNLocalDofWithGhost(const size_type proc, const size_type n, bool inWorld=true);
    void setFirstDof(const size_type proc, const size_type df, bool inWorld=true);
    void setLastDof(const size_type proc, const size_type df, bool inWorld=true);
    void setFirstDofGlobalCluster(const size_type proc, const size_type df, bool inWorld=true);
    void setLastDofGlobalCluster(const size_type proc, const size_type df, bool inWorld=true);

    void setMapGlobalProcessToGlobalCluster( std::vector<size_type> const& map);
    void setMapGlobalClusterToGlobalProcess( std::vector<size_type> const& map);
    void setMapGlobalProcessToGlobalCluster( size_type i, size_type j);
    void setMapGlobalClusterToGlobalProcess( size_type i, size_type j);
    void resizeMapGlobalProcessToGlobalCluster( size_type n);
    void resizeMapGlobalClusterToGlobalProcess( size_type n);

    void updateDataInWorld();

    //! \return true if DataMap is close, false otherwise
    bool closed() const { return M_closed; }

    void showMeMapGlobalProcessToGlobalCluster( std::ostream& __out = std::cout ) const;

    /**
     * \return the communicator
     */
    //mpi::communicator const& comm() const { return M_comm; }
    WorldComm const& worldComm() const { return M_worldComm; }
    // warning (vincent!!)
    WorldComm const& comm() const { return M_worldComm; }

    //@}



    /** @name  Mutators
     */
    //@{



    //@}

    /** @name  Methods
     */
    //@{

    void close() const;

    //@}



protected:

    mutable bool M_closed;

    /**
     * Total number of degrees of freedom.
     */
    size_type _M_n_dofs;

    /**
     * Number of degrees of freedom for each processor without ghosts.
     */
    std::vector<size_type> _M_n_localWithoutGhost_df;

    /**
     * Number of degrees of freedom for each processor with ghosts.
     */
    std::vector<size_type> _M_n_localWithGhost_df;

    /**
     * First DOF index on processor \p p.
     */
    std::vector<size_type> _M_first_df;

    /**
     * Last DOF index on processor \p p.
     */
    std::vector<size_type> _M_last_df;

    /**
     * First globalcluster DOF index on processor \p p.
     */
    std::vector<size_type> _M_first_df_globalcluster;

    /**
     * Last globalcluster DOF index on processor \p p.
     */
    std::vector<size_type> _M_last_df_globalcluster;

    // ??
    mutable std::vector<size_type> M_myglobalelements;


    /**
     * Map between Global Process To Global Cluster.
     */
    std::vector<size_type> M_mapGlobalProcessToGlobalCluster;

    /**
     * Map between Global Cluster To Global Process.
     */
    std::vector<size_type> M_mapGlobalClusterToGlobalProcess;

    /**
     * Communicator
     */
    //mpi::communicator M_comm;
    WorldComm M_worldComm;
private:

};

} // Feel
#endif /* __DataMap_H */
