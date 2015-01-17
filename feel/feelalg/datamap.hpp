/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-28
 */
#ifndef __DataMap_H
#define __DataMap_H 1

#include <vector>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/worldcomm.hpp>

namespace Feel
{

class DataMap;



/**
 *
 */
class IndexSplit : public std::vector<std::vector<size_type> >
{
    typedef IndexSplit self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
    typedef std::vector<std::vector<size_type> > super_type;
    typedef super_type container_type;
    typedef std::vector<size_type> subcontainer_type;

    public:

    IndexSplit()
        :
        super_type(),
        M_firstIndex(0),
        M_lastIndex(0),
        M_nIndex(0),
        M_nIndexForSmallerRankId(0)
    {}

    IndexSplit( int s )
        :
        super_type( s ),
        M_firstIndex( s, invalid_size_type_value ),
        M_lastIndex( s, invalid_size_type_value ),
        M_nIndex( s, invalid_size_type_value ),
        M_nIndexForSmallerRankId( s, invalid_size_type_value )
    {}

    IndexSplit( IndexSplit const& is )
        :
        super_type( is ),
        M_firstIndex( is.M_firstIndex ),
        M_lastIndex( is.M_lastIndex ),
        M_nIndex( is.M_nIndex ),
        M_nIndexForSmallerRankId( is.M_nIndexForSmallerRankId )
    {}

    subcontainer_type const& split( int k ) const { return this->operator[](k); }

    void resize( int s );
    void addSplit( size_type startSplit, self_ptrtype const& addedIndexSplit /*DataMap const& dm*/ );

    void showMe() const;

    size_type firstIndex( int i ) const { return M_firstIndex[i]; }
    void setFirstIndex( int i,size_type s ) { M_firstIndex[i] = s; }
    size_type lastIndex( int i ) const { return M_lastIndex[i]; }
    void setLastIndex( int i,size_type s ) { M_lastIndex[i] = s; }
    size_type nIndex( int i ) const { return M_nIndex[i]; }
    void setNIndex( int i,size_type s ) { M_nIndex[i] = s; }
    size_type nIndexForSmallerRankId( int i ) const { return M_nIndexForSmallerRankId[i]; }
    void setNIndexForSmallerRankId( int i,size_type s ) { M_nIndexForSmallerRankId[i] = s; }

    struct FieldsDef : public std::map<int,std::set<int> >
    {
        typedef std::map<int,std::set<int> > super_type;
        FieldsDef() : super_type() {}
        FieldsDef( super_type const& s ) : super_type( s ) {}
        void showMe() const;
    };

    static FieldsDef parseFieldsDef( std::string s );

    self_ptrtype applyFieldsDef( FieldsDef const& fieldsDef ) const;


private :

    std::vector<size_type> M_firstIndex, M_lastIndex, M_nIndex;
    std::vector<size_type> M_nIndexForSmallerRankId;

};

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
    typedef IndexSplit indexsplit_type;
    typedef boost::shared_ptr<indexsplit_type> indexsplit_ptrtype;

    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    DataMap( WorldComm const& _worldComm = Environment::worldComm() );

    /**
     * \param n total size of the vector
     * \param n_local local size of the vector on the curent processor
     */
    DataMap( size_type n, size_type n_local, WorldComm const& _worldComm = Environment::worldComm() );

    /**
     * \param n total size of the vector
     * \param firstdof array of size n_processors containing the first index on each processor
     * \param lastdof array of size n_processors containing the last index on each processor
     */
    DataMap( size_type n, std::vector<int> const& firstdof, std::vector<int> const& lastdof );

    DataMap( DataMap const & dm ) = default;
    DataMap( DataMap&& dm ) = default;
    virtual ~DataMap();

    //@}

    /** @name Operator overloads
     */
    //@{

    DataMap& operator=( DataMap const& dm ) = default;
    DataMap& operator=( DataMap && dm ) = default;

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @return the total number of degrees of freedom in the problem.
     */
    size_type nDof() const
    {
        return M_n_dofs;
    }

    /**
     * @return the number of degrees of freedom on this processor(with ghosts).
     */
    size_type nLocalDof () const
    {
        return this->nLocalDofWithGhost()/*nDofOnProcessor (M_comm.rank())*/;
    }

    /**
     * @return the number of degrees of freedom on this processor without ghosts.
     */
    size_type nLocalDofWithoutGhost() const
    {
        return nLocalDofWithoutGhost(this->worldComm().rank());
    }

    /**
     * @return the number of degrees of freedom on this processor without ghosts.
     */
    size_type nLocalDofWithoutGhost( const rank_type proc ) const
    {
        return M_n_localWithoutGhost_df[proc];
    }

    /**
     * @return the number of degrees of freedom on this processor with ghosts.
     */
    size_type nLocalDofWithGhost() const
    {
        return this->nLocalDofWithGhost(this->worldComm().rank());
    }

    /**
     * @return the number of degrees of freedom on this processor with ghosts.
     */
    size_type nLocalDofWithGhost( const rank_type proc ) const
    {
        return M_n_localWithGhost_df[proc];
    }

    /**
     * @return the number of degrees of freedom on this processor.
     */
    size_type nMyDof () const
    {
        return this->nDofOnProcessor ( M_worldComm.rank() );
    }

    /**
     * @return the number of degrees of freedom on subdomain \p proc.
     */
    size_type nDofOnProcessor( const rank_type proc ) const
    {
        DCHECK( proc < M_first_df.size() ) << "invalid proc id or dof table , proc: "
                                           <<  proc << ", first dof : " <<  M_first_df.size();
        return M_n_localWithoutGhost_df[proc];
        //return ( M_last_df[proc] - M_first_df[proc]+1);
    }

    rank_type nProcessors() const
    {
        return M_worldComm.size();
    }

    /**
     * @return the first dof index that is in  local subdomain
     */
    size_type firstDof() const
    {
        return this->firstDof( M_worldComm.rank() );
    }
    /**
     * @return the first dof index that is local to subdomain \p proc.
     */
    size_type firstDof( const rank_type proc ) const
    {
        DCHECK( proc < M_first_df.size() ) << "invalid proc id or dof table , proc: "
                                           <<  proc << ", first dof : " <<  M_first_df.size();
        return M_first_df[proc];
    }

    size_type firstDofGlobalCluster() const
    {
        return this->firstDofGlobalCluster( M_worldComm.rank() );
    }

    size_type firstDofGlobalCluster( rank_type proc ) const
    {
        DCHECK( proc < M_first_df_globalcluster.size() ) << "invalid proc id or dof table , proc: "
                                                         <<  proc << ", first dof global cluster : " <<  M_first_df_globalcluster.size();
        return M_first_df_globalcluster[proc];
    }

    std::vector<size_type> const& firstDofGlobalClusterWorld() const
    {
        return M_first_df_globalcluster;
    }

    /**
     * Returns the last dof index that is in local  subdomain
     */
    size_type lastDof() const
    {
        return this->lastDof( M_worldComm.rank() );
    }
    /**
     * Returns the last dof index that is local to subdomain \p proc.
     */
    size_type lastDof( const rank_type proc ) const
    {
        DCHECK( proc < M_last_df.size() ) << "invalid proc id or dof table , proc: "
                                          <<  proc << ", last dof : " <<  M_last_df.size();
        return M_last_df[proc];
    }

    /**
     * Returns the last dof index that is in local  subdomain
     */
    size_type lastDofGlobalCluster() const
    {
        return this->lastDofGlobalCluster( M_worldComm.rank() );
    }

    size_type lastDofGlobalCluster( rank_type proc ) const
    {
        DCHECK( proc < M_last_df_globalcluster.size() ) << "invalid proc id or dof table , proc: "
                                                         <<  proc << ", last dof global cluster : " <<  M_last_df_globalcluster.size();
        return M_last_df_globalcluster[proc];
    }

    std::vector<size_type> const& lastDofGlobalClusterWorld() const
    {
        return M_last_df_globalcluster;
    }

    rank_type procOnGlobalCluster( size_type globDof ) const;

    bool dofGlobalClusterIsOnProc( size_type globDof ) const
    {
        return this->dofGlobalClusterIsOnProc( globDof, this->worldComm().globalRank() );
    }

    bool dofGlobalClusterIsOnProc( size_type globDof, rank_type proc ) const
    {
        return ( ( this->nLocalDofWithoutGhost(proc) > 0 ) && ( globDof <= M_last_df_globalcluster[proc] ) && ( globDof >= M_first_df_globalcluster[proc] ) );
    }

    bool dofGlobalProcessIsGhost( size_type dof) const
    {
        return !this->dofGlobalClusterIsOnProc(this->mapGlobalProcessToGlobalCluster( dof ));
    }

    //! return global process id from global cluster id if available
    boost::tuple<bool,size_type> searchGlobalProcessDof( size_type gpdof ) const;

    //! Returns local ID of global ID, return invalid_size_type_value if not found on this processor.
    size_type  lid( size_type GID ) const
    {
        uint16_type pid = M_worldComm.rank();

        if ( GID >= firstDof( pid ) &&
                GID <= lastDof( pid ) )
            return GID - firstDof( pid );

        return invalid_size_type_value;
    }

    //! Returns global ID of local ID, return -1 if not found on this processor.
    size_type gid( size_type LID ) const
    {
        uint16_type pid = M_worldComm.rank();

        if ( LID < ( lastDof( pid )-firstDof( pid ) + 1 ) )
            return firstDof( pid ) + LID;

        return invalid_size_type_value;
    }

    //! Returns true if the GID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  myGID( size_type GID ) const
    {
        return( lid( GID )!=invalid_size_type_value );
    }

    //! Returns true if the LID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  myLID( size_type LID ) const
    {
        return( gid( LID )!=invalid_size_type_value );
    }

    //!Returns the minimum global ID across the entire map.
    size_type  minAllGID() const
    {
        return( firstDof( 0 ) );
    }

    //! Returns the maximum global ID across the entire map.
    size_type  maxAllGID() const
    {
        return( lastDof( M_worldComm.size()-1 ) );
    }

    //! Returns the maximum global ID owned by this processor.
    size_type  minMyGID() const
    {
        return firstDof( M_worldComm.rank() );
    }

    //! Returns the maximum global ID owned by this processor.
    size_type  maxMyGID() const
    {
        return lastDof( M_worldComm.rank() );
    };

    //!  The minimum local index value on the calling processor.
    size_type  minLID() const
    {
        return 0;
    };

    //! The maximum local index value on the calling processor.
    size_type  maxLID() const
    {
        return lastDof( M_worldComm.rank() )-firstDof( M_worldComm.rank() );
    };

    //! number of elements across all processors.
    size_type nGlobalElements() const
    {
        return M_n_dofs;
    };

    //! number of elements on the calling processor.
    //size_type nMyElements() const {return nLocalDofWithGhost();};
    size_type nMyElements() const
    {
        return nLocalDof();
    };

    //! Puts list of global elements on this processor size_typeo the user-provided array.
    std::vector<size_type> const& myGlobalElements() const;

    std::vector<size_type> const& mapGlobalProcessToGlobalCluster() const
    {
        return M_mapGlobalProcessToGlobalCluster;
    }
    std::vector<size_type> const& mapGlobalClusterToGlobalProcess() const
    {
        return M_mapGlobalClusterToGlobalProcess;
    }
    size_type mapGlobalProcessToGlobalCluster( size_type i ) const
    {
        return M_mapGlobalProcessToGlobalCluster[i];
    }
    size_type mapGlobalClusterToGlobalProcess( size_type i ) const
    {
        return M_mapGlobalClusterToGlobalProcess[i];
    }

    void setNDof( size_type ndof );

    void setNLocalDofWithoutGhost( const rank_type proc, const size_type n, bool inWorld=true );
    void setNLocalDofWithGhost( const rank_type proc, const size_type n, bool inWorld=true );
    void setFirstDof( const rank_type proc, const size_type df, bool inWorld=true );
    void setLastDof( const rank_type proc, const size_type df, bool inWorld=true );
    void setFirstDofGlobalCluster( const rank_type proc, const size_type df, bool inWorld=true );
    void setLastDofGlobalCluster( const rank_type proc, const size_type df, bool inWorld=true );

    void setMapGlobalProcessToGlobalCluster( std::vector<size_type> const& map );
    void setMapGlobalClusterToGlobalProcess( std::vector<size_type> const& map );
    void setMapGlobalProcessToGlobalCluster( size_type i, size_type j );
    void setMapGlobalClusterToGlobalProcess( size_type i, size_type j );
    void resizeMapGlobalProcessToGlobalCluster( size_type n );
    void resizeMapGlobalClusterToGlobalProcess( size_type n );

    void updateDataInWorld();


    typename std::set<rank_type>::const_iterator beginNeighborSubdomains() const { return M_neighbor_processors.begin(); }
    typename std::set<rank_type>::const_iterator endNeighborSubdomains() const { return M_neighbor_processors.end(); }
    std::set<rank_type> const& neighborSubdomains() const { return M_neighbor_processors; }
    void setNeighborSubdomains( std::set<rank_type> const& neigh) { M_neighbor_processors=neigh; }
    void addNeighborSubdomain( rank_type p ) { M_neighbor_processors.insert( p ); }
    void addNeighborSubdomains( std::set<rank_type> const& neigh )
    {
        for ( rank_type procIdNeigh : neigh )
            M_neighbor_processors.insert( procIdNeigh );
    }

    std::map<size_type, std::set<rank_type> > const& activeDofSharedOnCluster() const { return M_activeDofSharedOnCluster; }
    void setActiveDofSharedOnCluster(size_type j, std::set<rank_type> const& sharedRank ) { M_activeDofSharedOnCluster[j]=sharedRank; }


    //! \return true if DataMap is close, false otherwise
    bool closed() const
    {
        return M_closed;
    }

    void showMeMapGlobalProcessToGlobalCluster( bool showAll=false, std::ostream& __out = std::cout ) const;

    /**
     * \return the communicator
     */
    //mpi::communicator const& comm() const { return M_comm; }
    WorldComm const& worldComm() const
    {
        return M_worldComm;
    }
    // warning (vincent!!)
    WorldComm const& comm() const
    {
        return M_worldComm;
    }



    indexsplit_ptrtype const& indexSplit() const { return M_indexSplit; }
    void setIndexSplit( indexsplit_ptrtype const& is ) { M_indexSplit = is; }

    void buildIndexSplit();
    //@}



    /** @name  Mutators
     */
    //@{



    //@}

    /** @name  Methods
     */
    //@{

    void close() const;

    // add missing dof entries in // ( typically a ghost dof present in index set but not active dof associated )
    void updateIndexSetWithParallelMissingDof( std::vector<size_type> & _indexSet ) const;
    std::vector<size_type> buildIndexSetWithParallelMissingDof( std::vector<size_type> const& _indexSet ) const;

    // build sub data map from an index set
    boost::shared_ptr<DataMap> createSubDataMap( std::vector<size_type> const& idExtract,
                                                 bool _checkAndFixInputRange=true ) const;

    //@}



protected:

    mutable bool M_closed;

    /**
     * Total number of degrees of freedom.
     */
    size_type M_n_dofs;

    /**
     * Number of degrees of freedom for each processor without ghosts.
     */
    std::vector<size_type> M_n_localWithoutGhost_df;

    /**
     * Number of degrees of freedom for each processor with ghosts.
     */
    std::vector<size_type> M_n_localWithGhost_df;

    /**
     * First DOF index on processor \p p.
     */
    std::vector<size_type> M_first_df;

    /**
     * Last DOF index on processor \p p.
     */
    std::vector<size_type> M_last_df;

    /**
     * First globalcluster DOF index on processor \p p.
     */
    std::vector<size_type> M_first_df_globalcluster;

    /**
     * Last globalcluster DOF index on processor \p p.
     */
    std::vector<size_type> M_last_df_globalcluster;

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
     *The processors who neighbor the current processor
     */
    std::set<rank_type> M_neighbor_processors;

    /**
     * get set of rank which use an active dof (dofid on process, not cluster)
     */
    std::map<size_type, std::set<rank_type> > M_activeDofSharedOnCluster;

    /**
     * Communicator
     */
    WorldComm M_worldComm;

    /**
     * Index split ( differentiate multiphysic )
     */
    indexsplit_ptrtype M_indexSplit;


private:

};

} // Feel
#endif /* __DataMap_H */
