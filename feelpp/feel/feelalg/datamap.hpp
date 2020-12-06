/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
#include <feel/feelcore/commobject.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>

namespace Feel
{
template<typename SizeT>
class DataMap;

/**
 *
 */
class FEELPP_EXPORT IndexSplit : public std::vector<std::vector<size_type> >
{
    typedef IndexSplit self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;
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
        M_nIndexForSmallerRankId(0),
        M_tag(0)
    {}

    IndexSplit( int s )
        :
        super_type( s ),
        M_firstIndex( s, invalid_v<size_type> ),
        M_lastIndex( s, invalid_v<size_type> ),
        M_nIndex( s, invalid_v<size_type> ),
        M_nIndexForSmallerRankId( s, invalid_v<size_type> ),
        M_tag(s,0)
    {}

    IndexSplit( IndexSplit const& is )
        :
        super_type( is ),
        M_firstIndex( is.M_firstIndex ),
        M_lastIndex( is.M_lastIndex ),
        M_nIndex( is.M_nIndex ),
        M_nIndexForSmallerRankId( is.M_nIndexForSmallerRankId ),
        M_tag( is.M_tag )
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
    int tag( int i ) const { return M_tag[i]; }
    void setTag( int i,int tag ) { M_tag[i]=tag; }

    struct FieldsDef : public std::map<int,std::set<int> >
    {
        typedef std::map<int,std::set<int> > super_type;
        FieldsDef() : super_type() {}
        FieldsDef( super_type const& s ) : super_type( s ) {}
        void showMe() const;
    };

    static FieldsDef parseFieldsDef( std::string s );

    self_ptrtype applyFieldsDef( FieldsDef const& fieldsDef ) const;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(super_type);
        ar & BOOST_SERIALIZATION_NVP(M_firstIndex);
        ar & BOOST_SERIALIZATION_NVP(M_lastIndex);
        ar & BOOST_SERIALIZATION_NVP(M_nIndex);
        ar & BOOST_SERIALIZATION_NVP(M_nIndexForSmallerRankId);
        ar & BOOST_SERIALIZATION_NVP(M_tag);
    }

private :

    std::vector<size_type> M_firstIndex, M_lastIndex, M_nIndex;
    std::vector<size_type> M_nIndexForSmallerRankId;
    std::vector<int> M_tag;

};

/**
 * \class DataMap
 *  \brief data layout in a multi-processor environnement
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename SizeT=uint32_type>
class FEELPP_EXPORT DataMap : public CommObject
{

public:
    using super = CommObject;
    using size_type=SizeT;
    typedef IndexSplit indexsplit_type;
    typedef std::shared_ptr<indexsplit_type> indexsplit_ptrtype;

    typedef Eigen::Matrix<int, Eigen::Dynamic, 1>  localglobal_indices_type;
    typedef std::vector<localglobal_indices_type,Eigen::aligned_allocator<localglobal_indices_type> > vector_indices_type;

    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit DataMap( worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr() );

    /**
     * \param n total size of the vector
     * \param n_local local size of the vector on the curent processor
     */
    DataMap( size_type n, size_type n_local, worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr() );

    /**
     * \param n total size of the vector
     * \param firstdof array of size n_processors containing the first index on each processor
     * \param lastdof array of size n_processors containing the last index on each processor
     */
    DataMap( size_type n, std::vector<int> const& firstdof, std::vector<int> const& lastdof );

    /**
     *
     */
    DataMap( std::vector<std::shared_ptr<DataMap> > const& listofdm, worldcomm_ptr_t const& _worldComm );

    DataMap( DataMap const & dm ) = default;
    DataMap( DataMap&& dm ) = default;
    ~DataMap() override;

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
     * @return the number of local ghosts
     */
    size_type nLocalGhosts( const rank_type pid = invalid_rank_type_value ) const
    {
        rank_type thepid = ( pid == invalid_rank_type_value )? this->worldComm().rank() : pid;
        return nLocalDofWithGhost( thepid )-nLocalDofWithoutGhost( thepid );
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
        return this->nDofOnProcessor ( this->worldComm().rank() );
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
        return this->worldComm().size();
    }

    /**
     * @return the first dof index that is in  local subdomain
     */
    size_type firstDof() const
    {
        return this->firstDof( this->worldComm().rank() );
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
        return this->firstDofGlobalCluster( this->worldComm().rank() );
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
        return this->lastDof( this->worldComm().rank() );
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
        return this->lastDofGlobalCluster( this->worldComm().rank() );
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

    //! Returns local ID of global ID, return invalid_v<size_type> if not found on this processor.
    size_type  lid( size_type GID ) const
    {
        uint16_type pid = this->worldComm().rank();

        if ( GID >= firstDof( pid ) &&
                GID <= lastDof( pid ) )
            return GID - firstDof( pid );

        return invalid_v<size_type>;
    }

    //! Returns global ID of local ID, return -1 if not found on this processor.
    size_type gid( size_type LID ) const
    {
        uint16_type pid = this->worldComm().rank();

        if ( LID < ( lastDof( pid )-firstDof( pid ) + 1 ) )
            return firstDof( pid ) + LID;

        return invalid_v<size_type>;
    }

    //! Returns true if the GID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  myGID( size_type GID ) const
    {
        return( lid( GID )!=invalid_v<size_type> );
    }

    //! Returns true if the LID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  myLID( size_type LID ) const
    {
        return( gid( LID )!=invalid_v<size_type> );
    }

    //!Returns the minimum global ID across the entire map.
    size_type  minAllGID() const
    {
        return( firstDof( 0 ) );
    }

    //! Returns the maximum global ID across the entire map.
    size_type  maxAllGID() const
    {
        return( lastDof( this->worldComm().size()-1 ) );
    }

    //! Returns the maximum global ID owned by this processor.
    size_type  minMyGID() const
    {
        return firstDof( this->worldComm().rank() );
    }

    //! Returns the maximum global ID owned by this processor.
    size_type  maxMyGID() const
    {
        return lastDof( this->worldComm().rank() );
    };

    //!  The minimum local index value on the calling processor.
    size_type  minLID() const
    {
        return 0;
    };

    //! The maximum local index value on the calling processor.
    size_type  maxLID() const
    {
        return lastDof( this->worldComm().rank() )-firstDof( this->worldComm().rank() );
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

    //! processor numbering to world numbering
    std::vector<size_type> const& mapGlobalProcessToGlobalCluster() const
    {
        return M_mapGlobalProcessToGlobalCluster;
    }

    //! processor numbering to world numbering
    size_type mapGlobalProcessToGlobalCluster( size_type i ) const
    {
        return M_mapGlobalProcessToGlobalCluster[i];
    }

    FEELPP_DEPRECATED
    size_type mapGlobalClusterToGlobalProcess( size_type i ) const
    {
        return i;
    }

    void setNDof( size_type ndof );

    void setNLocalDofWithoutGhost( const rank_type proc, const size_type n, bool inWorld=true );
    void setNLocalDofWithGhost( const rank_type proc, const size_type n, bool inWorld=true );
    void setFirstDof( const rank_type proc, const size_type df, bool inWorld=true );
    void setLastDof( const rank_type proc, const size_type df, bool inWorld=true );
    void setFirstDofGlobalCluster( const rank_type proc, const size_type df, bool inWorld=true );
    void setLastDofGlobalCluster( const rank_type proc, const size_type df, bool inWorld=true );

    void setMapGlobalProcessToGlobalCluster( std::vector<size_type> const& map );
    void setMapGlobalProcessToGlobalCluster( size_type i, size_type j );
    void resizeMapGlobalProcessToGlobalCluster( size_type n );

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
    void setIsClosed( bool b )
    {
        M_closed = b;
    }

    /**
     * print some information
     * showAll : print more (can be large)
     */
    void showMe( bool showAll=false, std::ostream& __out = std::cout ) const;

    

    /**
     * \return the number of mapping (from functionspace id to container id with global process numbering)
     */
    int numberOfDofIdToContainerId() const { return M_dofIdToContainerId.size(); }
    /**
     * \return the mapping index which contains this global process id in container view
     */
    size_type databaseIndexFromContainerId( size_type gpdof ) const
        {
#if 0
            size_type currentStartId = 0;
            for ( int tag=0;tag<this->numberOfDofIdToContainerId();++tag )
            {
                size_type nGpDof = this->dofIdToContainerId( tag ).size();
                if ( gpdof >= currentStartId && gpdof < ( currentStartId + nGpDof ) )
                    return tag;
                currentStartId+=nGpDof;
            }
            return 0;
#endif
            for ( int tag=0;tag<this->numberOfDofIdToContainerId();++tag )
            {
                auto it = std::find( M_dofIdToContainerId[tag].begin(), M_dofIdToContainerId[tag].end(), gpdof );
                if ( it !=M_dofIdToContainerId[tag].end() )
                    return tag;
            }
            return invalid_v<size_type>;
        }
    /**
     * initialize the number of dofIdToContainerId mapping
     */
    void initNumberOfDofIdToContainerId( int nTag ) { M_dofIdToContainerId.resize(nTag); }
    /**
     * initialize the number of dof id in the mapping
     */
    void initDofIdToContainerId( int tag, int nDof ) { M_dofIdToContainerId[tag].resize(nDof,invalid_v<size_type>); }
    /**
     * initialize a mapping as identity
     */
    void initDofIdToContainerIdIdentity( int tag, int nDof )
        {
            M_dofIdToContainerId[tag].resize(nDof);
            std::iota( M_dofIdToContainerId[tag].begin(),M_dofIdToContainerId[tag].end(),0 );
        }
    /**
     * \return a reference of dofIdToContainerId mapping (from functionspace id to container id with global process numbering)
     */
    std::vector<size_type>& dofIdToContainerIdRef( int tag ) { return M_dofIdToContainerId[tag]; }
    void setDofIdToContainerId( int tag, std::vector<size_type> const& vec ) { M_dofIdToContainerId[tag] = vec; }

    /**
     * \return the dofIdToContainerId mapping (from functionspace id to container id with global process numbering)
     */
    std::vector<size_type> const& dofIdToContainerId( int tag ) const { return M_dofIdToContainerId[tag]; }
    /**
     * \return global process index in container from dof id
     */
    size_type dofIdToContainerId( int tag,size_type gpdof ) const { return M_dofIdToContainerId[tag][gpdof]; }

    //! update global process dof id in container from global process dof id in a doftable
    void dofIdToContainerId( int tag, std::set<size_type> const& gpdofs, std::set<size_type> & gpcont ) const
    {
        for ( size_type gp : gpdofs )
            gpcont.insert( this->dofIdToContainerId( tag, gp ) );
    }

    size_type containerIdToDofId( int tag, size_type gpdof ) const
    {
        size_type dof_id = invalid_v<size_type>;
        auto it = std::find( M_dofIdToContainerId[tag].begin(), M_dofIdToContainerId[tag].end(), gpdof );
        if ( it !=M_dofIdToContainerId[tag].end() )
            dof_id = std::distance( M_dofIdToContainerId[tag].begin(), it );
        return dof_id;
    }
    /**
     * \return the indexsplit description
     */
    indexsplit_ptrtype const& indexSplit() const { return M_indexSplit; }
    void setIndexSplit( indexsplit_ptrtype const& is ) { M_indexSplit = is; }
    void buildIndexSplit();

    bool hasIndexSplitWithComponents() const { return (M_indexSplitWithComponents.get() != 0); }
    indexsplit_ptrtype const& indexSplitWithComponents() const { return (M_indexSplitWithComponents)? M_indexSplitWithComponents : this->indexSplit(); }
    void setIndexSplitWithComponents( indexsplit_ptrtype const& is ) { M_indexSplitWithComponents = is; }
    void buildIndexSplitWithComponents( uint16_type nComp );
    //@}



    /** @name  Mutators
     */
    //@{



    //@}

    /** @name  Methods
     */
    //@{

    /**
     * \return true if this object is compatible with datamap dm
     */
    bool isCompatible( DataMap const& dm ) const;

    void close() const;

    // add missing dof entries in // ( typically a ghost dof present in index set but not active dof associated )
    void updateIndexSetWithParallelMissingDof( std::set<size_type> & _indexSet ) const;
    std::vector<size_type> buildIndexSetWithParallelMissingDof( std::vector<size_type> const& _indexSet ) const;
    // get process ids of active dof index (in cluster view) used in the input index set
    std::map<size_type, std::set<rank_type> >
    activeDofClusterUsedByProc( std::set<size_type> const& dofGlobalProcessPresent ) const;


    // build sub data map from an index set
    std::shared_ptr<DataMap> createSubDataMap( std::vector<size_type> const& idExtract,
                                                 bool _checkAndFixInputRange=true ) const;

    //@}

    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & BOOST_SERIALIZATION_NVP( M_closed );
        ar & BOOST_SERIALIZATION_NVP( M_n_dofs );
        ar & BOOST_SERIALIZATION_NVP( M_n_localWithoutGhost_df );
        ar & BOOST_SERIALIZATION_NVP( M_n_localWithGhost_df );
        ar & BOOST_SERIALIZATION_NVP( M_first_df );
        ar & BOOST_SERIALIZATION_NVP( M_last_df );
        ar & BOOST_SERIALIZATION_NVP( M_first_df_globalcluster );
        ar & BOOST_SERIALIZATION_NVP( M_last_df_globalcluster );
        ar & BOOST_SERIALIZATION_NVP( M_myglobalelements );
        ar & BOOST_SERIALIZATION_NVP( M_mapGlobalProcessToGlobalCluster );
        ar & BOOST_SERIALIZATION_NVP( M_neighbor_processors );
        ar & BOOST_SERIALIZATION_NVP( M_activeDofSharedOnCluster );
        ar & BOOST_SERIALIZATION_NVP( M_indexSplit );
        ar & BOOST_SERIALIZATION_NVP( M_indexSplitWithComponents );
        ar & BOOST_SERIALIZATION_NVP( M_dofIdToContainerId );
    }

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
     *The processors who neighbor the current processor
     */
    std::set<rank_type> M_neighbor_processors;

    /**
     * get set of rank which use an active dof (dofid on process, not cluster)
     */
    std::map<size_type, std::set<rank_type> > M_activeDofSharedOnCluster;

    
    /**
     * Index split ( differentiate multiphysic )
     */
    indexsplit_ptrtype M_indexSplit, M_indexSplitWithComponents;

    //! dof global process id (space def) to container global process id (identity with non composite space)
    std::vector<std::vector<size_type> > M_dofIdToContainerId;
private:

};

template<typename SizeT=uint32_type>
using datamap_t = DataMap<SizeT>;

template<typename SizeT=uint32_type>
using datamap_ptrtype = std::shared_ptr<datamap_t<SizeT>>;
template<typename SizeT=uint32_type>
using datamap_ptr_t = datamap_ptrtype<SizeT>;

} // Feel
#endif /* __DataMap_H */
