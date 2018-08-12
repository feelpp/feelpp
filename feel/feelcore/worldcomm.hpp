/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2012-02-14

  @copyright (C) 2012 Université Joseph Fourier (Grenoble I)
  @copyright (C) 2012-2017 Feel++ Consortium
  
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
//!
//! @file
//! @author Vincent Chabannes <vincent.chabannes@imag.fr> 
//!
#ifndef FEELPP_WORLDCOMM_HPP
#define FEELPP_WORLDCOMM_HPP 1

#include <boost/mpi.hpp>
#if defined(FEELPP_HAS_MPI_H)
#include <mpi.h>
#endif /* FEELPP_HAS_MPI_H */
#include <boost/smart_ptr/enable_shared_from_this.hpp>

#include <feel/feelcore/feel.hpp>

namespace Feel
{

//!
//! @ingroup Core
//! @brief Provides interface to MPI Communicators
//!
//! @todo define global, local and god WorldComm
//! @todo define sub WorldComm
//! 
//! @code
//! auto world = Environment::worldComm();
//! std::vector<int> mapColorWorld;
//! for( int rank : irange( world.size() ) )
//!   mapColorWorld.push_back( rank % 2 );
//! // Creates a WorldComm where we have split the processes into two sets
//! WorldComm worldColored( mapColorWorld );
//! @endcode
//!
class FEELPP_EXPORT WorldComm : public boost::mpi::communicator, public std::enable_shared_from_this<WorldComm>
{

    typedef boost::mpi::communicator super;

public:

    //! self type
    typedef WorldComm self_type;
    //! `std::shared_ptr` type
    typedef std::shared_ptr<WorldComm> self_ptrtype;
    //! underlying MPI communicator type
    typedef boost::mpi::communicator communicator_type;

    using worldcomm_t = WorldComm;
    using worldcomm_ptr_t = std::shared_ptr<WorldComm>;
    using worldscomm_ptr_t = std::vector<worldcomm_ptr_t>;

    //! Default constructor
    WorldComm();

    //! constructor from an `boost::mpi::communicator`
    WorldComm( super const& );

    //! build a sub WorldComm from a \c color
    WorldComm( int _color );

    //! Creates a sub WorldComm from a set of colors
    WorldComm( std::vector<int> const& _colorWorld );
    
    //! 
    WorldComm( std::vector<int> const& _colorWorld, int localRank,
               communicator_type const& _globalComm=communicator_type(),
               communicator_type const& _godComm=communicator_type()  );

    WorldComm( WorldComm const& _wc );

    WorldComm( communicator_type const& _globalComm,
               int _color,
               bool _isActive );

    WorldComm( int _colorLocal,int localRank,int _colorGlobal,int globalRank,
               communicator_type const& _godComm,
               bool _isActive,
               bool _doInitActiveMap=true );

    // A supp
    WorldComm( communicator_type const& _globalComm,
               communicator_type const& _godComm,
               int _localColor, int localRank,
               std::vector<int> const& isActive );
    WorldComm( communicator_type const& _globalComm,
               communicator_type const& _localComm,
               communicator_type const& _godComm,
               int _localColor,// int localRank,
               std::vector<int> const& isActive );

    //! copy a worldcomm
    WorldComm& operator=( WorldComm const& wc );

    //! comparison of worldcomm
    bool operator==( WorldComm const& wc ) const;

    static self_ptrtype New() { return self_ptrtype(new self_type); }
    static self_ptrtype New( super const& s ) { return self_ptrtype(new self_type( s )); }

    worldcomm_ptr_t clone() const { return std::make_shared<worldcomm_t>( *this ); }

    void init( int color = 0, bool colormap = false );

    //! Returns the current WorldComm
    communicator_type const& globalComm() const
    {
        return *this;
    }

    //! Returns the local WorldComm
    communicator_type const& localComm() const
    {
        return M_localComm;
    }
    communicator_type const& godComm() const
    {
        return M_godComm;
    }

    //! Returns the local WorldComm, @see localComm
    communicator_type const& comm() const
    {
        return M_localComm;
    }
    communicator_type const& selfComm() const
    {
        return this->subWorldCommSeq();
    }

    rank_type globalSize() const
    {
        return this->globalComm().size();
    }
    rank_type localSize() const
    {
        return this->localComm().size();
    }
    rank_type godSize() const
    {
        return this->godComm().size();
    }

    //! Returns the rank in the global WorldComm
    rank_type globalRank() const
    {
        return this->globalComm().rank();
    }

    //! Returns the rank in the local WorldComm
    rank_type localRank() const
    {
        return this->localComm().rank();
    }

    //! Returns the rank in the God WorldComm
    rank_type godRank() const
    {
        return this->godComm().rank();
    }

    bool hasSubWorlds( int n ) const;
    worldscomm_ptr_t & subWorlds( int n );
    worldscomm_ptr_t const& subWorlds( int n ) const;
    worldscomm_ptr_t & subWorldsGroupBySubspace( int n );
    worldcomm_ptr_t subWorld( int n );
    worldcomm_ptr_t subWorld( int n ) const;
    int subWorldId( int n ) const;

    std::vector<int> const& mapColorWorld() const
    {
        return M_mapColorWorld;
    }
    std::vector<rank_type> const& mapLocalRankToGlobalRank() const
    {
        return M_mapLocalRankToGlobalRank;
    }
    std::vector<rank_type> const& mapGlobalRankToGodRank() const
    {
        return M_mapGlobalRankToGodRank;
    }

    int mapColorWorld(int k) const
    {
        return M_mapColorWorld[k];
    }
    rank_type mapLocalRankToGlobalRank(int k) const
    {
        return M_mapLocalRankToGlobalRank[k];
    }
    rank_type mapGlobalRankToGodRank(int k) const
    {
        return M_mapGlobalRankToGodRank[k];
    }

    //! Returns the master rank
    rank_type masterRank() const
    {
        return M_masterRank;
    }

    //! Returns \c true if process has master rank, \c false otherwise
    bool isMasterRank() const
    {
        return ( this->globalRank() ==  this->masterRank() );
    }


    worldcomm_t & subWorldComm();
    worldcomm_t const & subWorldComm() const;
    worldcomm_ptr_t subWorldComm( int color ) const;
    worldcomm_t & subWorldComm( std::vector<int> const& colormap ) ;
    worldcomm_ptr_t  subWorldComm( int color, std::vector<int> const& colormap ) ;
    worldcomm_t & masterWorld( int n );
    int numberOfSubWorlds() const;

    //! Returns the sequential sub worldcomm
    WorldComm & subWorldCommSeq();
    //! Returns the sequential sub worldcomm
    worldcomm_ptr_t & subWorldCommSeqPtr();
    
    //! Returns the sequential sub worldcomm
    WorldComm const& subWorldCommSeq() const;
    //! Returns the sequential sub worldcomm
    worldcomm_ptr_t const& subWorldCommSeqPtr() const;

    //! Returns \c true if current worldcomm is active in God worldcomm
    bool isActive() const
    {
        return this->isActive(this->godRank());
    }

    //! Returns \c true if worldcomm \p rank is active in God worldcomm
    bool isActive( int rank ) const
    {
        return M_isActive[rank];
    }

    //! Returns the vector of activity in this world
    std::vector<int> const& activityOnWorld() const
    {
        return M_isActive;
    }

    rank_type localColorToGlobalRank( int _color,int _localRank ) const;

    void setColorMap( std::vector<int> const& colormap );

    //! prints information on \p out on the WorldComm
    void showMe( std::ostream& __out = std::cout ) const;

    WorldComm operator+( WorldComm const & _worldComm ) const;

    void setIsActive( std::vector<int> const& _isActive ) const { M_isActive=_isActive; }

    void upMasterRank();

    void applyActivityOnlyOn(int _localColor) const;

    boost::tuple<bool,std::set<int> > hasMultiLocalActivity() const;

    /**
     * register sub worlds associated to \p worldmap
     */
    void registerSubWorlds( int n ) const;
    void registerSubWorldsGroupBySubspace( int n );

private :

    FEELPP_NO_EXPORT void initSubWorldCommSeq();

private :

    communicator_type M_localComm;
    communicator_type M_godComm;
    std::shared_ptr<WorldComm> M_subWorldCommSeq;

    std::vector<int> M_mapColorWorld;
    std::vector<rank_type> M_mapLocalRankToGlobalRank;
    std::vector<rank_type> M_mapGlobalRankToGodRank;
    mutable std::map<int, std::pair<worldcomm_ptr_t,worldscomm_ptr_t > > M_subworlds;

    int M_masterRank;
    mutable std::vector<int/*bool*/> M_isActive;

};

using worldcomm_t = WorldComm;
using worldcomm_ptr_t = std::shared_ptr<WorldComm>;
using worldscomm_ptr_t = std::vector<worldcomm_ptr_t>;


inline worldscomm_ptr_t makeWorldsComm( int n, WorldComm & wc )
{
    worldscomm_ptr_t r( n );
    for( auto & w : r )
        w = wc.clone();
    return r;
}
inline worldscomm_ptr_t makeWorldsComm( int n, worldcomm_ptr_t const& wc )
{
    worldscomm_ptr_t r( n );
    for( auto & w : r )
        w = wc->clone();
    return r;
}

} //namespace Feel

#endif // __worldcomm_H
