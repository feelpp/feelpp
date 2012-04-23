/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2012-02-14

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file worldcomm.cpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2012-02-14
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/debug.hpp>
#include <feel/feelcore/worldcomm.hpp>

//#include <boost/mpi/communicator.hpp>

namespace Feel
{

//-------------------------------------------------------------------------------

WorldComm::WorldComm()
    :
    super(),
    M_localComm( super::split( 0 ) ),
    M_godComm(),
    M_mapColorWorld( this->globalSize() ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() ),
    M_isActive( this->godSize(),true )
{
    Debug() << "\n WorldComm : warning constructor empty!! on godRank " << this->godRank() << "\n";
    std::vector<int> globalRanks( this->globalSize() );
    mpi::all_gather( this->globalComm(),
                     this->globalRank(),
                     globalRanks );
    Debug() << "gather1\n";
    mpi::all_gather( this->globalComm(),
                     0,
                     M_mapColorWorld );
    Debug() << "gather2\n";
    mpi::all_gather( this->localComm(),
                     this->globalRank(),
                     M_mapLocalRankToGlobalRank );
    Debug() << "gather3\n";
    mpi::all_gather( this->globalComm(),
                     this->godRank(),
                     M_mapGlobalRankToGodRank );
    Debug() << "gather4\n";
    // choice : the smallest rank
    M_masterRank = *std::min_element( globalRanks.begin(),globalRanks.end() );
}

//-------------------------------------------------------------------------------

WorldComm::WorldComm( int color )
    :
    super(),
    M_localComm( super::split( color ) ),
    M_godComm(),
    M_mapColorWorld( this->globalSize() ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() ),
    M_isActive( this->godSize(),true )
{
    std::vector<int> globalRanks( this->globalSize() );
    mpi::all_gather( this->godComm(),
                     this->globalRank(),
                     globalRanks );

    mpi::all_gather( this->globalComm(),
                     color,
                     M_mapColorWorld );

    mpi::all_gather( this->localComm(),
                     this->globalRank(),
                     M_mapLocalRankToGlobalRank );

    mpi::all_gather( this->globalComm(),
                     this->godRank(),
                     M_mapGlobalRankToGodRank );

    // choice : the smallest rank
    M_masterRank = *std::min_element( globalRanks.begin(),globalRanks.end() );
}

//-------------------------------------------------------------------------------

WorldComm::WorldComm( std::vector<int> const& colorWorld )
    :
    super(),
    M_localComm( super::split( colorWorld[this->globalRank()] ) ),
    M_godComm(),
    M_mapColorWorld( colorWorld ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() ),
    M_isActive( this->godSize(),true )
{
    std::vector<int> globalRanks( this->globalSize() );
    mpi::all_gather( this->globalComm(),
                     this->globalRank(),
                     globalRanks );

    mpi::all_gather( this->localComm(),
                     this->globalRank(),
                     M_mapLocalRankToGlobalRank );

    mpi::all_gather( this->globalComm(),
                     this->godRank(),
                     M_mapGlobalRankToGodRank );

    // choice : the smallest rank
    M_masterRank = *std::min_element( globalRanks.begin(),globalRanks.end() );
}

//-------------------------------------------------------------------------------

WorldComm::WorldComm( WorldComm const& wc )
    :
    super( wc ),
    M_localComm( wc.M_localComm ),
    M_godComm( wc.M_godComm ),
    M_mapColorWorld( wc.M_mapColorWorld ),
    M_mapLocalRankToGlobalRank( wc.M_mapLocalRankToGlobalRank ),
    M_mapGlobalRankToGodRank( wc.M_mapGlobalRankToGodRank ),
    M_masterRank( wc.M_masterRank ),
    M_isActive( wc.M_isActive )
{}

//-------------------------------------------------------------------------------

WorldComm::WorldComm( communicator_type const& _globalComm, int _color, bool _isActive )
    :
    super( _globalComm ),
    M_localComm( super::split( _color ) ),
    M_godComm(),
    M_mapColorWorld( this->globalSize() ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() )
{
    //std::cout << "\n WorldComm COPY CONSTRUCTOR " << std::endl;

    mpi::all_gather( this->godComm(),
                     ( int )_isActive,
                     M_isActive );

    mpi::all_gather( this->globalComm(),
                     _color,
                     M_mapColorWorld );

    mpi::all_gather( this->localComm(),
                     this->globalRank(),
                     M_mapLocalRankToGlobalRank );

    mpi::all_gather( this->globalComm(),
                     this->godRank(),
                     M_mapGlobalRankToGodRank );

    // choice : the smallest rank
    M_masterRank = INT_MAX;
#if 0
    int indexLoc=0;

    //for (auto it=M_mapRealRank.begin(),en=M_mapRealRank.end();it!=en;++it)
    for ( auto it=this->mapLocalRankToGlobalRank().begin(),
            en=this->mapLocalRankToGlobalRank().end(); it!=en; ++it,++indexLoc )
    {
        if  ( _isActive )
        {
            if ( M_isActive[this->mapGlobalRankToGodRank()[*it]] )
            {
                if ( M_masterRank>*it ) M_masterRank=*it;
            }
        }

        else
        {
            if ( !M_isActive[this->mapGlobalRankToGodRank()[*it]] )
            {
                if ( M_masterRank>*it ) M_masterRank=*it;
            }

        }
    }

#else

    for ( int p=0; p<this->globalSize(); ++p )
    {
        if  ( _isActive )
        {
            if ( M_isActive[this->mapGlobalRankToGodRank()[p]] )
            {
                if ( M_masterRank>p ) M_masterRank=p;
            }
        }

        else
        {
            if ( !M_isActive[this->mapGlobalRankToGodRank()[p]] )
            {
                if ( M_masterRank>p ) M_masterRank=p;
            }

        }
    }

#endif

}

//-------------------------------------------------------------------------------

WorldComm::WorldComm( communicator_type const& _globalComm,
                      communicator_type const& _godComm,
                      int _color, std::vector<int> const& isActive )
    :
    super( _globalComm ),
    M_localComm( super::split( _color ) ),
    M_godComm(_godComm ),
    M_mapColorWorld( this->globalSize() ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() ),
    M_isActive(isActive)
{

    mpi::all_gather( this->globalComm(),
                     _color,
                     M_mapColorWorld );

    mpi::all_gather( this->localComm(),
                     this->globalRank(),
                     M_mapLocalRankToGlobalRank );

    mpi::all_gather( this->globalComm(),
                     this->godRank(),
                     M_mapGlobalRankToGodRank );

    // choice : the smallest rank
    M_masterRank = INT_MAX;

    for ( int p=0; p<this->globalSize(); ++p )
    {
        if  ( this->isActive() )
        {
            if ( M_isActive[this->mapGlobalRankToGodRank()[p]] )
            {
                if ( M_masterRank>p ) M_masterRank=p;
            }
        }

        else
        {
            if ( !M_isActive[this->mapGlobalRankToGodRank()[p]] )
            {
                if ( M_masterRank>p ) M_masterRank=p;
            }

        }
    }


}


//-------------------------------------------------------------------------------

WorldComm::self_type
WorldComm::subWorldComm() const
{
    return this->subWorldComm(this->mapColorWorld()[this->globalRank()]);
}

//-------------------------------------------------------------------------------

WorldComm::self_type
WorldComm::subWorldComm( int _color ) const
{
#if 0
    bool isActive;
    int myColor = this->mapColorWorld()[this->globalRank()];

    if ( myColor==_color )
        isActive=true;

    else
        isActive=false;

    return self_type( this->localComm(), myColor, isActive );
#else
    bool isActive;
    int myColor = this->mapColorWorld()[this->globalRank()];

    if ( myColor==_color )
        isActive=true;

    else
        isActive=false;

    std::vector<int> newIsActive(this->godSize(),false);
    if (isActive)
        {
            for ( int p=0;p<this->localSize();++p)
                {
                    newIsActive[ mapGlobalRankToGodRank()[mapLocalRankToGlobalRank()[p]] ]=true;
                }
        }

    return self_type( this->localComm(),
                      this->godComm(),
                      myColor,
                      newIsActive );
#endif
}

//-------------------------------------------------------------------------------

WorldComm::self_type
WorldComm::subWorldCommSeq() const
{
    std::vector<int> newIsActive(this->godSize(),false);
    newIsActive[this->godRank()]=true;

    return self_type( communicator_type(MPI_COMM_SELF,boost::mpi::comm_attach),
                      this->godComm(),
                      this->globalRank(),
                      newIsActive/*this->M_isActive*/ );
}

//-------------------------------------------------------------------------------

void
WorldComm::showMe( std::ostream& __out ) const
{
    this->globalComm().barrier();

    for ( int proc = 0 ; proc < this->globalSize(); ++proc )
    {
        if ( this->globalRank()==proc && this->isActive() )
        {
            std::cout << "\n"
                      << "----------------------------------------------------------------------\n"
                      << "-------------WorldComm[proc "<<proc<<"]------------------------\n"
                      << " godrank " << this->godRank() << "\n"
                      << " globalrank " << this->globalRank() << "\n"
                      << " localrank " << this->localRank() << "\n"
                      << " godsize " << this->godSize() << "\n"
                      << " globalsize " << this->globalSize() << "\n"
                      << " localsize " << this->localSize() << "\n"
                      << " masterRank " << this->masterRank() << "\n"
                      << " isActive " << this->isActive() << "\n"
                      << "------------------------------------------------\n";

            for ( int k=0; k<( int )this->mapLocalRankToGlobalRank().size(); ++k )
                std::cout << k << " " << this->mapLocalRankToGlobalRank()[k] << std::endl;

            std::cout << "------------------------------------------------"<< std::endl;

            for ( int k=0; k<( int )this->mapColorWorld().size(); ++k )
                std::cout << k << " " << this->mapColorWorld()[k] << std::endl;

            std::cout << "------------------------------------------------"<< std::endl;

            for ( int k=0; k<( int )this->M_isActive.size(); ++k )
                std::cout << k << " " << this->M_isActive[k] << std::endl;

            std::cout << "------------------------------------------------"<< std::endl;

            std::cout << "\n----------------------------------------------------------------------\n";
            std::cout << std::endl;
        }

        this->globalComm().barrier();
    }

    this->godComm().barrier();
}

//-------------------------------------------------------------------------------

WorldComm::self_type
WorldComm::operator+( WorldComm const & _worldComm ) const
{
    int color;
    bool active;

    //int myColor = this->mapColorWorld()[this->globalRank()];
    //int otherColor = _worldComm.mapColorWorld()[this->globalRank()];
    if ( this->isActive() || _worldComm.isActive() )
    {
        color=1;
        active=true;
    }

    else
    {
        color=0;
        active=false;
    }

    //auto fusionComm = super::split(color);
    auto fusionComm = this->godComm().split( color );

    int colorOnProc;

    if ( this->isActive() ) colorOnProc=this->mapColorWorld()[this->globalRank()];

    else if ( _worldComm.isActive() ) colorOnProc=_worldComm.mapColorWorld()[_worldComm.globalRank()];

    else colorOnProc=INT_MAX;// maybe todo

    return WorldComm( fusionComm,colorOnProc,active );
}

//-------------------------------------------------------------------------------

int
WorldComm::localColorToGlobalRank( int _color,int _localRank ) const
{
    int res=0,cptLoc=0;
    bool find=false;

    while ( !find )
    {
        if ( this->mapColorWorld()[res]==_color )
            if ( cptLoc == _localRank ) find=true;

            else
            {
                ++cptLoc;
                ++res;
            }

        else ++res;
    }

    return res;
}

//-------------------------------------------------------------------------------


} //namespace Feel

