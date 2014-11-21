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
#include <feel/feelcore/environment.hpp>

//#include <boost/mpi/communicator.hpp>

namespace Feel
{

//-------------------------------------------------------------------------------

WorldComm::WorldComm()
    :
    super(),
    M_localComm( this->globalComm() ),
    M_godComm( this->globalComm() ),
    M_subWorldCommSeq(),
    M_mapColorWorld( this->globalSize(), 0 ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() ),
    M_masterRank( 0 ),
    M_isActive( this->godSize(),true )
{
    for (int p=0 ; p<this->globalSize() ; ++p )
    {
        M_mapLocalRankToGlobalRank[p] = p;
        M_mapGlobalRankToGodRank[p] = p;
    }

    this->initSubWorldCommSeq();
}

WorldComm::WorldComm( super const& s )
    :
    super(),
    M_localComm( super::split( 0, this->globalRank() ) ),
    M_godComm(s),
    M_subWorldCommSeq(),
    M_mapColorWorld( this->globalSize() ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() ),
    M_isActive( this->godSize(),true )
{
    //LOG(INFO) << "\n WorldComm : warning constructor empty!! on godRank " << this->godRank() << "\n";
    init( 0, false );
}

//-------------------------------------------------------------------------------

WorldComm::WorldComm( int color )
    :
    super(),
    //M_selfComm( MPI_COMM_SELF, mpi::comm_attach ),
    M_localComm( super::split( color, this->globalRank() ) ),
    M_godComm(),
    M_subWorldCommSeq(),
    M_mapColorWorld( this->globalSize() ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() ),
    M_isActive( this->godSize(),true )
{
    init( color, false );
}

//-------------------------------------------------------------------------------

    WorldComm::WorldComm( std::vector<int> const& colorWorld )
    :
    super(),
    M_localComm( super::split( colorWorld[this->globalRank()] ) ),
    M_godComm(),
    M_subWorldCommSeq(),
    M_mapColorWorld( colorWorld ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() ),
    M_isActive( this->godSize(),true )
{
    init( 0, false );
}

//-------------------------------------------------------------------------------

WorldComm::WorldComm( std::vector<int> const& colorWorld,
                      int _localRank,
                      communicator_type const& _globalComm,
                      communicator_type const& _godComm )
    :
    super(_globalComm),
    M_localComm( super::split( colorWorld[this->globalRank()],_localRank ) ),
    M_godComm(_godComm),
    M_subWorldCommSeq(),
    M_mapColorWorld( colorWorld ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() ),
    M_isActive( this->godSize(),true )
{
    init( 0, false );
}

//-------------------------------------------------------------------------------

void
WorldComm::init( int color, bool colormap )
{
    mpi::all_gather( this->localComm(),
                     this->globalRank(),
                     M_mapLocalRankToGlobalRank );

    std::vector<int> globalRanks( this->globalSize() );
    if ( !colormap )
    {
        auto dataSendToGather = boost::make_tuple( this->globalRank(),this->godRank() );
        std::vector<boost::tuple<rank_type,rank_type> > dataRecvToGather( this->globalSize() );
        mpi::all_gather( this->globalComm(),
                         dataSendToGather,
                         dataRecvToGather);
        for (int p=0 ; p<this->globalSize() ; ++p )
        {
            globalRanks[p] = dataRecvToGather[p].get<0>();
            M_mapGlobalRankToGodRank[p] = dataRecvToGather[p].get<1>();
        }
    }
    else
    {
        auto dataSendToGather = boost::make_tuple( this->globalRank(),this->godRank(), color );
        std::vector<boost::tuple<rank_type,rank_type,int> > dataRecvToGather( this->globalSize() );

        mpi::all_gather( this->globalComm(),
                         dataSendToGather,
                         dataRecvToGather);
        for (int p=0 ; p<this->globalSize() ; ++p )
        {
            globalRanks[p] = dataRecvToGather[p].get<0>();
            M_mapGlobalRankToGodRank[p] = dataRecvToGather[p].get<1>();
            M_mapColorWorld[p] = dataRecvToGather[p].get<2>();
        }

    }

    // choice : the smallest rank
    M_masterRank = *std::min_element( globalRanks.begin(),globalRanks.end() );

    this->initSubWorldCommSeq();

}

//-------------------------------------------------------------------------------

WorldComm::WorldComm( WorldComm const& wc )
    :
    super( wc ),
    M_localComm( wc.M_localComm ),
    M_godComm( wc.M_godComm ),
    M_subWorldCommSeq( wc.M_subWorldCommSeq ),
    M_mapColorWorld( wc.M_mapColorWorld ),
    M_mapLocalRankToGlobalRank( wc.M_mapLocalRankToGlobalRank ),
    M_mapGlobalRankToGodRank( wc.M_mapGlobalRankToGodRank ),
    M_subworlds( wc.M_subworlds ),
    M_masterRank( wc.M_masterRank ),
    M_isActive( wc.M_isActive )
{}

//-------------------------------------------------------------------------------

WorldComm::WorldComm( communicator_type const& _globalComm, int _color, bool _isActive )
    :
    super( _globalComm ),
    M_localComm( super::split( _color ) ),
    M_godComm(),
    M_subWorldCommSeq(),
    M_mapColorWorld( this->globalSize() ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() )
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

    mpi::all_gather( this->godComm(),
                     ( int )_isActive,
                     M_isActive );


    this->upMasterRank();

    this->initSubWorldCommSeq();

}

//-------------------------------------------------------------------------------

    WorldComm::WorldComm( int _colorLocal, int localRank,
                          int _colorGlobal, int globalRank,
                          communicator_type const& _godComm, bool _isActive, bool _doInitActiveMap )
    :
    super( _godComm.split( _colorGlobal, globalRank ) ),
    M_localComm( super::split( _colorLocal, localRank ) ),
    M_godComm(_godComm),
    M_subWorldCommSeq(),
    M_mapColorWorld( this->globalSize() ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() )
{
    mpi::all_gather( this->globalComm(),
                     _colorLocal,
                     M_mapColorWorld );
    mpi::all_gather( this->localComm(),
                     this->globalRank(),
                     M_mapLocalRankToGlobalRank );
    mpi::all_gather( this->globalComm(),
                     this->godRank(),
                     M_mapGlobalRankToGodRank );
    M_masterRank = 0;

    if (_doInitActiveMap)
    {
        mpi::all_gather( this->godComm(),
                         ( int )_isActive,
                         M_isActive );

        this->upMasterRank();
    }

    this->initSubWorldCommSeq();
}


//-------------------------------------------------------------------------------

WorldComm::WorldComm( communicator_type const& _globalComm,
                      communicator_type const& _godComm,
                      int _color, int localRank,
                      std::vector<int> const& isActive )
    :
    super( _globalComm ),
    M_localComm( _globalComm.split( _color, localRank ) ),
    M_godComm(_godComm ),
    M_subWorldCommSeq(),
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


    this->upMasterRank();

    this->initSubWorldCommSeq();
}

//-------------------------------------------------------------------------------

WorldComm::WorldComm( communicator_type const& _globalComm,
                      communicator_type const& _localComm,
                      communicator_type const& _godComm,
                      int _color,
                      std::vector<int> const& isActive )
    :
    super( _globalComm ),
    M_localComm(_localComm ),
    M_godComm(_godComm ),
    M_subWorldCommSeq(),
    M_mapColorWorld( this->globalSize() ),
    M_mapLocalRankToGlobalRank( this->localSize() ),
    M_mapGlobalRankToGodRank( this->globalSize() ),
    M_isActive(isActive)
{
    mpi::all_gather( this->localComm(),
                     this->globalRank(),
                     M_mapLocalRankToGlobalRank );

    auto dataSendToGather = boost::make_tuple( _color,this->godRank() );
    std::vector<boost::tuple<int,rank_type> > dataRecvToGather( this->globalSize() );
    mpi::all_gather( this->globalComm(),
                     dataSendToGather,
                     dataRecvToGather);
    for (int p=0 ; p<this->globalSize() ; ++p )
    {
        M_mapColorWorld[p] = dataRecvToGather[p].get<0>();
        M_mapGlobalRankToGodRank[p] = dataRecvToGather[p].get<1>();
    }

    this->upMasterRank();

    if ( this->M_isActive[this->godRank()] && this->globalSize() == 1 && this->localSize() == 1 &&
         std::accumulate( this->M_isActive.begin(),this->M_isActive.end(), 0, std::plus<int>() ) == 1 )
    {
        M_subWorldCommSeq.reset( new self_type(*this) );
    }
    else
    {
        this->initSubWorldCommSeq();
    }

}

WorldComm&
WorldComm::operator=( WorldComm const& wc )
{
    if ( this != &wc )
    {
        super::operator=( wc );
        M_localComm = wc.M_localComm;
        M_godComm = wc.M_godComm;
        M_subWorldCommSeq = wc.M_subWorldCommSeq;
        M_mapColorWorld = wc.M_mapColorWorld;
        M_mapLocalRankToGlobalRank = wc.M_mapLocalRankToGlobalRank;
        M_mapGlobalRankToGodRank = wc.M_mapGlobalRankToGodRank;
        M_subworlds = wc.M_subworlds;
        M_masterRank = wc.M_masterRank;
        M_isActive = wc.M_isActive;
    }
    return *this;
}


bool
WorldComm::operator==( WorldComm const& wc ) const
{
    if ( this->localSize() != wc.localSize() || this->globalSize() != wc.globalSize() )
        return false;

    for ( int p = 0 ; p< this->mapColorWorld().size() ; ++p )
        if ( this->mapColorWorld()[p] != wc.mapColorWorld()[p] )
            return false;

    for ( int p = 0 ; p< this->mapLocalRankToGlobalRank().size() ; ++p )
        if ( this->mapLocalRankToGlobalRank()[p] != wc.mapLocalRankToGlobalRank()[p] )
            return false;

    for ( int p = 0 ; p< this->mapGlobalRankToGodRank().size() ; ++p )
        if ( this->mapGlobalRankToGodRank()[p] != wc.mapGlobalRankToGodRank()[p] )
            return false;

    for ( int p = 0 ; p< this->activityOnWorld().size() ; ++p )
        if ( this->activityOnWorld()[p] != wc.activityOnWorld()[p] )
            return false;

    if ( this->masterRank() != wc.masterRank() )
        return false;

    return true;
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
    //std::cout << " WorldComm::subWorldComm( int _color ) " << std::endl;

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
                      this->localComm(),
                      this->godComm(),
                      myColor,
                      newIsActive );

}

//-------------------------------------------------------------------------------
WorldComm::self_type
WorldComm::subWorldComm( std::vector<int> const& colormap )
{
    M_mapColorWorld = colormap;
    return this->subWorldComm();
}

//-------------------------------------------------------------------------------

WorldComm::self_type
WorldComm::subWorldComm( int _color, std::vector<int> const& colormap )
{
    M_mapColorWorld = colormap;
    return this->subWorldComm( _color );
}

//-------------------------------------------------------------------------------

WorldComm::self_type const&
WorldComm::subWorldCommSeq() const
{
    CHECK(M_subWorldCommSeq) << "M_subWorldCommSeq is not initialized()\n";
    return *M_subWorldCommSeq;
}

//-------------------------------------------------------------------------------

void
WorldComm::initSubWorldCommSeq()
{
    std::vector<int> newIsActive(this->godSize(),false);
    newIsActive[this->godRank()]=true;

    M_subWorldCommSeq.reset( new WorldComm( communicator_type(MPI_COMM_SELF,boost::mpi::comm_attach),
                                            communicator_type(MPI_COMM_SELF,boost::mpi::comm_attach),
                                            this->godComm(),
                                            this->godRank(), // local color
                                            newIsActive ) );
}

//-------------------------------------------------------------------------------

void
WorldComm::showMe( std::ostream& __out ) const
{
    this->globalComm().barrier();

    for ( int proc = 0 ; proc < this->globalSize(); ++proc )
    {
        if ( this->globalRank()==proc  && this->isActive() )
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
            std::cout << "mapLocalRankToGlobalRank :\n";
            for ( int k=0; k<( int )this->mapLocalRankToGlobalRank().size(); ++k )
                std::cout << k << " " << this->mapLocalRankToGlobalRank()[k] << std::endl;
            std::cout << "------------------------------------------------"<< std::endl;
            std::cout << "mapGlobalRankToGodRank :\n";
            for ( int k=0; k<( int )this->mapGlobalRankToGodRank().size(); ++k )
                std::cout << k << " " << this->mapGlobalRankToGodRank()[k] << std::endl;
            std::cout << "------------------------------------------------"<< std::endl;
            std::cout << "mapColorWorld :\n";
            for ( int k=0; k<( int )this->mapColorWorld().size(); ++k )
                std::cout << k << " " << this->mapColorWorld()[k] << std::endl;
            std::cout << "------------------------------------------------"<< std::endl;
            std::cout << "isActive :\n";
            for ( int k=0; k<( int )this->M_isActive.size(); ++k )
                std::cout << k << " " << this->M_isActive[k] << std::endl;
            std::cout << "------------------------------------------------"<< std::endl;
            std::cout << "----------------------------------------------------------------------\n";
            std::cout << std::endl;
        }

        this->globalComm().barrier();
    }

    this->globalComm().barrier();
}

//-------------------------------------------------------------------------------

WorldComm::self_type
WorldComm::operator+( WorldComm const & _worldComm ) const
{
    int color;
    bool active;
#if 0
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
#else
    int godRankActif=0;
    int colorSplit=0;
    int colorOnProc=0;
    int newLocRank=0;
    int newGlobRank=0;
    if ( this->isActive()  )
    {
        colorSplit=1;
        active=true;
        godRankActif=this->godRank();
        colorOnProc=this->mapColorWorld()[this->globalRank()];
        newLocRank=this->localRank();
    }
    else if (_worldComm.isActive() )
    {
        colorSplit=1;
        active=true;
        godRankActif=_worldComm.godRank();
        colorOnProc=_worldComm.mapColorWorld()[_worldComm.globalRank()];
        newLocRank= _worldComm.localRank();
    }
    else
    {
        colorSplit=0;
        active=false;
    }
    /*std::cout << "OP+ ---1----"
              << this->localRank()<< " " << this->globalRank() << " " << this->godRank()
              <<"     "  << _worldComm.localRank()<< " " << _worldComm.globalRank() << " " << _worldComm.godRank()
              <<std::endl;this->godComm().barrier();*/

    self_type res(colorOnProc,newLocRank,//this->localRank(),
                  colorSplit,this->godRank(),//this->globalRank(),
                  this->godComm(),
                  active,
                  false );

    std::vector<int> godRankActivefusion(res.globalSize());
    mpi::all_gather( res.globalComm(),
                     godRankActif,
                     godRankActivefusion );
    std::vector<int> newIsActive(this->godSize(),false);
    for (int k=0;k<godRankActivefusion.size();++k)
    {
        newIsActive[godRankActivefusion[k]]=true;
    }
    res.setIsActive(newIsActive);
    res.upMasterRank();

    return res;

#endif
}

//-------------------------------------------------------------------------------

rank_type
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

void
WorldComm::upMasterRank()
{
    // choice : the smallest rank
    M_masterRank = INT_MAX;
    for ( int p=0; p<this->globalSize(); ++p )
        {
            if  (this->isActive() )
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

// Warning : god communication
boost::tuple<bool,std::set<int> >
WorldComm::hasMultiLocalActivity() const
{
    const int _nProc = this->godSize();
    std::set<int> colorWhichIsActive;
    std::vector<int> activities_god(this->godSize());
    std::vector<int> colors_god(this->godSize());
    for (int p=0;p<_nProc;++p)
        {
            mpi::all_gather( this->godComm(),
                             this->activityOnWorld()[p],
                             activities_god );
            mpi::all_gather( this->godComm(),
                             this->mapColorWorld(this->globalRank()),
                             colors_god );

            bool atLesatOneActiveProc  = *std::max_element( activities_god.begin(),activities_god.end() );
            if (atLesatOneActiveProc)
                {
                    const int myactivity=activities_god[0];
                    for (int p2=0;p2<_nProc;++p2)
                        {
                            if (activities_god[p2]) colorWhichIsActive.insert(colors_god[p2]);
                        }
                }
        }
#if 0
    std::for_each(colorWhichIsActive.begin(),
                  colorWhichIsActive.end(),
                  []( int e ) { std::cout << e << std::endl;} );
#endif
    if (colorWhichIsActive.size()>1)
        return boost::make_tuple(true,colorWhichIsActive);
    else
        return boost::make_tuple(false,colorWhichIsActive);
}

//-------------------------------------------------------------------------------

void
WorldComm::applyActivityOnlyOn(int _localColor) const
{
    std::vector<int> _isActive(this->godSize(),0);
    const int _nProc = this->globalSize();
    for (int p=0;p<_nProc;++p)
        {
            if ( this->mapColorWorld()[p]== _localColor )
                _isActive[this->mapGlobalRankToGodRank(p)]=1;
        }
    this->setIsActive(_isActive);
}

//-------------------------------------------------------------------------------

WorldComm const&
WorldComm::subWorld( int n ) const
{
    if ( this->globalSize() == 1 )
        return *this;
    if ( !hasSubWorlds( n ) )
        registerSubWorlds( n );
    return M_subworlds.find(n)->second.second[subWorldId(n)];
}
int
WorldComm::subWorldId( int n ) const
{
    if ( this->globalSize() == 1 )
        return 0;

    if ( !hasSubWorlds( n ) )
        registerSubWorlds( n );
    return M_subworlds.find(n)->second.first.M_mapColorWorld[this->globalRank()];
}
bool
WorldComm::hasSubWorlds( int n ) const
{
    return M_subworlds.find( n ) != M_subworlds.end();
}

std::vector<WorldComm> const&
WorldComm::subWorlds( int n ) const
{
    if ( !hasSubWorlds( n ) )
        registerSubWorlds( n );
    return M_subworlds[n].second;
}

std::vector<WorldComm> const&
WorldComm::subWorldsGroupBySubspace( int n )
{
    if ( !hasSubWorlds( n ) )
        registerSubWorldsGroupBySubspace( n );
    return M_subworlds[n].second;
}

int
WorldComm::numberOfSubWorlds() const
{
    if ( this->globalSize() == 1 || M_subworlds.size() == 0 )
        return 1;
    return M_subworlds.begin()->first;
}

WorldComm const&
WorldComm::masterWorld( int n )
{
    if ( !hasSubWorlds( n ) )
        registerSubWorlds( n );
    return M_subworlds[n].first;
}

void
WorldComm::registerSubWorlds( int n ) const
{
    std::vector<WorldComm> subworlds( n, *this );
    M_subworlds.insert( std::make_pair( n, std::make_pair( *this, subworlds ) ) );
}
void
WorldComm::registerSubWorldsGroupBySubspace( int n )
{
    if ( ( n > 1 ) && ( this->globalSize() > 1 ) )
    {
        // if necessary register a world with n subworlds
        if ( !WorldComm::hasSubWorlds( n ) )
        {
            std::vector<int> MapWorld(this->globalSize());

            // be careful the number of processors must be a multiple of nSpaces for now
            int nprocs_per_space = 0;
            if ( globalSize() <= n )
                nprocs_per_space = 1;
            else
                nprocs_per_space = globalSize()/n;
            for (int proc = 0 ; proc < nprocs_per_space; ++proc)
            {
                for( int s = 0; s < n; ++s )
                    MapWorld[nprocs_per_space*s+proc] = s;
            }
            WorldComm wc(*this);
            wc.setColorMap( MapWorld );

            std::vector<WorldComm> subworlds( n, wc );
            for( int s = 0; s < n; ++s )
            {
                if (this->globalSize()>1)
                    subworlds[s] = wc.subWorldComm( s, MapWorld );
                else
                    subworlds[s] = wc;
            }
            M_subworlds.insert( std::make_pair( n, std::make_pair( wc,  subworlds ) ) );

        }
    }
    else if ( ( n == 1 ) || ( this->globalSize() == 1 ) )
    {
        std::vector<WorldComm> subworlds( n, *this );
        M_subworlds.insert( std::make_pair( n, std::make_pair( *this, subworlds ) ) );
    }
}

void
WorldComm::setColorMap( std::vector<int> const& colormap )
{
    M_mapColorWorld = colormap;
    M_localComm = super::split( M_mapColorWorld[this->globalRank()] );
    LOG(INFO) << "WorldComm::setColorMap: localSize = " << M_localComm.size() << "\n";
    mpi::all_gather( this->localComm(),
                     this->globalRank(),
                     M_mapLocalRankToGlobalRank );
}

} //namespace Feel
