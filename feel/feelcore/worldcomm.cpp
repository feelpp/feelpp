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

#include <feel/feelcore/worldcomm.hpp>

namespace Feel
{

using namespace boost;


WorldComm::WorldComm()
    :
    super(),
    M_localComm(super::split(0)),
    M_godComm(),
    M_mapColorWorld(this->globalSize()),
    M_mapLocalRankToGlobalRank(this->localSize()),
    M_mapGlobalRankToGodRank(this->globalSize()),
    M_isActive(this->godSize(),true)
{
    std::vector<int> globalRanks(this->globalSize());
    mpi::all_gather( this->globalComm(),
                     this->globalRank(),
                     globalRanks );

    mpi::all_gather( this->globalComm(),
                     0,
                     M_mapColorWorld );

    mpi::all_gather( this->localComm(),
                     this->globalRank(),
                     M_mapLocalRankToGlobalRank );

   mpi::all_gather( this->globalComm(),
                    this->godRank(),
                    M_mapGlobalRankToGodRank );

    // choice : the smallest rank
    M_masterRank = *std::min_element(globalRanks.begin(),globalRanks.end());
}

WorldComm::WorldComm(int color)
    :
    super(),
    M_localComm(super::split(color)),
    M_godComm(),
    M_mapColorWorld(this->globalSize()),
    M_mapLocalRankToGlobalRank(this->localSize()),
    M_mapGlobalRankToGodRank(this->globalSize()),
    M_isActive(this->godSize(),true)
{
    std::vector<int> globalRanks(this->globalSize());
    mpi::all_gather( this->godComm(),
                     this->globalRank(),
                     globalRanks);

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
    M_masterRank = *std::min_element(globalRanks.begin(),globalRanks.end());
}

WorldComm::WorldComm(std::vector<int> const& colorWorld)
    :
    super(),
    M_localComm(super::split(colorWorld[this->globalRank()])),
    M_godComm(),
    M_mapColorWorld(colorWorld),
    M_mapLocalRankToGlobalRank(this->localSize()),
    M_mapGlobalRankToGodRank(this->globalSize()),
    M_isActive(this->godSize(),true)
{
    std::vector<int> globalRanks(this->globalSize());
    mpi::all_gather( this->globalComm(),
                     this->globalRank(),
                     globalRanks);

    mpi::all_gather( this->localComm(),
                     this->globalRank(),
                     M_mapLocalRankToGlobalRank );

    mpi::all_gather( this->globalComm(),
                     this->godRank(),
                     M_mapGlobalRankToGodRank );

    // choice : the smallest rank
    M_masterRank = *std::min_element(globalRanks.begin(),globalRanks.end());
}

WorldComm::WorldComm( WorldComm const& wc )
    :
    super(wc),
    M_localComm(wc.M_localComm),
    M_godComm(wc.M_godComm),
    M_mapColorWorld(wc.M_mapColorWorld),
    M_mapLocalRankToGlobalRank(wc.M_mapLocalRankToGlobalRank),
    M_mapGlobalRankToGodRank(wc.M_mapGlobalRankToGodRank),
    M_masterRank(wc.M_masterRank),
    M_isActive(wc.M_isActive)
{}

    WorldComm::WorldComm( communicator_type const& _globalComm, int _color, bool _isActive)
    :
    super(_globalComm),
    M_localComm(super::split(_color)),
    M_godComm(),
    M_mapColorWorld(this->globalSize()),
    M_mapLocalRankToGlobalRank(this->localSize()),
    M_mapGlobalRankToGodRank(this->globalSize())
{
    mpi::all_gather( this->godComm(),
                     (int)_isActive,
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
    //M_masterRank = *std::min_element(M_mapRealRank.begin(),M_mapRealRank.end());
    M_masterRank = INT_MAX;
    int indexLoc=0;
    //for (auto it=M_mapRealRank.begin(),en=M_mapRealRank.end();it!=en;++it)
    for (auto it=this->mapLocalRankToGlobalRank().begin(),
             en=this->mapLocalRankToGlobalRank().end();it!=en;++it,++indexLoc)
        {
            if  (_isActive)
                {
                    if (M_isActive[this->mapGlobalRankToGodRank()[*it]])
                        {
                            if (M_masterRank>*it) M_masterRank=*it;
                        }
                }
            else
                {
                    if (!M_isActive[this->mapGlobalRankToGodRank()[*it]])
                        {
                            if (M_masterRank>*it) M_masterRank=*it;
                        }

                }
        }

}

WorldComm::self_type
WorldComm::subWorldComm(int _color)
{
    bool isActive;
    int myColor = this->mapColorWorld()[this->globalRank()];
    if (myColor==_color)
        isActive=true;
    else
        isActive=false;

    return self_type(this->localComm(), myColor, isActive);
}

} //namespace Feel

