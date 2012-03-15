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
   \file worldcomm.hpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2012-02-14
 */
#ifndef __worldcomm_H
#define __worldcomm_H 1

#include <boost/mpi.hpp>
#if defined(FEELPP_HAS_MPI_H)
#include <mpi.h>
#endif /* FEELPP_HAS_MPI_H */

namespace Feel
{


class WorldComm : public boost::mpi::communicator
{

    typedef boost::mpi::communicator super;

public:

    typedef WorldComm self_type;
    typedef boost::mpi::communicator communicator_type;

    WorldComm();

    WorldComm(int _color);

    WorldComm(std::vector<int> const& _colorWorld);

    WorldComm( WorldComm const& _wc );

    WorldComm( communicator_type const& _globalComm,
               int _color,
               bool _isActive);

    communicator_type const& globalComm() const { return *this; }
    communicator_type const& localComm() const { return M_localComm; }
    communicator_type const& godComm() const { return M_godComm; }
    communicator_type const& comm() const { return M_localComm; }

    int globalSize() const { return this->globalComm().size(); }
    int localSize() const { return this->localComm().size(); }
    int godSize() const { return this->godComm().size(); }

    int globalRank() const { return this->globalComm().rank(); }
    int localRank() const { return this->localComm().rank(); }
    int godRank() const { return this->godComm().rank(); }


    std::vector<int> const& mapColorWorld() const { return M_mapColorWorld; }
    std::vector<int> const& mapLocalRankToGlobalRank() const { return M_mapLocalRankToGlobalRank; }
    std::vector<int> const& mapGlobalRankToGodRank() const { return M_mapGlobalRankToGodRank; }

    int masterRank() const { return M_masterRank; }

    WorldComm subWorldComm(int color) const;

    bool isActive() const { return M_isActive[this->godRank()]; }

    int localColorToGlobalRank(int _color,int _localRank) const;

    /**
     * showMe
     */
    void showMe( std::ostream& __out = std::cout ) const;

    WorldComm operator+(WorldComm const & _worldComm) const;

private :

    communicator_type M_localComm;
    communicator_type M_godComm;

    std::vector<int> M_mapColorWorld;
    std::vector<int> M_mapLocalRankToGlobalRank;
    std::vector<int> M_mapGlobalRankToGodRank;

    int M_masterRank;
    std::vector<int/*bool*/> M_isActive;

};

} //namespace Feel

#endif // __worldcomm_H
