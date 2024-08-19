//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 10 Aug 2018
//! @copyright 2018 Feel++ Consortium
//!
#ifndef FEELPP_COMMOBJECT_HPP
#define FEELPP_COMMOBJECT_HPP 1



#include <feel/feelcore/worldcomm.hpp>

namespace Feel {

class CommObject 

{
public:

    CommObject() = default;
    explicit CommObject( worldcomm_ptr_t& w )
        :
        M_worldComm( w )
        {}
    explicit CommObject( worldcomm_ptr_t const& w )
        :
        M_worldComm( w )
        {}
    CommObject( CommObject const& ) = default;
    CommObject( CommObject && ) = default;
    virtual ~CommObject() = default;

    CommObject& operator=( CommObject const& ) = default;
    CommObject& operator=( CommObject && ) = default;

    //!
    //! returns true if worldcomm is allocated, false othewise
    //!
    bool hasWorldComm() const { return bool{M_worldComm}; }
    
    /**
     * \return the mpi communicator
     */
    worldcomm_ptr_t const& worldCommPtr() const
        {
            return M_worldComm;
        }
    worldcomm_ptr_t & worldCommPtr() 
        {
            return M_worldComm;
        }
    /**
     * \return the sequential mpi communicator
     */
    worldcomm_ptr_t const& subWorldCommSeqPtr() const
        {
            return M_worldComm->subWorldCommSeqPtr();
        }
    worldcomm_ptr_t & subWorldCommSeqPtr() 
        {
            return M_worldComm->subWorldCommSeqPtr();
        }

    /**
     * \return the mpi communicator
     */
    WorldComm & worldComm() 
        {
            return *M_worldComm;
        }
    /**
     * \return the mpi communicator
     */
    WorldComm const& worldComm() const
        {
            return *M_worldComm;
        }
    /**
     * \return the mpi communicator
     */
    WorldComm const& comm() const
        {
            return *M_worldComm;
        }

    virtual void setWorldComm( worldcomm_ptr_t const& w )
        {
            M_worldComm = w;
        }
    virtual void setWorldComm( worldcomm_t & w )
        {
            M_worldComm = w.shared_from_this();
        }
private:

    worldcomm_ptr_t M_worldComm;
    
};

/**
 * @fn get worldComm
 * 
 * @param commObj 
 * @return worldcomm_t const& 
 */
inline worldcomm_t const&
worldComm( CommObject const& commObj )
{
    return commObj.worldComm();
}


//!
//! @fn get worldComm 
//! @return the worldcomm associated to the communication object
//!
inline worldcomm_ptr_t const&
worldCommPtr( CommObject const& commObj )
{
    return commObj.worldCommPtr();
}


} // Feel
#endif
