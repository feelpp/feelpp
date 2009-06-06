/* -*- mode: c++ -*-

  This file is part of the Life library

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

#include <life/lifecore/application.hpp>


namespace Life
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

    DataMap();

    /**
     * \param n total size of the vector
     * \param n_local local size of the vector on the curent processor
     */
    DataMap( size_type n, size_type n_local );

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
     * @return the number of degrees of freedom on this processor.
     */
    size_type nLocalDof () const
    { return this->nDofOnProcessor (Application::processId()); }

    /**
     * @return the number of degrees of freedom on this processor.
     */
    size_type nMyDof () const
    { return this->nDofOnProcessor (Application::processId()); }

    /**
     * @return the number of degrees of freedom on subdomain \p proc.
     */
    size_type nDofOnProcessor(const size_type proc = Application::processId() ) const
    {
        LIFE_ASSERT(proc < _M_first_df.size())( proc )( _M_first_df.size() ).error( "invalid proc id or dof table" );
        return ( _M_last_df[proc] - _M_first_df[proc]+1);
    }

    size_type nProcessors() const
    {
        return ( Application::nProcess() );
    }

    /**
     * @return the first dof index that is local to subdomain \p proc.
     */
    size_type firstDof(const size_type proc = Application::processId()) const
    {
        LIFE_ASSERT(proc < _M_first_df.size())( proc )( _M_first_df.size() ).error( "invalid proc id or dof table" );
        return _M_first_df[proc];
    }

    /**
     * Returns the last dof index that is local to subdomain \p proc.
     */
    size_type lastDof(const unsigned int proc = Application::processId()) const
    {
        LIFE_ASSERT(proc < _M_last_df.size())( proc )( _M_last_df.size() ).error( "invalid proc id or dof table" );
        return _M_last_df[proc];
    }


    //! Returns local ID of global ID, return invalid_size_type_value if not found on this processor.
    size_type  lid(size_type GID) const
    {
        uint16_type pid = Application::processId();
        if ( GID >= firstDof( pid ) &&
             GID <= lastDof( pid ) )
            return GID - firstDof( pid );
        return invalid_size_type_value;
    }

    //! Returns global ID of local ID, return -1 if not found on this processor.
    size_type gid( size_type LID) const
    {
        uint16_type pid = Application::processId();
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
    size_type  maxAllGID() const {return(lastDof( Application::nProcess()-1 ) );}

    //! Returns the maximum global ID owned by this processor.
    size_type  minMyGID() const {return firstDof( Application::processId() );}

    //! Returns the maximum global ID owned by this processor.
    size_type  maxMyGID() const {return lastDof( Application::processId() );};

    //!  The minimum local index value on the calling processor.
    size_type  minLID() const {return 0;};

    //! The maximum local index value on the calling processor.
    size_type  maxLID() const {return lastDof( Application::processId() )-firstDof( Application::processId() );};

    //! number of elements across all processors.
    size_type nGlobalElements() const {return _M_n_dofs;};

    //! number of elements on the calling processor.
    size_type nMyElements() const {return nLocalDof();};

    //! Puts list of global elements on this processor size_typeo the user-provided array.
    std::vector<size_type> const& myGlobalElements() const;

    //! \return true if DataMap is close, false otherwise
    bool closed() const { return M_closed; }

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
     * First DOF index on processor \p p.
     */
    std::vector<size_type> _M_first_df;

    /**
     * Last DOF index (plus 1) on processor \p p.
     */
    std::vector<size_type> _M_last_df;

    mutable std::vector<size_type> M_myglobalelements;

private:

};

} // Life
#endif /* __DataMap_H */
