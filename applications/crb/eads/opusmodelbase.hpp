/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-07-20

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file opusmodelbase.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-07-20
 */
#ifndef __OpusModelBase_H
#define __OpusModelBase_H 1

#include <opusdata.hpp>
#include <eads.hpp>

namespace Feel
{
/**
 * \addtogroup Models
 * \\@{
 */

/**
 * \class OpusModelBase
 * \brief base class for all Opus Models
 *
 * @author Christophe Prud'homme
 * @see
 */
class OpusModelBase
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{

    typedef OpusData opusdata_type;
    typedef boost::shared_ptr<opusdata_type> opusdata_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    OpusModelBase();

    //! constructor frm po::variables_map
    OpusModelBase( po::variables_map const& vm );

    //! copy constructor
    OpusModelBase( OpusModelBase const & );

    //! destructor
    ~OpusModelBase();

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    OpusModelBase& operator=( OpusModelBase const & o )
    {
        if ( this != &o )
        {
            M_data = o.M_data;
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    //! \return the data associated with the opus model
    opusdata_ptrtype data() const
    {
        return M_data;
    }

    //! return the variables_map
    po::variables_map vm()
    {
        return M_data->vm();
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //! set the data for the model
    void setData( opusdata_ptrtype data )
    {
        M_data = data;
    }

    //@}

    /** @name  Methods
     */
    //@{
    virtual void run( const double * X, unsigned long N, double * Y, unsigned long P ) {};
    virtual void run()
    {
        M_data->print();
    }

    //@}



protected:

private:
    opusdata_ptrtype M_data;
};
}
/**
 * \\@}
 */
#endif /* __OpusModelBase_H */
