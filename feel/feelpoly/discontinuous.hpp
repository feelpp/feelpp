/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-02-17

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2013 Feel++ Consortium

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
   \file discontinuous.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-02-17
 */
#ifndef __Discontinuous_H
#define __Discontinuous_H 1

#include <boost/fusion/sequence.hpp>
#include <feel/feelpoly/continuity.hpp>

namespace Feel
{
// import fusion namespace in Feel
namespace fusion = boost::fusion;

/**
 * \class Continuous
 * \brief describe continuous functions
 *
 * @author Christophe Prud'homme
 * @see
 */
class Discontinuous
    :
    // necessary for boost.parameters
public Feel::detail::continuity_base
{
public:


    /** @name Constants
     */
    //@{

    static const bool is_continuous = false;
    static const bool is_discontinuous_locally = false;
    static const bool is_discontinuous_totally = true;


    //@}

    /** @name Typedefs
     */
    //@{
    typedef fusion::vector<> discontinuity_markers_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    Discontinuous();
    //! copy constructor
    Discontinuous( Continuous const & );
    //! destructor
    ~Discontinuous();

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    Discontinuous& operator=( Discontinuous const & o )
    {
        if ( this != &o )
        {
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}

    template<typename MeshType, typename DofType>
    class apply
    {
    public:
        typedef size_type result_type;
        typedef MeshType mesh_type;
        typedef DofType dof_type;
        typedef typename dof_type::fe_type fe_type;

        apply( MeshType& M, DofType& D )
            :
            M_mesh( M ),
            M_dof( D )
        {}
        template<typename T>
        result_type operator()( const T& t, const size_type& start ) const
        {
            //return build( T::value, start );
            return start;
        }
    private:
    private:
        MeshType& M_mesh;
        DofType M_dof;
    };

protected:

private:

};

} // Feel
#endif /* __Discontinuous_H */
