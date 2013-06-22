/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-29

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
   \file dofcomposite.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-29
 */
#ifndef __DofComposite_H
#define __DofComposite_H 1

#include <feel/feelalg/datamap.hpp>


namespace Feel
{
/**
 * \class DofComposite
 * \brief Compositing of degree of freedom table
 *
 *  @author Christophe Prud'homme
 *  @see
 */
class DofComposite : public DataMap
{
    typedef DataMap super;
public:


    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    DofComposite( WorldComm const& _worldComm = Environment::worldComm() ): super( _worldComm ) {}
    DofComposite( size_type n, size_type n_local, WorldComm const& _worldComm = Environment::worldComm() ) : super( n, n_local, _worldComm ) {}
    DofComposite( DofComposite const & dc ) : super( dc ) {}
    ~DofComposite() {}

    //@}

    /** @name Operator overloads
     */
    //@{


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
    
    std::pair<std::map<size_type,size_type>,std::map<size_type,size_type> >
    pointIdToDofRelation(std::string fname="") const
    {
        return std::pair<std::map<size_type,size_type>,std::map<size_type,size_type> >();
    }

    //@}



protected:

private:

};
} // Feel
#endif /* __DofComposite_H */
