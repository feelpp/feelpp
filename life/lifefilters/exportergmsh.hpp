/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-03-30

  Copyright (C) 2005-2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file ExporterGmsh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-03-30
 */
#ifndef __ExporterGmsh_H
#define __ExporterGmsh_H 1

#include <iostream>
#include <fstream>


#include <boost/lambda/lambda.hpp>

#include <life/lifecore/debug.hpp>

#include <life/lifefilters/exporter.hpp>

namespace Life
{
/**
 * \class ExporterGmsh
 * \brief Exporter to GMSH format
 *
 * \ingroup Exporter
 * @author Christophe Prud'homme
 */
template<typename MeshType>
class ExporterGmsh
    :
        public Exporter<MeshType>
{
public:


    /** @name Typedefs
     */
    //@{

    typedef MeshType mesh_type;

    typedef Exporter<MeshType> super;
    typedef typename mesh_type::value_type value_type;
    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_iterator timeset_iterator;
    typedef typename super::timeset_const_iterator timeset_const_iterator;

    typedef typename matrix_node<value_type>::type matrix_node_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    ExporterGmsh( std::string const& __p = "default", int freq = 1 );

    ExporterGmsh( po::variables_map const& vm, std::string const& exp_prefix = "" );

    ExporterGmsh( ExporterGmsh const & __ex );

    ~ExporterGmsh();

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

    Exporter<MeshType>* setOptions( po::variables_map const& vm, std::string const& exp_prefix = "" )
    {
        super::setOptions( vm, exp_prefix );

        return this;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
       save the timeset
     */
    void save() const;

    /**
     * export mesh
     */
    void visit( mesh_type* mesh );

    //@}



protected:

private:
};

} // Life

#if !defined( LIFE_INSTANTIATION_MODE )
# include <life/lifefilters/exportergmsh.cpp>
#endif // LIFE_INSTANTIATION_MODE

#endif /* __ExporterGmsh_H */

