/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-04-25

  Copyright (C) 2013 Université de Strasbourg

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
   \file exporterexodus.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-04-25
 */
#ifndef __ExporterExodus_H
#define __ExporterExodus_H 1

#include <iostream>
#include <fstream>


#include <boost/lambda/lambda.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <feel/feelmesh/filters.hpp>

namespace Feel
{
namespace fs = boost::filesystem;

/**
 * \class ExporterExodus
 * \brief exporter to Exodus format
 *
 * \ingroup Exporter
 * @author Christophe Prud'homme
 */
template<typename MeshType, int N>
class ExporterExodus
    :
public Exporter<MeshType, N>
{
    typedef Exporter<MeshType, N> super;
public:


    /** @name Typedefs
     */
    //@{

    typedef MeshType mesh_type;

    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_iterator timeset_iterator;
    typedef typename super::timeset_const_iterator timeset_const_iterator;

    //@}

    /** @name Constructors, destructor
     */
    //@{
    /**
     *
     The elements that are supported by the EnSight6 format are:

     \htmlonly
     <pre>
     1                 1------------------2        1----------2--------3
     point                   two node bar                three node bar


     7
     4-------------3          4-------------3
     3                 |             |          |             |
     3                        /\                |             |          |             |
     /\                      /  \               |             |        8 |             | 6
     /  \               6    /    \  5           |             |          |             |
     /    \                  /      \             |             |          |             |
     /      \                /        \            |             |          |             |
     /        \              /          \           |             |          |      5      |
     /          \            /    4       \          1-------------2          1-------------2
     1------------2           1------------2
     three node triangle       six node triangle       four node quadrangle     eight node quadrangle


     /\
     / |\
     /  |4\
     /   |  \
     /    |   \
     /     |    \
     1------|-----\
     \     |    3/
     \    |    /
     \  2|   /
     \  |  /
     \ | /
     \\2/

     four node tetrahedron
     </pre>
     \endhtmlonly

    */
    ExporterExodus( WorldComm const& worldComm = Environment::worldComm() );
    ExporterExodus( std::string const& __p = "default", int freq = 1, WorldComm const& worldComm = Environment::worldComm() );
    ExporterExodus( po::variables_map const& vm=Environment::vm(), std::string const& exp_prefix = "", WorldComm const& worldComm = Environment::worldComm() ) FEELPP_DEPRECATED;

    ExporterExodus( ExporterExodus const & __ex );

    ~ExporterExodus();


    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the ensight element type
     */
    std::string const& elementType() const
    {
        return M_element_type;
    }


    //@}

    /** @name  Mutators
     */
    //@{

    Exporter<MeshType,N>* setOptions( po::variables_map const& vm, std::string const& exp_prefix = "" ) FEELPP_DEPRECATED
    {
        super::setOptions( exp_prefix );

        return this;
    }

    Exporter<MeshType,N>* setOptions( std::string const& exp_prefix = "" )
    {
        super::setOptions( exp_prefix );

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

    void visit( mesh_type* mesh );

    //@}



protected:

private:

    /**
     * init the ensight exporter
     */
    void init();


private:

    mutable std::string M_filename;
    std::string M_element_type;
};

template<typename MeshType, int N>
ExporterExodus<MeshType,N>::ExporterExodus( WorldComm const& worldComm )
:
super( worldComm ),
M_element_type()

{
    init();
}
template<typename MeshType, int N>
ExporterExodus<MeshType,N>::ExporterExodus( std::string const& __p, int freq, WorldComm const& worldComm )
    :
    super( "exodus", __p, freq, worldComm ),
    M_element_type()
{
    init();
}
template<typename MeshType, int N>
ExporterExodus<MeshType,N>::ExporterExodus( po::variables_map const& vm, std::string const& exp_prefix, WorldComm const& worldComm )
    :
    super( vm, exp_prefix, worldComm )
{
    init();
}

template<typename MeshType, int N>
ExporterExodus<MeshType,N>::ExporterExodus( ExporterExodus const & __ex )
    :
    super( __ex ),
    M_element_type( __ex.M_element_type )
{
}

template<typename MeshType, int N>
ExporterExodus<MeshType,N>::~ExporterExodus()
{}

template<typename MeshType, int N>
void
ExporterExodus<MeshType,N>::init()
{
    if ( mesh_type::nDim == 1 )
        if ( mesh_type::Shape == SHAPE_LINE )
            M_element_type = ( mesh_type::nOrder == 1 )?"bar2":"bar3";

    if ( mesh_type::nDim == 2 )
    {
        if ( mesh_type::Shape == SHAPE_TRIANGLE )
            M_element_type = ( mesh_type::nOrder == 1 )?"tria3":"tria6";

        else if ( mesh_type::Shape == SHAPE_QUAD )
            M_element_type = ( mesh_type::nOrder == 1 )?"quad4":"quad8";
    }

    if ( mesh_type::nDim == 3 )
    {
        if ( mesh_type::Shape == SHAPE_TETRA )
            M_element_type = ( mesh_type::nOrder == 1 )?"tetra4":"tetra10";

        else if ( mesh_type::Shape == SHAPE_HEXA )
            M_element_type = ( mesh_type::nOrder == 1 )?"hexa8":"hexa20";
    }
}
template<typename MeshType, int N>
void
ExporterExodus<MeshType,N>::save() const
{
    if ( !this->worldComm().isActive() ) return;

}


template<typename MeshType, int N>
void
ExporterExodus<MeshType,N>::visit( mesh_type* __mesh )
{
}

} // Feel


#endif /* __ExporterExodus_H */
