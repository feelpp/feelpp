/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2004-11-06

  Copyright (C) 2004 EPFL

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
   \file importer.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2004-11-06
 */
#ifndef __importer_H
#define __importer_H 1

#include <feel/feelcore/visitor.hpp>

namespace Feel
{

/**
  \enum MeshFormat
*/
enum MeshFormat {
    MESHPP,
    INRIA,
    GMSH,
    NETGEN,
    GAMBIT
};

/**
 * \class Importer
 *
 * import mesh data formats into Feel mesh data structure.
 *
 * \ingroup Importer
 * \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
 */
template<typename MeshType>
class Importer
    :
        public VisitorBase,
        public Visitor<MeshType>
{
public:

    typedef MeshType mesh_type;
    typedef typename mesh_type::point_type point_type;
    typedef typename point_type::node_type node_type;
    typedef typename mesh_type::edge_type edge_type;
    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_type::element_type element_type;

    /**
     * default constructor. use GMSH as default mesh format
     */
    Importer( MeshFormat const& format = GMSH )
        :
        _M_filename(),
        _M_format( format )
        {}

    /**
     * constructor
     * @param filename mesh filename to import
     * @param format format of the file
     */
    Importer( std::string const& _filename,  MeshFormat const& _format = GMSH )
        :
        _M_filename( _filename ),
        _M_format( _format )
        {}

    virtual ~Importer()
    {}

    /**
     * set the file name
     * @param __filename
     */
    void setFilename( std::string const& __filename )
    {
        _M_filename = __filename;
    }

    /**
     * set the format of the mesh file
     * @param __format  format
     */
    void setFormat( MeshFormat const& __format )
    {
        _M_format = __format;
    }

    /**
     * \return the filename
     */
    std::string const& filename() const { return _M_filename; }


    /**
     * \return the mesh format
     */
    MeshFormat format() const { return _M_format; }

private:

    //! name of the file to import
    std::string _M_filename;

    //! format of the file to import
    MeshFormat _M_format;
};
}

#endif /* __Importer_H */
