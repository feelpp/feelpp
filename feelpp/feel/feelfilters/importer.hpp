/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-11-06
 */
#ifndef __importer_H
#define __importer_H 1

#include <feel/feelcore/visitor.hpp>
#include <feel/feelcore/commobject.hpp>
#include <feel/feelfilters/enums.hpp>


namespace Feel
{


/**
 * \class Importer
 *
 * import mesh data formats into Feel mesh data structure.
 *
 * \ingroup Importer
 * \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
template<typename MeshType>
class FEELPP_EXPORT Importer
    :
        public CommObject,
        public VisitorBase,
        public Visitor<MeshType>
{
  public:

    using super = CommObject;
    typedef MeshType mesh_type;
    typedef typename mesh_type::point_type point_type;
    typedef typename point_type::node_type node_type;
    typedef typename mesh_type::edge_type edge_type;
    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_type::element_type element_type;

    /**
     * default constructor. use GMSH as default mesh format
     */
    Importer( MeshFormat const& _format = GMSH, 
              FileType t = ASCII,
              worldcomm_ptr_t const& _worldcomm = Environment::worldCommPtr() )
        :
        super( _worldcomm ),
        M_filename(),
        M_format( _format ),
        M_file_type( t )
        {}

    /**
     * constructor
     * @param filename mesh filename to import
     * @param format format of the file
     */
    Importer( std::string const& _filename,  
              MeshFormat const& _format = GMSH, 
              FileType const& _file_type = ASCII, 
              worldcomm_ptr_t const& _worldcomm = Environment::worldCommPtr() )
        :
        super( _worldcomm ),
        M_filename( _filename ),
        M_format( _format ),
        M_file_type( _file_type )
        {}

    ~Importer() override
    {}

    /**
     * set the file name
     * @param __filename
     */
    void setFilename( std::string const& __filename )
    {
        M_filename = __filename;
    }

    /**
     * set the format of the mesh file
     * @param __format  format
     */
    void setFormat( MeshFormat const& __format )
    {
        M_format = __format;
    }

    /**
     * set file type see \p FileType
     */
    void setFileType( FileType t ) 
    {
        M_file_type = t;
    }

    /**
     * \return the filename
     */
    std::string const& filename() const
    {
        return M_filename;
    }


    /**
     * \return the mesh format
     */
    MeshFormat format() const
    {
        return M_format;
    }

    /**
     * file type see \p FileType
     */
    FileType fileType() const { return M_file_type; }

  private:

    //! name of the file to import
    std::string M_filename;

    //! format of the file to import
    MeshFormat M_format;

    //! file type
    FileType M_file_type;
};
}

#endif /* __Importer_H */
