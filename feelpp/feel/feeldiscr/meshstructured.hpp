/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 01 Jun 2016

 Copyright (C) 2016 Feel++ Consortium

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
#ifndef FEELPP_MESHSTRUCTURED_HPP
#define FEELPP_MESHSTRUCTURED_HPP 1

#include <feel/feeldiscr/mesh.hpp>

namespace Feel
{

template <typename T>
using holo3_image = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;


//! Structured mesh class
//!
//! A structured mesh is such that points and elements can
//! be located with a simple index pair (i,j)
//!
//! \code
//! TODO add code example here
//! \endcode
//!
class MeshStructured : public Mesh<Hypercube<2>>
{

  public:
    using super = Mesh<Hypercube<2, 1, 2>>;
    using point_type = super::point_type;
    using element_type = super::element_type;
    using face_type = super::face_type;
    using node_type = super::node_type;

    MeshStructured(): super() { this->setStructureProperty( "00010" ); }
    MeshStructured( MeshStructured const& ) = delete;
    MeshStructured( MeshStructured&& ) = delete;
    MeshStructured& operator=( MeshStructured const& ) = delete;
    MeshStructured& operator=( MeshStructured&& ) = delete;

    //!
    //!
    //!
    //MeshStructured( int nx, int ny, double pixelsize, WorldComm const& );
    MeshStructured( int nx, int ny, double pixelsize,
                    std::optional<holo3_image<float>> const& cx, std::optional<holo3_image<float>> const& cy,
                    worldcomm_ptr_t const& = Environment::worldCommPtr(), bool withCoord = false, std::string pathPoly = "", bool withPoly = false );
    //MeshStructured( int nx, int ny,holo3_image<float> cx,holo3_image<float> cy, WorldComm const& );

    void updateGhostCellInfoByUsingNonBlockingComm(
        std::unordered_map<int, int> const& idStructuredMeshToFeelMesh,
        std::unordered_map<int, boost::tuple<int, rank_type>> const& mapGhostElt,
        std::vector<int> const& nbMsgToRecv );

    /*
     * Create an element
     * g_i = global x index
     * g_j = global y index
     */
    element_type newElement( int g_i, int g_j )
    {
        element_type e;
        return e;
    }

  private:
    size_type M_nx;                         // Global X number of elements
    size_type M_ny;                         // Global Y number of elements
    std::optional<holo3_image<float>> M_cx; // X-coordinates for nodes
    std::optional<holo3_image<float>> M_cy; // Y-coordinates for nodes
    int M_l_nx;                             // local X number of elements (ghost excluded!)
    int M_l_ny;                             // local Y number of elements
    int M_s_x;                              // local first x index (0 for first element)
    int M_s_y;                              // local first y index (0 for first element)
    double M_pixelsize;
    // std::map<int,boost::tuple<int,rank_type> > mapGhostElt;
    // std::vector<rank_type> ghosts;
    // std::map<int,int> __idGmshToFeel;
#if 0
    int localToGlobal( int ii, int jj, int rank )
    {
        return ii * M_ny + ( jj + ( rank * M_ny / this->worldComm().godSize() ) );
    }
    int globalToLocal( int i, int j, int rank )
    {
        int LM_ny = ( M_ny / this->worldComm().godSize() ) * ( rank + 1 ) - ( M_ny / this->worldComm().godSize() ) * rank;
        return i * LM_ny + j - rank * M_ny / this->worldComm().godSize();
    }
#endif    
};

} // namespace Feel

#endif
