/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 18 Dec 2014
 
 Copyright (C) 2014-2016 Feel++ Consortium
 
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
#define FEELPP_INSTANTIATE_STRAIGHTENMESH
#include <feel/feelfilters/straightenmesh_impl.hpp>

namespace Feel
{
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 3 )
template boost::shared_ptr<Mesh<Simplex<3, 3>>>
    straightenMesh<Mesh<Simplex<3, 3>>>( boost::shared_ptr<Mesh<Simplex<3, 3>>>,
                                         WorldComm const&, bool, bool );
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 4 )
template boost::shared_ptr<Mesh<Simplex<3, 4>>>
    straightenMesh<Mesh<Simplex<3, 4>>>( boost::shared_ptr<Mesh<Simplex<3, 4>>>,
                                         WorldComm const&, bool, bool );
#endif
}
