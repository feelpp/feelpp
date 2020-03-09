/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013-2016 Feel++ Consortium

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
   \file straightenmesh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_STRAIGHTENMESH_HPP)
#define FEELPP_STRAIGHTENMESH_HPP 1

#include <feel/feeldiscr/mesh.hpp>
#define FEELPP_INSTANTIATE_STRAIGHTENMESH 
namespace Feel 
{ 
template<typename MeshType>
std::shared_ptr<MeshType>
straightenMesh( std::shared_ptr<MeshType> m, 
                worldcomm_ptr_t const& comm = Environment::worldCommPtr(),
                bool refine = false,
                bool save = false );



#if !defined(FEELPP_INSTANTIATE_STRAIGHTENMESH)
// 1D
extern template std::shared_ptr<Mesh<Simplex<1,1>>>
straightenMesh<Mesh<Simplex<1,1>>>( std::shared_ptr<Mesh<Simplex<1,1>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );
extern template std::shared_ptr<Mesh<Hypercube<1,1>>>
straightenMesh<Mesh<Hypercube<1,1>>>( std::shared_ptr<Mesh<Hypercube<1,1>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );
extern template std::shared_ptr<Mesh<Simplex<1,1,2>>>
straightenMesh<Mesh<Simplex<1,1,2>>>( std::shared_ptr<Mesh<Simplex<1,1,2>>>, 
                                      worldcomm_ptr_t const& , bool, bool  );
extern template std::shared_ptr<Mesh<Hypercube<1,1,2>>>
straightenMesh<Mesh<Hypercube<1,1,2>>>( std::shared_ptr<Mesh<Hypercube<1,1,2>>>, 
                                      worldcomm_ptr_t const& , bool, bool  );
extern template std::shared_ptr<Mesh<Simplex<1,2>>>
straightenMesh<Mesh<Simplex<1,2>>>( std::shared_ptr<Mesh<Simplex<1,2>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );

// 2D
extern template std::shared_ptr<Mesh<Simplex<2,1>>>
straightenMesh<Mesh<Simplex<2,1>>>( std::shared_ptr<Mesh<Simplex<2,1>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );
extern template std::shared_ptr<Mesh<Hypercube<2,1>>>
straightenMesh<Mesh<Hypercube<2,1>>>( std::shared_ptr<Mesh<Hypercube<2,1>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );
extern template std::shared_ptr<Mesh<Simplex<2,2>>>
straightenMesh<Mesh<Simplex<2,2>>>( std::shared_ptr<Mesh<Simplex<2,2>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );

extern template std::shared_ptr<Mesh<Simplex<2,3>>>
straightenMesh<Mesh<Simplex<2,3>>>( std::shared_ptr<Mesh<Simplex<2,3>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );
extern template std::shared_ptr<Mesh<Simplex<2,4>>>
straightenMesh<Mesh<Simplex<2,4>>>( std::shared_ptr<Mesh<Simplex<2,4>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );

// 3D
extern template std::shared_ptr<Mesh<Simplex<3,1>>>
straightenMesh<Mesh<Simplex<3,1>>>( std::shared_ptr<Mesh<Simplex<3,1>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );
extern template std::shared_ptr<Mesh<Hypercube<3,1>>>
straightenMesh<Mesh<Hypercube<3,1>>>( std::shared_ptr<Mesh<Hypercube<3,1>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );
extern template std::shared_ptr<Mesh<Simplex<3,2>>>
straightenMesh<Mesh<Simplex<3,2>>>( std::shared_ptr<Mesh<Simplex<3,2>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );

extern template std::shared_ptr<Mesh<Simplex<3,3>>>
straightenMesh<Mesh<Simplex<3,3>>>( std::shared_ptr<Mesh<Simplex<3,3>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );
extern template std::shared_ptr<Mesh<Simplex<3,4>>>
straightenMesh<Mesh<Simplex<3,4>>>( std::shared_ptr<Mesh<Simplex<3,4>>>, 
                                    worldcomm_ptr_t const& , bool, bool  );
#endif
}


#include <feel/feelfilters/straightenmesh_impl.hpp>

#endif /* FEELPP_STRAIGHTENMESH_HPP */
