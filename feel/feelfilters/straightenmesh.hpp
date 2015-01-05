/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013-2015 Feel++ Consortium

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


namespace Feel 
{ 
template<typename MeshType>
boost::shared_ptr<MeshType>
straightenMesh( boost::shared_ptr<MeshType> m, 
                WorldComm const& comm = Environment::worldComm(),
                bool refine = false,
                bool save = false );

#if !defined(FEELPP_INSTANTIATE_STRAIGHTENMESH)
// 1D
extern template boost::shared_ptr<Mesh<Simplex<1,1>>>
straightenMesh<Mesh<Simplex<1,1>>>( boost::shared_ptr<Mesh<Simplex<1,1>>>, 
                                    WorldComm const& , bool, bool  );
extern template boost::shared_ptr<Mesh<Simplex<1,1,2>>>
straightenMesh<Mesh<Simplex<1,1,2>>>( boost::shared_ptr<Mesh<Simplex<1,1,2>>>, 
                                      WorldComm const& , bool, bool  );
extern template boost::shared_ptr<Mesh<Hypercube<1,1,2>>>
straightenMesh<Mesh<Hypercube<1,1,2>>>( boost::shared_ptr<Mesh<Hypercube<1,1,2>>>, 
                                      WorldComm const& , bool, bool  );
extern template boost::shared_ptr<Mesh<Simplex<1,2>>>
straightenMesh<Mesh<Simplex<1,2>>>( boost::shared_ptr<Mesh<Simplex<1,2>>>, 
                                    WorldComm const& , bool, bool  );

// 2D
extern template boost::shared_ptr<Mesh<Simplex<2,1>>>
straightenMesh<Mesh<Simplex<2,1>>>( boost::shared_ptr<Mesh<Simplex<2,1>>>, 
                                    WorldComm const& , bool, bool  );
extern template boost::shared_ptr<Mesh<Hypercube<2,1>>>
straightenMesh<Mesh<Hypercube<2,1>>>( boost::shared_ptr<Mesh<Hypercube<2,1>>>, 
                                    WorldComm const& , bool, bool  );
extern template boost::shared_ptr<Mesh<Simplex<2,2>>>
straightenMesh<Mesh<Simplex<2,2>>>( boost::shared_ptr<Mesh<Simplex<2,2>>>, 
                                    WorldComm const& , bool, bool  );

extern template boost::shared_ptr<Mesh<Simplex<2,3>>>
straightenMesh<Mesh<Simplex<2,3>>>( boost::shared_ptr<Mesh<Simplex<2,3>>>, 
                                    WorldComm const& , bool, bool  );
extern template boost::shared_ptr<Mesh<Simplex<2,4>>>
straightenMesh<Mesh<Simplex<2,4>>>( boost::shared_ptr<Mesh<Simplex<2,4>>>, 
                                    WorldComm const& , bool, bool  );

// 3D
extern template boost::shared_ptr<Mesh<Simplex<3,1>>>
straightenMesh<Mesh<Simplex<3,1>>>( boost::shared_ptr<Mesh<Simplex<3,1>>>, 
                                    WorldComm const& , bool, bool  );
extern template boost::shared_ptr<Mesh<Simplex<3,2>>>
straightenMesh<Mesh<Simplex<3,2>>>( boost::shared_ptr<Mesh<Simplex<3,2>>>, 
                                    WorldComm const& , bool, bool  );
#endif
}
#endif /* FEELPP_STRAIGHTENMESH_HPP */
