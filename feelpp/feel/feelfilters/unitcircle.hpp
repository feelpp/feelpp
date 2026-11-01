/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Alexandre Ancel <alexandre.ancel@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2013-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file unitcircle.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_UNITCIRCLE_HPP)
#define FEELPP_UNITCIRCLE_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#ifdef FEELPP_HAS_GMSH
#include <feel/feelfilters/creategmshmesh.hpp>
#endif
#include <feel/feelfilters/domain.hpp>

namespace Feel {

/**
 * build a mesh of the unit circle using triangles
 */
template<int Ngeo=1>
inline
std::shared_ptr<Mesh<Simplex<2,Ngeo> > >
unitCircle( double h = doption(_name="gmsh.hsize") )
{
#ifdef FEELPP_HAS_GMSH
    return createGMSHMesh(_mesh=new Mesh<Simplex<2,Ngeo> >,
                          _desc=domain( _name="circle",
                                        _shape="ellipsoid",
                                        _dim=2,
                                        _xmin=-1,
                                        _ymin=-1,
                                        _h=h ) );
#else
    LOG(WARNING) << "unitCircle: Feel++ was not built with Gmsh. This function will return a empty mesh.";
    return std::make_shared<Mesh<Simplex<2, Ngeo> > >();
#endif
}

}

#endif /* FEELPP_UNITCIRCLE_HPP */
