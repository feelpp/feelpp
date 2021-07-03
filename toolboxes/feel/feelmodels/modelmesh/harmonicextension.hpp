/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Copyright (C) 2014 Université Joseph Fourier (Grenoble I)

  This library Is free software; you can redistribute it and/or
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
   \file harmonicextension.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-05-16
*/

#ifndef FEELPP_MODELS_HARMONICEXTENSION_H
#define FEELPP_MODELS_HARMONICEXTENSION_H 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
//#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelvf/operators.hpp>

#include <feel/feelmodels/modelcore/modelalgebraic.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename MeshType, int Order >
class HarmonicExtension : public ModelAlgebraic,
                          public std::enable_shared_from_this< HarmonicExtension<MeshType,Order> >
{
public :
    typedef HarmonicExtension<MeshType,Order> self_type;
    typedef ModelAlgebraic super_type;

    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef bases<Lagrange<Order,Vectorial> > basis_type;
    typedef FunctionSpace<mesh_type,basis_type> space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;

    // typedef ModelAlgebraicFactory model_algebraic_factory_type;
    // typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    typedef super_type::backend_ptrtype backend_ptrtype;
    typedef super_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef super_type::vector_ptrtype vector_ptrtype;
    //typedef super_type::indexsplit_type indexsplit_type;
    //typedef super_type::indexsplit_ptrtype indexsplit_ptrtype;

    typedef std::map< std::string, std::vector<flag_type> > flagSet_type;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Scalar,Discontinuous> > > space_P0_type;
    typedef std::shared_ptr<space_P0_type> space_P0_ptrtype;

    HarmonicExtension(mesh_ptrtype mesh, backend_ptrtype const& backend,
                      std::string prefix="",
                      worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(),
                      bool useGhostEltFromExtendedStencil=false,
                      ModelBaseRepository const& modelRep = ModelBaseRepository() );

    HarmonicExtension(space_ptrtype space, backend_ptrtype const& backend,
                      std::string prefix="",
                      ModelBaseRepository const& modelRep = ModelBaseRepository() );

    void init();

    std::shared_ptr<std::ostringstream> getInfo() const;

    void updateLinearPDE( DataUpdateLinear & data ) const;

    void solve();

    mesh_ptrtype const& mesh() const;
    space_ptrtype const& functionSpace() const;
    element_ptrtype const& displacement() const;
    element_ptrtype const& dispImposedOnBoundary() const;

    flagSet_type const& flagSet() const;
    std::vector<flag_type> const& flagSet(std::string key) const;
    flag_type flagSet(std::string key, int k) const;
    void setflagSet( flagSet_type const & fl );


    template < typename elem_type >
    void
    generateALEMap( elem_type const & dispOnBoundary )
    {
        bool useGhostEltFromExtendedStencil = this->functionSpace()->dof()->buildDofTableMPIExtended() && this->mesh()->worldComm().localSize()>1;
        EntityProcessType entityProcess = (useGhostEltFromExtendedStencil)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
        *M_dispImposedOnBoundary = vf::project(_space=this->functionSpace(),
                                               _range=elements(this->mesh(),entityProcess),
                                               _expr=vf::idv(dispOnBoundary) );
        this->solve();
    }

private :

    mesh_ptrtype M_mesh;
    space_ptrtype M_Xh;
    element_ptrtype M_displacement;
    element_ptrtype M_dispImposedOnBoundary;
    vector_ptrtype M_vectorSolution;

    flagSet_type M_flagSet;

    space_P0_ptrtype M_XhP0;

    bool M_useAdaptPenal;
};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_HARMONICEXTENSION_H
