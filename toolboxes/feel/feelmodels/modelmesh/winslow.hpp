/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2018-05-21

  Copyright (C) 2018 Feel++ Consortium

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

#ifndef FEELPP_TOOLBOXES_MESH_WINSLOW_H
#define FEELPP_TOOLBOXES_MESH_WINSLOW_H 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/projector.hpp>
//#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feelmodels/modelcore/modelalgebraic.hpp>
#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename MeshType, int Order >
class Winslow : public ModelAlgebraic,
                public std::enable_shared_from_this< Winslow<MeshType,Order> >
{
public :
    typedef Winslow<MeshType,Order> self_type;
    typedef ModelAlgebraic super_type;

    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    typedef typename super_type::backend_ptrtype backend_ptrtype;
    typedef typename super_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super_type::vector_ptrtype vector_ptrtype;

    typedef Preconditioner<double> preconditioner_type;
    typedef std::shared_ptr<preconditioner_type> preconditioner_ptrtype;


    typedef bases<Lagrange<Order,Vectorial> > basis_type;
    typedef FunctionSpace<mesh_type,basis_type> space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;

    typedef bases<Lagrange<Order-1,Scalar> > basis_scal_m1_type;
    typedef FunctionSpace<mesh_type,basis_scal_m1_type> space_scal_m1_type;
    typedef std::shared_ptr<space_scal_m1_type> space_scal_m1_ptrtype;
    typedef typename space_scal_m1_type::element_type element_scal_m1_type;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Scalar,Discontinuous> > > space_p0_type;
    typedef std::shared_ptr<space_p0_type> space_p0_ptrtype;
    typedef typename space_p0_type::element_type element_p0_type;
    typedef std::shared_ptr<element_p0_type> element_p0_ptrtype;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Tensor2,Discontinuous> > > space_p0_tensor2_type;
    typedef std::shared_ptr<space_p0_tensor2_type> space_p0_tensor2_ptrtype;
    typedef typename space_p0_tensor2_type::element_type element_p0_tensor2_type;
    typedef std::shared_ptr<element_p0_tensor2_type> element_p0_tensor2_ptrtype;

    typedef Projector<space_scal_m1_type,space_scal_m1_type> projector_scal_m1_type;
    typedef std::shared_ptr<projector_scal_m1_type> projector_scal_m1_ptrtype;

    typedef std::map< std::string, std::vector<flag_type> > flagSet_type;

    typedef Exporter<mesh_type,mesh_type::nOrder> exporter_type;
    typedef std::shared_ptr<exporter_type> exporter_ptrtype;

    Winslow( mesh_ptrtype mesh, std::string const& prefix="",
             worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(), bool useGhostEltFromExtendedStencil=false,
             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    Winslow( space_ptrtype const& space, std::string const& prefix="",
             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    void init();

    std::shared_ptr<std::ostringstream> getInfo() const;


    backend_ptrtype const& backend() const { return M_backend; }
    space_ptrtype const& functionSpace() const { return M_Xh; }
    element_ptrtype const& displacement() const { return M_displacement; }
    element_ptrtype const& dispImposedOnBoundary() const { return M_dispImposedOnBoundary; }
    element_ptrtype const& identity() const { return M_identity; }
    space_scal_m1_ptrtype const& functionSpaceScalM1() const { return M_XhScalM1; }

    flagSet_type const& flagSet() const { return M_flagSet; }
    bool hasFlagSet( std::string const& key ) const { return ( M_flagSet.find(key) != M_flagSet.end() ); }
    std::vector<flag_type> const& flagSet( std::string const& key ) const
        {
            CHECK( hasFlagSet( key ) ) << "the flag type " << key << " is unknown \n";
            return M_flagSet.find(key)->second;
        }
    void setflagSet( flagSet_type const & fl ) { M_flagSet=fl; }

    projector_scal_m1_ptrtype /*const&*/ l2projector() const { return M_l2projector; }

    template < typename elem_type, typename elem2_type >
    void generateALEMap( elem_type const & elem, elem2_type const & elem2 );


    void solve();

    void updateLinearPDE( DataUpdateLinear & data ) const;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const;

    void updateNewtonInitialGuess( vector_ptrtype& U ) const;
    void updateJacobian( DataUpdateJacobian & data ) const;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const;
    void updateResidual( DataUpdateResidual & data ) const;
    void updateResidualDofElimination( DataUpdateResidual & data ) const;

private :
    void updateMeshAdaptation();

    void updateLinearPDE( DataUpdateLinear & data, mpl::int_<2> /**/ ) const;
    void updateLinearPDE( DataUpdateLinear & data, mpl::int_<3> /**/ ) const;
    void updateJacobian( DataUpdateJacobian & data, mpl::int_<2> /**/ ) const;
    void updateJacobian( DataUpdateJacobian & data, mpl::int_<3> /**/ ) const;
    void updateResidual( DataUpdateResidual & data, mpl::int_<2> /**/ ) const;
    void updateResidual( DataUpdateResidual & data, mpl::int_<3> /**/ ) const;


private :

    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    vector_ptrtype M_vectorSolution;
    std::set<size_type> M_dofsWithValueImposed;

    std::string M_solverType;

    mesh_ptrtype M_mesh;
    flagSet_type M_flagSet;
    space_ptrtype M_Xh;
    element_ptrtype M_displacement;
    element_ptrtype M_displacementOld;

    element_ptrtype M_dispImposedOnBoundary;
    element_ptrtype M_identity;

    space_scal_m1_ptrtype M_XhScalM1;
    space_p0_ptrtype M_XhScalP0Disc;

    space_p0_tensor2_ptrtype M_XhTensor2P0Disc;

    projector_scal_m1_ptrtype M_l2projector;

    backend_ptrtype M_backendMetricDerivative;
    sparse_matrix_ptrtype M_matrixMetricDerivative;
    vector_ptrtype M_vectorMetricDerivative;
    element_ptrtype M_fieldMetricDerivative;

    bool M_useMeshAdapation, M_useMeshAdapationScalar;
    element_p0_ptrtype M_weightFunctionScalar;
    element_p0_tensor2_ptrtype M_weightFunctionTensor2;
    element_p0_ptrtype M_hMinRadius;


}; // class Winslow



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

template< typename MeshType, int Order >
template < typename elem_type,typename elem2_type >
void
Winslow<MeshType,Order>::generateALEMap( elem_type const & elem,elem2_type const & elem2 )
{
    bool useGhostEltFromExtendedStencil = M_displacement->functionSpace()->dof()->buildDofTableMPIExtended() && M_displacement->mesh()->worldComm().localSize()>1;
    EntityProcessType entityProcess = (useGhostEltFromExtendedStencil)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    M_dispImposedOnBoundary->on( _range=elements(M_displacement->mesh(),entityProcess),
                                 _expr=idv(elem) );


    // initial value with a previous solution (before solve)
    *M_displacement = *M_identity;
    *M_displacement += vf::project(M_displacement->functionSpace(),elements(M_displacement->mesh(),entityProcess),vf::idv(elem2) );
    *M_displacementOld = *M_identity;
    *M_displacementOld += vf::project(M_displacement->functionSpace(),elements(M_displacement->mesh(),entityProcess),vf::idv(elem2) );

    // solve winslow model
    this->solve();
}

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_TOOLBOXES_MESH_WINSLOW_H
