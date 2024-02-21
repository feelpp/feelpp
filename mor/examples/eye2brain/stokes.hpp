//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Thomas Saigre @thomas-saigre
//! @date 20 Feb 2024
//! @copyright 2024 Feel++ Consortium
//!
//!

#ifndef FEELPP_EYE2BRAIN_STOKES_HPP
#define FEELPP_EYE2BRAIN_STOKES_HPP

#include <feel/feelmor/modelcrbbase.hpp>
#include <feel/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>

namespace Feel
{

class ParameterDefinition
{
public:
    using parameterspace_type = ParameterSpaceX;
};

class FunctionSpaceDefinition
{
public:
    using value_type = double;
    static const uint16_type Order = 1;
    using convex_type =  Simplex<3>;
    using mesh_type = Mesh<convex_type>;
    using fct_base_type = Lagrange<Order, Scalar>;
    using basis_type = bases<fct_base_type, fct_base_type>;
    using space_type = FunctionSpace<mesh_type, basis_type, value_type>;
    using element_type = typename space_type::element_type;

    using V_basis_type = bases<fct_base_type>;
    using V_space_type = FunctionSpace<mesh_type, V_basis_type, value_type>;

    using J_basis_type = bases<Lagrange<Order+2, Scalar,Discontinuous> >;
    using J_space_type = FunctionSpace<mesh_type, J_basis_type, Discontinuous, value_type>;
};

template<int Order, int Dim>
class FEELPP_EXPORT Eye2BrainStokes : public ModelCrbBase<ParameterDefinition,
                                                         FunctionSpaceDefinition >
{
public:
    using super_type = ModelCrbBase<ParameterDefinition, FunctionSpaceDefinition>;
    using parameter_space_type = ParameterDefinition;
    using function_space_type = FunctionSpaceDefinition;

    using value_type = double;

    using mesh_type = super_type::mesh_type;
    using mesh_ptrtype = super_type::mesh_ptrtype;

    using space_type = super_type::space_type;
    using element_type = super_type::element_type;
    using element_ptrtype = super_type::element_ptrtype;

    using J_space_type = FunctionSpaceDefinition::J_space_type;
    using J_space_ptrtype = std::shared_ptr<J_space_type>;
    using q_sigma_space_type = space_type::template sub_functionspace<0>::type;
    using q_sigma_element_type = q_sigma_space_type::element_type;
    using V_view_type = typename element_type::template sub_element_type<0>;
    using V_view_ptrtype = typename element_type::template sub_element_ptrtype<0>;
    using T_view_ptrtype = typename element_type::template sub_element_ptrtype<1>;

    using current_space_type = FunctionSpace<mesh_type, bases<Lagrange<1, Vectorial> > >;
    using current_element_type = typename current_space_type::element_type;

    using prop_type = ModelProperties;
    using prop_ptrtype = std::shared_ptr<prop_type>;

    using parameter_type = super_type::parameter_type;
    using vectorN_type = super_type::vectorN_type;
    using beta_vector_type = typename super_type::beta_vector_type;
    using beta_type = boost::tuple<beta_vector_type,  std::vector<beta_vector_type> >;
    using affine_decomposition_type = typename super_type::affine_decomposition_type;

    using sparse_matrix_ptrtype = typename super_type::sparse_matrix_ptrtype;

private:
    mesh_ptrtype M_mesh;
    prop_ptrtype M_modelProps;
    std::vector< std::vector< element_ptrtype > > M_InitialGuess;

    element_ptrtype M_VT;
    V_view_ptrtype M_V;
    T_view_ptrtype M_T;
    parameter_type M_mu;

    J_space_ptrtype M_Jh;

public:
    // Constructors
    Eye2BrainStokes();
    Eye2BrainStokes( mesh_ptrtype mesh );

    // Helpers
    int Qa() const;                               // ???
    int Nl() const;                               // ???
    int Ql( int l ) const;                        // ???
    int mQA( int q ) const;                       // ???
    int mLQF( int l, int q ) const;               // ???
    int mCompliantQ( int q ) const;               // ???
    int mIntensityQ( int q ) const;               // ???
    int mAverageTempQ( int q ) const;             // ???
    void resizeQm( bool resizeMatrix = true );
    parameter_type newParameter() { return Dmu->element(); }
    // parameter_type paramFromProperties() const;

    void initModel() override;
    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir ) override;

    void decomposition();

    beta_vector_type computeBetaInitialGuess( parameter_type const& mu ) override;
    beta_type computeBetaQm( element_type const& T, parameter_type const& mu ) override;
    beta_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu ) override;
    beta_type computeBetaQm( parameter_type const& mu ) override;
    beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu, double time=1e30 ) override;
    void fillBetaQm( parameter_type const& mu );

    affine_decomposition_type computeAffineDecomposition() override;
    std::vector<std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA() override;
    std::vector<std::vector<element_ptrtype> > computeInitialGuessAffineDecomposition() override;

    element_type solve( parameter_type const& mu ) override;
    value_type
    output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false) override;

    int mMaxSigma();
    q_sigma_element_type eimSigmaQ(int m);
    vectorN_type eimSigmaBeta( parameter_type const& mu );
    void computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu );
    void computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu, element_type& VT );
}; // Eye2BrainStokes class

} // namespace Feel

#endif
