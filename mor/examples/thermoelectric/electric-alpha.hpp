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
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!

#ifndef FEELPP_ALPHAELECTRIC_HPP
#define FEELPP_ALPHAELECTRIC_HPP

#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelcrb/deim.hpp>

namespace Feel
{

FEELPP_EXPORT AboutData
makeAbout( std::string const& str = "electriccrbmodel" )
{
    AboutData about( str.c_str() );
    return about;
}

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

    using fct_base_type = Lagrange<Order,Scalar>;
    using basis_type = bases<fct_base_type>;
    using space_type = FunctionSpace<mesh_type, basis_type, value_type>;
    using element_type = typename space_type::element_type;

    // using V_basis_type = bases<fct_base_type>;
    // using V_space_type = FunctionSpace<mesh_type, V_basis_type, value_type>;

    using J_basis_type = bases<Lagrange<Order+2, Scalar,Discontinuous> >;
    using J_space_type = FunctionSpace<mesh_type, J_basis_type, Discontinuous, value_type>;
};

template <typename ParameterDefinition, typename FunctionSpaceDefinition >
class EimDefinition
{
public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename FunctionSpaceDefinition::space_type space_type;
    typedef typename FunctionSpaceDefinition::J_space_type J_space_type;

    /* EIM */
    // Scalar continuous //
    typedef EIMFunctionBase<space_type, space_type , parameterspace_type> fun_type;
    // Scalar Discontinuous //
    typedef EIMFunctionBase<J_space_type, space_type , parameterspace_type> fund_type;

};

class FEELPP_EXPORT AlphaElectric : public ModelCrbBase<ParameterDefinition,
                                                        FunctionSpaceDefinition,
                                                        Linear ,
                                                        EimDefinition<ParameterDefinition,
                                                                      FunctionSpaceDefinition> >
{
public:
    using super_type = ModelCrbBase<ParameterDefinition,
                                    FunctionSpaceDefinition,
                                    Linear,
                                    EimDefinition<ParameterDefinition,
                                                  FunctionSpaceDefinition> >;
    using self_type = AlphaElectric;
    using parameterspace_type = ParameterDefinition::parameterspace_type;
    using function_space_type = FunctionSpaceDefinition;
    using eim_definition_type = EimDefinition<ParameterDefinition, FunctionSpaceDefinition>;

    using value_type = double;

    using mesh_type = super_type::mesh_type;
    using mesh_ptrtype = super_type::mesh_ptrtype;

    using space_type = super_type::space_type;
    using element_type = super_type::element_type;
    using element_ptrtype = super_type::element_ptrtype;

    using J_space_type = FunctionSpaceDefinition::J_space_type;
    // using V_space_type = FunctionSpaceDefinition::V_space_type;
    // using Vh_element_type = typename V_space_type::element_type;
    using q_sigma_space_type = space_type;
    using q_sigma_element_type = q_sigma_space_type::element_type;
    // using V_view_type = typename element_type::template sub_element_type<0>;
    // using V_view_ptrtype = typename element_type::template sub_element_ptrtype<0>;
    // using T_view_ptrtype = typename element_type::template sub_element_ptrtype<1>;

    using current_space_type = FunctionSpace<mesh_type, bases<Lagrange<1, Vectorial> > >;
    using current_element_type = typename current_space_type::element_type;

    using prop_type = ModelProperties;
    using prop_ptrtype = std::shared_ptr<prop_type>;
    using mat_type = ModelMaterial;
    using map_mat_type = std::map<std::string,mat_type>;

    using parameter_type = super_type::parameter_type;
    using vectorN_type = super_type::vectorN_type;
    using beta_vector_type = typename super_type::beta_vector_type;
    using beta_type = boost::tuple<beta_vector_type,  std::vector<beta_vector_type> >;
    using affine_decomposition_type = typename super_type::affine_decomposition_type;

    using sparse_matrix_ptrtype = typename super_type::sparse_matrix_ptrtype;

    using mdeim_type = MDEIM<self_type>;
    using mdeim_ptrtype = std::shared_ptr<mdeim_type>;
    using deim_type = DEIM<self_type>;
    using deim_ptrtype = std::shared_ptr<deim_type>;

private:
    mesh_ptrtype M_mesh;
    prop_ptrtype M_modelProps;
    map_mat_type M_materials;
    map_mat_type M_materialsWithGeo;
    map_mat_type M_materialsWithoutGeo;

    std::vector< std::vector< element_ptrtype > > M_InitialGuess;

    element_ptrtype pT;
    parameter_type M_mu;

    double M_sigma;

    int M_trainsetDeimSize;
    int M_trainsetMdeimSize;
    double M_penalDir;
    std::string M_propertyPath;
    bool M_testDeim;
    std::string M_dbBasename;
    int M_verbose;

public:
    static po::options_description makeOptions( std::string const& prefix="" );
    // Constructors
    explicit AlphaElectric(std::string const& prefix = "");
    AlphaElectric( mesh_ptrtype mesh, std::string const& prefix = "" );

    // Helpers
    int Qa();
    int Nl();
    int Ql( int l );
    int mQA( int q );
    int mLQF( int l, int q );
    int mCompliantQ( int q );
    void resizeQm( bool resizeMat = true );
    parameter_type newParameter() { return Dmu->element(); }
    parameter_type paramFromVec( std::vector<double> const& x );
    parameter_type paramFromProperties() const;
    parameter_type param0();

    std::string alpha( parameter_type const& mu, ModelMaterial const& mat );
    std::string alphaPrime( parameter_type const& mu, ModelMaterial const& mat );
    std::string alphaRef( ModelMaterial const& mat );
    std::string alphaPrimeRef( ModelMaterial const& mat );
    sparse_matrix_ptrtype assembleForMDEIM( parameter_type const& mu, int const& tag ) override;
    vector_ptrtype assembleForDEIM( parameter_type const& mu, int const& tag ) override;
    void initModel() override;
    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir ) override;

    void decomposition();

    beta_vector_type computeBetaInitialGuess( parameter_type const& mu ) override;
    beta_type computeBetaQm( element_type const& T, parameter_type const& mu ) override;
    beta_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu ) override;
    beta_type computeBetaQm( parameter_type const& mu ) override;
    beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu, double time=1e30 ) override;
    void fillBetaQm( parameter_type const& mu, vectorN_type betaA, vectorN_type betaF );

    affine_decomposition_type computeAffineDecomposition() override;
    std::vector<std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA() override;
    std::vector<std::vector<element_ptrtype> > computeInitialGuessAffineDecomposition() override;

    element_type solve( parameter_type const& mu ) override;
    value_type
    output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false) override;

    double sigma(std::string mat);
    void computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu );
    map_mat_type const& materials() const { return M_materials; }
    map_mat_type const& materialsWithGeo() const { return M_materialsWithGeo; }
    map_mat_type const& materialsWithoutGeo() const { return M_materialsWithoutGeo; }
    bool isInMaterials(ExpressionStringAtMarker const& ex) const;
    std::string const& dbBasename() const { return M_dbBasename; }
}; // AlphaElectric class

} // namespace Feel

#endif
