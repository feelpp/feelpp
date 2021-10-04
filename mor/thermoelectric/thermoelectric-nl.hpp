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

#ifndef FEELPP_CRB_THERMOELECTRIC_NL_HPP
#define FEELPP_CRB_THERMOELECTRIC_NL_HPP

#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>

namespace Feel
{

class ParameterDefinition
{
public:
    using parameterspace_type = ParameterSpaceX;
};

template<int Dim, int Order, int G_Order>
class FunctionSpaceDefinition
{
public:
    using value_type = double;
    using convex_type = Simplex<Dim, G_Order>;
    using mesh_type = Mesh<convex_type>;

    using fct_base_type = Lagrange<Order,Scalar>;
    using basis1_type = bases<fct_base_type>;
    using space1_type = FunctionSpace<mesh_type, basis1_type, value_type>;
    using basis_type = bases<fct_base_type, fct_base_type>;
    using space_type = FunctionSpace<mesh_type, basis_type, value_type>;
    using element_type = typename space_type::element_type;

    using basisd_type = bases<Lagrange<Order-1, Scalar, Discontinuous> >;
    using spaced_type = FunctionSpace<mesh_type, basisd_type, Discontinuous, value_type>;
};

template <typename ParameterDefinition, typename FunctionSpaceDefinition >
class EimDefinition
{
public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename FunctionSpaceDefinition::space_type space_type;
    typedef typename FunctionSpaceDefinition::space1_type space1_type;
    typedef typename FunctionSpaceDefinition::spaced_type spaced_type;

    /* EIM */
    // Scalar continuous //
    typedef EIMFunctionBase<space1_type, space_type , parameterspace_type> fun_type;
    // Scalar Discontinuous //
    typedef EIMFunctionBase<spaced_type, space_type , parameterspace_type> fund_type;

};

class ThermoElectricNLBase
{
public:
    static po::options_description makeOptions();
    static AboutData makeAbout( std::string const& str = "thermoelectricnlcrbmodel" );
};

template<int Dim, int Order, int G_Order>
class FEELPP_EXPORT ThermoElectricNL
    : public ModelCrbBase<ParameterDefinition,
                          FunctionSpaceDefinition<Dim, Order, G_Order>,
                          NonLinear,
                          EimDefinition<ParameterDefinition,
                                        FunctionSpaceDefinition<Dim, Order, G_Order>> >,
      public ThermoElectricNLBase
{
public:
    using super_type = ModelCrbBase<ParameterDefinition,
                                    FunctionSpaceDefinition<Dim, Order, G_Order>,
                                    NonLinear,
                                    EimDefinition<ParameterDefinition,
                                                  FunctionSpaceDefinition<Dim, Order, G_Order>> >;
    using self_type = ThermoElectricNL<Dim, Order, G_Order>;

    using parameter_space_type = ParameterDefinition;
    using function_space_type = FunctionSpaceDefinition<Dim, Order, G_Order>;
    using eim_definition_type = EimDefinition<ParameterDefinition, FunctionSpaceDefinition<Dim, Order, G_Order>>;

    using eim_space_type = typename eim_definition_type::space1_type;
    using eim_space_ptrtype = std::shared_ptr<eim_space_type>;
    using q_sigma_element_type = typename eim_space_type::element_type;
    using eimd_space_type = typename eim_definition_type::spaced_type;
    using eimd_space_ptrtype = std::shared_ptr<eimd_space_type>;
    using map_eim_type = std::map<std::string, int>;
    using value_type = double;
    using element_type = typename super_type::element_type;
    using element_ptrtype = std::shared_ptr<element_type>;
    using parameter_type = typename super_type::parameter_type;
    using mesh_type = typename super_type::mesh_type;
    using mesh_ptrtype = typename super_type::mesh_ptrtype;
    using prop_type = ModelProperties;
    using prop_ptrtype = std::shared_ptr<prop_type>;
    using materials_type = ModelMaterials;
    using mat_type = ModelMaterial;
    using map_mat_type = std::map<std::string, mat_type>;

    using functionspace_type = typename super_type::functionspace_type;
    using betaqm_type = typename super_type::betaqm_type;
    using beta_vector_type = typename super_type::beta_vector_type;
    using vectorN_type = typename super_type::vectorN_type;
    using sparse_matrix_ptrtype = typename super_type::sparse_matrix_ptrtype;
    using vector_ptrtype = typename super_type::vector_ptrtype;

    using rb_space_type = typename super_type::rbfunctionspace_type;
    using ctx_electro_type = typename rb_space_type::template sub_rbfunctionspace<0>::type::ctxrbset_type;
    using ctx_thermo_type = typename rb_space_type::template sub_rbfunctionspace<1>::type::ctxrbset_type;

    using mdeim_type = MDEIM<self_type>;
    using mdeim_ptrtype = std::shared_ptr<mdeim_type>;
    using deim_type = DEIM<self_type>;
    using deim_ptrtype = std::shared_ptr<deim_type>;


private:
    mesh_ptrtype M_mesh;
    prop_ptrtype M_modelProps;
    materials_type M_materials;
    map_mat_type M_elecMaterials;
    map_mat_type M_therMaterials;
    int M_nbElecMat;
    int M_nbTherMat;
    int M_nbPotDir;
    int M_nbTempRobin;
    map_eim_type M_elecEimIndex;
    map_eim_type M_therEimIndex;
    int M_eimGradSize;

    parameter_type M_mu;
    element_type M_u;

    double M_sigmaMax;

    std::vector< std::vector< element_ptrtype > > M_InitialGuess;

    ctx_electro_type M_ctxElectro;
    ctx_thermo_type M_ctxThermo;

    std::string M_propertyPath;
    double M_gamma;
    double M_tolerance;
    int M_maxit;
    int M_trainsetEimSize;
    int M_verbose;
    bool M_useDEIM;

public:
    explicit ThermoElectricNL( std::string prefix = "" );
    prop_ptrtype modelProperties() const { return M_modelProps; }
    parameter_type parameterProperties() const;
    int Qa( bool isLinear = false );
    int mMaxA(int q );
    int Nl();
    int Ql( int l );
    int mMaxF( int l, int q);
    int QIntensity( ModelOutput const& out ) const;
    int QAverageTemp( ModelOutput const& out ) const;
    int mMaxIntensity( int q, ModelOutput const& out ) const;
    int mMaxAverageTemp( int q, ModelOutput const& out ) const;
    int QInitialGuess() const override { return 2; }
    int mMaxInitialGuess(int q) const override { return 1; }
    void resize();
    map_mat_type const& elecMaterials() const { return M_elecMaterials; }
    map_mat_type const& therMaterials() const { return M_therMaterials; }

    typename functionspace_type::mesh_support_vector_type functionspaceMeshSupport( mesh_ptrtype const& mesh ) const override;
    void initModel() override;
    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir ) override;
    void initOutputsPoints();
    vectorN_type computeOutputsPointsElectro( vectorN_type const& urb );
    vectorN_type computeOutputsPointsThermo( vectorN_type const& urb );

    void assemble() override;
    void assembleWithEIM();
    void assembleWithDEIM();
    sparse_matrix_ptrtype assembleForMDEIMnl( parameter_type const& mu, element_type const& u, int const& tag ) override;
    vector_ptrtype assembleForDEIMnl( parameter_type const& mu, element_type const& u, int const& tag ) override;
    betaqm_type computeBetaQm( parameter_type const& mu ) override;
    betaqm_type computeBetaQm( element_type const& u, parameter_type const& mu) override;
    betaqm_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu) override;
    void fillBetaQm( parameter_type const& mu, vectorN_type const& betaGrad, std::vector<vectorN_type> const& betaK, std::vector<vectorN_type> const& betaSigma );
    void fillBetaQmWithDEIM( parameter_type const& mu, vectorN_type const& betaDEIM, vectorN_type const& betaMDEIM );
    std::vector<std::vector<element_ptrtype> > computeInitialGuessAffineDecomposition() override;
    beta_vector_type computeBetaInitialGuess( parameter_type const& mu) override;
    std::vector< std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA() override;
    beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu, double time=0 ) override;

    element_type solve(parameter_type const& mu) override;
    element_type solveLinear(parameter_type const& mu);
    value_type
    output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false) override;

    int mMaxSigma( std::string mat);
    q_sigma_element_type eimSigmaQ(std::string mat, int m);
    vectorN_type eimSigmaBeta(std::string mat, parameter_type const& mu, vectorN_type const& vtN);

}; // ThermoElectricNL class

} // Feel namespace

#endif
