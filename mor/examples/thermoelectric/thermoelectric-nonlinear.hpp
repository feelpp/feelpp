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

#ifndef FEELPP_THERMOELECTRIC_NONLINEAR_HPP
#define FEELPP_THERMOELECTRIC_NONLINEAR_HPP

#include <feel/feelmor/modelcrbbase.hpp>
#include <feel/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>

namespace Feel
{

FEELPP_EXPORT po::options_description
makeThermoElectricOptions()
{
    po::options_description options( "Thermoelectric" );
    options.add_options()
        ( "thermoelectric.basename", Feel::po::value<std::string>()->default_value("thermoelectric_linear"),
          "name of the database" )
        ( "thermoelectric.filename", Feel::po::value<std::string>()->default_value("thermoelectric.json"),
          "json file containing application parameters and boundary conditions")
        ( "thermoelectric.penal-dir", po::value<double>()->default_value( 1e4 ), "penalisation term" )
        ( "thermoelectric.trainset-eim-size", po::value<int>()->default_value(40), "size of the eim trainset" )
        ( "thermoelectric.export-FE", po::value<bool>()->default_value(true), "export FE solution" )
        ( "thermoelectric.picard.maxit", po::value<int>()->default_value(5), "maximum number of iterations for Picard" )
        ( "thermoelectric.picard.tol", po::value<double>()->default_value(1e-8), "tolerance for Picard" )
        ;
    options.add(backend_options("thermo-electro") );
    return options;
}

FEELPP_EXPORT AboutData
makeThermoElectricAbout( std::string const& str = "thermoelectriccrbmodel" )
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
    using basis_type = bases<fct_base_type, fct_base_type>;
    using space_type = FunctionSpace<mesh_type, basis_type, value_type>;
    using element_type = typename space_type::element_type;

    using V_basis_type = bases<fct_base_type>;
    using V_space_type = FunctionSpace<mesh_type, V_basis_type, value_type>;

    // using fct_eim_base_type = Lagrange<Order+2, Scalar>;
    // using basis_eim_type = bases<fct_eim_base_type>;
    // using space_eim_type = FunctionSpace<mesh_type, basis_eim_type, value_type>;

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
    // using space_eim_type = typename FunctionSpaceDefinition::space_eim_type;
    // using element_eim_type = typename space_eim_type::element_type;
    using V_space_type = typename FunctionSpaceDefinition::V_space_type;

    /* EIM */
    // Scalar continuous //
    typedef EIMFunctionBase<V_space_type, space_type , parameterspace_type> fun_type;
    // typedef EIMFunctionBase<space_eim_type, space_type , parameterspace_type> fun_type;
    // Scalar Discontinuous //
    typedef EIMFunctionBase<J_space_type, space_type , parameterspace_type> fund_type;

};

class FEELPP_EXPORT ThermoElectric : public ModelCrbBase<ParameterDefinition,
                                                         FunctionSpaceDefinition,
                                                         NonLinear,
                                                         EimDefinition<ParameterDefinition,
                                                                       FunctionSpaceDefinition> >
{
public:
    using super_type = ModelCrbBase<ParameterDefinition,
                                    FunctionSpaceDefinition,
                                    NonLinear,
                                    EimDefinition<ParameterDefinition,
                                                  FunctionSpaceDefinition> >;
    using parameter_space_type = ParameterDefinition;
    using function_space_type = FunctionSpaceDefinition;
    using eim_definition_type = EimDefinition<ParameterDefinition, FunctionSpaceDefinition>;

    using value_type = double;

    using mesh_type = super_type::mesh_type;
    using mesh_ptrtype = super_type::mesh_ptrtype;

    using space_type = super_type::space_type;
    using element_type = super_type::element_type;
    using element_ptrtype = super_type::element_ptrtype;

    using J_space_type = function_space_type::J_space_type;
    using J_space_ptrtype = std::shared_ptr<J_space_type>;
    using q_sigma_space_type = eim_definition_type::V_space_type;
    using q_sigma_space_ptrtype = std::shared_ptr<q_sigma_space_type>;
    using q_sigma_element_type = q_sigma_space_type::element_type;
    using V_view_type = typename element_type::template sub_element_type<0>;
    using V_view_ptrtype = typename element_type::template sub_element_ptrtype<0>;
    using T_view_ptrtype = typename element_type::template sub_element_ptrtype<1>;

    using current_space_type = FunctionSpace<mesh_type, bases<Lagrange<1, Vectorial> > >;
    using current_element_type = typename current_space_type::element_type;

    using prop_type = ModelProperties;
    using prop_ptrtype = std::shared_ptr<prop_type>;
    using mat_type = ModelMaterial;
    using map_mat_type = std::map<std::string, mat_type>;

    using parameter_type = super_type::parameter_type;
    using vectorN_type = super_type::vectorN_type;
    using beta_vector_type = typename super_type::beta_vector_type;
    using beta_type = boost::tuple<beta_vector_type, std::vector<beta_vector_type> >;
    using affine_decomposition_type = typename super_type::affine_decomposition_type;

    using sparse_matrix_ptrtype = typename super_type::sparse_matrix_ptrtype;

private:
    mesh_ptrtype M_mesh;
    prop_ptrtype M_modelProps;
    map_mat_type M_materials;
    map_mat_type M_elecMaterials;
    map_mat_type M_therMaterials;

    std::vector< std::vector< element_ptrtype > > M_InitialGuess;

    element_ptrtype M_VT;
    V_view_ptrtype M_V;
    T_view_ptrtype M_T;
    parameter_type M_mu;

    J_space_ptrtype M_Jh;
    // q_sigma_space_ptrtype M_Th;

    sparse_matrix_ptrtype M_M;

    std::string M_propertyPath;
    double M_penalDir;
    double M_picardTol;
    int M_picardMaxit;
    int M_trainsetEimSize;
    bool M_exportFE;

public:
    static po::options_description makeOptions( std::string const& prefix="thermoelectric" );
    // Constructors
    explicit ThermoElectric( std::string const& prefix = "thermoelectric" );
    ThermoElectric( mesh_ptrtype mesh, std::string const& prefix = "thermoelectric" );

    int indexOfMatV(std::string mat ) const { return std::distance(M_elecMaterials.begin(),M_elecMaterials.find(mat)); }
    int indexOfMatT(std::string mat ) const { return std::distance(M_therMaterials.begin(),M_therMaterials.find(mat)); }
    // Size of the decomposition
    int Qa() const;
    int Nl() const;
    int Ql( int l ) const;
    int mMaxA( int q ) const;
    int mMaxL( int l, int q ) const;
    int mMaxCompliant( int q ) const;
    int QIntensity( ModelOutput const& out ) const;
    int QAverageTemp( ModelOutput const& out ) const;
    int mMaxIntensity( int q, ModelOutput const& out ) const;
    int mMaxAverageTemp( int q, ModelOutput const& out ) const;
    int QInitialGuess() const override;
    int mMaxInitialGuess( int q ) const override;
    void resizeQm( bool resizeMatrix = true );

    // Parameters
    parameter_type newParameter() { return Dmu->element(); }
    parameter_type paramFromProperties() const;

    // Initialization
    void initModel() override;
    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir ) override;

    // mesh support of functionspace
    functionspace_type::mesh_support_vector_type functionspaceMeshSupport( mesh_ptrtype const& mesh ) const override;

    // Decomposition
    void assemble() override;
    affine_decomposition_type computeAffineDecomposition() override;
    std::vector<std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA() override;
    std::vector<std::vector<element_ptrtype> > computeInitialGuessAffineDecomposition() override;

    // Beta coefficients
    beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu , double time=1e30 ) override;
    beta_vector_type computeBetaInitialGuess( parameter_type const& mu ) override;
    beta_type computeBetaQm( parameter_type const& mu ,  double time , bool only_terms_time_dependent=false) override;
    beta_type computeBetaQm( element_type const& T, parameter_type const& mu ) override;
    beta_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu ) override;
    beta_type computeBetaQm( parameter_type const& mu ) override;
    void fillBetaQm( parameter_type const& mu, vectorN_type betaEimGrad, std::vector<vectorN_type> betaEimsSigma, std::vector<vectorN_type> betaEimsK );

    // FE resolution
    element_type solve( parameter_type const& mu ) override;

    // Scalar product
    double scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y );
    double scalarProduct( vector_type const& x, vector_type const& y );
    sparse_matrix_ptrtype energyMatrix() override;

    // Output
    value_type
    output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false) override;

    // BiotSavart API
    int mMaxJoule();
    int mMaxSigma( std::string const& mat );
    q_sigma_element_type eimSigmaQ( std::string const& mat, int m );
    vectorN_type eimSigmaBeta( std::string const& mat, parameter_type const& mu );
    vectorN_type eimSigmaBeta( std::string const& mat, parameter_type const& mu, element_type const& U );
    vectorN_type eimSigmaBeta( std::string const& mat, parameter_type const& mu, vectorN_type const& Urb );
    void computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu );
    void computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu, element_type& VT );
}; // Thermoelectric class

} // namespace Feel

#endif
