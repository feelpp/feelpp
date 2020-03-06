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

#ifndef FEELPP_CRB_MINIMAL_NL_HPP
#define FEELPP_CRB_MINIMAL_NL_HPP

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

class FunctionSpaceDefinition
{
public:
    using value_type = double;
    static const uint16_type Order = 1;
    using convex_type =  Simplex<3>;
    using mesh_type = Mesh<convex_type>;

    using fct_base_type = Lagrange<Order,Scalar>;
    using basis1_type = bases<fct_base_type>;
    using space1_type = FunctionSpace<mesh_type, basis1_type, value_type>;
    using basis_type = bases<fct_base_type, fct_base_type>;
    using space_type = FunctionSpace<mesh_type, basis_type, value_type>;
    using element_type = typename space_type::element_type;

    using basisd_type = bases<Lagrange<0, Scalar, Discontinuous> >;
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

class FEELPP_EXPORT PoissonNL : public ModelCrbBase<ParameterDefinition,
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

    using eim_space_type = typename eim_definition_type::space1_type;
    using eim_space_ptrtype = std::shared_ptr<eim_space_type>;
    using eimd_space_type = typename eim_definition_type::spaced_type;
    using eimd_space_ptrtype = std::shared_ptr<eimd_space_type>;
    using map_eim_type = std::map<std::string, int>;
    using value_type = double;
    using element_type = super_type::element_type;
    using element_ptrtype = std::shared_ptr<element_type>;
    using parameter_type = super_type::parameter_type;
    using mesh_type = super_type::mesh_type;
    using mesh_ptrtype = super_type::mesh_ptrtype;
    using prop_type = ModelProperties;
    using prop_ptrtype = std::shared_ptr<prop_type>;
    using materials_type = ModelMaterials;
    using mat_type = ModelMaterial;
    using map_mat_type = std::map<std::string, mat_type>;

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

    parameter_type M_mu;
    element_type M_u;

    double M_sigmaMax;

    std::vector< std::vector< element_ptrtype > > M_InitialGuess;

public:
    static po::options_description makeOptions();
    static AboutData makeAbout( std::string const& str = "poissonnlcrbmodel" );
    PoissonNL();
    int Qa();
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
    int indexOfElecMat(std::string const& mat ) const;
    int indexOfTherMat(std::string const& mat ) const;

    functionspace_type::mesh_support_vector_type functionspaceMeshSupport( mesh_ptrtype const& mesh ) const override;
    void initModel() override;
    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir ) override;
    value_type
    output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false) override;
    void assemble() override;
    betaqm_type computeBetaQm( parameter_type const& mu ) override;
    betaqm_type computeBetaQm( element_type const& u, parameter_type const& mu) override;
    betaqm_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu) override;
    void fillBetaQm( parameter_type const& mu, vectorN_type const& betaGrad, std::vector<vectorN_type> const& betaK, std::vector<vectorN_type> const& betaSigma );
    std::vector<std::vector<element_ptrtype> > computeInitialGuessAffineDecomposition() override;
    beta_vector_type computeBetaInitialGuess( parameter_type const& mu) override;
    std::vector< std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA() override;
    beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu, double time=0 ) override;
    element_type solve(parameter_type const& mu) override;
    element_type solveLinear(parameter_type const& mu);

}; // PoissonNL class

} // Feel namespace

#endif
