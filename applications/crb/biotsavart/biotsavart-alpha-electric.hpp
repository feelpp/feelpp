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

#ifndef FEELPP_BIOTSAVART_ALPHA_ELECTRIC_HPP
#define FEELPP_BIOTSAVART_ALPHA_ELECTRIC_HPP

#include <feel/feel.hpp>
#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/eim.hpp>
#include <feel/feelcrb/deim.hpp>
#include <feel/feelcrb/crbmodel.hpp>
#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <electric-alpha.hpp>

namespace Feel
{

FEELPP_EXPORT po::options_description biotsavartOptions()
{
    po::options_description opt("BiotSavart options");
    opt.add_options()
        ( "biotsavart.filename", po::value<std::string>()->default_value("biotsavart.json"),
          "json file")
        ( "biotsavart.repart", po::value<bool>()->default_value( false ),
          "repartition the mesh" )
        ( "biotsavart.compute-fe", po::value<bool>()->default_value(false),
          "compute the finite element version of Biot Savart" )
        ( "biotsavart.compute-offline", po::value<bool>()->default_value(true),
          "compute the offline part of Biot Savart" )
        ( "biotsavart.compute-online", po::value<bool>()->default_value(true),
          "compute the online part of Biot Savart" )
        ( "biotsavart.rebuild-database", po::value<bool>()->default_value(false),
          "rebuild the integrals or not" )
        ( "biotsavart.path-to-database", po::value<std::string>()->default_value("BiotSavart"),
          "path to the database" )
        ( "biotsavart.trainset-deim-size", po::value<int>()->default_value(10),
          "size of the trainset for DEIM of BiotSavart" )
        ;
    opt.add(deimOptions("bs"));
    return opt;
}

class BSParameterDefinition
{
public:
    using parameterspace_type = ParameterSpaceX;
};

class BSFunctionSpaceDefinition
{
public:
    using value_type = double;
    static const uint16_type Order = 1;
    using convex_type =  Simplex<3>;
    using mesh_type = Mesh<convex_type>;

    using fct_base_type = Lagrange<Order,Vectorial>;
    using basis_type = bases<fct_base_type>;
    using space_type = FunctionSpace<mesh_type, basis_type, value_type>;
    using element_type = typename space_type::element_type;
};

template <typename ParameterDefinition, typename FunctionSpaceDefinition >
class BSEimDefinition
{
public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename FunctionSpaceDefinition::space_type space_type;
    // typedef typename FunctionSpaceDefinition::J_space_type J_space_type;

    /* EIM */
    // Scalar continuous //
    typedef EIMFunctionBase<space_type, space_type , parameterspace_type> fun_type;
    // Scalar Discontinuous //
    typedef EIMFunctionBase<space_type, space_type , parameterspace_type> fund_type;

};

template<typename te_rb_model_type>
class FEELPP_EXPORT BiotSavartAlphaElectricCRB
    : public ModelCrbBase<BSParameterDefinition,
                          BSFunctionSpaceDefinition,
                          Linear,
                          BSEimDefinition<BSParameterDefinition,BSFunctionSpaceDefinition> >
{
  public:
    using super_type = ModelCrbBase<BSParameterDefinition,
                                    BSFunctionSpaceDefinition,
                                    Linear,
                                    BSEimDefinition<BSParameterDefinition,
                                                    BSFunctionSpaceDefinition> >;

    using self_type = BiotSavartAlphaElectricCRB<te_rb_model_type>;
    using self_ptrtype = std::shared_ptr<self_type>;

    using te_rb_model_ptrtype = std::shared_ptr<te_rb_model_type>;
    using crb_model_type = CRBModel<te_rb_model_type>;
    using crb_model_ptrtype = std::shared_ptr<crb_model_type>;
    using crb_type = CRB<crb_model_type>;
    using crb_ptrtype = std::shared_ptr<crb_type>;

    using value_type = typename te_rb_model_type::value_type;
    using param_space_type = typename te_rb_model_type::parameterspace_type;
    using parameter_type = typename param_space_type::element_type;
    using vectorN_type = Eigen::VectorXd;
    using eigen_vector_type = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using beta_vector_type = typename te_rb_model_type::beta_vector_type;

    using mesh_type = typename super_type::mesh_type;
    using mesh_ptrtype = typename super_type::mesh_ptrtype;
    using exporter_type = Exporter<mesh_type>;
    using exporter_ptrtype = std::shared_ptr<exporter_type>;

    using space_type = typename super_type::space_type;
    using space_ptrtype = typename super_type::space_ptrtype;
    using element_type = typename super_type::element_type;
    using element_ptrtype = typename super_type::element_ptrtype;

    using cond_space_type = typename te_rb_model_type::space_type;
    using cond_space_ptrtype = typename te_rb_model_type::space_ptrtype;
    using cond_element_type = typename te_rb_model_type::element_type;
    using cond_element_ptrtype = typename te_rb_model_type::element_ptrtype;

    using current_element_type = typename te_rb_model_type::current_element_type;
    using current_space_type = typename te_rb_model_type::current_space_type;

    using dof_point_type = boost::tuple<node_type, size_type, uint16_type >;
    using dof_points_type = typename std::vector<dof_point_type>;

    using deim_type = DEIM<self_type>;
    using deim_ptrtype = std::shared_ptr<deim_type>;
    using vector_ptrtype = typename super_type::vector_ptrtype;

    using prop_type = ModelProperties;
    using prop_ptrtype = std::shared_ptr<prop_type>;
    using mat_type = ModelMaterial;
    using map_mat_type = std::map<std::string,mat_type>;

  public:
    static self_ptrtype New(crb::stage stage = crb::stage::online) { return std::make_shared<self_type>(stage); }
    BiotSavartAlphaElectricCRB(crb::stage stage = crb::stage::online);

    void initModel();
    // setupSpecificityModel
    void setupCRB(crb_ptrtype crb);

    //!
    //!
    //!
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);

    parameter_type paramFromOption();
    parameter_type param0();
    parameter_type paramFromProperties() const;
    void runBS();
    void setupCommunicatorsBS();
    vector_ptrtype assembleForDEIM( parameter_type const& mu, int const& tag );
    vectorN_type beta( parameter_type& mu, int M = -1 ) { return M_deim->beta(mu, M); }
    std::vector<vector_ptrtype> q() { return M_deim->q(); }
    void online( parameter_type & mu, int M = -1 );
    void expand( int N = -1 );
    void computeFE( parameter_type & mu );
    std::vector<double> computeErrors();
    void exportResults( parameter_type const& mu );
    double homogeneity( element_type& B );

    parameter_type newParameter() const { return M_teCrbModel->newParameter(); }
    void setParameter( parameter_type& mu ) { M_mu = mu; }
    int nbParameters() const { return M_crbModel->parameterSpace()->dimension();  }
    auto parameterSpace() const { return M_crbModel->parameterSpace(); }
    int dimension() const { return M_deim->size(); }
    int crbDimension() const { return M_crb->dimension(); }
    cond_element_type potential() const { return M_V; }
    cond_element_type potentialFE() const { return M_VFe; }
    element_type magneticFlux() const { return M_B; }
    element_type magneticFluxFE() const { return M_BFe; }
    mesh_ptrtype mesh() const { return M_mesh; }
    cond_space_ptrtype spaceCond() const { return M_XhCond; }
    space_ptrtype spaceMgn() const { return this->Xh; }
    cond_element_type alpha( parameter_type const& mu );

protected:
    te_rb_model_ptrtype M_teCrbModel;
    crb_model_ptrtype M_crbModel;
    crb_ptrtype M_crb;

    deim_ptrtype M_deim;

    std::vector<Eigen::MatrixXd> M_intMND;

    prop_ptrtype M_modelProps;
    map_mat_type M_materials;

    mesh_ptrtype M_mesh;
    cond_space_ptrtype M_XhCond;
    cond_element_type M_V;
    cond_element_type M_VFe;
    element_type M_B;
    element_type M_BFe;
    current_element_type M_j;
    vectorN_type M_uN;
    eigen_vector_type M_betaMu;
    parameter_type M_mu;
    int M_N;
    int M_M;

    std::vector< mpi::communicator > M_commsC1M;
    std::map<int, dof_points_type> M_dofMgn;

    std::string M_propertyPath;
    bool M_repart;
    bool M_computeFe;
    bool M_computeOffline;
    bool M_computeOnline;
    bool M_rebuildDb;
    std::string M_pathToDb;
    int M_trainsetDeimSize;
}; // class BiotSavartAlphaElectroCRB

#if !defined(FEELPP_INSTANTIATE_BIOTSAVARTALPHAELECTRIC_ELECTRIC)
extern template class FEELPP_EXPORT BiotSavartAlphaElectricCRB<AlphaElectric>;
#endif
} // namespace Feel


#endif
