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

#ifndef FEELPP_BIOTSAVART_HPP
#define FEELPP_BIOTSAVART_HPP

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/eim.hpp>
#include <feel/feelcrb/crbmodel.hpp>
#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <thermoelectric-linear.hpp>

namespace Feel
{

#if 0
/*
 * base class for thermoelectric model for biotsavart
 * need also:
 * - a constructor taking a mesh
 * - a function makeThermoElectricCRBOptions() returning the options
 */
class ThermoElectricBase
{
public:
    virtual int mMaxSigma() = 0;
    virtual auto eimSigmaQ(int m) = 0;
    virtual Eigen::VectorXd eimSigmaBeta( ParameterSpace<5> mu ) = 0;
    virtual template<typename space_type> computeTruthCurrentDensity( ParameterSpace<5> mu ) = 0;
};
#endif

FEELPP_EXPORT po::options_description biotsavartOptions()
{
    po::options_description opt("BiotSavart options");
    opt.add_options()
        ( "biotsavart.conductor", po::value<std::vector<std::string> >()->multitoken(),
          "marker for the conductor" )
        ( "biotsavart.mgn", po::value<std::vector<std::string> >()->multitoken(),
          "marker for the magnetic part" )
        ( "biotsavart.repart", po::value<bool>()->default_value( false ),
          "repartition the mesh" )
        ( "biotsavart.run-mode", po::value<int>()->default_value( 2 ),
          "execution mode: pfem=0, scm=1, crb=2, scm_online=3, crb_online=4" )
        ( "biotsavart.param", po::value<std::vector<double> >()->multitoken(),
          "parameter to evaluate" )
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
        ( "biotsavart.crb-dimension", po::value<int>()->default_value(-1),
          "number of reduced basis to use" )
        ( "biotsavart.eim-dimension", po::value<int>()->default_value(-1),
          "number of eim basis to use" )
        ;
    return opt;
}

template<typename te_rb_model_type>
class FEELPP_EXPORT BiotSavartCRB
    : public ModelCrbBase<typename te_rb_model_type::parameter_space_type,
                          typename te_rb_model_type::function_space_type,
                          NonLinear, // BiotSavart ??
                          typename te_rb_model_type::eim_definition_type >
{
  public:
    using super_type = ModelCrbBase<typename te_rb_model_type::parameter_space_type,
                                    typename te_rb_model_type::function_space_type,
                                    NonLinear,
                                    typename te_rb_model_type::eim_definition_type >;

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

    using mesh_type = typename te_rb_model_type::mesh_type;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using exporter_type = Exporter<mesh_type>;
    using exporter_ptrtype = std::shared_ptr<exporter_type>;

    using space_type = typename te_rb_model_type::space_type;
    using space_ptrtype = std::shared_ptr<space_type>;
    using element_type = typename space_type::element_type;
    using element_ptrtype = std::shared_ptr<element_type>;
    using V_space_type = typename space_type::template sub_functionspace<0>::type;
    using V_space_ptrtype = std::shared_ptr<V_space_type>;
    using V_element_type = typename V_space_type::element_type;
    using V_element_ptrtype = std::shared_ptr<V_element_type>;
    using T_space_type = typename space_type::template sub_functionspace<1>::type;
    using T_space_ptrtype = std::shared_ptr<T_space_type>;
    using T_element_type = typename T_space_type::element_type;
    using T_element_ptrtype = std::shared_ptr<T_element_type>;

    using vec_fct_type = Lagrange<1, Vectorial>;
    using vec_basis_type = bases<vec_fct_type>;
    using vec_space_type = FunctionSpace<mesh_type, vec_basis_type>;
    using vec_space_ptrtype = std::shared_ptr<vec_space_type>;
    using vec_element_type = typename vec_space_type::element_type;
    using vec_element_ptrtype = std::shared_ptr<vec_element_type>;

    // using disc_vec_fct_type = Lagrange<0, Vectorial, Discontinuous>;
    // using disc_vec_basis_type = bases<disc_vec_fct_type>;
    // using disc_vec_space_type = FunctionSpace<mesh_type, disc_vec_basis_type>;
    // using disc_vec_space_ptrtype = std::shared_ptr<disc_vec_space_type>;
    // using disc_vec_element_type = typename disc_vec_space_type::element_type;
    // using disc_vec_element_ptrtype = std::shared_ptr<disc_vec_element_type>;
    using current_element_type = typename te_rb_model_type::current_element_type;
    using current_space_type = typename te_rb_model_type::current_space_type;

    using dof_point_type = boost::tuple<node_type, size_type, uint16_type >;
    using dof_points_type = typename std::vector<dof_point_type>;


public:
    BiotSavartCRB();

    void initModel();
    // setupSpecificityModel

    //!
    //! 
    //!
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);

    parameter_type paramFromOption();
    void runBS();
    void offline();
    void setupCommunicatorsBS();
    void setDimensions();
    void loadIntegrals();
    void saveIntegrals();
    void computeIntegrals(int M = 0, int N = 0);
    void online( parameter_type & mu );
    void computeUn( parameter_type &mu, int N);
    void computeB( vectorN_type & uN );
    void computeFE( parameter_type & mu );
    void exportResults();

    parameter_type newParameter() { return M_teCrbModel->newParameter(); }
    int nbParameters() const { return M_crbModel->parameterSpace()->dimension();  }
    auto parameterSpace() const { return M_crbModel->parameterSpace(); }
    int nMax() const { return M_N; }
    vectorN_type uN() const { return M_uN; }
    element_type potentialTemperature() const { return M_VT; }
    element_type potentialTemperatureFE() const { return M_VTFe; }
    vec_element_type magneticFlux() const { return M_B; }
    vec_element_type magneticFluxFE() const { return M_BFe; }
    mesh_ptrtype meshCond() const { return M_meshCond; }
    mesh_ptrtype meshMgn() const { return M_meshMgn; }

protected:
    CRBModelMode M_mode;
    te_rb_model_ptrtype M_teCrbModel;
    crb_model_ptrtype M_crbModel;
    crb_ptrtype M_crb;

    std::vector<Eigen::MatrixXd> M_intMND;

    mesh_ptrtype M_meshCond;
    mesh_ptrtype M_meshMgn;
    space_ptrtype M_Xh;
    vec_space_ptrtype M_XhMgn;
    element_type M_VT;
    element_type M_VTFe;
    vec_element_type M_B;
    current_element_type M_j;
    vec_element_type M_BFe;
    vectorN_type M_uN;
    eigen_vector_type M_betaMu;
    parameter_type M_mu;
    int M_N;
    int M_M;

    std::vector< mpi::communicator > M_commsC1M;
    std::map<int, dof_points_type> M_dofMgn;

}; // class BiotSavartCRB

#if !defined(FEELPP_INSTANTIATE_BIOTSAVART_THERMOELECTRIC)
extern template class FEELPP_EXPORT BiotSavartCRB<ThermoElectric>;
#endif
} // namespace Feel


#endif
