/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Samuel Quinodoz
             Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-02-25

  Copyright (C) 2007 Samuel Quinodoz
  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
#ifndef __Convection_crb_H
#define __Convection_crb_H 1

/**
   \file convection_crb.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \author Elisa Schenone
 \date 2012-08-13
 */
#include <feel/options.hpp>
//#include <feel/feelcore/applicationxml.hpp>
#include <feel/feelcore/application.hpp>

// (non)linear algebra backend
#include <feel/feelalg/backend.hpp>

// quadrature rules
#include <feel/feelpoly/im.hpp>

// function space
#include <feel/feeldiscr/functionspace.hpp>

// linear operators
#include <feel/feeldiscr/operatorlinear.hpp>

// exporter
#include <feel/feelfilters/exporter.hpp>

#include <feel/feelcrb/parameterspace.hpp>


// use the Feel namespace
using namespace Feel;
using namespace Feel::vf;

#if !defined( CONVECTION_DIM )
#define CONVECTION_DIM 2
#endif
#if !defined( CONVECTION_ORDER_U )
#define CONVECTION_ORDER_U 2
#endif
#if !defined( CONVECTION_ORDER_P )
#define CONVECTION_ORDER_P 1
#endif
#if !defined( CONVECTION_ORDER_T )
#define CONVECTION_ORDER_T 2
#endif
#if !defined( CRB_SOLVER )
#define CRB_SOLVER 0
#endif

/**
 * \class Convection_crb
 * The class derives from the Application class
 * the template arguments are :
 * \tparam Order_s velocity polynomial order
 * \tparam Order_t temperature polynomial order
 * \tparam Order_p pressure polynomial order
 */
//template< int Order_s, int Order_p, int Order_t >
class Convection_crb 
{
public:

    static const uint16_type Order = 1;
    static const uint16_type ParameterSpaceDimension = 2;
    static const bool is_time_dependent = false;

    //typedef Convection<Order_s, Order_p, Order_t> self_type;
    static const int Order_s = CONVECTION_ORDER_U;
    static const int Order_p = CONVECTION_ORDER_P;
    static const int Order_t = CONVECTION_ORDER_T;
    typedef Convection_crb self_type;

    // Definitions pour mesh
    typedef Simplex<CONVECTION_DIM> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef backend_type::vector_type vector_type;

    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;
    
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef boost::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;


    // space and associated elements definitions
    typedef Lagrange<Order_s, Vectorial,Continuous,PointSetFekete> basis_u_type; // velocity space
    typedef Lagrange<Order_p, Scalar,Continuous,PointSetFekete> basis_p_type; // pressure space
    typedef Lagrange<Order_t, Scalar,Continuous,PointSetFekete> basis_t_type; // temperature space
    
    /* parameter space */
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;
    typedef parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef parameterspace_type::sampling_type sampling_type;
    typedef parameterspace_type::sampling_ptrtype sampling_ptrtype;
    typedef Eigen::VectorXd theta_vector_type;

#if defined( FEELPP_USE_LM )
    typedef Lagrange<0, Scalar> basis_l_type; // multipliers for pressure space
    typedef fusion::vector< basis_u_type , basis_p_type , basis_t_type,basis_l_type> basis_type;
#else
    typedef fusion::vector< basis_u_type , basis_p_type , basis_t_type> basis_type;
#endif

    //! numerical type is double
    typedef double value_type;

    typedef FunctionSpace<mesh_type, basis_type> space_type;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename element_type:: sub_element<0>::type element_0_type;
    typedef typename element_type:: sub_element<1>::type element_1_type;
    typedef typename element_type:: sub_element<2>::type element_2_type;
#if defined( FEELPP_USE_LM )
    typedef typename element_type:: sub_element<3>::type element_3_type;
#endif

    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;

    typedef boost::shared_ptr<element_type> element_ptrtype;
    
    typedef OperatorLinear<space_type,space_type> oplin_type;
    typedef boost::shared_ptr<oplin_type> oplin_ptrtype;
    typedef FsFunctionalLinear<space_type> funlin_type;
    typedef boost::shared_ptr<funlin_type> funlin_ptrtype;

    // Definition pour les exportations
    typedef Exporter<mesh_type> export_type;

    typedef boost::tuple<std::vector<sparse_matrix_ptrtype>, std::vector<std::vector<vector_ptrtype>>> affine_decomposition_type;
    
    // Constructeur
    Convection_crb( );
    Convection_crb( po::variables_map const& vm );

    // generate the mesh
    Feel::gmsh_ptrtype createMesh();
    
    // Functions usefull for crb resolution :
    
    void init();
    
    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
    {
        return M_Dmu;
    };

    affine_decomposition_type computeAffineDecomposition()
    {
        return boost::make_tuple( M_Aq, M_Fq );
    }

    
    // \return the number of terms in affine decomposition of left hand
    // side bilinear form
    int Qa() const;
    
    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     *
     * \return number of outputs associated to the model
     */
    int Nl() const;
    
    /**
     * \param l the index of output
     * \return number of terms  in affine decomposition of the \p q th output term
     */
    int Ql( int l ) const;

    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<theta_vector_type, std::vector<theta_vector_type> >
    computeThetaq( parameter_type const& mu, double time=0  ) ;
    
    /**
     * \brief return the coefficient vector
     */
    theta_vector_type const& thetaAq() const
    {
        return M_thetaAq;
    };
    
    /**
     * \brief return the coefficient vector
     */
    std::vector<theta_vector_type> const& thetaFq() const
    {
        return M_thetaFq;
    };
    
    /**
     * \brief return the coefficient vector \p q component
     *
     */
    value_type thetaAq( int q ) const
    {
        return M_thetaAq( q );
    };
    
    /**
     * \return the \p q -th term of the \p l -th output
     */
    value_type thetaL( int l, int q ) const
    {
        return M_thetaFq[l]( q );
    };

    void update( parameter_type const& mu );
    
    
    
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );
    
    /**
     * \brief solve the model for parameter \p mu
     * \param mu the model parameter
     * \param T the temperature field
     */
    void solve( parameter_type const& mu, element_ptrtype& T );
    
    /**
     * solve for a given parameter \p mu
     */
    void solve( parameter_type const& mu );
    
    /**
     * solve \f$ M u = f \f$
     */
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f );

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u );
    void exportResults( element_ptrtype& U, int t );
    void exportResults( element_type& U, double t );
        
    /**
     * returns the scalar product of the boost::shared_ptr vector x and
     * boost::shared_ptr vector y
     */
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y );
    
    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_type const& x, vector_type const& y );
    
    /**
     * specific interface for OpenTURNS
     *
     * \param X input vector of size N
     * \param N size of input vector X
     * \param Y input vector of size P
     * \param P size of input vector Y
     */
    void run( const double * X, unsigned long N, double * Y, unsigned long P );
    
    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu );

    sparse_matrix_ptrtype newMatrix() const
    {
        return M_backend->newMatrix( Xh, Xh );
    };
    
    space_ptrtype functionSpace()
    {
        return Xh;
    };
    
    void setMeshSize( double s )
    {
        meshSize = s;
    };

    
    po::options_description const& optionsDescription() const
    {
        return _M_desc;
    }
    
    /**
     * get the variable map
     *
     *
     * @return the variable map
     */
    po::variables_map const& vm() const
    {
        return M_vm;
    }
    

    
    sparse_matrix_ptrtype innerProduct()
    {
        return M;
    }

    
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J );
//    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R, parameter_type const& mu );
    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R );

private:

    
    po::options_description _M_desc;

    po::variables_map M_vm;
    backend_ptrtype M_backend;

    space_ptrtype Xh;

    oplin_ptrtype M_oplin;
    funlin_ptrtype M_lf;

    sparse_matrix_ptrtype M_L;
    sparse_matrix_ptrtype M_D;
    vector_ptrtype F;

    sparse_matrix_ptrtype D,M;

    //pas de temps
    //value_type dt;

    // Exporters
    boost::shared_ptr<export_type> exporter;

    // Timers
    std::map<std::string,std::pair<boost::timer,double> > timers;

    std::vector <double> Grashofs;
    double M_current_Grashofs;
    double M_current_Prandtl;
    
    // Variables usefull for crb resolution :
    
    double meshSize;
    
    element_ptrtype pT;

    std::vector<sparse_matrix_ptrtype> M_Aq;
    std::vector<std::vector<vector_ptrtype> > M_Fq;
    
    parameterspace_ptrtype M_Dmu;
    theta_vector_type M_thetaAq;
    std::vector<theta_vector_type> M_thetaFq;
    
    

};
#endif /* __Convection_crb_H */
