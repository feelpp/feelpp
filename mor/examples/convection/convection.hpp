/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#ifndef CONVECTIONCRB_H
#define CONVECTIONCRB_H 1

#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmor/parameterspace.hpp>
#include <feel/feelmor/eim.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <feel/feelmor/modelcrbbase.hpp>
#include <feel/feelmor/reducedbasisspace.hpp>
#include <feel/feelmor/crbtrilinearplugin.hpp>

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

#if !defined( CONVECTION_PARAMETER_SPACE_DIMENSION )
#define CONVECTION_PARAMETER_SPACE_DIMENSION 2
#endif

#define A_OUT if(Environment::isMasterRank())std::cout


class FEELPP_EXPORT FunctionSpaceDefinition
{
public :
    static const bool is_time_dependent = false;
    static const bool is_linear = false;
    static const uint16_type Order = 1;

    static const int Order_s = CONVECTION_ORDER_U;
    static const int Order_p = CONVECTION_ORDER_P;
    static const int Order_t = CONVECTION_ORDER_T;

    // Definitions pour mesh
    typedef Simplex<CONVECTION_DIM> entity_type;
    typedef Mesh<entity_type> mesh_type;

    // space and associated elements definitions
    typedef Lagrange<Order_s, Vectorial,Continuous,PointSetFekete> basis_u_type; // velocity space
    typedef Lagrange<Order_p, Scalar,Continuous,PointSetFekete> basis_p_type; // pressure space
    typedef Lagrange<0, Scalar> basis_l_type; // multipliers for pressure space
    typedef Lagrange<Order_t, Scalar,Continuous,PointSetFekete> basis_t_type; // temperature space

    typedef bases< basis_u_type , basis_p_type , basis_l_type, basis_t_type> basis_type;

    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef typename space_type::element_type element_type;

    typedef bases< Lagrange<Order_p, Scalar> > single_basis_type;
    typedef FunctionSpace<mesh_type, single_basis_type> mono_space_type;
};



class FEELPP_EXPORT ConvectionCrb :
    public ModelCrbBase< ParameterSpace<>,
                         FunctionSpaceDefinition >
{
public:
    typedef double value_type;

    //@{ /// Model types
    typedef ModelCrbBase< ParameterSpace<>,
                          FunctionSpaceDefinition > super_type;

    static const uint16_type Order = 1;
    static const bool is_time_dependent = false;
    typedef ConvectionCrb self_type; //@}

    //@{ /// Mesh Setting
    typedef Simplex<CONVECTION_DIM> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype; //@}

    //@{ /// Backend and Matrix
    typedef Backend<double> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef std::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;//@}

    typedef typename FunctionSpaceDefinition::space_type space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;

    typedef typename space_type::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;
    typedef typename element_type:: sub_element<0>::type element_u_type;
    typedef typename element_type:: sub_element<1>::type element_p_type;
    typedef typename element_type:: sub_element<2>::type element_l_type;
    typedef typename element_type:: sub_element<3>::type element_t_type; //@}


    //@{ /// Parameters space
    typedef ParameterSpace<> parameterspace_type;
    typedef std::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;
    typedef parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef std::vector< std::vector< double > > beta_vector_type;
    typedef parameterspace_type::sampling_type sampling_type;
    typedef parameterspace_type::sampling_ptrtype sampling_ptrtype; //@}

    typedef std::shared_ptr< Exporter<mesh_type> > exporter_ptrtype;

    typedef boost::tuple<
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<std::vector<vector_ptrtype> > >
        > affine_decomposition_type;

    typedef Eigen::MatrixXd matrixN_type;


    /// Constructors
    ConvectionCrb();

    /// Initialize the model
    void initModel();

    int Qa() const; ///< number of terms in AD
    int QaTri() const; ///< number of terms in tri AD
    int Ql( int l ) const;
    int mMaxA( int q );
    int mMaxF( int output_index, int q );
    int Nl() const; ///< number of outputs

    static std::string name();

    affine_decomposition_type computeAffineDecomposition()
    {
        return boost::make_tuple( M_Aqm, M_Fqm );
    }

    using super_type::computeBetaQm;
    betaqm_type computeBetaQm( parameter_type const& mu ) ;
    betaqm_type computeBetaQm(element_type const& T,  parameter_type const& mu )
    {
        return computeBetaQm( mu );
    }

    void update( parameter_type const& mu );

    /// \brief solve the model for parameter \p mu
    void solve( parameter_type const& mu, element_ptrtype& T );
    element_type solve( parameter_type const& mu );

    /// solve \f$ M u = f \f$
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f );

    /// Scalar product : <X,Y>
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y );
    double scalarProduct( vector_type const& x, vector_type const& y );

    /// run using specific interface for OpenTURNS
    void run( const double * X, unsigned long N, double * Y, unsigned long P );

    /// Evaluate the value of a specific output for a given parameter
    value_type output( int output_index,
                       parameter_type const& mu ,
                       element_type& unknown, bool need_to_solve=false);

    /// \name Accessors
    //@{
    /// \return the energy matrix
    sparse_matrix_ptrtype energyMatrix() { return M; }
    /// \return the fem function space
    space_ptrtype functionSpace() { return Xh; }
    //@}

    void updateJ( const vector_ptrtype& X, sparse_matrix_ptrtype& J );
    void updateR( const vector_ptrtype& X, vector_ptrtype& R );
    sparse_matrix_ptrtype computeTrilinearForm( const element_type& X );
    sparse_matrix_ptrtype jacobian( const element_type& X );
    vector_ptrtype residual( const element_type& X );

private:
    bool M_psiT;
    double M_delta, M_rez;
    element_ptrtype pT;

    // Timers
    std::map<std::string,std::pair<boost::timer,double> > timers;

    backend_ptrtype M_backend;

    beta_vector_type M_betaAqm;
    std::vector<beta_vector_type> M_betaFqm;

    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm, M_Mqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;

    sparse_matrix_ptrtype D, M, M_A_tril, M_L, M_D;
    vector_ptrtype F;

};






#endif // CONVECTIONCRB_H
