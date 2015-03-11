#ifndef __CRBMODELBASE_H
#define __CRBMODELBASE_H 1

#include <boost/shared_ptr.hpp>
#include<boost/regex.hpp>
#include <vector>

#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feeldiscr/operatorlinearfree.hpp>
#include <feel/feeldiscr/operatorlinearcomposite.hpp>
#include <feel/feeldiscr/fsfunctionallinearfree.hpp>
#include <feel/feeldiscr/fsfunctionallinearcomposite.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>


#include <feel/feelcrb/eim.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcore/pslogger.hpp>

#include "affinedecomposition.hpp"

namespace Feel
{

enum class CRBModelMode
{
    PFEM = 0, SCM = 1, CRB = 2, SCM_ONLINE=3, CRB_ONLINE=4
};

enum {
    /** TimeIndependent */
    TimeIndependent=0,
    /** TimeDependent */
    TimeDependent = 0x1,
    /**  */
    Linear = 0,
    /**  */
    NonLinear = 0x2,
    /** Coercive PDE */
    Coercive = 0,
    /** Inf-Sup PDE */
    InfSup = 0x4,
    /** Use block Structure */
    BlockStruct = 0x8
};

class FunctionSpaceDefinitionBase
{
public :
    /*mesh*/
    typedef Simplex<1> entity_type ;
    typedef Mesh<entity_type > mesh_type ;

    /*basis*/
    typedef bases< Lagrange<1,Scalar> > basis_type ;

    /*space*/
    typedef FunctionSpace<mesh_type , basis_type > space_type ;

    typedef typename space_type::element_type element_type;

    static const bool is_time_dependent=false;
    static const bool is_linear=true;
};

template <typename ParameterDefinition, typename FunctionSpaceDefinition>
class EimDefinitionBase
{

public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename mpl::if_<is_shared_ptr<FunctionSpaceDefinition>,
                              mpl::identity<typename FunctionSpaceDefinition::element_type>,
                              mpl::identity<FunctionSpaceDefinition>>::type::type::space_type space_type;
    //typedef typename FunctionSpaceDefinition::space_type space_type;

    /* EIM */
    typedef EIMFunctionBase<space_type , space_type  , parameterspace_type > fun_type ;
    typedef EIMFunctionBase<space_type , space_type  , parameterspace_type > fund_type ;

};


/**
 * \class CRBModel
 *
 */
template < typename ParameterDefinition=ParameterSpace<1>,
           typename FunctionSpaceDefinition=FunctionSpaceDefinitionBase,
           int _Options=0,
           typename EimDefinition=EimDefinitionBase<ParameterDefinition,FunctionSpaceDefinition>
           >
class CRBModelBase :
        public ModelCrbBaseBase,
        public boost::enable_shared_from_this<CRBModelBase<ParameterDefinition,
                                                           FunctionSpaceDefinition,
                                                           _Options,
                                                           EimDefinition> >
{
public:
    static const uint16_type ParameterSpaceDimension = ParameterDefinition::ParameterSpaceDimension ;
    static const bool is_time_dependent = ((_Options&TimeDependent)==TimeDependent);
    static const bool is_linear = !((_Options&NonLinear)==NonLinear);
    static const bool use_block_structure = (_Options&BlockStruct)==BlockStruct;
    static const int Options = _Options;


    typedef CRBModelBase self_type;
    typedef double value_type;

    typedef ParameterDefinition parameter_definition;
    typedef FunctionSpaceDefinition function_space_definition;
    typedef EimDefinition eim_definition;
    // Functions Space
    typedef typename mpl::if_<is_shared_ptr<FunctionSpaceDefinition>,
                              mpl::identity<typename FunctionSpaceDefinition::element_type>,
                              mpl::identity<FunctionSpaceDefinition>>::type::type::space_type space_type;
    typedef space_type functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef functionspace_ptrtype space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    template <int Row>
    using subelement_type = mpl::if_< mpl::bool_<use_block_structure>,
                                      typename space_type::template sub_functionspace<Row>::type::element_type::element_type,
                                      typename space_type::element_type >;

    // RB Function Space
    typedef ReducedBasisSpace<functionspace_type> rbfunctionspace_type;
    typedef boost::shared_ptr< rbfunctionspace_type > rbfunctionspace_ptrtype;

    // Parameters Space
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;

    // Mesh
    typedef typename functionspace_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    // Backend and Tensors
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    // Eim
    typedef typename EimDefinition::fun_type fun_type;
    typedef typename EimDefinition::fund_type fund_type;
    typedef boost::shared_ptr<fun_type> fun_ptrtype;
    typedef std::vector<fun_ptrtype> funs_type;
    typedef boost::shared_ptr<fund_type> fund_ptrtype;
    typedef std::vector<fund_ptrtype> funsd_type;
    typedef Eigen::VectorXd vectorN_type;

    typedef typename boost::tuple<sparse_matrix_ptrtype,
                                  sparse_matrix_ptrtype,
                                  std::vector<vector_ptrtype>
                                  > offline_merge_type;

    typedef typename boost::tuple< sparse_matrix_ptrtype,
                                   sparse_matrix_ptrtype,
                                   std::vector<vector_ptrtype>
                                   > monolithic_type;
    typedef typename boost::tuple< std::map<int,double>,
                                   std::map<int,double>,
                                   std::vector< std::map<int,double> >
                                   > eim_interpolation_error_type ;

    //! time discretization
    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    typedef OperatorLinearComposite< space_type , space_type > operatorcomposite_type;
    typedef boost::shared_ptr<operatorcomposite_type> operatorcomposite_ptrtype;
    typedef FsFunctionalLinearComposite< space_type > functionalcomposite_type;
    typedef boost::shared_ptr<functionalcomposite_type> functionalcomposite_ptrtype;


    typedef vf::detail::BilinearForm<space_type, space_type, VectorUblas<value_type>> form2_type;


    typedef typename std::vector< std::vector < element_ptrtype > > initial_guess_type;

    typedef std::vector< std::vector< value_type > > beta_vector_type;


    typedef std::vector< std::vector< sparse_matrix_ptrtype >> affine_matrix_type;
    typedef std::vector< std::vector< affine_matrix_type >> block_affine_matrix_type;
    typedef std::vector< std::vector < vector_ptrtype >> affine_vector_type;
    typedef std::vector< affine_vector_type > affine_output_type;
    typedef std::vector< affine_output_type > block_affine_output_type;

    typedef AffineDecompositionVector<self_type> vector_ad_type;
    typedef AffineDecompositionMatrix<self_type> matrix_ad_type;

    typedef boost::tuple<
        std::vector< std::vector<sparse_matrix_ptrtype> >,//Mq
        std::vector< std::vector<sparse_matrix_ptrtype> >,//Aq
        std::vector< std::vector<std::vector<vector_ptrtype> > >//Fq
        > affine_decomposition_type;

    typedef typename boost::tuple< beta_vector_type,
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   > betaqm_type;
    typedef typename boost::tuple< beta_vector_type,
                                   std::vector<beta_vector_type>
                                   > steady_betaqm_type;

    static const int nb_spaces = functionspace_type::nSpaces;

    /*typedef typename mpl::if_< boost::is_same< mpl::int_<nb_spaces>, mpl::int_<2> >,
                               fusion::vector< mpl::int_<0>, mpl::int_<1> >,
                               typename mpl::if_ < boost::is_same< mpl::int_<nb_spaces> , mpl::int_<3> >,
                                                   fusion::vector < mpl::int_<0>, mpl::int_<1>, mpl::int_<2> > ,
                                                   typename mpl::if_< boost::is_same< mpl::int_<nb_spaces>, mpl::int_<4> >,
                                                                      fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3> >,
                                                                      fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3>, mpl::int_<4> >
                                                                      >::type
                                                   >::type
     >::type index_vector_type;*/

    typedef typename mpl::range_c< int, 0, functionspace_type::nSpaces > index_vector_type;




    CRBModelBase() :
        M_is_initialized( false ),
        M_backend( backend() ),
        M_backend_primal( backend( _name="backend-primal") ),
        M_backend_dual( backend( _name="backend-dual") ),
        Dmu( new parameterspace_type ),
        M_Nl( 0 ),
        M_has_eim( false )
        {
            M_backend_l2.resize(1);
            M_backend_l2[0] = backend( _name="backend-l2" );
        }

    virtual void initModel()=0;

    void init( CRBModelMode mode=CRBModelMode::PFEM )
        {
            if ( M_is_initialized )
                return;

            M_mode = mode;
            M_inner_product_matrix.resize(1);

            M_A = matrix_ad_type( this->shared_from_this() );
            M_M = matrix_ad_type( this->shared_from_this() );
            M_F.push_back( vector_ad_type(this->shared_from_this() ) );

            LOG(INFO)<< "Model Initialization";
            initModel();
            initDerived();
            if ( !is_time_dependent )
                initializeMassMatrix();

            countAffineDecompositionTerms();

            checkBdf();

            u = Xh->element();
            v = Xh->element();
            M_is_initialized = true;
        }

    virtual void initDerived()
        {
            bool is_symmetric = boption(_name="crb.use-symmetric-matrix");

            double norm=0;
            if ( hasEim() || !is_symmetric )
            {
                M_inner_product_matrix[0] = mergeLinearA();
                norm = M_inner_product_matrix[0]->l1Norm();
            }
            else if ( norm==0 && is_symmetric )
                M_inner_product_matrix[0] = M_energy_matrix;
        }

    template< int Row=0, int Col=0, typename OpeType, typename BetaType >
    void addMass( OpeType ope, BetaType beta )
        {
            M_M.template createBlock<Row,Col>();
            std::string filename = ( boost::format( "M%1%%2%-" ) %Row %Col ) .str();
            M_M.template get<Row,Col>->addOpe( ope, beta, filename );
        }
    template< int Row=0, int Col=0, typename OpeType, typename BetaType >
    void addLhs( OpeType ope, BetaType beta )
        {
            M_A.template createBlock<Row,Col>();
            std::string filename = ( boost::format( "A%1%%2%-" ) %Row %Col ) .str();
            M_A.template get<Row,Col>()->addOpe( ope, beta, filename );
        }
    template< int Row=0, typename FunType, typename BetaType >
    void addRhs( FunType fun, BetaType beta )
        {
            int output = 0;
            M_Nl = std::max( M_Nl, output+1 );
            M_F[output].template createBlock<Row>();
            std::string filename = ( boost::format( "F%1%-%2%-" ) %Row %output ) .str();
            M_F[output].template get<Row>()->addFun( fun, beta, filename );
        }
    template< int Row=0, typename FunType, typename BetaType >
    void addOutput( FunType fun, BetaType beta, int output=1 )
        {
            M_Nl = std::max( M_Nl, output+1 );
            if ( M_F.size()<=output )
            {
                M_F.resize( output+1 );
                M_F[output] = vector_ad_type(this->shared_from_this() );
            }
            M_F[output].template createBlock<Row>();
            std::string filename = ( boost::format( "F%1%-%2%-" ) %Row %output ) .str();
            M_F[output].template get<Row>()->addFun( fun, beta, filename );
        }

    matrix_ad_type getA()
        {
            return M_A;
        }
    std::vector< vector_ad_type > getF()
        {
            return M_F;
        }



    virtual std::string modelName()
        {
            return "generic-model-name";
        }



    virtual void run( const double * X, unsigned long N, double * Y, unsigned long P )
        {
            if( Environment::worldComm().isMasterRank() )
            {
                std::cout<<"**************************************************"<<std::endl;
                std::cout<<"** You are using the function run whereas       **"<<std::endl;
                std::cout<<"** your model has not implemented this function **"<<std::endl;
                std::cout<<"**************************************************"<<std::endl;
            }
        }

    virtual void setFunctionSpaces( functionspace_ptrtype Vh )
        {
            Xh = Vh;
            XN = rbfunctionspace_type::New( _space=Xh );
            if (Environment::isMasterRank() )
                std::cout << "Number of dof : " << Xh->nDof() << std::endl
                          << "Number of local dof : " << Xh->nLocalDof() << std::endl;
        }

    void initializeMassMatrix()
        {
            M_M.initializeMassMatrix();
        }

    parameterspace_ptrtype parameterSpace()
        {
            return Dmu;
        }
    functionspace_ptrtype functionSpace() const
        {
            return Xh;
        }
    rbfunctionspace_ptrtype rBFunctionSpace()
        {
            return XN;
        }

    bool isSteady()
        {
            return !is_time_dependent || boption(_name="crb.is-model-executed-in-steady-mode") ;
        }
    bool isLinear()
        {
            return is_linear;
        }

    virtual parameter_type refParameter()
        {
            LOG(INFO)<< "You did not specified reference parameter"
                     << "You can do it with the implementation of the function refParameter()"
                     << "Default : reference parameter = minimal value";
            return Dmu->min();
        }


    virtual steady_betaqm_type computeBetaQm( parameter_type const& mu )
        {
            return boost::make_tuple( computeBetaAqm(mu), computeBetaFqm(mu) );
        }

    virtual betaqm_type computeBetaQm( parameter_type const& mu, double time,
                               bool only_time_dependent_terms=false )
        {
            betaqm_type beta_coeff;
            steady_betaqm_type steady_beta_coeff;

            if ( !is_time_dependent )
            {
                steady_beta_coeff = computeBetaQm( mu );
                beta_coeff = boost::make_tuple( computeBetaMqm( mu, time ),
                                                steady_beta_coeff.template get<0>(),
                                                steady_beta_coeff.template get<1>() );
            }
            else
                beta_coeff = boost::make_tuple( computeBetaMqm( mu, time ),
                                                computeBetaAqm( mu, time ),
                                                computeBetaFqm( mu, time ) );
            return beta_coeff;
        }

    virtual steady_betaqm_type computeBetaQm( element_type const& T, parameter_type const& mu )
        {
            return computeBetaQm( T, mu, mpl::bool_<use_block_structure>() );
        }

    steady_betaqm_type computeBetaQm( element_type const& T, parameter_type const& mu, mpl::bool_<false> )
        {
            return boost::make_tuple( computeBetaAqm( T,mu ), computeBetaFqm( T,mu ) );
        }
    steady_betaqm_type computeBetaQm( element_type const& T, parameter_type const& mu, mpl::bool_<true> )
        {
            if (Environment::isMasterRank())
                std::cout << "Warning : call to computeBetaQm with block structure !!!\n";
            steady_betaqm_type dummy;
            return dummy;
        }

    betaqm_type computeBetaQm( vector_ptrtype const& T, parameter_type const& mu ,
                        double time , bool only_time_dependent_terms=false )
        {
            auto solution = functionSpace()->element();
            solution = *T;
            return computeBetaQm( solution , mu , time , only_time_dependent_terms );
        }

    virtual betaqm_type computeBetaQm( element_type const& T, parameter_type const& mu ,
                                       double time , bool only_time_dependent_terms=false )
        {
            betaqm_type beta_coeff;
            steady_betaqm_type steady_beta_coeff;
            beta_vector_type betaM;

            if ( !is_time_dependent )
            {
                steady_beta_coeff = computeBetaQm( T, mu );
                betaM.resize( Qm() );
                for ( int q=0; q<betaM.size(); q++ )
                    betaM[q].push_back( 1 );
                beta_coeff = boost::make_tuple( computeBetaMqm(mu, time ),
                                                steady_beta_coeff.template get<0>(),
                                                steady_beta_coeff.template get<1>() );
            }
            else
                throw std::logic_error( "Since youre model is not linear (and time dependent), you have to implement the function computeBetaQm( u, mu, time ) ");

            return beta_coeff;
        }

    template<int Row=0, int Col=0>
    double& setBetaM( int const& q, int const& m=0 )
        {
            M_M.template createBlock<Row,Col>();
            return M_M.template get<Row,Col>()->setBeta(q,m);
        }
    template<int Row=0, int Col=0>
    beta_vector_type computeBetaMqm( parameter_type const& mu, double time=0 )
        {
            return M_M.template get<Row,Col>()->computeBetaQm( mu, time );
        }
    template<int Row=0, int Col=0>
    beta_vector_type computeBetaMqm( element_type const& T, parameter_type const& mu, double time=0 )
        {
            return M_M.template get<Row,Col>()->computeBetaQm( T, mu, time );
        }
    template<int Row=0, int Col=0>
    double betaM( int const& q, int const& m=0 )
        {
            return M_M.template get<Row,Col>()->beta( q, m );
        }
    template<int Row=0, int Col=0>
    beta_vector_type betaMqm()
        {
            return M_M.template get<Row,Col>()->betaQm();
        }

    template<int Row=0, int Col=0>
    double& setBetaA( int const& q, int const& m=0 )
        {
            M_A.template createBlock<Row,Col>();
            return M_A.template get<Row,Col>()->setBeta(q,m);
        }
    template<int Row=0, int Col=0>
    beta_vector_type computeBetaAqm( parameter_type const& mu, double time=0 )
        {
            return M_A.template get<Row,Col>()->computeBetaQm( mu, time );
        }
    template<int Row=0, int Col=0>
    beta_vector_type computeBetaAqm( element_type const& T, parameter_type const& mu, double time=0 )
        {
            return M_A.template get<Row,Col>()->computeBetaQm( T, mu, time );
        }
    template<int Row=0, int Col=0>
    double betaA( int const& q, int const& m=0 )
        {
            return M_A.template get<Row,Col>()->beta( q, m );
        }
    template<int Row=0, int Col=0>
    beta_vector_type betaAqm()
        {
            return M_A.template get<Row,Col>()->betaQm();
        }

    template<int Row=0>
    double& setBetaF( int const& output, int const& q, int const m )
        {
            M_Nl = std::max( M_Nl, output+1 );
            if ( M_F.size()<=output )
            {
                M_F.resize( output+1 );
                M_F[output] = vector_ad_type(this->shared_from_this() );
            }
            M_F[output].template createBlock<Row>();
            return M_F[output].template get<Row>()->setBeta( q, m );
        }
    template<int Row=0>
    beta_vector_type computeBetaFqm( int output, parameter_type const& mu,
                                     double time=0 )
        {
            return M_F[output].template get<Row>()->computeBetaQm( mu, time );
        }
    template<int Row=0>
    std::vector< beta_vector_type > computeBetaFqm( parameter_type const& mu, double time=0 )
        {
            std::vector< beta_vector_type > betaF;
            for ( int output=0; output<M_Nl; output++ )
                betaF.push_back( M_F[output].template get<Row>()->computeBetaQm( mu, time ) );
            return betaF;
        }
    template<int Row=0>
    std::vector< beta_vector_type > computeBetaFqm( element_type const& T,
                                                    parameter_type const& mu, double time=0 )
        {
            std::vector< beta_vector_type > betaF;
            for ( int output=0; output<M_Nl; output++ )
                betaF.push_back( M_F[output].template get<Row>()->computeBetaQm( T, mu, time ) );
            return betaF;
        }
    template<int Row=0>
    double betaF( int const& output, int const& q, int const& m=0 )
        {
            return M_F[output].template get<Row>()->beta(q,m);
        }
    template<int Row=0>
    std::vector< beta_vector_type > betaFqm()
        {
            std::vector< beta_vector_type > betaF;
            for ( int output=0; output<M_Nl; output++ )
                betaF.push_back( M_F[output].template get<Row>()->betaQm() );
            return betaF;
        }

    virtual affine_decomposition_type computeAffineDecomposition()
        {
            affine_output_type Fqm;
            for ( int output=0; output<M_Nl; output++)
                Fqm.push_back( M_F[output].template get<0>()->compute() );

            return boost::make_tuple( M_M.template get<0,0>()->compute(),
                                      M_A.template get<0,0>()->compute(),
                                      Fqm );
        }
    template<int Row=0, int Col=0>
    sparse_matrix_ptrtype ptrM( int const& q, int const& m=0 )
        {
            M_M.template createBlock<Row,Col>();
            return M_M.template get<Row,Col>()->ptrMat( q, m );
        }

    template<int Row=0, int Col=0>
    sparse_matrix_ptrtype ptrA( int const& q, int const& m=0 )
        {
            M_A.template createBlock<Row,Col>();
            return M_A.template get<Row,Col>()->ptrMat( q, m );
        }
    template<int Row=0 >
    vector_ptrtype ptrF( int const& output, int const& q, int const& m=0 )
        {
            M_Nl = std::max( M_Nl, output+1 );
            if ( M_F.size()<=output )
            {
                M_F.resize( output+1 );
                M_F[output] = vector_ad_type(this->shared_from_this() );
            }
            M_F[output].template createBlock<Row>();
            return M_F[output].template get<Row>()->ptrVec(q,m);
       }



    void countAffineDecompositionTerms()
        {
            M_M.check();
            M_A.check();
            for (int output=0; output<M_Nl; output++ )
                M_F[output].check();
        }
    template<int Row=0,int Col=0>
    size_type Qm() const
        {
            auto M=M_M.template get<Row,Col>();
            if(M)
                return M->Q();
            return 0;
        }
    template<int Row=0,int Col=0>
    size_type Qa() const
        {
            auto A = M_A.template get<Row,Col>();
            if ( A )
                return A->Q();
            return 0;
        }
    template<int Row=0>
    size_type Ql( int output) const
        {
            auto F=M_F[output].template get<Row>();
            if (F)
                return F->Q();
            return 0;
        }
    size_type Nl() const
        {
            return M_Nl;
        }

    template<int Row=0,int Col=0>
    int mMaxM( int q )
        {
            auto M=M_M.template get<Row,Col>();
            if(M)
                return M->mMax(q);
            return 0;
        }
    template<int Row=0,int Col=0>
    int mMaxA( int q )
        {
            auto A=M_A.template get<Row,Col>();
            if(A)
                return A->mMax(q);
            return 0;
        }
    template<int Row=0>
    int mMaxF( int output, int q )
        {
            auto F=M_F[output].template get<Row>();
            if(F)
                return F->mMax(q);
            return 0;
        }

    template<int Row=0,int Col=0>
    sparse_matrix_ptrtype Mqm( uint16_type q, uint16_type m, bool transpose = false ) const
        {
            return M_M.template get<Row,Col>()->compute( q, m, transpose );
        }
    template<int Row=0,int Col=0>
    value_type Mqm( uint16_type q, uint16_type m,
                    element_type const& xi_i, element_type const& xi_j, bool transpose = false ) const
        {
            return Mqm<Row,Col>( q, m, transpose )->energy( xi_j, xi_i );
        }
    template<int Row=0,int Col=0>
    sparse_matrix_ptrtype Aqm( uint16_type q, uint16_type m, bool transpose = false ) const
        {
            return M_A.template get<Row,Col>()->compute( q, m, transpose );
        }
    template<int Row=0,int Col=0, typename Type1, typename Type2>
    value_type Aqm( uint16_type q, uint16_type m,
                    Type1 const& xi_i, Type2 const& xi_j,
                    bool transpose = false ) const
        {
            auto A=M_A.template get<Row,Col>();
            if(A)
                return A->compute( q, m, transpose )->energy( xi_j, xi_i );
            else if ( boption("crb.saddlepoint.transpose") )
                return M_A.template get<Col,Row>()->compute( q, m, !transpose )->energy(xi_j, xi_i );
            return 0;
        }
    template<int Row=0>
    vector_ptrtype Fqm( uint16_type output, uint16_type q, int m=0 ) const
        {
            return M_F[output].template get<Row>()->compute(q,m);
        }
    template<int Row=0, typename Type>
    value_type Fqm( uint16_type output, uint16_type q,  uint16_type m, Type const& xi )
        {
            return inner_product( *Fqm<Row>( output, q, m ), xi );
        }


    void checkBdf()
        {
            if ( is_time_dependent )
                CHECK( M_bdf )<< "You have implemented a transient problem but you did not initialize M_bdf";
        }
    double timeStep()
        {
            if ( isSteady() )
                return 1e30;
            else
                return M_bdf->timeStep();
        }
    double timeFinal()
        {
            if ( isSteady() )
                return 1e30;
            else
                return M_bdf->timeFinal();
        }
    double timeInitial()
        {
            if( isSteady() )
                return 0;
            else
                return M_bdf->timeInitial();
        }
    double timeOrder()
        {
            if( is_time_dependent )
                return M_bdf->timeOrder();
            else
                return 0;
        }


    virtual void adaptMesh( parameter_type const& mu )
        {}
    void partitionMesh( std::string mshfile , std::string target , int dimension, int order )
    {
        int N = Environment::worldComm().globalSize();

        if( Environment::worldComm().isMasterRank() )
            std::cout<<"[ModelCrbBase] generate target file : "<<target<<" from "<<mshfile<<std::endl;

        Gmsh gmsh( dimension,
                   order,
                   Environment::worldComm() );
        gmsh.setNumberOfPartitions( N );
        gmsh.rebuildPartitionMsh( mshfile /*mesh with 1 partition*/, target /*mesh with N partitions*/ );
    }

    void setHasEim( bool has_eim )
        {
            M_has_eim = has_eim;
        }
    bool hasEim()
        {
            return M_has_eim;
        }
    bool hasEimError()
        {
            auto all_errors = eimInterpolationErrorEstimation();
            auto errorsM = all_errors.template get<0>();
            auto errorsA = all_errors.template get<1>();
            auto errorsF = all_errors.template get<2>();
            int sizeM = errorsM.size();
            int sizeA = errorsA.size();
            int sizeF = errorsF.size();
            bool b=false;
            return sizeM > 0 || sizeA > 0 || sizeF > 0 ;
        }
    virtual eim_interpolation_error_type eimInterpolationErrorEstimation( parameter_type const& mu ,
                                                                          vectorN_type const& uN )
        {
            return boost::make_tuple( M_eim_error_mq , M_eim_error_aq, M_eim_error_fq);
        }
    virtual eim_interpolation_error_type eimInterpolationErrorEstimation()
        {
            return boost::make_tuple( M_eim_error_mq , M_eim_error_aq, M_eim_error_fq);
        }
    template<int Row=0,int Col=0>
    sparse_matrix_ptrtype mergeLinearA()
        {
            auto linear_Aqm = M_A.template get<Row,Col>()->linearMqm();
            auto muref = refParameter();
            auto linear_beta = computeBetaLinearDecompositionA( muref );
            int linear_Q = M_A.template get<Row,Col>()->QLinear();

            auto mat = backend()->newMatrix( _trial=functionSpace(),
                                             _test=functionSpace() );
            for ( int q=0; q<linear_Q; q++ )
            {
                int linear_mMax = M_A.template get<Row,Col>()->mMaxLinear(q);
                for( int m=0; m<linear_mMax; m++)
                    mat->addMatrix( linear_beta[q][m], linear_Aqm[q][m] );
            }

            return mat;
        }
    virtual std::vector< std::vector< sparse_matrix_ptrtype >>
    computeLinearDecompositionA()
        {
            std::vector< std::vector< sparse_matrix_ptrtype >> empty_one;
            LOG(INFO) << "Warning : you did not write computeLinearDecompositionA(), the linear part of A is automatically evaluate.\n";
            return empty_one;
        }

    virtual beta_vector_type
    computeBetaLinearDecompositionA ( parameter_type const& mu, double time=0 )
        {
            LOG(INFO) << "Warning : you did not write computeBetaLinearDecompositionA(), the linear coefficients of A are automatically evaluate.\n";
            return M_A.template get()->computeBetaLinear( mu );
        }
    template<int Row=0,int Col=0>
    int QLinearDecompositionA() const
        {
            return M_A.template get<Row,Col>()->QLinear();
        }
    template<int Row=0,int Col=0>
    int mMaxLinearDecompositionA( int q ) const
        {
            return M_A.template get<Row,Col>()->mMaxLinear(q);
        }



    vectorN_type computeStatistics ( Eigen::VectorXd vector , std::string name )
        {
            double min=0,max=0,mean=0,mean1=0,mean2=0,standard_deviation=0,variance=0;
            Eigen::MatrixXf::Index index;
            Eigen::VectorXd square;

            std::ofstream file;
            std::string filename="OnlineStatistics-"+name;

            file.open( filename , std::ios::app );

            if( vector.size() > 0 )
            {
                bool force = boption("eim.use-dimension-max-functions");
                int Neim=0;
                if( force )
                    Neim = ioption("eim.dimension-max");

                int N = vector.size();

                min = vector.minCoeff(&index);
                max = vector.maxCoeff(&index);
                mean = vector.mean();
                mean1 = mean * mean;
                square  = vector.array().pow(2);
                mean2 = square.mean();
                standard_deviation = math::sqrt( mean2 - mean1 );

                if( Environment::worldComm().isMasterRank() )
                {
                    if( force )
                    {
                        std::cout<<"statistics  for "<<name<<" (  was called "<< N << " times with "<<Neim<<" basis functions )"<<std::endl;
                        file <<" statistics  for "<<name<<" (  was called "<< N << " times with "<<Neim<<" basis functions )\n";
                    }
                    else
                    {
                        std::cout<<" statistics  for "<<name<<" (  was called "<< N << " times )"<<std::endl;
                        file <<" statistics  for "<<name<<" (  was called "<< N << " times )\n";
                    }
                    std::cout<<"min : "<<min<<" - max : "<<max<<" mean : "<<mean<<" standard deviation : "<<standard_deviation<<"  (see "<<filename<<")"<<std::endl;
                    file<<"min : "<<min<<" - max : "<<max<<" mean : "<<mean<<" standard deviation : "<<standard_deviation<<"\n";
                }
            }
            vectorN_type result(4);
            result(0)=min;
            result(1)=max;
            result(2)=mean;
            result(3)=standard_deviation;
            return result;
        }
    void readConvergenceDataFromFile( std::vector< vectorN_type > & vector, std::string filename )
        {
            if( Environment::worldComm().isMasterRank() )
            {
                std::vector< std::vector< double > > tmpvector;
                std::ifstream file ( filename );

                std::string str;
                int N,i;
                int Nmax=0;
                double value;
                boost::regex re("[+-]?([[:digit:]]*\\.)?[[:digit:]]+([eE][+-]?[[:digit:]]+)?");
                if( file )
                {
                    //first, determine max elements of the RB
                    //because we could have performed several runs
                    //and the size of the RB can vary between two runs
                    std::string _s;
                    while( std::getline(file, _s) )
                    {
                        std::vector< double > _v;
                        boost::sregex_iterator it(_s.begin(),_s.end(),re); //parse the list of numbers

                        std::string first=it->str();
                        N=std::atoi(first.c_str());
                        it++;
                        if( N > Nmax )
                            Nmax = N;

                        boost::sregex_iterator j;
                        i=0;
                        _v.resize(N);
                        for(; it!=j; ++it)
                        {
                            if( (*it)[0].matched )
                            {
                                std::string s = it->str();
                                auto value = std::atof(s.c_str());
                                if(value != std::numeric_limits<double>::min())
                                {
                                    _v[i]=value;
                                }
                            }
                            else
                                throw std::logic_error( "[ModelCrbBase::fillVectorFromFile] ERROR : Cannot read value" );
                            i++;
                        }
                        //CHECK(i == N);
                        tmpvector.push_back(_v);
                    }
                    vector.resize( Nmax );
                }
                else
                {
                    std::cout<<"The file "<<filename<<" was not found "<<std::endl;
                    throw std::logic_error( "[ModelCrbBase::fillVectorFromFile] ERROR loading the file " );
                }
                file.close();

                //now copy std::vector into eigen vector
                int nbvalues = tmpvector.size();
                for(int n=0; n<Nmax; n++)
                {
                    vector[n].resize(nbvalues);
                    for(int i=0; i<nbvalues; i++)
                    {
                        vector[n](i) = tmpvector[i][n];
                    }
                }

            }//master proc
        }//end of function
    void writeConvergenceStatistics( std::vector< vectorN_type > const& vector,
                                     std::string filename , std::string extra="")
        {
            if( Environment::worldComm().isMasterRank() )
            {
                Eigen::MatrixXf::Index index_max;
                Eigen::MatrixXf::Index index_min;

                double totaltime=0;
                double globaltotaltime=0;
                std::ofstream file;
                file.open( filename );
                if( extra == "totaltime" )
                    file << "NbBasis" << "\t" << "Min" << "\t" << "Max" << "\t" << "Mean" << "\t" << "Variance" << "\t"<< "Total time" << "\n";
                else
                    file << "NbBasis" << "\t" << "Min" << "\t" << "Max" << "\t" << "Mean" << "\t" << "Variance" << "\n";
                int Nmax = vector.size();
                std::vector<double> nbruns( Nmax );
                for(int n=0; n<Nmax; n++)
                {
                    totaltime=0;
                    double variance=0;
                    int sampling_size = vector[n].size();
                    double mean=vector[n].mean();
                    double max = vector[n].maxCoeff(&index_max);
                    double min = vector[n].minCoeff(&index_min);
                    for(int i=0; i<sampling_size; i++)
                    {
                        variance += (vector[n](i) - mean)*(vector[n](i) - mean)/sampling_size;
                        totaltime+= vector[n](i);
                    }
                    globaltotaltime += totaltime;
                    if( extra == "totaltime" )
                        file <<n+1<<"\t"<<min<<"\t"<<max<<"\t"<<mean<<"\t"<<variance<<"\t"<<totaltime<<"\n";
                    else
                        file <<n+1<<"\t"<<min<<"\t"<<max<<"\t"<<mean<<"\t"<<variance<<"\n";
                    nbruns[n]=vector[n].size();
                }

                if( extra == "totaltime" )
                {
                    file << "#global total time : "<<globaltotaltime<<"\n";
                }

                //write information about number of runs
                for(int n=0; n<Nmax; n++)
                {
                    file <<"#N = "<<n<<" -- number of runs : "<<nbruns[n]<<"\n";
                }

            }
        }
    void writeVectorsExtremumsRatio( std::vector< vectorN_type > const& error,
                                     std::vector< vectorN_type > const& estimated,
                                     std::string filename )
        {
            if( Environment::worldComm().isMasterRank() )
            {

                Eigen::MatrixXf::Index index_max;
                Eigen::MatrixXf::Index index_min;

                std::ofstream file;
                file.open( filename );
                file << "NbBasis" << "\t" << "Min" << "\t" << "Max" << "\n";
                int Nmax = error.size();
                int estimatedsize = estimated.size();
                CHECK( Nmax == estimatedsize ) << "vectors have not the same size ! v1.size() : "<<error.size()<<" and v2.size() : "<<estimated.size()<<"\n";
                for(int n=0; n<Nmax; n++)
                {
                    double min_error = error[n].maxCoeff(&index_min);
                    double max_error = error[n].minCoeff(&index_max);
                    double min_estimated = estimated[n](index_min);
                    double max_estimated = estimated[n](index_max);
                    double min_ratio = min_estimated / min_error;
                    double max_ratio = max_estimated / max_error;
                    file << n+1 <<"\t"<< min_ratio <<" \t" << max_ratio<< "\n";
                }
            }
        }


    sparse_matrix_ptrtype newMatrix() const
        {
            return M_backend->newMatrix( _test=functionSpace(), _trial=functionSpace() );
        }
    vector_ptrtype newVector() const
        {
            return M_backend->newVector( functionSpace() );
        }
    void addEnergyMatrix( form2_type const & f )
        {
            M_energy_matrix = f.matrixPtr() ;
        }
    void addEnergyMatrix( sparse_matrix_ptrtype const & matrix )
        {
            M_energy_matrix = matrix ;
        }

    sparse_matrix_ptrtype energyMatrix ()
        {
            double norm = M_inner_product_matrix[0]->l1Norm();
            CHECK( norm > 0 )<<"The energy matrix has not be filled !\n";
            return M_inner_product_matrix[0];
        }
    sparse_matrix_ptrtype energyMatrix () const
        {
            double norm = M_inner_product_matrix[0]->l1Norm();
            CHECK( norm > 0 )<<"The energy matrix has not be filled ! \n";
            return M_inner_product_matrix[0];
        }

    void addMassMatrix( form2_type const & f )
        {
            M_mass_matrix = f.matrixPtr() ;
        }
    void addMassMatrix( sparse_matrix_ptrtype const & matrix )
        {
            M_mass_matrix = matrix ;
        }
    virtual sparse_matrix_ptrtype const& massMatrix () const
        {
            double norm = M_mass_matrix->l1Norm();
            CHECK( norm > 0 )<<"The mass matrix has not be filled !\n";
            return M_mass_matrix;
        }
    virtual sparse_matrix_ptrtype massMatrix ()
        {
            double norm = M_mass_matrix->l1Norm();
            CHECK( norm > 0 )<<"The mass matrix has not be filled !\n";
            return M_mass_matrix;
        }

    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y, int n_space=0 )
        {
            double norm = M_inner_product_matrix[n_space]->l1Norm();
            CHECK( norm > 0 )<<"The energy matrix has not be filled ! (space #)"<<n_space<<std::endl;
            return M_inner_product_matrix[n_space]->energy( X, Y );
        }
    double scalarProduct( vector_type const& X, vector_type const& Y, int n_space=0 )
        {
            double norm = M_inner_product_matrix[n_space]->l1Norm();
            CHECK( norm > 0 )<<"The energy matrix has not be filled ! (space #)"<<n_space<<std::endl;
            return M_inner_product_matrix[n_space]->energy( X, Y );
        }

    /**
     * returns the scalar product used for the mass matrix of the vector x and vector y
     */
    double scalarProductForMassMatrix( vector_type const& X, vector_type const& Y )
        {
            auto M = massMatrix();
            return M->energy( X, Y );
        }

    /**
     * returns the scalar product used for the mass matrix of the vector x and vector y
     */
    double scalarProductForMassMatrix( vector_ptrtype const& X, vector_ptrtype const& Y )
        {
            auto M = massMatrix();
            return M->energy( X, Y );
        }
    /**
     * returns the scalar product used to assemble POD matrix of the vector x and vector y
     */
    double scalarProductForPod( vector_type const& X, vector_type const& Y, int n_space=0 )
        {
            if ( is_time_dependent )
                return M_inner_product_matrix[n_space]->energy( X, Y );
            else
                return 0;
        }
    /**
     * returns the scalar product used to assemble POD matrix of the vector x and vector y
     */
    double scalarProductForPod( vector_ptrtype const& X, vector_ptrtype const& Y )
        {
            return scalarProductForPod( *X, *Y );
        }
    double computeNormL2( element_type u1 , element_type u2 )
        {
            static const bool is_composite = functionspace_type::is_composite;
            return computeNormL2( u1, u2 , mpl::bool_< is_composite >() );
        }
    double computeNormL2( element_type u1 , element_type u2 , mpl::bool_<false>)
        {
            double norm=normL2(_range=elements( u2.mesh() ),_expr=( idv(u1)-idv(u2) ) );
            return norm;
        }
    double computeNormL2( element_type u1 , element_type u2 , mpl::bool_<true>)
        {
            //in that case we accumulate the norm of every components
            ComputeNormL2InCompositeCase compute_normL2_in_composite_case( u1, u2 );
            index_vector_type index_vector;
            fusion::for_each( index_vector, compute_normL2_in_composite_case );
            return compute_normL2_in_composite_case.norm();
        }
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f, int const& n_space=0 )
        {
            M_backend_l2[n_space]->solve( _matrix=M_inner_product_matrix[n_space],  _solution=u, _rhs=f );
        }




    void generateGeoFileForOutputPlot(  vectorN_type outputs, vectorN_type parameter, vectorN_type estimated_error )
        {
            bool use_estimated_error=true;
            if( estimated_error(0) < 0 )
                use_estimated_error=false;

            Eigen::MatrixXf::Index index;
            double min_output = outputs.minCoeff(&index);
            //double min_scale=std::floor(min_output);
            double min_scale=min_output;
            double x=0;
            double output=0;
            double estimated_down=0;
            double estimated_up=0;
            double delta=0;

            std::string plotFile = "GMSH-outputs.geo";
            std::ofstream file_outputs_geo_gmsh ( plotFile, std::ios::out );
            file_outputs_geo_gmsh << "View \"outputs\" {\n";
            for(int i=0; i<outputs.size(); i++)
            {
                if( use_estimated_error )
                {
                    delta=estimated_error(i);
                    x=parameter(i);
                    estimated_down=outputs(i);
                    output=outputs(i)+delta/2.0;
                    estimated_up=estimated_down+delta;
                    file_outputs_geo_gmsh <<std::setprecision(14)<< "SP("<<x<<",0,0){"<<output<<", "<<min_scale<<"};\n";
                    file_outputs_geo_gmsh <<std::setprecision(14)<< "SP("<<x<<",0,0){"<<estimated_down<<", "<<min_scale<<"};\n";
                    file_outputs_geo_gmsh <<std::setprecision(14)<< "SP("<<x<<",0,0){"<<estimated_up<<", "<<min_scale<<"};\n";
                    file_outputs_geo_gmsh <<std::setprecision(14)<< "SP("<<x<<",0,0){"<<output<<", "<<min_scale<<"};\n";
                }
                else
                {
                    x=parameter(i);
                    output=outputs(i);
                    file_outputs_geo_gmsh << "SP("<<x<<",0,0){"<<output<<", "<<min_scale<<"};\n";
                }
            }

            std::string conclude=" }; \n ";
            conclude += "vid = PostProcessing.NbViews-1;\n";
            conclude += "View[vid].Axes = 1;\n";
            conclude += "View[vid].Type = 2;\n\n";
            conclude += "For i In {0:vid-1}\n";
            conclude += "  View[i].Visible=0;\n";
            conclude += "EndFor\n";

            file_outputs_geo_gmsh<<conclude;
            file_outputs_geo_gmsh.close();

            /* Adds the generated file for automatic loading in Gmsh */
            Environment::olLoadInGmsh(plotFile);
        }
    gmsh_ptrtype createStructuredGrid( std::vector<int> components_vary,
                                       std::vector<parameter_type> extremums,
                                       std::vector<int> cuttings,
                                       std::vector<double> time_cuttings, bool time_vary )
        {
            auto min=extremums[0];
            auto max=extremums[1];
            double min0 = min(components_vary[0]);
            double min1 = min(components_vary[1]);
            double max0 = max(components_vary[0]);
            double max1 = max(components_vary[1]);
            double Ti=time_cuttings[0];
            double Tf=time_cuttings[1];
            double dt=time_cuttings[2];
            int nb0=cuttings[0];
            int nb1=cuttings[1];
            if( time_vary )
            {
                nb0=(Tf-Ti)/dt;
                min0=Ti;
                max0=Tf;
            }
            gmsh_ptrtype gmshp( new Gmsh );
            std::ostringstream ostr;

            //we want that each cell created here will be a cell of the mesh
            //so we take a large hsize
            int p=1;
            ostr <<"min0 = "<<min0<<";\n"
                 <<"max0 = "<<max0<<";\n"
                 <<"min1 = "<<min1<<";\n"
                 <<"max1 = "<<max1<<";\n"
                 <<"Point (1) = { min0, min1, 1, hsize };\n"
                 <<"Point (2) = { max0, min1, 1, hsize };\n"
                 <<"Point (3) = { max0, max1, 1, hsize };\n"
                 <<"Point (4) = { min0, max1, 1, hsize };\n"
                 <<"Line (1) = { 1 , 2 };\n"
                 <<"Line (2) = { 2 , 3 };\n"
                 <<"Line (3) = { 3 , 4 };\n"
                 <<"Line (4) = { 4 , 1 };\n"
                 <<"Line Loop (1) = {1,2,3,4};\n"
                 <<"Plane Surface (100) = {1};\n"
                 <<"Transfinite Line{1,-3} = "<<nb0<<";\n"
                 <<"Transfinite Line{2,-4} = "<<nb1<<";\n"
                 <<"Transfinite Surface{100} = {1,2,3,4};\n"
                 <<"Physical Surface (\"Omega\") = {100};\n"
                ;

            std::ostringstream nameStr;
            nameStr.precision( 3 );
            nameStr << "StructuredGrid";
            gmshp->setPrefix( nameStr.str() );
            gmshp->setDescription( ostr.str() );
            return gmshp;
        }





    struct ComputeNormL2InCompositeCase
    {

        ComputeNormL2InCompositeCase( element_type const composite_u1 ,
                                      element_type const composite_u2 )
            :
            M_composite_u1 ( composite_u1 ),
            M_composite_u2 ( composite_u2 )
            {}

        template< typename T >
        void
        operator()( const T& t ) const
            {
                int i = T::value;
                if( i == 0 )
                    M_vec.resize( 1 );
                else
                    M_vec.conservativeResize( i+1 );

                auto u1 = M_composite_u1.template element< T::value >();
                auto u2 = M_composite_u2.template element< T::value >();
                mesh_ptrtype mesh = u1.functionSpace()->mesh();
                double norm  = normL2(_range=elements( mesh ),_expr=( idv(u1)-idv(u2) ) );
                M_vec(i)= norm ;
            }

        double norm()
            {
                return M_vec.sum();
            }

        mutable vectorN_type M_vec;
        element_type M_composite_u1;
        element_type M_composite_u2;
    };

    //for linear steady models, mass matrix does not exist
    //non-linear steady models need mass matrix for the initial guess
    //this class can provide the function operatorCompositeM() to ensure compilation
    //but we need to know if the model can provide an operator composite for mass matrix
    //because non-linear problems have to do these operators
    //because CRBModel is not able to construct such operators in the general case.
    //Elements of function space are members of CRBModel
    //then for composite spaces, we need a view of these elements
    //BUT we can't have a view of an element of a non-composite space
    //this function returns true if the model provides an operator composite for mass matrix
    //false by default
    virtual bool constructOperatorCompositeM() const
        {
            return false;
        }
    virtual operatorcomposite_ptrtype operatorCompositeM() const
        {
            CHECK( false ) << "You chose to use operatorCompositeM but you did not write the function operatorCompositeM().\nYou have to do it or set constructOperatorCompositeM() to false.\n";
                operatorcomposite_ptrtype dummy;
            return dummy;
        }

    virtual beta_vector_type computeBetaInitialGuess( parameter_type const& mu ) const
        {
            beta_vector_type beta;
            beta.resize(1); //q=1
            beta[0].resize(1); //m=1
            beta[0][0]=0;
            return beta;
        }
    element_ptrtype assembleInitialGuess( parameter_type const& mu )
        {
            element_ptrtype initial_guess = Xh->elementPtr();
            initial_guess_type vector_initial_guess;
            beta_vector_type beta;
            vector_initial_guess =computeInitialGuessAffineDecomposition();
            beta = computeBetaInitialGuess( mu );
            int q_max = vector_initial_guess.size();
            for ( size_type q = 0; q < q_max; ++q )
            {
                int m_max=vector_initial_guess[q].size();
                for ( size_type m = 0; m < m_max ; ++m )
                {
                    element_type temp = Xh->element();
                    temp = *vector_initial_guess[q][m];
                    temp.scale( beta[q][m] );
                    *initial_guess += temp;
                }
            }
            return initial_guess;
        }
    virtual funs_type scalarContinuousEim()
        {
            return M_funs;
        }
    virtual funsd_type scalarDiscontinuousEim () const
        {
            return M_funs_d;
        }

    virtual int mMaxInitialGuess( int q ) const
        {
            if( q == 0 )
                return 1;
            else
                throw std::logic_error( "[ModelCrbBase::mMaxInitialGuess(q)] ERROR wrong index q, should be 0" );
        }
    virtual int QInitialGuess() const
        {
            if ( is_linear )
                return 0;
            else
                return 1;
        }
    virtual void initializationField( element_ptrtype& initial_field,parameter_type const& mu )
        {
            if ( is_time_dependent )
                initial_field->zero();
        }
    template< int Row=0, int Col=0>
    value_type linearDecompositionAqm( uint16_type q, uint16_type m, element_type const& xi_i,
                                       element_type const& xi_j, bool transpose = false ) const
        {
            bool stock = boption(_name="crb.stock-matrices");
            if( stock )
            {
                //in this case matrices have already been stocked
                return M_A.template get<Row,Col>()->computeLinear(q,m,transpose)->energy( xi_i, xi_j );
            }
            CHECK( stock )<<" This check is here because nothing is done to deal with linear decomposition of a (model using EIM) when using "
                          <<" operators free and option crb.stock-matrices=false.\nFor now we suppose that crb.stock-matrices=true.\n "
                          <<" So make sure that all developments have been done to deal with crb.stock-matrices=false before delete this CHECK.\n";
            return 0;
        }


    virtual monolithic_type computeMonolithicFormulation( parameter_type const& mu )
        {
            return boost::make_tuple( M_monoM, M_monoA, M_monoF );
        }

    virtual monolithic_type computeMonolithicFormulationU( parameter_type const& mu ,
                                                           element_type const& u )
        {
            return boost::make_tuple( M_monoM, M_monoA, M_monoF );
        }


    virtual element_type solve( parameter_type const& mu )
        {
            element_type solution;// = M_model->functionSpace()->element();
            if( is_linear )
                solution = this->solveFemUsingAffineDecompositionFixedPoint( mu );

            return solution;
        }
    element_type solveFemMonolithicFormulation( parameter_type const& mu )
        {
            bdf_ptrtype mybdf;
            mybdf = bdf( _space=Xh, _name="mybdf" );
            sparse_matrix_ptrtype A;
            sparse_matrix_ptrtype M;
            std::vector<vector_ptrtype> F;
            element_ptrtype InitialGuess = Xh->elementPtr();

            auto u = Xh->element();

            double time_initial;
            double time_step;
            double time_final;

            if ( this->isSteady() )
            {
                time_initial=0;
                time_step = 1e30;
                time_final = 1e30;
                // !!
                //some stuff needs to be done here
                //to deal with non linear problems
                //and have an initial guess
                // !!
            }
            else
            {
                time_initial=this->timeInitial();
                time_step=this->timeStep();
                time_final=this->timeFinal();
                this->initializationField( InitialGuess, mu );
            }

            mybdf->setTimeInitial( time_initial );
            mybdf->setTimeStep( time_step );
            mybdf->setTimeFinal( time_final );

            double bdf_coeff ;
            auto vec_bdf_poly = M_backend->newVector( Xh );

            auto uold = Xh->element();
            u=uold;

            int max_fixedpoint_iterations  = ioption("crb.max-fixedpoint-iterations");
            double increment_fixedpoint_tol  = doption("crb.increment-fixedpoint-tol");
            int iter=0;
            double norm=0;

            for( mybdf->start(*InitialGuess); !mybdf->isFinished(); mybdf->next() )
            {
                bdf_coeff = mybdf->polyDerivCoefficient( 0 );
                auto bdf_poly = mybdf->polyDeriv();
                *vec_bdf_poly = bdf_poly;
                iter=0;
                do {

                    if( is_linear )
                        boost::tie(M, A, F) = this->computeMonolithicFormulation( mu );
                    else
                        boost::tie(M, A, F) = this->computeMonolithicFormulationU( mu , u );

                    if( !isSteady() )
                    {
                        A->addMatrix( bdf_coeff, M );
                        F[0]->addVector( *vec_bdf_poly, *M );
                    }
                    uold=u;
                    M_backend_primal->solve( _matrix=A , _solution=u, _rhs=F[0] );

                    mybdf->shiftRight(u);

                    if( is_linear )
                        norm = 0;
                    else
                        norm = this->computeNormL2( uold , u );
                    iter++;
                } while( norm > increment_fixedpoint_tol && iter<max_fixedpoint_iterations );
            }

            return u;
        }
    element_type solveFemDualMonolithicFormulation( parameter_type const& mu )
        {
            int output_index = ioption(_name="crb.output-index");

            bdf_ptrtype mybdf;
            mybdf = bdf( _space=Xh , _name="mybdf" );
            sparse_matrix_ptrtype A,Adu;
            sparse_matrix_ptrtype M;
            std::vector<vector_ptrtype> F;
            vector_ptrtype Rhs( M_backend->newVector( Xh ) );
            element_ptrtype InitialGuess = Xh->elementPtr();
            auto dual_initial_field = Xh->elementPtr();

            auto udu = Xh->element();

            double time_initial;
            double time_step;
            double time_final;

            if ( this->isSteady() )
            {
                time_initial=0;
                time_step = 1e30;
                time_final = 1e30;
                // !!
                //some stuff needs to be done here
                //to deal with non linear problems
                //and have an initial guess
                // !!
            }
            else
            {
                time_initial=this->timeFinal()+this->timeStep();
                time_step=-this->timeStep();
                time_final=this->timeInitial()+this->timeStep();
            }

            mybdf->setTimeInitial( time_initial );
            mybdf->setTimeStep( time_step );
            mybdf->setTimeFinal( time_final );

            double bdf_coeff ;
            auto vec_bdf_poly = M_backend->newVector( Xh );

            if ( this->isSteady() )
                udu.zero() ;
            else
            {
                boost::tie(M, A, F) = this->computeMonolithicFormulation( mu );
                *Rhs=*F[output_index];
                M_backend_dual->solve( _matrix=M, _solution=dual_initial_field, _rhs=Rhs );
                udu=*dual_initial_field;
            }

            for( mybdf->start(*InitialGuess); !mybdf->isFinished(); mybdf->next() )
            {
                bdf_coeff = mybdf->polyDerivCoefficient( 0 );
                auto bdf_poly = mybdf->polyDeriv();
                *vec_bdf_poly = bdf_poly;

                boost::tie(M, A, F) = this->computeMonolithicFormulation( mu );

                if( ! isSteady() )
                {
                    A->addMatrix( bdf_coeff, M );
                    Rhs->zero();
                    *vec_bdf_poly = bdf_poly;
                    Rhs->addVector( *vec_bdf_poly, *M );
                }
                else
                {
                    *Rhs = *F[output_index];
                    Rhs->close();
                    Rhs->scale( -1 );
                }

                if( boption("crb.use-symmetric-matrix") )
                    Adu = A;
                else
                    A->transpose( Adu );

                M_backend_dual->solve( _matrix=Adu , _solution=udu, _rhs=Rhs );

                mybdf->shiftRight(udu);
            }

            return udu;

        }
    element_type solveFemDualUsingAffineDecompositionFixedPoint( parameter_type const& mu )
        {
            int output_index = ioption(_name="crb.output-index");

            bdf_ptrtype mybdf;
            mybdf = bdf( _space=Xh , _name="mybdf" );
            sparse_matrix_ptrtype A,Adu;
            sparse_matrix_ptrtype M;
            std::vector<vector_ptrtype> F;
            auto udu = Xh->element();
            auto uold = Xh->element();
            vector_ptrtype Rhs( M_backend->newVector( Xh ) );
            auto dual_initial_field = Xh->elementPtr();

            double time_initial;
            double time_step;
            double time_final;

            if ( this->isSteady() )
            {
                time_initial=0;
                time_step = 1e30;
                time_final = 1e30;
                //InitialGuess = this->assembleInitialGuess( mu ) ;
            }
            else
            {
                time_initial=this->timeFinal()+this->timeStep();
                time_step=-this->timeStep();
                time_final=this->timeInitial()+this->timeStep();
            }

            mybdf->setTimeInitial( time_initial );
            mybdf->setTimeStep( time_step );
            mybdf->setTimeFinal( time_final );

            double norm=0;
            int iter=0;

            double bdf_coeff ;
            auto vec_bdf_poly = M_backend->newVector( Xh );

            if ( this->isSteady() )
                udu.zero() ;
            else
            {
                boost::tie( M, A, F) = this->update( mu , mybdf->timeInitial() );
                *Rhs=*F[output_index];
                M_backend_dual->solve( _matrix=M, _solution=dual_initial_field, _rhs=Rhs);
                udu=*dual_initial_field;
            }


            int max_fixedpoint_iterations  = ioption(_name="crb.max-fixedpoint-iterations");
            double increment_fixedpoint_tol  = doption(_name="crb.increment-fixedpoint-tol");
            for( mybdf->start(udu); !mybdf->isFinished(); mybdf->next() )
            {
                iter=0;
                bdf_coeff = mybdf->polyDerivCoefficient( 0 );
                auto bdf_poly = mybdf->polyDeriv();
                *vec_bdf_poly = bdf_poly;
                do {
                    if( is_linear )
                        boost::tie(M, A, F) = this->update( mu , mybdf->time() );
                    else
                        boost::tie(M, A, F) = this->update( mu , udu , mybdf->time() );

                    if( ! isSteady() )
                    {
                        A->addMatrix( bdf_coeff, M );
                        Rhs->zero();
                        *vec_bdf_poly = bdf_poly;
                        Rhs->addVector( *vec_bdf_poly, *M );
                    }
                    else
                    {
                        *Rhs = *F[output_index];
                        Rhs->close();
                        Rhs->scale( -1 );
                    }

                    if( boption("crb.use-symmetric-matrix") )
                        Adu = A;
                    else
                        A->transpose( Adu );

                    uold = udu;
                    M_backend_dual->solve( _matrix=Adu , _solution=udu, _rhs=Rhs);

                    if( boption(_name="crb.use-linear-model") )
                        norm = 0;
                    else
                        norm = this->computeNormL2( uold , udu );
                    iter++;
                } while( norm > increment_fixedpoint_tol && iter<max_fixedpoint_iterations );
                mybdf->shiftRight(udu);
            }
            return udu;
        }
    element_type solveFemUsingOfflineEim( parameter_type const& mu )
        {
            bdf_ptrtype mybdf;
            mybdf = bdf( _space=Xh, _name="mybdf" );
            sparse_matrix_ptrtype A;
            sparse_matrix_ptrtype M;
            std::vector<vector_ptrtype> F;
            element_ptrtype InitialGuess = Xh->elementPtr();
            vector_ptrtype Rhs( M_backend->newVector( Xh ) );
            auto u = Xh->element();

            double time_initial;
            double time_step;
            double time_final;

            if ( this->isSteady() )
            {
                time_initial=0;
                time_step = 1e30;
                time_final = 1e30;
            }
            else
            {
                time_initial=this->timeInitial();
                time_step=this->timeStep();
                time_final=this->timeFinal();
                this->initializationField( InitialGuess, mu );
            }

            mybdf->setTimeInitial( time_initial );
            mybdf->setTimeStep( time_step );
            mybdf->setTimeFinal( time_final );

            double bdf_coeff ;
            auto vec_bdf_poly = M_backend->newVector( Xh );

            for( mybdf->start(*InitialGuess); !mybdf->isFinished(); mybdf->next() )
            {
                bdf_coeff = mybdf->polyDerivCoefficient( 0 );
                auto bdf_poly = mybdf->polyDeriv();
                *vec_bdf_poly = bdf_poly;
                boost::tie(M, A, F) = this->update( mu , mybdf->time() );
                *Rhs = *F[0];
                if( !isSteady() )
                {
                    A->addMatrix( bdf_coeff, M );
                    Rhs->addVector( *vec_bdf_poly, *M );
                }
                M_backend_primal->solve( _matrix=A , _solution=u, _rhs=Rhs);
                mybdf->shiftRight(u);
            }

            return u;

        }
    element_type solveFemUsingAffineDecompositionFixedPoint( parameter_type const& mu )
        {
            bdf_ptrtype mybdf;
            mybdf = bdf( _space=Xh, _name="mybdf" );
            sparse_matrix_ptrtype A;
            sparse_matrix_ptrtype M;
            std::vector<vector_ptrtype> F;
            element_ptrtype InitialGuess = Xh->elementPtr();
            auto u = Xh->element("u");
            auto uold = Xh->element("u_old");
            vector_ptrtype Rhs( M_backend->newVector( Xh ) );

            double time_initial;
            double time_step;
            double time_final;

            if ( this->isSteady() )
            {
                time_initial=0;
                time_step = 1e30;
                time_final = 1e30;
                //we want to have the initial guess given by function update
                InitialGuess = this->assembleInitialGuess( mu ) ;
            }
            else
            {
                time_initial=this->timeInitial();
                time_step=this->timeStep();
                time_final=this->timeFinal();
                this->initializationField( InitialGuess, mu );
            }

            mybdf->setTimeInitial( time_initial );
            mybdf->setTimeStep( time_step );
            mybdf->setTimeFinal( time_final );

            u=*InitialGuess;
            double norm=0;
            int iter=0;

            double bdf_coeff ;
            auto vec_bdf_poly = M_backend->newVector( Xh );

            int max_fixedpoint_iterations  = ioption("crb.max-fixedpoint-iterations");
            double increment_fixedpoint_tol  = doption("crb.increment-fixedpoint-tol");

            for( mybdf->start(*InitialGuess); !mybdf->isFinished(); mybdf->next() )
            {
                iter=0;
                bdf_coeff = mybdf->polyDerivCoefficient( 0 );
                auto bdf_poly = mybdf->polyDeriv();
                *vec_bdf_poly = bdf_poly;
                do {
                    if( is_linear )
                        boost::tie(M, A, F) = update( mu , mybdf->time() );
                    else
                        boost::tie(M, A, F) = update( mu , u , mybdf->time() );
                    *Rhs = *F[0];

                    if( !isSteady() )
                    {
                        A->addMatrix( bdf_coeff, M );
                        Rhs->addVector( *vec_bdf_poly, *M );
                    }
                    uold = u;
                    M_backend_primal->solve( _matrix=A , _solution=u, _rhs=Rhs);

                    if( is_linear )
                        norm = 0;
                    else
                        norm = this->computeNormL2( uold , u );
                    iter++;
                } while( norm > increment_fixedpoint_tol && iter<max_fixedpoint_iterations );
                mybdf->shiftRight(u);
            }
            return u;
        }
    offline_merge_type update( parameter_type const& mu,  double time=0,
                               bool only_time_dependent_terms=false )
        {
            auto all_beta = this->computeBetaQm( mu , time , only_time_dependent_terms );
            return offlineMerge( all_beta, only_time_dependent_terms);
        }
    offline_merge_type update( parameter_type const& mu, element_type const& T, double time=0,
                               bool only_time_dependent_terms=false )
        {
#if !defined(NDEBUG)
            mu.check();
#endif
            auto all_beta = this->computeBetaQm(  T , mu , time , only_time_dependent_terms );
            return offlineMerge( all_beta, only_time_dependent_terms);
        }

    template<int Row=0,int Col=0>
    offline_merge_type offlineMerge( betaqm_type const& all_beta, bool only_time_dependent_terms )
        {
            sparse_matrix_ptrtype A;
            sparse_matrix_ptrtype M;
            std::vector<vector_ptrtype> F( M_Nl );

            if ( !only_time_dependent_terms )
            {
                M = M_M.template get<Row,Col>()->merge( all_beta.template get<0>() );
                A = M_A.template get<Row,Col>()->merge( all_beta.template get<1>() );
            }

            for ( int output =0; output<M_Nl; output++ )
                F[output] = M_F[output].template get<Row>()->merge( all_beta.template get<2>()[output] );

            return boost::make_tuple( M, A, F );
        }


    std::vector< std::vector<element_ptrtype> > computeInitialGuessVAffineDecomposition( )
        {
            initial_guess_type initial_guess_v;
            initial_guess_v = computeInitialGuessAffineDecomposition();
            assembleInitialGuessV( initial_guess_v);
            for(int q=0; q<M_InitialGuessVector.size(); q++)
            {
                for(int m=0; m<M_InitialGuessVector[q].size(); m++)
                    *M_InitialGuessV[q][m] = *M_InitialGuessVector[q][m];
            }
            return M_InitialGuessV;
        }
    virtual std::vector< std::vector< element_ptrtype > >
    computeInitialGuessAffineDecomposition()
        {
            std::vector< std::vector<element_ptrtype> > q;
            q.resize(1);
            q[0].resize(1);
            auto mesh = mesh_type::New();
            auto Xh=space_type::New( mesh );
            element_ptrtype elt ( new element_type ( Xh ) );
            q[0][0] = elt;
            return q;
        }
    void assembleInitialGuessV( initial_guess_type & initial_guess )
        {
            static const bool is_composite = functionspace_type::is_composite;
            return assembleInitialGuessV( initial_guess, mpl::bool_< is_composite >() );
        }
    void assembleInitialGuessV( initial_guess_type & initial_guess, mpl::bool_<true> )
        {
            auto mesh = Xh->mesh();

            int q_max= this->QInitialGuess();
            M_InitialGuessV.resize( q_max );
            M_InitialGuessVector.resize( q_max );
            for(int q = 0; q < q_max; q++ )
            {
                int m_max= this->mMaxInitialGuess(q);
                M_InitialGuessV[q].resize( m_max );
                M_InitialGuessVector[q].resize( m_max );
                for(int m = 0; m < m_max; m++ )
                {
                    M_InitialGuessV[q][m] = Xh->elementPtr();
                    M_InitialGuessVector[q][m] = this->newVector();
                }
            }

            index_vector_type index_vector;
            AssembleInitialGuessVInCompositeCase assemble_initial_guess_v_in_composite_case ( v , initial_guess , this->shared_from_this());
            fusion::for_each( index_vector, assemble_initial_guess_v_in_composite_case );

            for(int q = 0; q < q_max; q++ )
            {
                int m_max = this->mMaxInitialGuess(q) ;
                for(int m = 0; m < m_max; m++ )
                    M_InitialGuessVector[q][m]->close();
            }
        }
    void assembleInitialGuessV( initial_guess_type & initial_guess, mpl::bool_<false> )
        {
            auto mesh = Xh->mesh();

            int q_max= this->QInitialGuess();
            M_InitialGuessV.resize( q_max );
            M_InitialGuessVector.resize( q_max );
            for(int q = 0; q < q_max; q++ )
            {
                int m_max= this->mMaxInitialGuess(q);
                M_InitialGuessV[q].resize( m_max );
                M_InitialGuessVector[q].resize( m_max );
                for(int m = 0; m < m_max; m++ )
                {
                    M_InitialGuessV[q][m] = Xh->elementPtr();
                    M_InitialGuessVector[q][m] = this->newVector();
                    form1( _test=Xh, _vector=M_InitialGuessVector[q][m]) =
                        integrate( _range=elements( mesh ), _expr=idv( initial_guess[q][m] )*id( v )  );
                    M_InitialGuessVector[q][m]->close();
                }
            }
        }

    struct AssembleInitialGuessVInCompositeCase
    {
        AssembleInitialGuessVInCompositeCase( element_type  const v ,
                                              initial_guess_type  const initial_guess ,
                                              boost::shared_ptr<self_type > crb_model)
            :
            M_composite_v ( v ),
            M_composite_initial_guess ( initial_guess ),
            M_crb_model ( crb_model )
            {}

        template< typename T >
        void
        operator()( const T& t ) const
            {
                auto v = M_composite_v.template element< T::value >();
                auto Xh = M_composite_v.functionSpace();
                mesh_ptrtype mesh = Xh->mesh();
                int q_max = M_crb_model->QInitialGuess();
                for(int q = 0; q < q_max; q++)
                {
                    int m_max = M_crb_model->mMaxInitialGuess(q);
                    for( int m = 0; m < m_max ; m++)
                    {
                        auto initial_guess_qm = M_crb_model->InitialGuessVector(q,m);
                        auto view = M_composite_initial_guess[q][m]->template element< T::value >();
                        form1( _test=Xh, _vector=initial_guess_qm ) +=
                            integrate ( _range=elements( mesh ), _expr=trans( Feel::vf::idv( view ) )*Feel::vf::id( v ) );
                    }
                }

            }

        element_type  M_composite_v;
        initial_guess_type  M_composite_initial_guess;
        mutable boost::shared_ptr<self_type> M_crb_model;
    };

protected:
    bool M_is_initialized;
    CRBModelMode M_mode;
    backend_ptrtype M_backend;
    backend_ptrtype M_backend_primal;
    backend_ptrtype M_backend_dual;
    std::vector< backend_ptrtype > M_backend_l2;

    parameterspace_ptrtype Dmu;
    functionspace_ptrtype Xh;
    rbfunctionspace_ptrtype XN;

    std::vector< sparse_matrix_ptrtype > M_inner_product_matrix;
    sparse_matrix_ptrtype M_energy_matrix;;
    sparse_matrix_ptrtype M_mass_matrix;

    bdf_ptrtype M_bdf;

    funs_type M_funs;
    funsd_type M_funs_d;

    int M_QLinearDecompositionA;

    sparse_matrix_ptrtype M_monoA;
    sparse_matrix_ptrtype M_monoM;
    std::vector<vector_ptrtype> M_monoF;

    int M_Nl;
    matrix_ad_type M_M;
    matrix_ad_type M_A;
    std::vector< vector_ad_type > M_F;

    mutable std::vector< std::vector<element_ptrtype> > M_InitialGuessV;
    mutable std::vector< std::vector<vector_ptrtype> > M_InitialGuessVector;

    bool M_has_eim;
    std::vector<int> M_mMaxLinearDecompositionA;
    std::vector< std::vector<sparse_matrix_ptrtype> > M_linearAqm;
    std::map<int,double> M_eim_error_mq;
    std::map<int,double> M_eim_error_aq;
    std::vector< std::map<int,double> > M_eim_error_fq;

    element_type u;
    element_type v;

}; //class CRBModel

} // namespace Feel

#endif //__CRBMODELBASE_H
