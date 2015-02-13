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
    InfSup = 0x4
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
    static const int Options = _Options;

    typedef CRBModelBase self_type;
    typedef double value_type;

    // Functions Space
    typedef typename mpl::if_<is_shared_ptr<FunctionSpaceDefinition>,
                              mpl::identity<typename FunctionSpaceDefinition::element_type>,
                              mpl::identity<FunctionSpaceDefinition>>::type::type::space_type space_type;
    typedef space_type functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef functionspace_ptrtype space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    // RB Function Space
    typedef ReducedBasisSpace<self_type> rbfunctionspace_type;
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

    typedef FsFunctionalLinear< space_type > functional_type;
    typedef boost::shared_ptr<functional_type> functional_ptrtype;

    typedef OperatorLinear< space_type , space_type > operator_type;
    typedef boost::shared_ptr<operator_type> operator_ptrtype;

    typedef OperatorLinearComposite< space_type , space_type > operatorcomposite_type;
    typedef boost::shared_ptr<operatorcomposite_type> operatorcomposite_ptrtype;
    typedef FsFunctionalLinearComposite< space_type > functionalcomposite_type;
    typedef boost::shared_ptr<functionalcomposite_type> functionalcomposite_ptrtype;

    typedef vf::detail::BilinearForm<functionspace_type, functionspace_type,VectorUblas<value_type>> form2_type;
    typedef vf::detail::LinearForm<functionspace_type,vector_type,vector_type> form1_type;


    typedef typename std::vector< std::vector < element_ptrtype > > initial_guess_type;

    typedef std::vector< std::vector< value_type > > beta_vector_type;
    typedef AffineDecomposition<EimDefinition,self_type> affinedecomposition_type;
    typedef boost::shared_ptr<affinedecomposition_type> affinedecomposition_ptrtype;

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
    typedef typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<2> > , fusion::vector< mpl::int_<0>, mpl::int_<1> >  ,
                               typename mpl::if_ < boost::is_same< mpl::int_<nb_spaces> , mpl::int_<3> > , fusion::vector < mpl::int_<0> , mpl::int_<1> , mpl::int_<2> > ,
                                                   typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<4> >, fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3> >,
                                                                      fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3>, mpl::int_<4> >
                                                                      >::type >::type >::type index_vector_type;



    CRBModelBase() :
        M_is_initialized( false ),
        M_backend( backend() ),
        M_backend_primal( backend( _name="backend-primal") ),
        M_backend_dual( backend( _name="backend-dual") ),
        M_backend_l2( backend( _name="backend-l2" ) ),
        Dmu( new parameterspace_type ),
        M_has_eim( false )
        {}

    virtual void initModel()=0;

    void init( CRBModelMode mode=CRBModelMode::PFEM )
        {
            if ( M_is_initialized )
                return;


            M_mode = mode;
            M_AD = affinedecomposition_ptrtype(
                new affinedecomposition_type( this->shared_from_this() ) );

            LOG(INFO)<< "Model Initialization";
            initModel();
            //M_AD->check();
            checkBdf();

            u = Xh->element();
            v = Xh->element();
            M_is_initialized = true;
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

    void setFunctionSpaces( functionspace_ptrtype Vh, int row=1 )
        {
            Xh = Vh;
            XN = rbfunctionspace_type::New( _model=this->shared_from_this() );
            if (Environment::isMasterRank() )
                std::cout << "Number of dof : " << Xh->nDof() << std::endl
                          << "Number of local dof : " << Xh->nLocalDof() << std::endl;
        }

    template<typename ExprType>
    void addMass( ExprType const& expr, std::string const& symbol,
                  int const row=1, int const col=1 )
        {
            auto ope = opLinearFree( _domainSpace=functionSpace(col-1),
                                     _imageSpace=functionSpace(row-1),
                                     _expr=expr );
            M_AD->addMass( ope, symbol, row, col) ;
        }
    void addMass( sparse_matrix_ptrtype const& mat, std::string const& symbol,
                  int row=1, int col=1 )
        {
            M_AD->addMass( mat, symbol, row, col );
        }
    void addMass( form2_type const& form2, std::string const& symbol,
                  int row=1, int col=1 )
        {
            M_AD->addMass( form2.matrixPtr(), symbol,row, col );
        }

    template<typename ExprType>
    void addLhs( ExprType const& expr, std::string const & symbol,
                 int const row=1, int const col=1 )
        {
            auto ope = opLinearFree( _domainSpace=functionSpace(col-1),
                                     _imageSpace=functionSpace(row-1),
                                     _expr=expr );
            M_AD->addLhs( ope, symbol, row, col) ;
        }
    void addLhs( sparse_matrix_ptrtype const& mat, std::string const& symbol,
                 int row=1, int col=1 )
        {
            M_AD->addLhs( mat, symbol, row, col );
        }
    void addLhs( form2_type const& form2, std::string const& symbol,
                 int row=1, int col=1 )
        {
            M_AD->addLhs( form2.matrixPtr(), symbol, row, col );
        }

    template<typename ExprType>
    void addRhs( ExprType const& expr, std::string const& symbol,
                 int const row=1 )
        {
            auto fun = functionalLinearFree( _space=functionSpace(row-1),
                                             _expr=expr );
            M_AD->addOutput( fun, symbol, 0 , row) ;
        }
    void addRhs( vector_ptrtype const& vec, std::string const& symbol,
                 int row=1 )
        {
            M_AD->addOutput( vec, symbol, 0, row );
        }

    void addRhs( form1_type const& form1, std::string const& symbol, int row=1  )
        {
            M_AD->addOutput( form1.vectorPtr(), symbol, 0, row );
        }

    template<typename ExprType>
    void addOutput( ExprType const& expr, std::string const& symbol, int output=1, int const row=1 )
        {
            auto fun = functionalLinearFree( _space=functionSpace(row-1),
                                             _expr=expr );
            M_AD->addOutput( fun, symbol, output, row) ;
        }
    void addOutput( vector_ptrtype const& vec, std::string const& symbol, int output=1, int row=1 )
        {
            M_AD->addOutput( vec, symbol, output, row );
        }
    void addOutput( form1_type const& form1, std::string const& symbol, int output=1, int row=1 )
        {
            M_AD->addOutput( form1.vectorPtr(), symbol, output, row );
        }



    parameterspace_ptrtype parameterSpace()
        {
            return Dmu;
        }

    functionspace_ptrtype functionSpace( int num=1 ) const
        {
            return Xh;
        }
    rbfunctionspace_ptrtype rBFunctionSpace( int num=1 )
        {
            return XN;
        }

    affinedecomposition_ptrtype AD()
        {
            return M_AD;
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

            return boost::make_tuple( M_AD->betaAqm(mu), M_AD->betaFqm(mu) );
        }

    virtual betaqm_type computeBetaQm( parameter_type const& mu, double time,
                               bool only_time_dependent_terms=false )
        {
            betaqm_type beta_coeff;
            steady_betaqm_type steady_beta_coeff;
            beta_vector_type betaM;

            if ( !is_time_dependent )
            {
                steady_beta_coeff = computeBetaQm( mu );
                betaM.resize(1);
                betaM[0].resize(1);
                betaM[0][0]=1;
                beta_coeff = boost::make_tuple( betaM, steady_beta_coeff.template get<0>(),
                                                steady_beta_coeff.template get<1>() );
            }
            else
                beta_coeff = boost::make_tuple( M_AD->betaMqm( mu, time ),
                                                M_AD->betaAqm( mu, time ),
                                                M_AD->betaFqm( mu, time ) );
            return beta_coeff;
        }

    virtual steady_betaqm_type BetaQm( element_type const& T, parameter_type const& mu )
        {
            steady_betaqm_type dummy;
            throw std::logic_error( "Since youre model is not linear (neither time dependent) you have to implement the function BetaQm( u, mu ) ");
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
                                       double time=0 , bool only_time_dependent_terms=false )
        {
        betaqm_type beta_coeff;
        steady_betaqm_type steady_beta_coeff;
        beta_vector_type betaM;

        if ( !is_time_dependent )
        {
            steady_beta_coeff = BetaQm( T, mu );
            betaM.resize(1);
            betaM[0].resize(1);
            betaM[0][0]=1;
            beta_coeff = boost::make_tuple( betaM, steady_beta_coeff.template get<0>(),
                                            steady_beta_coeff.template get<1>() );
        }
        else
            throw std::logic_error( "Since youre model is not linear (and time dependent), you have to implement the function computeBetaQm( u, mu, time ) ");

        return beta_coeff;
    }

    virtual affine_decomposition_type computeAffineDecomposition()
        {
            return M_AD->compute();
        }
    void countAffineDecompositionTerms()
        {
            M_AD->count();
        }
    virtual size_type Qm() const
        {
            if (is_time_dependent)
                return M_AD->Qm();
            else
            {
                if ( constructOperatorCompositeM() )
                    return functionspace_type::nSpaces;
                else
                    return 0;
            }
        }
    virtual size_type Ql( int l ) const
        {
            return M_AD->Ql(l);
        }
    virtual size_type Qa() const
        {
            return M_AD->Qa();
        }
    size_type Nl() const
        {
            return M_AD->Nl();
        }


    int mMaxF( int output, int q )
        {
            return M_AD->mMaxF( output, q );
        }
    int mMaxM( int q )
        {
            if ( is_time_dependent )
                return M_AD->mMaxM(q);
            else
                return 1;
        }
    int mMaxA( int q )
        {
            return M_AD->mMaxA(q);
        }

    sparse_matrix_ptrtype Mqm( uint16_type q, uint16_type m, bool transpose = false,
                               int row=1, int col=1) const
        {
            return M_AD->M( q, m, transpose, row, col );
        }
    value_type Mqm( uint16_type q, uint16_type m,
                    element_type const& xi_i, element_type const& xi_j, bool transpose = false,
                    int row=1, int col=1) const
        {
            return M_AD->M( q, m, false, row, col)->energy( xi_j, xi_i, transpose );
        }
    sparse_matrix_ptrtype Aqm( uint16_type q, uint16_type m, bool transpose = false,
                                     int row=1, int col=1) const
        {
            return M_AD->A( q, m, transpose, row, col );
        }

    value_type Aqm( uint16_type q, uint16_type m,
                    element_type const& xi_i, element_type const& xi_j, bool transpose = false,
                    int row=1, int col=1 ) const
        {
            return M_AD->A( q, m, false, row, col )->energy( xi_j, xi_i, transpose);
        }

    vector_ptrtype Fqm( uint16_type l, uint16_type q, int m=0,
                        int row=1 ) const
        {
            return M_AD->F( l, q, m, row );
        }
    value_type Fqm( uint16_type l, uint16_type q,  uint16_type m, element_type const& xi,
                    int row=1 )
        {
            return inner_product( *(M_AD->F( l, q, m, row )), xi );
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


    sparse_matrix_ptrtype newMatrix( int row=1, int col=1 ) const
        {
            return M_backend->newMatrix( _test=functionSpace(row-1), _trial=functionSpace(col-1) );
        }
    vector_ptrtype newVector( int row=1 ) const
        {
            return M_backend->newVector( functionSpace(row-1) );
        }
    void addEnergyMatrix( form2_type const & f )
        {
            M_inner_product_matrix = f.matrixPtr() ;
        }
    void addEnergyMatrix( sparse_matrix_ptrtype const & matrix )
        {
            M_inner_product_matrix = matrix ;
        }

    virtual sparse_matrix_ptrtype energyMatrix ()
        {
            double norm = M_inner_product_matrix->l1Norm();
            CHECK( norm > 0 )<<"The energy matrix has not be filled !\n";
            return M_inner_product_matrix;
        }
    virtual sparse_matrix_ptrtype energyMatrix () const
        {
            double norm = M_inner_product_matrix->l1Norm();
            CHECK( norm > 0 )<<"The energy matrix has not be filled ! \n";
            return M_inner_product_matrix;
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
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y )
        {
            return M_inner_product_matrix->energy( X, Y );
        }
    double scalarProduct( vector_type const& X, vector_type const& Y )
        {
            return M_inner_product_matrix->energy( X, Y );
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
    double scalarProductForPod( vector_type const& X, vector_type const& Y )
        {
            if ( is_time_dependent )
                return M_inner_product_matrix->energy( X, Y );
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
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f )
        {
            M_backend_l2->solve( _matrix=M_inner_product_matrix,  _solution=u, _rhs=f );
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
    virtual operatorcomposite_ptrtype operatorCompositeA() const
        {
            return M_compositeA;
        }
    virtual operatorcomposite_ptrtype operatorCompositeM()
        {
            return M_compositeM;
        }
    virtual std::vector< functionalcomposite_ptrtype > functionalCompositeF() const
        {
            return M_compositeF;
        }

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
    virtual beta_vector_type computeBetaLinearDecompositionA ( parameter_type const& mu ,  double time=0 )
        {
            auto tuple = computeBetaQm( mu, time );
            return tuple.template get<1>();
        }
    int QLinearDecompositionA() const
        {
            return M_QLinearDecompositionA;
        }
    int mMaxLinearDecompositionA( int q ) const
        {
            return M_mMaxLinearDecompositionA[q];
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
    value_type linearDecompositionAqm( uint16_type q, uint16_type m, element_type const& xi_i,
                                       element_type const& xi_j, bool transpose = false ) const
        {
            bool stock = boption(_name="crb.stock-matrices");
            if( stock )
            {
                //in this case matrices have already been stocked
                return M_linearAqm[q][m]->energy( xi_j, xi_i, transpose );
            }
            CHECK( stock )<<" This check is here because nothing is done to deal with linear decomposition of a (model using EIM) when using "
                          <<" operators free and option crb.stock-matrices=false.\nFor now we suppose that crb.stock-matrices=true.\n "
                          <<" So make sure that all developments have been done to deal with crb.stock-matrices=false before delete this CHECK.\n";
            return 0;
        }
    virtual operatorcomposite_ptrtype operatorCompositeLightA() const
        {
            return M_compositeA;
        }
    virtual operatorcomposite_ptrtype operatorCompositeLightM() const
        {
            if ( is_time_dependent )
                return M_compositeM;
            else
            {    if ( constructOperatorCompositeM() )
                    return M_compositeM;
                else
                    return preAssembleMassMatrix( true );
            }
        }
    operatorcomposite_ptrtype preAssembleMassMatrix( bool light_version ) const
        {
            static const bool is_composite = functionspace_type::is_composite;
            return preAssembleMassMatrix( mpl::bool_< is_composite >() , light_version );
        }
    operatorcomposite_ptrtype preAssembleMassMatrix( mpl::bool_<false> , bool light_version ) const
        {
            auto mesh = Xh->mesh();

            auto expr=integrate( _range=elements( mesh ) , _expr=idt( u )*id( v ) );
            auto op_mass = opLinearComposite( _domainSpace=Xh , _imageSpace=Xh  );
            auto opfree = opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr );
            opfree->setName("mass operator (automatically created)");
            //in this case, the affine decompositon has only one element
            if( light_version )
            {
                //note that only time independent problems need to
                //inform if we use light version or not of the affine decomposition
                //i.e. if we use EIM expansion or not
                //because transient models provide already a decomposition of the mass matrix
                //so it is not automatically created
                op_mass->addElement( 0 , opfree );
            }
            else
            {
                op_mass->addElement( boost::make_tuple(0,0) , opfree );
            }
            return op_mass;
        }
    operatorcomposite_ptrtype preAssembleMassMatrix( mpl::bool_<true> , bool light_version ) const
        {
            index_vector_type index_vector;
            PreAssembleMassMatrixInCompositeCase<self_type> preassemble_mass_matrix_in_composite_case ( u , v );
            fusion::for_each( index_vector, preassemble_mass_matrix_in_composite_case );

            auto op_mass = preassemble_mass_matrix_in_composite_case.opmass();
            return op_mass;
        }

    virtual std::vector< functionalcomposite_ptrtype > functionalCompositeLightF() const
        {
            return M_compositeF;
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
            offline_merge_type offline_merge;
            if( boption(_name="crb.stock-matrices") )
                offline_merge = offlineMerge( all_beta , only_time_dependent_terms );
            else
                offline_merge = offlineMergeOnFly( all_beta, only_time_dependent_terms );

            return offline_merge;
        }
    offline_merge_type update( parameter_type const& mu, element_type const& T, double time=0,
                               bool only_time_dependent_terms=false )
        {
#if !defined(NDEBUG)
            mu.check();
#endif
            auto all_beta = this->computeBetaQm(  T , mu , time , only_time_dependent_terms );

            offline_merge_type offline_merge;

            if( boption(_name="crb.stock-matrices") )
                offline_merge = offlineMerge( all_beta , only_time_dependent_terms );
            else
                offline_merge = offlineMergeOnFly( all_beta, only_time_dependent_terms );

            return offline_merge;

        }


    offline_merge_type offlineMergeOnFly( betaqm_type const& all_beta,
                                          bool only_time_dependent_terms )
        {

            //recovery from the model of free operators Aqm in A = \sum_q \sum_m \beta_qm( u ; mu ) A_qm( u, v )
            //idem with other operator / functional

            operatorcomposite_ptrtype compositeM;
            operatorcomposite_ptrtype compositeA;
            std::vector< functionalcomposite_ptrtype > vector_compositeF;

            if( ! only_time_dependent_terms )
            {

                auto compositeAlight=operatorCompositeLightA();
                if( compositeAlight )
                {
                    compositeA=compositeAlight;
                    compositeM = this->operatorCompositeLightM();
                    vector_compositeF = this->functionalCompositeLightF();
                }
                else
                {
                    compositeA = this->operatorCompositeA();
                    compositeM = this->operatorCompositeM();
                    vector_compositeF = this->functionalCompositeF();
                }
                //acces to beta coefficients
                auto beta_M = all_beta.template get<0>();
                auto beta_A = all_beta.template get<1>();

                //associate beta coefficients to operators
                compositeA->setScalars( beta_A );
                compositeM->setScalars( beta_M );
            }

            //merge
            auto A = this->newMatrix();
            auto M = this->newMatrix();
            if( ! only_time_dependent_terms )
            {
                compositeA->sumAllMatrices( A );
                compositeM->sumAllMatrices( M );
            }

            //warning : beta_F is a vector of beta_coefficients
            auto beta_F = all_beta.template get<2>();
            std::vector<vector_ptrtype> F( Nl() );
            for(int output=0; output<Nl(); output++)
            {
                auto compositeF = vector_compositeF[output];
                compositeF->setScalars( beta_F[output] );
                F[output] = this->newVector();
                compositeF->sumAllVectors( F[output] );
            }

            return boost::make_tuple( M, A, F );

        }
    offline_merge_type offlineMerge( betaqm_type const& all_beta , bool only_time_dependent_terms )
        {
            auto A = this->newMatrix();
            auto M = this->newMatrix();

            if( ! only_time_dependent_terms )
            {
                //acces to beta coefficients
                auto beta_M = all_beta.template get<0>();
                auto beta_A = all_beta.template get<1>();

                A->zero();
                for ( size_type q = 0; q < Qa(); ++q )
                {
                    for(size_type m = 0; m < mMaxA(q); ++m )
                    {
                        A->addMatrix( beta_A[q][m], M_AD->Aqm()[q][m] );
                    }
                }

                if( Qm() > 0 )
                {
                    for ( size_type q = 0; q < Qm(); ++q )
                    {
                        for(size_type m = 0; m < mMaxM(q) ; ++m )
                            M->addMatrix( beta_M[q][m] , M_AD->Mqm()[q][m] );
                    }
                }
            }

            std::vector<vector_ptrtype> F( Nl() );

            //warning : beta_F is a vector of beta_coefficients
            auto beta_F = all_beta.template get<2>();

            for ( size_type l = 0; l < Nl(); ++l )
            {
                F[l] = M_backend->newVector( Xh );
                F[l]->zero();

                for ( size_type q = 0; q < Ql( l ); ++q )
                {
                    for ( size_type m = 0; m < mMaxF(l,q); ++m )
                    {
                        F[l]->add( beta_F[l][q][m] , M_AD->Fqm(l)[q][m] );
                    }
                }
                F[l]->close();
            }

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
            AssembleInitialGuessVInCompositeCase<self_type> assemble_initial_guess_v_in_composite_case ( v , initial_guess , this->shared_from_this());
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
            using namespace Feel::vf;
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



protected:
    bool M_is_initialized;
    CRBModelMode M_mode;
    backend_ptrtype M_backend;
    backend_ptrtype M_backend_primal;
    backend_ptrtype M_backend_dual;
    backend_ptrtype M_backend_l2;

    parameterspace_ptrtype Dmu;
    functionspace_ptrtype Xh;
    rbfunctionspace_ptrtype XN;

    affinedecomposition_ptrtype M_AD;

    sparse_matrix_ptrtype M_inner_product_matrix;
    sparse_matrix_ptrtype M_mass_matrix;

    bdf_ptrtype M_bdf;

    funs_type M_funs;
    funsd_type M_funs_d;

    int M_QLinearDecompositionA;
    operatorcomposite_ptrtype M_compositeA;
    operatorcomposite_ptrtype M_compositeM;
    std::vector< functionalcomposite_ptrtype > M_compositeF;

    sparse_matrix_ptrtype M_monoA;
    sparse_matrix_ptrtype M_monoM;
    std::vector<vector_ptrtype> M_monoF;


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
