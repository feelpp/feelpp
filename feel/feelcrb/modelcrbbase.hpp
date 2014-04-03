/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Stephane Veys <stephane.veys@imag.fr>
 Date: 2013-02-22

 Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
/**
   \file model.hpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-02-22
*/
#ifndef ModelCrbBase_H
#define ModelCrbBase_H 1

#include <feel/feel.hpp>
#include <feel/feelcrb/eim.hpp>

namespace Feel
{



class ParameterDefinitionBase
{
public :
    typedef ParameterSpace<1> parameterspace_type ;
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

    static const bool is_time_dependent=false;
    static const bool is_linear=true;
};

template <typename ParameterDefinition, typename FunctionSpaceDefinition>
class EimDefinitionBase
{

public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename FunctionSpaceDefinition::space_type space_type;

    /* EIM */
    typedef EIMFunctionBase<space_type , space_type  , parameterspace_type > fun_type ;
    typedef EIMFunctionBase<space_type , space_type  , parameterspace_type > fund_type ;

};



template <typename ParameterDefinition=ParameterDefinitionBase, typename FunctionSpaceDefinition=FunctionSpaceDefinitionBase, typename EimDefinition = EimDefinitionBase<ParameterDefinition,FunctionSpaceDefinition> >
class ModelCrbBase : public ModelCrbBaseBase
{

public :

    typedef ModelCrbBase self_type;
    typedef typename EimDefinition::fun_type fun_type;
    typedef typename EimDefinition::fund_type fund_type;

    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename parameterspace_type::element_type parameter_type;

    typedef typename FunctionSpaceDefinition::space_type space_type;
    typedef space_type functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    typedef boost::shared_ptr<fun_type> fun_ptrtype;
    typedef std::vector<fun_ptrtype> funs_type;

    typedef boost::shared_ptr<fund_type> fund_ptrtype;
    typedef std::vector<fund_ptrtype> funsd_type;

    typedef Eigen::VectorXd vectorN_type;

    typedef OperatorLinear< space_type , space_type > operator_type;
    typedef boost::shared_ptr<operator_type> operator_ptrtype;

    typedef OperatorLinearComposite< space_type , space_type > operatorcomposite_type;
    typedef boost::shared_ptr<operatorcomposite_type> operatorcomposite_ptrtype;

    typedef FsFunctionalLinearComposite< space_type > functionalcomposite_type;
    typedef boost::shared_ptr<functionalcomposite_type> functionalcomposite_ptrtype;

    typedef FsFunctionalLinear< space_type > functional_type;
    typedef boost::shared_ptr<functional_type> functional_ptrtype;

    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    static const uint16_type ParameterSpaceDimension = ParameterDefinition::ParameterSpaceDimension ;

    typedef std::vector< std::vector< double > > beta_vector_type;
    typedef std::vector< double > beta_vector_light_type;


    static const bool is_time_dependent = FunctionSpaceDefinition::is_time_dependent;
    static const bool is_linear = FunctionSpaceDefinition::is_linear;

    typedef double value_type;

    typedef typename FunctionSpaceDefinition::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   std::vector< std::vector<sparse_matrix_ptrtype> >,//Mq
                                   std::vector< std::vector<sparse_matrix_ptrtype> >,//Aq
                                   std::vector< std::vector<std::vector<vector_ptrtype> > >//Fq
                                   >,
                               boost::tuple<
                                   std::vector< std::vector<sparse_matrix_ptrtype> >,//Aq
                                   std::vector< std::vector<std::vector<vector_ptrtype> > >//Fq
                                   >
                               >::type affine_decomposition_type;

    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   std::vector< sparse_matrix_ptrtype >,//Mq
                                   std::vector< sparse_matrix_ptrtype >,//Aq
                                   std::vector< std::vector< vector_ptrtype > >//Fq
                                   >,
                               boost::tuple<
                                   std::vector< sparse_matrix_ptrtype >,//Aq
                                   std::vector< std::vector< vector_ptrtype > >//Fq
                                   >
                               >::type affine_decomposition_light_type;

    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   sparse_matrix_ptrtype, sparse_matrix_ptrtype, std::vector<vector_ptrtype>
                                   >,
                               boost::tuple<
                                   sparse_matrix_ptrtype, std::vector<vector_ptrtype>
                                   >
                               >::type monolithic_type;


    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   beta_vector_type,
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   >,
                               boost::tuple<
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   >
                               >::type betaqm_type;

    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   beta_vector_light_type,
                                   beta_vector_light_type,
                                   std::vector<beta_vector_light_type>
                                   >,
                               boost::tuple<
                                   beta_vector_light_type,
                                   std::vector<beta_vector_light_type>
                                   >
                               >::type betaq_type;

    typedef std::vector< boost::tuple< sparse_matrix_ptrtype , std::string > > lhs_light_type;
    typedef boost::shared_ptr<lhs_light_type> lhs_light_ptrtype;
    typedef std::vector< boost::tuple<  vector_ptrtype , std::string > > rhs_light_type;
    typedef boost::shared_ptr<rhs_light_type> rhs_light_ptrtype ;
    typedef std::vector< boost::tuple<  std::vector< vector_ptrtype > , std::string > > outputs_light_type;
    typedef boost::shared_ptr<outputs_light_type> outputs_light_ptrtype ;

    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    ModelCrbBase()
        :
        M_is_initialized( false )
    {
    }




    /*
     * return the left hand side terms from the affine decomposition
     * that is to say bilinear forms Aq and beta coefficients associated
     */
    void addLhs( boost::tuple< sparse_matrix_ptrtype , std::string > const & tuple )
    {
        M_Aq.push_back( tuple.template get<0>() );
        M_betaAqString.push_back( tuple.template get<1>() );
    }

    /*
     * return the terms from the affine decomposition linked to the mass matrix
     * that is to say bilinear forms Mq and beta coefficients associated
     */
    void addMass( boost::tuple< sparse_matrix_ptrtype , std::string > const & tuple )
    {
        M_Mq.push_back( tuple.template get<0>() );
        M_betaMqString.push_back( tuple.template get<1>() );
    }

    /*
     * return the right hand side terms from the affine decomposition
     * that is to say linear forms Fq[0] and beta coefficients associated
     */
    void addRhs( boost::tuple<  vector_ptrtype , std::string > const & tuple )
    {
        if( M_Fq.size() == 0 )
        {
            M_Fq.resize(2);
            M_betaFqString.resize(2);
        }
        M_Fq[0].push_back( tuple.template get<0>() );
        M_betaFqString[0].push_back( tuple.template get<1>() );
    }

    /*
     * return the terms from the affine decomposition linked to output
     * that is to say linear forms and beta coefficients associated
     */
    void addOutput( boost::tuple<  vector_ptrtype , std::string > const & tuple )
    {
        int size = M_Fq.size();
        if( size == 0 )
        {
            //there is no rhs or output yet
            M_Fq.resize(2);
            M_betaFqString.resize(2);
        }
        M_Fq[1].push_back( tuple.template get<0>() );
        M_betaFqString[1].push_back( tuple.template get<1>() );
    }

    void addEnergyMatrix( sparse_matrix_ptrtype const & matrix )
    {
        M_energy_matrix = matrix ;
    }

    void addMassMatrix( sparse_matrix_ptrtype const & matrix )
    {
        M_mass_matrix = matrix ;
    }

    void setInitialized( const bool & b )
    {
        M_is_initialized = b ;
    }

    bool isInitialized()
    {
        return M_is_initialized;
    }

    virtual funs_type scalarContinuousEim () const
    {
        return M_funs;
    }

    virtual funsd_type scalarDiscontinuousEim () const
    {
        return M_funs_d;
    }

    virtual void initModel() = 0;

    /*
     * the user has to provide the affine decomposition
     * he has two choices : use vectors/matrices or operator free
     * if he uses operators free he must implement :
     * - operatorCompositeM (only if the model is time dependent)
     * - operatorCompositeA
     * - functionalCompositeF
     * else, if he uses vectors/matrices (classical way) he must implement
     * - computeAffineDecomposition
     */

    virtual operatorcomposite_ptrtype operatorCompositeA()
    {
        return M_compositeA;
    }
    virtual operatorcomposite_ptrtype operatorCompositeM()
    {
        return M_compositeM;
    }
    virtual std::vector< functionalcomposite_ptrtype > functionalCompositeF()
    {
        return M_compositeF;
    }

    virtual operatorcomposite_ptrtype operatorCompositeLightA()
    {
        return M_compositeA;
    }
    virtual operatorcomposite_ptrtype operatorCompositeLightM()
    {
        return M_compositeM;
    }
    virtual std::vector< functionalcomposite_ptrtype > functionalCompositeLightF()
    {
        return M_compositeF;
    }

    virtual affine_decomposition_type computeAffineDecomposition()
    {
        return computeAffineDecomposition( mpl::bool_< is_time_dependent >() );
    }
    affine_decomposition_type computeAffineDecomposition( mpl::bool_<true> )
    {
        return boost::make_tuple( M_Mqm , M_Aqm , M_Fqm );
    }
    affine_decomposition_type computeAffineDecomposition( mpl::bool_<false> )
    {
        return boost::make_tuple( M_Aqm , M_Fqm );
    }

    virtual affine_decomposition_light_type computeAffineDecompositionLight()
    {
        return computeAffineDecompositionLight( mpl::bool_< is_time_dependent >() );
    }
    affine_decomposition_light_type computeAffineDecompositionLight( mpl::bool_<true> )
    {
        return boost::make_tuple( M_Mq , M_Aq , M_Fq );
    }
    affine_decomposition_light_type computeAffineDecompositionLight( mpl::bool_<false> )
    {
        return boost::make_tuple( M_Aq , M_Fq );
    }

    virtual monolithic_type computeMonolithicFormulation( parameter_type const& mu )
    {
        return computeMonolithicFormulation( mu , mpl::bool_< is_time_dependent >() );
    }
    monolithic_type computeMonolithicFormulation( parameter_type const& mu, mpl::bool_<true> )
    {
        return boost::make_tuple( M_monoM , M_monoA , M_monoF );
    }
    monolithic_type computeMonolithicFormulation( parameter_type const& mu, mpl::bool_<false> )
    {
        return boost::make_tuple( M_monoA , M_monoF );
    }


    virtual betaq_type computeBetaQ( parameter_type const& mu ,  double time=0 )
    {
        return computeBetaQ( mu, mpl::bool_< is_time_dependent >(), time );
    }
    betaq_type computeBetaQ( parameter_type const& mu , mpl::bool_<true>, double time=0 )
    {
        int sizeA=M_ginacAq.size();
        if( sizeA > 0 )
        {
            //if the user doesn't implement computeBetaQ function
            //and so gives affne decomposition terms using addLhs ect...
            //then we consider we have 2 outputs
            int nboutputs=2;
            int sizeM=M_ginacMq.size();
            int sizeF0=M_ginacFq[0].size();
            int sizeF1=M_ginacFq[1].size();
            int musize=mu.size();
            std::string symbol;

            std::map<std::string,double> map_symbols;
            for(int i=0; i<musize; i++)
            {
                symbol = ( boost::format("mu%1%") %i ).str();
                map_symbols.insert( std::pair< std::string, double > (symbol,mu(i)) );
            }

            M_betaAq.resize( sizeA );
            for(int q=0; q<sizeA; q++)
            {
                double coeff = M_ginacAq[q].evaluate(map_symbols);
                M_betaAq[q]=coeff;
            }
            M_betaMq.resize( sizeM );
            for(int q=0; q<sizeM; q++)
            {
                double coeff = M_ginacMq[q].evaluate( map_symbols );
                M_betaMq[q]=coeff;
            }
            M_betaFq.resize( nboutputs );
            for(int output=0; output<nboutputs;output++)
            {
                int size=M_betaFqString[output].size();
                M_betaFq[output].resize( size );
                for(int q=0; q<size; q++)
                {
                    double coeff = M_ginacFq[output][q].evaluate( map_symbols );
                    M_betaFq[output][q]=coeff;
                }
            }
        }

        return boost::make_tuple( M_betaMq, M_betaAq, M_betaFq );
    }
    betaq_type computeBetaQ( parameter_type const& mu , mpl::bool_<false>, double time=0 )
    {
        int sizeA=M_ginacAq.size();

        if( sizeA > 0 )
        {
            //if the user doesn't implement computeBetaQ function
            //and so gives affne decomposition terms using addLhs ect...
            //then we consider we have 2 outputs
            int nboutputs=2;
            int sizeM=M_ginacMq.size();
            int sizeF0=M_ginacFq[0].size();
            int sizeF1=M_ginacFq[1].size();
            int musize=mu.size();
            std::string symbol;

            std::map<std::string,double> map_symbols;
            for(int i=0; i<musize; i++)
            {
                symbol = ( boost::format("mu%1%") %i ).str();
                map_symbols.insert( std::pair< std::string, double > (symbol,mu(i)) );
            }

            M_betaAq.resize( sizeA );
            for(int q=0; q<sizeA; q++)
            {
                double coeff = M_ginacAq[q].evaluate(map_symbols);
                M_betaAq[q]=coeff;
            }
            M_betaFq.resize( nboutputs );
            for(int output=0; output<nboutputs;output++)
            {
                int size=M_betaFqString[output].size();
                M_betaFq[output].resize( size );
                for(int q=0; q<size; q++)
                {
                    double coeff = M_ginacFq[output][q].evaluate( map_symbols );
                    M_betaFq[output][q]=coeff;
                }
            }
        }

        return boost::make_tuple( M_betaAq, M_betaFq );
    }

    virtual betaqm_type computeBetaQm( parameter_type const& mu ,  double time=0 )
    {
        return computeBetaQm( mu, mpl::bool_< is_time_dependent >(), time );
    }
    betaqm_type computeBetaQm( parameter_type const& mu , mpl::bool_<true>, double time=0 )
    {
        return boost::make_tuple( M_betaMqm, M_betaAqm, M_betaFqm );
    }
    betaqm_type computeBetaQm( parameter_type const& mu , mpl::bool_<false>, double time=0 )
    {
        return boost::make_tuple( M_betaAqm, M_betaFqm );
    }


    void buildGinacBetaExpressions( parameter_type const& mu )
    {
        //not that the parameter mu is here to indicates
        //the dimension of the parameter space
        //and to construct symbols
        int musize = mu.size();
        for(int i=0; i<musize; i++)
        {
            std::string symbol = ( boost::format("mu%1%") %i ).str();
            M_symbols_vec.push_back( symbol );
        }
        if( M_betaAqString.size() > 0 )
        {
            int size = M_betaAqString.size();
            for(int q=0; q<size; q++)
            {
                std::string filename = ( boost::format("GinacA%1%") %q ).str();
                M_ginacAq.push_back( expr( M_betaAqString[q], Symbols( M_symbols_vec ) , filename ) );
            }
        }
        if( M_betaMqString.size() > 0 )
        {
            int size = M_betaMqString.size();
            for(int q=0; q<size; q++)
            {
                std::string filename = ( boost::format("GinacM%1%") %q ).str();
                M_ginacMq.push_back( expr( M_betaMqString[q], Symbols( M_symbols_vec ) , filename ) );
            }
        }
        if( M_betaFqString.size() > 0 )
        {
            int nboutputs = M_betaFqString.size();
            M_ginacFq.resize( nboutputs );
            for(int output=0; output<nboutputs; output++)
            {
                int size = M_betaFqString[output].size();
                for(int q=0; q<size; q++)
                {
                    std::string filename = ( boost::format("GinacF%1%.%2%") %output %q ).str();
                    M_ginacFq[output].push_back( expr( M_betaFqString[output][q], Symbols( M_symbols_vec ) , filename ) );
                }
            }
        }
    }


    //this function is not called bdf() to not interfere with bdf constructor
    virtual bdf_ptrtype bdfModel()
    {
        if( is_time_dependent )
        {
            if( Environment::worldComm().isMasterRank() )
            {
                std::cout<<"*******************************************************************"<<std::endl;
                std::cout<<"** Error ! You have implemented a transient problem but you      **"<<std::endl;
                std::cout<<"** forgot to implement bdfModel() function that returns your bdf **"<<std::endl;
                std::cout<<"*******************************************************************"<<std::endl;
            }
            bool go=false;
            CHECK( go );
        }
        bdf_ptrtype dummy_bdf;
        return dummy_bdf;
    }

    //initialize the field for transient problems
    virtual void initializationField( element_ptrtype& initial_field,parameter_type const& mu )
    {
        initial_field->setZero();
    }


    //for linear models, beta coefficients don't depend on solution u
    //so the user doesn't have to specify this function
    virtual betaqm_type computeBetaQm( element_type const& u, parameter_type const& mu ,  double time=0 )
    {
        if( Environment::worldComm().isMasterRank() )
        {
            std::cout<<"****************************************************************"<<std::endl;
            std::cout<<"** You are using the function computeBetaQm( u , mu ) whereas **"<<std::endl;
            std::cout<<"** your model has only implemented computeBetaQm( mu )        **"<<std::endl;
            std::cout<<"****************************************************************"<<std::endl;
        }
        betaqm_type dummy_beta_coeff;
        return dummy_beta_coeff;
    }



    /**
     * By default reference parameters are min parameters
     * If the user needs to specify them
     * then this function should returns true
     */
    virtual bool referenceParametersGivenByUser()
    {
        return false;
    }
    virtual parameter_type refParameter()
    {
        if( Environment::worldComm().isMasterRank() )
        {
            std::cout<<"**************************************************************************************************"<<std::endl;
            std::cout<<"** You want to specify reference parameters because referenceParametersGivenByUser returns true **"<<std::endl;
            std::cout<<"** your must impelement the function refParameter() !                                           **"<<std::endl;
            std::cout<<"**************************************************************************************************"<<std::endl;
        }
        bool go=false;
        CHECK( go );
        parameter_type muref;
        return muref;
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
    virtual bool constructOperatorCompositeM()
    {
        return false;
    }



    /**
     * note about the initial guess :
     * by default, if the model doesn't give an initial guess
     * then the initial guess is zero
     */
    virtual int QInitialGuess() const
    {
        return 1;
    }

    virtual int mMaxInitialGuess( int q ) const
    {
        if( q == 0 )
            return 1;
        else
            throw std::logic_error( "[ModelCrbBase::mMaxInitialGuess(q)] ERROR wrong index q, should be 0" );
    }

    virtual beta_vector_type computeBetaInitialGuess( parameter_type const& mu )
    {
        beta_vector_type beta;
        beta.resize(1); //q=1
        beta[0].resize(1); //m=1
        beta[0][0]=0;
        return beta;
    }



    /**
     * solve the model for a given parameter \p mu
     * linear models don't need to provide this fuction
     * but non-linear : yes
     * and also those using CRBTrilinear
     */
    virtual element_type solve( parameter_type const& mu )
    {
        element_type solution;
        return solution;
    }


    /**
     * inner product
     */
    virtual sparse_matrix_ptrtype energyMatrix ()
    {
        double norm = M_energy_matrix->l1Norm();
        CHECK( norm > 0 )<<"The energy matrix has not be filled !\n";
        return M_energy_matrix;
    }

    virtual sparse_matrix_ptrtype energyMatrix () const
    {
        double norm = M_energy_matrix->l1Norm();
        CHECK( norm > 0 )<<"The energy matrix has not be filled !\n";
        return M_energy_matrix;
    }

    /**
     * inner product for mass matrix
     * Transient models need to implement these functions.
     */
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
#if 0
    virtual sparse_matrix_ptrtype const& innerProductForPod () const
    {
        throw std::logic_error("Your model is time-dependant so you MUST implement innerProductForPod function");
        return M;
    }
    virtual sparse_matrix_ptrtype innerProductForPod ()
    {
        throw std::logic_error("Your model is time-dependant so you MUST implement innerProductForPod function");
        return M;
    }
#endif

    /*
     * If the model is nonlinear and then need to implement an initial guess
     * it has to implement it !
     * Here is just a code to compile if model is a linear model.
     */
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


    /**
     * compute statistics on vectors
     * arguments : vectors of double and associated names
     * results : statistics on data
     */
    vectorN_type computeStatistics( Eigen::VectorXd vector , std::string name )
    {
        double min=0,max=0,mean=0,mean1=0,mean2=0,standard_deviation=0,variance=0;
        Eigen::MatrixXf::Index index;
        Eigen::VectorXd square;

        std::ofstream file;
        std::string filename="OnlineStatistics-"+name;

        file.open( filename , std::ios::app );

        if( vector.size() > 0 )
        {
            bool force = option("eim.use-dimension-max-functions").template as<bool>();
            int Neim=0;
            if( force )
                Neim = option("eim.dimension-max").template as<int>();

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


    void writeConvergenceStatistics( std::vector< vectorN_type > const& vector, std::string filename , std::string extra="")
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

    /**
    * \param : vector1
    * \param : vector2
    * look for the index minIdx for which we have min( vector1 )
    * then compute vector2( minIdx) / vector1( minIdx )
    * and same thing with max : i.e. vector2( maxIdx ) / vector1( maxIdx )
    * usefull to compute error estimation efficiency associated to
    * min/max errors
    **/
    void writeVectorsExtremumsRatio( std::vector< vectorN_type > const& error, std::vector< vectorN_type > const& estimated, std::string filename )
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

    void readConvergenceDataFromFile( std::vector< vectorN_type > & vector, std::string filename )
    {
        if( Environment::worldComm().isMasterRank() )
        {

            std::vector< std::vector< double > > tmpvector;
            std::ifstream file ( filename );
            std::string str;
            int N;
            int Nmax=0;
            double value;
            if( file )
            {
                //first, determine max elements of the RB
                //because we could have performed several runs
                //and the size of the RB can vary between two runs
                while( ! file.eof() )
                {
                    file >> N;
                    if( N > Nmax )
                        Nmax = N;
                    for(int n=0; n<N; n++)
                    {
                        file >> value;
                    }
                }

                vector.resize( Nmax );
                tmpvector.resize( Nmax );
                //go to the begining of the file
                file.clear();
                file.seekg(0,std::ios::beg);


                file >> N;
                while( ! file.eof() )
                {
                    for(int n=0; n<N; n++)
                    {
                        file >> value;
                        tmpvector[n].push_back( value );
                    }
                    file >> N;
                }

            }
            else
            {
                std::cout<<"The file "<<filename<<" was not found "<<std::endl;
                throw std::logic_error( "[ModelCrbBase::fillVectorFromFile] ERROR loading the file " );
            }

            file.close();

            //now copy std::vector into eigen vector
            for(int n=0; n<Nmax; n++)
            {
                int nbvalues = tmpvector[n].size();
                vector[n].resize(nbvalues);
                for(int i=0; i<nbvalues; i++)
                {
                    vector[n](i) = tmpvector[n][i];
                }
            }

        }//master proc
    }//end of function


    /*
     * if the model has geometric parameters then the mesh can be adapted to
     * the current parameter to visualize solution field
     */
    void adaptMesh( parameter_type const& mu ){ /*by default nothing to be done*/ ; }

    /*
     * \param components_vary : vector of indices for components vary
     * \param extremums : vector containing min and max parameters valuers
     * \param cuttings : vector containing the cutting in each direction + time initial and time final and time step used
     * \param time_vary : (bool) the time vary if true
     */
    gmsh_ptrtype createStructuredGrid( std::vector<int> components_vary, std::vector<parameter_type> extremums, std::vector<int> cuttings,
                                       std::vector<double> time_cuttings, bool time_vary)
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

    /*
     * \param filename : name of the file to be generated
     * \param outputs : vector containing outputs
     * \param parameter : vector containing parameter values
     * \param estimated_error : vector containing estimated error on outputs
     *
     */
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


    /**
     * specific interface for OpenTURNS
     *
     * \param X input vector of size N
     * \param N size of input vector X
     * \param Y input vector of size P
     * \param P size of input vector Y
     */
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

protected :

    funs_type M_funs;
    funsd_type M_funs_d;
    bool M_is_initialized;

    sparse_matrix_ptrtype M;
    sparse_matrix_ptrtype M_energy_matrix;
    sparse_matrix_ptrtype M_mass_matrix;

    operatorcomposite_ptrtype M_compositeA;
    operatorcomposite_ptrtype M_compositeM;
    std::vector< functionalcomposite_ptrtype > M_compositeF;

    std::vector<sparse_matrix_ptrtype> M_Aq;
    std::vector<sparse_matrix_ptrtype> M_Mq;
    std::vector<std::vector<vector_ptrtype> > M_Fq;

    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > M_Mqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;

    sparse_matrix_ptrtype M_monoA;
    sparse_matrix_ptrtype M_monoM;
    std::vector<vector_ptrtype> M_monoF;

    beta_vector_light_type M_betaAq;
    beta_vector_light_type M_betaMq;
    std::vector<beta_vector_light_type> M_betaFq;

    std::vector< std::string > M_betaAqString;
    std::vector< std::string > M_betaMqString;
    std::vector< std::vector< std::string > > M_betaFqString;
    std::vector< Expr<GinacEx<2> > > M_ginacAq;
    std::vector< Expr<GinacEx<2> > > M_ginacMq;
    std::vector< std::vector< Expr<GinacEx<2> > > > M_ginacFq;
    std::vector< std::string > M_symbols_vec;

    beta_vector_type M_betaAqm;
    beta_vector_type M_betaMqm;
    std::vector<beta_vector_type> M_betaFqm;

    lhs_light_type M_lhs;
    lhs_light_type M_mass;
    rhs_light_type M_rhs;
    rhs_light_type M_output;
};

}//Feel
#endif /* __Model_H */
