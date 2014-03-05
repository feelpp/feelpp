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

    static const bool is_time_dependent = FunctionSpaceDefinition::is_time_dependent;
    static const bool is_linear = FunctionSpaceDefinition::is_linear;

    typedef double value_type;

    typedef typename FunctionSpaceDefinition::mesh_type mesh_type;

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
                                   beta_vector_type,
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   >,
                               boost::tuple<
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   >
                               >::type betaqm_type;


    ModelCrbBase()
        :
        M_is_initialized( false )
    {
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

    virtual affine_decomposition_type computeAffineDecomposition()
    {
        return computeAffineDecomposition( mpl::bool_< is_time_dependent >() );
    }
    affine_decomposition_type computeAffineDecomposition( mpl::bool_<true> )
    {
        std::vector< std::vector<sparse_matrix_ptrtype> > M;
        std::vector< std::vector<sparse_matrix_ptrtype> > A;
        std::vector< std::vector<std::vector<vector_ptrtype> > > F;
        return boost::make_tuple( M , A , F );
    }
    affine_decomposition_type computeAffineDecomposition( mpl::bool_<false> )
    {
        std::vector< std::vector<sparse_matrix_ptrtype> > A;
        std::vector< std::vector<std::vector<vector_ptrtype> > > F;
        return boost::make_tuple( A , F );
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
     * inner product for mass matrix
     * Transient models need to implement these functions.
     */
    virtual sparse_matrix_ptrtype const& innerProductForMassMatrix () const
    {
        throw std::logic_error("Your model is time-dependant so you MUST implement innerProductForMassMatrix function");
        return M;
    }
    virtual sparse_matrix_ptrtype innerProductForMassMatrix ()
    {
        throw std::logic_error("Your model is time-dependant so you MUST implement innerProductForMassMatrix function");
        return M;
    }

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


    void writeConvergenceStatistics( std::vector< vectorN_type > const& vector, std::string filename )
    {
        if( Environment::worldComm().isMasterRank() )
        {
            Eigen::MatrixXf::Index index_max;
            Eigen::MatrixXf::Index index_min;

            std::ofstream file;
            file.open( filename );
            file << "NbBasis" << "\t" << "Min" << "\t" << "Max" << "\t" << "Mean" << "\t" << "Variance" << "\n";
            int Nmax = vector.size();
            std::vector<double> nbruns( Nmax );
            for(int n=0; n<Nmax; n++)
            {
                double variance=0;
                int sampling_size = vector[n].size();
                double mean=vector[n].mean();
                double max = vector[n].maxCoeff(&index_max);
                double min = vector[n].minCoeff(&index_min);
                for(int i=0; i<sampling_size; i++)
                    variance += (vector[n](i) - mean)*(vector[n](i) - mean)/sampling_size;
                file <<n+1<<"\t"<<min<<"\t"<<max<<"\t"<<mean<<"\t"<<variance<<"\n";
                nbruns[n]=vector[n].size();
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
        double min_scale=std::floor(min_output);
        double x=0;
        double output=0;
        double estimated_down=0;
        double estimated_up=0;
        double delta=0;

        std::ofstream file_outputs_geo_gmsh ( "GMSH-outputs.geo", std::ios::out );
        file_outputs_geo_gmsh << "View \" outputs \" {\n";
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

    operatorcomposite_ptrtype M_compositeA;
    operatorcomposite_ptrtype M_compositeM;
    std::vector< functionalcomposite_ptrtype > M_compositeF;

};

}//Feel
#endif /* __Model_H */
