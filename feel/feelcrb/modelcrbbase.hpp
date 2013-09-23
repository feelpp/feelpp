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
    typedef Lagrange<1,Scalar> basis_type ;

    /*space*/
    typedef FunctionSpace<mesh_type , basis_type > space_type ;
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



template <typename ParameterDefinition, typename FunctionSpaceDefinition, typename EimDefinition = EimDefinitionBase<ParameterDefinition,FunctionSpaceDefinition> >
class ModelCrbBase : public ModelCrbBaseBase
{

public :
    typedef typename EimDefinition::fun_type fun_type;
    typedef typename EimDefinition::fund_type fund_type;

    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename parameterspace_type::element_type parameter_type;

    typedef typename FunctionSpaceDefinition::space_type space_type;

    typedef boost::shared_ptr<fun_type> fun_ptrtype;
    typedef std::vector<fun_ptrtype> funs_type;

    typedef boost::shared_ptr<fund_type> fund_ptrtype;
    typedef std::vector<fund_ptrtype> funsd_type;

    typedef Eigen::VectorXd vectorN_type;

    typedef std::vector< std::vector< double > > beta_vector_type;

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
     * returns the scalar product used fior mass matrix ( to solve eigen values problem )
     * of the boost::shared_ptr vector x and boost::shared_ptr vector
     * Transient models need to implement these functions.
     */
    virtual double scalarProductForMassMatrix( vector_ptrtype const& X, vector_ptrtype const& Y )
    {
        throw std::logic_error("Your model is time-dependant so you MUST implement scalarProductForMassMatrix function");
        return 0;
    }
    virtual double scalarProductForMassMatrix( vector_type const& x, vector_type const& y )
    {
        throw std::logic_error("Your model is time-dependant so you MUST implement scalarProductForMassMatrix function");
        return 0;
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

        if( vector.size() > 0 )
        {
            bool force = option("eim.use-dimension-max-functions").template as<bool>();
            int Neim=0;
            if( force )
                Neim = option("eim.dimension-max").template as<int>();

            int N = vector.size();

            if( force )
                LOG( INFO ) <<" statistics  for "<<name<<" (  was called "<< N << " times with "<<Neim<<" basis functions )";
            else
                LOG( INFO ) <<" statistics  for "<<name<<" (  was called "<< N << " times )";

            min = vector.minCoeff(&index);
            max = vector.maxCoeff(&index);
            mean = vector.mean();
            mean1 = mean * mean;
            square  = vector.array().pow(2);
            mean2 = square.mean();
            standard_deviation = math::sqrt( mean2 - mean1 );
            LOG(INFO)<<"min : "<<min<<" - max : "<<max<<" mean : "<<mean<<" standard deviation : "<<standard_deviation;
        }
        vectorN_type result(4);
        result(0)=min;
        result(1)=max;
        result(2)=mean;
        result(3)=standard_deviation;
        return result;
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
