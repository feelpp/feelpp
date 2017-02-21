/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-05-02

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file deim.hpp
   \author JB Wahl
   \date 2017-01-27
 */
#ifndef _FEELPP_DEIM_HPP
#define _FEELPP_DEIM_HPP 1

#include <limits>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>

#include <boost/ref.hpp>
#include <boost/next_prior.hpp>
#include <boost/type_traits.hpp>
#include <boost/tuple/tuple.hpp>

#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <Eigen/Core>


using Feel::cout;

namespace Feel
{

template <typename ParameterType>
struct paramCompare
{
public :
    typedef ParameterType parameter_type;

    bool operator() ( parameter_type mu1, parameter_type mu2 )
    {
        return compare( mu1, mu2, 0 );
    }
private :
    bool compare( parameter_type mu1, parameter_type mu2, int i )
    {
        if ( i==mu1.size() )
            return false;
        else if ( mu1[i]==mu2[i] )
            return compare( mu1, mu2, i+1 );
        else
            return mu1[i]<mu2[i];
    }
};


template <typename ModelType, typename TensorType>
class DEIMBase
{
public :
    typedef ModelType model_type;
    typedef typename model_type::value_type value_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;

    typedef typename ModelType::parameterspace_type parameterspace_type;
    typedef typename ModelType::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename ModelType::parameter_type parameter_type;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef TensorType tensor_type;
    typedef boost::shared_ptr<tensor_type> tensor_ptrtype;

    typedef std::function<tensor_ptrtype ( const parameter_type& mu)> assemble_function_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef Eigen::MatrixXd matrixN_type;
    typedef Eigen::VectorXd vectorN_type;

    typedef boost::tuple<parameter_type,double> bestfit_type;

    static const bool is_matrix = std::is_same<tensor_type,sparse_matrix_type>::value;
    typedef typename mpl::if_< mpl::bool_<is_matrix>,
                               std::pair<int,int>, int >::type indice_type;

    typedef boost::tuple<double,int> vectormax_type;
    typedef boost::tuple<double,std::pair<int,int>> matrixmax_type;

    typedef std::map<parameter_type,tensor_ptrtype,paramCompare<parameter_type>> solutionsmap_type;

    DEIMBase()
    {}

    DEIMBase( parameterspace_ptrtype Dmu, sampling_ptrtype sampling ) :
        M_parameter_space( Dmu ),
        M_trainset( sampling ),
        M_M(0),
        M_tol(1e-8)
    {
        if ( !M_trainset )
            M_trainset = Dmu->sampling();
        if ( M_trainset->empty() )
        {
            int sampling_size = ioption(_name="eim.sampling-size");
            std::string file_name = ( boost::format("eim_trainset_%1%") % sampling_size ).str();
            std::string sampling_mode = "log-equidistribute";
            std::ifstream file ( file_name );
            if( ! file )
            {
                M_trainset->sample( sampling_size, sampling_mode, true, file_name );
            }
            else
            {
                M_trainset->clear();
                M_trainset->readFromFile(file_name);
            }
            cout << "DEIM sampling created with " << sampling_size << " points.\n";
        }
    }

    virtual ~DEIMBase()
    {}

    assemble_function_type assemble;

    void run()
    {
        tic();
        int mMax = ioption("eim.dimension-max");
        double error=0;
        auto mu = M_parameter_space->max();

        cout << "DEIM : Start algorithm with mu=["<<mu[0];
        for ( int i=1; i<mu.size(); i++ )
            cout << "," << mu[i];
        cout << "]\n";
        cout <<"===========================================\n";
        do{
            cout << "DEIM : Construction of basis "<<M_M+1<<"/"<<mMax<<", with mu=["<<mu[0];
            for ( int i=1; i<mu.size(); i++ )
                cout << "," << mu[i];
            cout << "]\n";

            tic();
            addNewVector(mu);
            toc("Add new vector in DEIM basis");

            if ( M_M<mMax )
            {
                auto best_fit = computeBestFit();
                error = best_fit.template get<1>();
                mu = best_fit.template get<0>();
                cout << "DEIM : Current error : "<< error
                           <<", tolerance : " << M_tol <<std::endl;
                cout <<"===========================================\n";
            }

        }while( M_M<mMax && error>M_tol );

        cout << "DEIM : Stopping greedy algorithm. Number of basis function : "<<M_M<<std::endl;

        toc("DEIM : Total Time");
    }

    vectorN_type beta( parameter_type const& mu )
    {
        return computeCoefficient( mu );
    }

    std::vector<tensor_ptrtype> q()
    {
        return M_bases;
    }

    int size()
    {
        return M_M;
    }


protected :
    void addNewVector( parameter_type const& mu )
    {
        tensor_ptrtype Phi = residual( mu );

        auto vec_max = vectorMaxAbs( Phi );
        auto i = vec_max.template get<1>();
        double max = vec_max.template get<0>();

        M_M++;

        Phi->scale( 1./max );
        M_bases.push_back( Phi );
        M_index.push_back(i);

        M_B.conservativeResize(M_M,M_M);
        // update last row of M_B
        for ( int j = 0; j<M_M; j++ )
            M_B(M_M-1, j) = evaluate( M_bases[j], M_index[M_M-1]);

        //update last col of M_B
        for ( int i=0; i<M_M-1; i++ )
            M_B(i, M_M-1) = evaluate( M_bases[M_M-1], M_index[i] );
    }

    double evaluate( vector_ptrtype V, int const& index )
    {
        double value=0;
        int proc_number = V->map().procOnGlobalCluster(index);

        if ( Environment::worldComm().globalRank()==proc_number )
            value = V->operator()( index - V->map().firstDofGlobalCluster() );

        boost::mpi::broadcast( Environment::worldComm(), value, proc_number );
        return value;
    }

    double evaluate( sparse_matrix_ptrtype M, std::pair<int,int> const&  idx )
    {
        int i = idx.first;
        int j =idx.second;

        double value=0;
        int proc_number = M->mapRow().procOnGlobalCluster(i);

        if ( !M->closed() )
            M->close();
        if ( Environment::worldComm().globalRank()==proc_number )
            value = M->operator() (i,j);

        boost::mpi::broadcast( Environment::worldComm(), value, proc_number );
        return value;
    }

    vectormax_type vectorMaxAbs( vector_ptrtype V )
    {
        auto newV = V->clone();
        *newV = *V;
        newV->abs();
        int index = 0;
        double max = newV->maxWithIndex( &index );
        double eval=evaluate( newV, index );

        return boost::make_tuple( max, index );
    }

    matrixmax_type vectorMaxAbs( sparse_matrix_ptrtype M )
    {
        DCHECK( M->isInitialized() ) << "MatrixPetsc<> not initialized";
        if ( !M->closed() )
            M->close();

        auto V = backend()->newVector( M->mapRowPtr() );
        PetscInt idx[M->mapRow().nLocalDof()];

        int ierr=MatGetRowMaxAbs( toPETSc(M)->mat(), toPETSc(V)->vec(), idx );
        CHKERRABORT( M->comm(),ierr );

        int i_row = 0;
        int i_col = 0;
        PetscReal val=0;
        VecMax(toPETSc(V)->vec(), &i_row, &val);
        double max = static_cast<Real>( val );

        int proc_number = V->map().procOnGlobalCluster(i_row);
        if ( Environment::worldComm().globalRank()==proc_number )
            i_col = idx[i_row - V->map().firstDofGlobalCluster()];
        boost::mpi::broadcast( Environment::worldComm(), i_col, proc_number );

        std::pair<int,int> index (i_row,i_col);

        return boost::make_tuple( max, index );
    }

    vectorN_type computeCoefficient( parameter_type const& mu )
    {
        tensor_ptrtype T = assemble( mu );
        return computeCoefficient( T );
    }

    vectorN_type computeCoefficient( tensor_ptrtype T )
    {
        vectorN_type rhs (M_M);
        vectorN_type coeff (M_M);
        if ( M_M>0 )
        {
            for ( int i=0; i<M_M; i++ )
                rhs(i) = evaluate( T, M_index[i] );
            coeff = M_B.lu().solve( rhs );
        }
        return coeff;
    }

    bestfit_type computeBestFit()
    {
        tic();
        double max=0;
        auto mu_max = M_parameter_space->element();
        for ( auto const& mu : *M_trainset )
        {
            tensor_ptrtype T = residual( mu );

            double norm = T->linftyNorm();
            if ( norm>max )
            {
                max = norm;
                mu_max = mu;
            }
        }
        toc("DEIM : compute best fit");

        return boost::make_tuple( mu_max, max );
    }

    tensor_ptrtype residual( parameter_type const& mu )
    {
        tensor_ptrtype T;
        if ( !M_solutions[mu] )
        {
            T = this->assemble( mu );
            M_solutions[mu] = copyTensor(T);
        }
        else
            T = M_solutions[mu];

        vectorN_type coeff = computeCoefficient( T );

        auto newT = copyTensor( T );

        for ( int i=0; i<M_M; i++ )
            add( newT, -coeff(i), M_bases[i] );

        return newT;
    }

    void add( vector_ptrtype V, double const& a, vector_ptrtype vec )
    {
        V->add( a, vec );
    }
    void add( sparse_matrix_ptrtype M, double const& a, sparse_matrix_ptrtype mat )
    {
        M->addMatrix( a, mat );
    }

    vector_ptrtype copyTensor( vector_ptrtype V )
    {
        vector_ptrtype newV = V->clone();
        *newV = *V;
        return newV;
    }

    sparse_matrix_ptrtype copyTensor( sparse_matrix_ptrtype M )
    {
        sparse_matrix_ptrtype newM = backend()->newMatrix( M->mapColPtr(),
                                                           M->mapRowPtr(),
                                                           M->graph() );
        *newM=*M;
        return newM;
    }


    parameterspace_ptrtype M_parameter_space;
    sampling_ptrtype M_trainset;
    int M_M;
    double M_tol;
    matrixN_type M_B;
    std::vector< tensor_ptrtype > M_bases;
    std::vector<indice_type> M_index;
    solutionsmap_type M_solutions;
};


template <typename ModelType>
class DEIM :
        public DEIMBase<ModelType, typename Backend<typename ModelType::value_type>::vector_type>
{
public :
    typedef DEIMBase<ModelType, typename Backend<typename ModelType::value_type>::vector_type> super_type;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename super_type::sampling_ptrtype sampling_ptrtype;
    DEIM() :
        super_type()
    {}

    DEIM( parameterspace_ptrtype Dmu, sampling_ptrtype sampling=NULL ) :
        super_type( Dmu, sampling )
    {}

    ~DEIM()
    {}


private :
};


template <typename ModelType>
class MDEIM :
        public DEIMBase<ModelType, typename Backend<typename ModelType::value_type>::sparse_matrix_type>
{
public :
    typedef DEIMBase<ModelType, typename Backend<typename ModelType::value_type>::sparse_matrix_type> super_type;

    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename super_type::sampling_ptrtype sampling_ptrtype;

    MDEIM() :
        super_type()
    {}

    MDEIM( parameterspace_ptrtype Dmu, sampling_ptrtype sampling=NULL ) :
        super_type( Dmu, sampling )
    {}

    ~MDEIM()
    {}

private :

};


    po::options_description eimOptions( std::string const& prefix ="");
}

#endif
