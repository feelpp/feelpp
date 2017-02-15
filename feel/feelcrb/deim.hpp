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

namespace Feel
{

template <typename ModelType, typename TensorType>
class DEIMBase
{
public :
    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;

    typedef typename ModelType::parameterspace_type parameterspace_type;
    typedef typename ModelType::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename ModelType::parameter_type parameter_type;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef TensorType tensor_type;
    typedef boost::shared_ptr<tensor_type> tensor_ptrtype;

    typedef std::function<tensor_ptrtype ( const parameter_type& mu)> assemble_function_type;

    typedef Eigen::MatrixXd matrixN_type;
    typedef Eigen::VectorXd vectorN_type;

    typedef boost::tuple<parameter_type,double> bestfit_type;

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
            Feel::cout << "DEIM sampling created with " << sampling_size << " points.\n";
        }
    }

    virtual ~DEIMBase()
    {}

    assemble_function_type assemble;

    void run()
    {
        int mMax = ioption("eim.dimension-max");
        auto mu = M_parameter_space->max();

        addNewVector(mu);

        tic();
        for ( ; M_M<=mMax; )
        {
            // choose next mu
            auto best_fit = computeBestFit();

            Feel::cout << "DEIM : Current error : "<< best_fit.template get<1>()
                       <<", tolerance : " << M_tol <<std::endl;
            if ( best_fit.template get<1>()<M_tol )
            {
                Feel::cout << "DEIM : Stopping greedy algorithm\n";
                break;
            }

            mu = best_fit.template get<0>();
            addNewVector(mu);

        }
        toc("DEIM : greedy");
        Feel::cout << "DEIM number of basis function vector : "<< M_M << std::endl;
    }

    vectorN_type beta( parameter_type mu )
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
    virtual void addNewVector( parameter_type mu )=0;
    virtual vectorN_type computeCoefficient( tensor_ptrtype V )=0;
    virtual tensor_ptrtype residual( parameter_type mu )=0;

    vectorN_type computeCoefficient( parameter_type mu )
    {
        tensor_ptrtype T = assemble( mu );
        return computeCoefficient( T );
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




    parameterspace_ptrtype M_parameter_space;
    sampling_ptrtype M_trainset;
    int M_M;
    double M_tol;
    matrixN_type M_B;
    std::vector< tensor_ptrtype > M_bases;
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
    typedef typename super_type::vectorN_type vectorN_type;
    typedef typename super_type::tensor_ptrtype tensor_ptrtype;

    typedef typename ModelType::value_type value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef boost::tuple<double,int> vectormax_type;

    DEIM() :
        super_type()
    {}

    DEIM( parameterspace_ptrtype Dmu, sampling_ptrtype sampling=NULL ) :
        super_type( Dmu, sampling )
    {}

    ~DEIM()
    {}


private :
    void addNewVector( parameter_type mu )
    {
        vector_ptrtype Phi = residual( mu );

        auto vec_max = vectorMaxAbs( Phi );
        int i = vec_max.template get<1>();
        double max = vec_max.template get<0>();

        M_M++;

        Phi->scale( 1./max );
        M_bases.push_back( Phi );
        M_index.push_back(i);

        M_B.conservativeResize(M_M,M_M);
        // update last row of M_B
        for ( int j = 0; j<M_M; j++ )
        {
            M_B(M_M-1, j) = M_bases[j]->operator() (M_index[M_M-1] );
        }
        //update last col of M_B
        for ( int i=0; i<M_M-1; i++ )
        {
            M_B(i, M_M-1) = M_bases[M_M-1]->operator() (M_index[i]);
        }
    }

    vectormax_type vectorMaxAbs( vector_ptrtype V )
    {
        auto newV = V->clone();
        *newV = *V;
        newV->abs();
        int index = 0;
        double max = newV->maxWithIndex( &index );
        return boost::make_tuple( max, index );
    }

    tensor_ptrtype residual( parameter_type mu )
    {
        vector_ptrtype V = this->assemble( mu );
        vectorN_type coeff = computeCoefficient( V );

        vector_ptrtype newV = V->clone();
        *newV = *V;

        for ( int i=0; i<M_M; i++ )
            newV->add( -coeff(i), M_bases[i] );

        return newV;
    }

    using super_type::computeCoefficient;

    vectorN_type computeCoefficient( tensor_ptrtype V )
    {
        vectorN_type rhs (M_M);
        vectorN_type coeff (M_M);
        if ( M_M>0 )
        {
            for ( int i=0; i<M_M; i++ )
                rhs(i) = V->operator() ( M_index[i] );
            coeff = M_B.lu().solve( rhs );
        }
        return coeff;
    }

    std::vector<int> M_index;

    using super_type::M_M;
    using super_type::M_B;
    using super_type::M_bases;
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
    typedef typename super_type::vectorN_type vectorN_type;
    typedef typename super_type::tensor_ptrtype tensor_ptrtype;


    typedef typename ModelType::value_type value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;


    typedef boost::tuple<double,std::pair<int,int>> matrixmax_type;

    MDEIM() :
        super_type()
    {}

    MDEIM( parameterspace_ptrtype Dmu, sampling_ptrtype sampling=NULL ) :
        super_type( Dmu, sampling )
    {}

    ~MDEIM()
    {}

private :
    void addNewVector( parameter_type mu )
    {
        sparse_matrix_ptrtype Phi = residual( mu );

        auto mat_max = matrixMaxAbs( Phi );
        auto idx = mat_max.template get<1>();
        double max = mat_max.template get<0>();

        M_M++;

        Phi->scale( 1./max );
        M_bases.push_back( Phi );
        M_index.push_back( idx );

        M_B.conservativeResize(M_M,M_M);
        for ( int j = 0; j<M_M; j++ )
        {
            M_B(M_M-1, j) = M_bases[j]->operator() ( M_index[M_M-1].first, M_index[M_M-1].second );
        }
        //update last col of M_B
        for ( int i=0; i<M_M-1; i++ )
        {
            M_B(i, M_M-1) = M_bases[M_M-1]->operator() ( M_index[i].first, M_index[i].second );
        }

    }

    matrixmax_type matrixMaxAbs( sparse_matrix_ptrtype M )
    {
        DCHECK( M->isInitialized() ) << "MatrixPetsc<> not initialized";
        if ( !M->closed() )
            M->close();

        vector_ptrtype V = M->diagonal();
        int ierr=0;
        int nrow = M->size1();
        int idx [nrow];

        ierr=MatGetRowMaxAbs( toPETSc(M)->mat(), toPETSc(V)->vec(), idx );
        CHKERRABORT( M->comm(),ierr );

        int i_row = 0;
        double max = V->maxWithIndex( &i_row );
        int i_col = idx[i_row];
        std::pair<int,int> index (i_row,i_col);

        return boost::make_tuple( max, index );
    }

    tensor_ptrtype residual( parameter_type mu )
    {
        sparse_matrix_ptrtype M = this->assemble( mu );
        vectorN_type coeff = computeCoefficient( M );

        sparse_matrix_ptrtype newM = backend()->newMatrix( M->mapRowPtr(), M->mapColPtr() );
        //Mat* m;
        MatConvert(toPETSc(M)->mat(), MATSAME, MAT_INITIAL_MATRIX, &toPETSc(newM)->mat() );
        //sparse_matrix_ptrtype newM ( new MatrixPetsc<value_type>( *m ) );
        //auto newM = boost::make_shared<decltype(*M)>(*M);

        for ( int i=0; i<M_M; i++ )
            newM->addMatrix( -coeff(i), M_bases[i] );

        return newM;
    }

    using super_type::computeCoefficient;

    vectorN_type computeCoefficient( tensor_ptrtype M )
    {
        vectorN_type rhs (M_M);
        vectorN_type coeff (M_M);
        if ( M_M>0 )
        {
            for ( int i=0; i<M_M; i++ )
                rhs(i) = M->operator() ( M_index[i].first, M_index[i].second );
            coeff = M_B.lu().solve( rhs );
        }
        return coeff;
    }


    std::vector<std::pair<int,int>> M_index;

    using super_type::M_M;
    using super_type::M_B;
    using super_type::M_bases;

};


po::options_description eimOptions( std::string const& prefix ="");
}

#endif
