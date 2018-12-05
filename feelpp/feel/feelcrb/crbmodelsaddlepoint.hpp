/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-08-09

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
 * @file   crbmodelsaddlepoint.hpp
 * @author wahl
 * @date   Tue Nov 25 16:24:55 2014
 */


#ifndef __CRBModelSaddlePoint_H
#define __CRBModelSaddlePoint_H 1

#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/crbmodel.hpp>

namespace Feel
{
template<typename ModelType>
class CRBModelSaddlePoint :
        public CRBModel<ModelType>
{
    typedef CRBModel<ModelType> super;
public :
    typedef ModelType model_type;
    typedef std::shared_ptr<ModelType> model_ptrtype;

    typedef typename model_type::value_type value_type;
    typedef typename model_type::parameter_type parameter_type;

    typedef typename ModelType::mesh_type mesh_type;
    typedef typename ModelType::mesh_ptrtype mesh_ptrtype;

    typedef typename model_type::space_type space_type;
    typedef typename model_type::element_type element_type;

    template <int T>
    using subspace_type = typename space_type::template sub_functionspace<T>::type;
    template <int T>
    using subspace_ptrtype = std::shared_ptr<subspace_type<T>>;
    template<int T>
    using subelement_type = typename subspace_type<T>::element_type;

    CRBModelSaddlePoint( crb::stage stage, int level = 0 ) :
        super ( stage, level ),
        M_block_initialized( false )
    {
        this->initModelSaddlePoint();
    }

    CRBModelSaddlePoint( model_ptrtype const& model , crb::stage stage, int level = 0 ) :
        super ( model, stage, level ),
        M_block_initialized( false )
    {
        this->initModelSaddlePoint();
    }

    void initModelSaddlePoint()
    {

        this->M_backend_l2_vec.resize(2);
        this->M_backend_l2_vec[0]=backend( _name="backend-Xh0" );
        this->M_backend_l2_vec[1]=backend( _name="backend-Xh1" );

        auto Xh0 = this->functionSpace()->template functionSpace<0>();
        auto Xh1 = this->functionSpace()->template functionSpace<1>();

        if ( this->M_inner_product_matrix_vec.size()==0 )
            this->M_inner_product_matrix_vec.resize(2);

        auto M = this->energyMatrix();
        if ( M )
        {
            Feel::cout << "Using energy matrix as inner product matrix ! (SP)\n";
            M_inner_product_matrix_vec[0] = M->createSubMatrix( M->mapRow().dofIdToContainerId(0),
                                                                M->mapCol().dofIdToContainerId(0) );
            M_inner_product_matrix_vec[1] = M->createSubMatrix( M->mapRow().dofIdToContainerId(1),
                                                                M->mapCol().dofIdToContainerId(1) );
        }
        else
        {
            auto u = Xh0->element();
            this->M_inner_product_matrix_vec[0] = backend()->newMatrix( _test=Xh0, _trial=Xh0 );
            auto m2 = form2( _trial=Xh0, _test=Xh0, _matrix=this->M_inner_product_matrix_vec[0] );
            m2 = integrate( elements(Xh0->mesh()), inner( gradt(u), grad(u) ) );

            this->M_inner_product_matrix_vec[0]->close();

            auto p = Xh1->element();
            this->M_inner_product_matrix_vec[1] = backend()->newMatrix( _test=Xh1, _trial=Xh1 );
            form2( _test=Xh1, _trial=Xh1, _matrix=this->M_inner_product_matrix_vec[1] )
                =integrate( elements(Xh1->mesh()),
                            inner( idt(p), id(p) ) );
            this->M_inner_product_matrix_vec[1]->close();
        }
    }

    virtual bool useMonolithicRbSpace() { return false; }

    void l2solveSP( vector_ptrtype& u, vector_ptrtype const& f, int n_space )
    {
        M_backend_l2_vec[n_space]->solve( _matrix=M_inner_product_matrix_vec[n_space],
                                          _solution=u, _rhs=f );
    }

    subelement_type<0> supremizer( parameter_type const& mu, element_type const& U )
    {
        auto vec = U.template elementPtr<1>();
        auto Xh0 = this->functionSpace()->template functionSpace<0>();
        auto Xh1 = this->functionSpace()->template functionSpace<1>();
        auto us = Xh0->element();

        vector_ptrtype rhs = this->M_backend_l2_vec[0]->newVector( Xh0 );
        rhs->zero();

        sparse_matrix_ptrtype A01 = this->M_backend_l2_vec[0]->newMatrix( _test=Xh0, _trial=Xh1 );
        A01 = offlineMergeForSupremizer(mu,U);

        rhs->addVector( vec, A01 );
        this->M_backend_l2_vec[0]->solve( _matrix=this->M_inner_product_matrix_vec[0],
                                          _rhs=rhs, _solution=us );
        return us;
    }

    sparse_matrix_ptrtype offlineMergeForSupremizer( parameter_type const& mu,
                                                     element_type const& U )
    {
        sparse_matrix_ptrtype A;
        if ( this->isLinear() )
            boost::tie( boost::tuples::ignore, A, boost::tuples::ignore ) = this->update(mu);
        else
            boost::tie( boost::tuples::ignore, A, boost::tuples::ignore ) = this->update(mu,U);

        auto const& Xh0_indices = A->mapRow().dofIdToContainerId( 0 );
        auto const& Xh1_indices = A->mapRow().dofIdToContainerId( 1 );
        auto A01 = A->createSubMatrix( Xh0_indices, Xh1_indices );
        A01->close();
        return A01;
    }

    template<typename EType1, typename EType2>
    value_type AqmBlock( uint16_type q, uint16_type m,
                         EType1 const& xi_i, EType2 const& xi_j,
                         uint16_type nSpace1, uint16_type nSpace2,
                         bool transpose = false ) const
    {
        auto A = this->M_Aqm_block[nSpace1][nSpace2][q][m];
        if ( A->linftyNorm()<1e-12 )
            return 0;
        else
            return A->energy( xi_i, xi_j, transpose );
    }

    std::vector< std::vector< sparse_matrix_ptrtype >>
    AqmBlock( uint16_type n_space1, uint16_type n_space2 )
    {
        return this->M_Aqm_block[n_space1][n_space2];
    }
    std::vector< std::vector< vector_ptrtype >>
    FqmBlock( uint16_type l, uint16_type n_space )
    {
        return this->M_Fqm_block[l][n_space];
    }

    template <typename EType>
    value_type FqmBlock( uint16_type l, uint16_type q, uint16_type m,
                         EType const& xi, uint16_type nSpace )
    {
        auto F = this->M_Fqm_block[l][nSpace][q][m];
        if ( F->linftyNorm()<1e-12 )
            return 0;
        else
            return inner_product( *F, xi );
    }

    void initBlockMatrix()
    {
        if ( !M_block_initialized )
        {
            M_block_initialized=true;
            M_Aqm_block.resize( 2 );

            for ( size_type r=0; r<2; r++ )
            {
                M_Aqm_block[r].resize( 2 );
                for ( size_type c=0; c<2; c++ )
                {
                    M_Aqm_block[r][c].resize( this->Qa() );
                    for ( size_type q=0; q<this->Qa(); q++ )
                    {
                        int mMax = this->mMaxA(q);
                        M_Aqm_block[r][c][q].resize( mMax );
                        for ( size_type m=0; m<mMax; m++ )
                        {
                            auto A = this->M_Aqm[q][m];
                            auto const& i_row = A->mapRow().dofIdToContainerId( r );
                            auto const& i_col = A->mapCol().dofIdToContainerId( c );
                            M_Aqm_block[r][c][q][m] = A->createSubMatrix( i_row, i_col );
                        }
                    }
                }
            }

            int n_output = this->M_Fqm.size();
            M_Fqm_block.resize( n_output );

            for ( size_type l=0; l<n_output; l++ )
            {
                M_Fqm_block[l].resize( 2 );
                for ( size_type r=0; r<2; r++ )
                {
                    int qMax = this->Ql(l);
                    M_Fqm_block[l][r].resize( qMax );
                    for ( size_type q=0; q<qMax; q++ )
                    {
                        int mMax = this->mMaxF( l, q );
                        M_Fqm_block[l][r][q].resize( mMax );
                        for ( size_type m=0; m<mMax; m++ )
                        {
                            auto V = this->M_Fqm[l][q][m];
                            auto const& i_row = V->map().dofIdToContainerId(r);
                            M_Fqm_block[l][r][q][m] = V->createSubVector( i_row );
                        }
                    }
                }
            }
        }
    }

    void clearBlockMatix()
    {
        M_Aqm_block.clear();
        M_Fqm_block.clear();
        M_block_initialized = false;
    }

    using super::scalarProduct;
    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_type const& X, vector_type const& Y, int n_space )
    {
        return M_inner_product_matrix_vec[n_space]->energy( X, Y );
    }
    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y, int n_space )
    {
        return M_inner_product_matrix_vec[n_space]->energy( X, Y );
    }



protected :
    std::vector<sparse_matrix_ptrtype> M_inner_product_matrix_vec;
    std::vector< backend_ptrtype > M_backend_l2_vec;

    std::vector< std::vector< std::vector< std::vector<sparse_matrix_ptrtype>>>> M_Aqm_block;
    std::vector< std::vector< std::vector< std::vector<vector_ptrtype>>>> M_Fqm_block, M_Lqm_block;
    bool M_block_initialized;

}; // class CRBModelSaddlepoint


} // namespace Feel

#endif // __CRBModelSaddlePoint_H
