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
#include <feel/feelcrb/crbmodelbase.hpp>
#include <feel/feelalg/vectorblock.hpp>

#define blockA(X,Y) this->M_A.template get<X,Y>()
#define blockF(X,Y) this->M_F[X].template get<Y>()

namespace Feel
{
class FunctionSpaceDefinitionSaddlePoint
{
public :
    /*mesh*/
    typedef Simplex<1> entity_type ;
    typedef Mesh<entity_type > mesh_type ;

    /*basis*/
    typedef Lagrange< 2, Vectorial, Continuous >  basis_u_type ;
    typedef Lagrange< 1, Scalar, Continuous > basis_p_type ;
    typedef bases< basis_u_type, basis_p_type > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type , basis_type > space_type ;
    typedef typename space_type::element_type element_type;

    static const bool is_time_dependent=false;
    static const bool is_linear=true;
};


template < typename ParameterDefinition=ParameterSpace<1>,
           typename FunctionSpaceDefinition=FunctionSpaceDefinitionSaddlePoint,
           int _Options=0,
           typename EimDefinition=EimDefinitionBase<ParameterDefinition,FunctionSpaceDefinition>
           >
class CRBModelSaddlePoint :
        public CRBModelBase< ParameterDefinition, FunctionSpaceDefinition, _Options, EimDefinition >
{
    typedef CRBModelBase< ParameterDefinition, FunctionSpaceDefinition, _Options, EimDefinition > super_type;
    typedef CRBModelSaddlePoint< ParameterDefinition, FunctionSpaceDefinition, _Options, EimDefinition > self_type;

public :
    typedef typename super_type::value_type value_type;

    //! space_type
    typedef typename super_type::space_type space_type;
    typedef typename super_type::functionspace_ptrtype functionspace_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::element_ptrtype element_ptrtype;
    template <int T>
    using subspace_type = typename space_type::template sub_functionspace<T>::type::element_type;
    template <int T>
    using subspace_ptrtype = boost::shared_ptr<subspace_type<T>>;

    template <int T>
    using rbspace_type = ReducedBasisSpace<subspace_type<T>>;
    template <int T>
    using rbspace_ptrtype = boost::shared_ptr<rbspace_type<T>>;

    typedef typename super_type::backend_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename super_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super_type::vector_ptrtype vector_ptrtype;
    typedef typename super_type::vector_type vector_type;

    typedef typename super_type::parameterspace_type parameterspace_type;
    typedef typename super_type::parameter_type parameter_type;

    typedef Eigen::VectorXd vectorN_type;
    typedef std::vector< std::vector< double > > beta_vector_type;

    using super_type::rBFunctionSpace;

    CRBModelSaddlePoint() :
        super_type()
        {}

    virtual element_type solve( parameter_type const& mu )
        {
            bool transpose = boption("crb.saddlepoint.transpose");

            mergeBlock( A00, blockA(0,0), mu );
            mergeBlock( A01, blockA(0,1), mu );
            mergeBlock( A10, blockA(1,0), mu );
            mergeBlock( A11, blockA(1,1), mu );

            mergeBlock( F0, blockF(0,0), mu );
            mergeBlock( F1, blockF(0,1), mu );

            if ( transpose && blockA(0,1) && !blockA(1,0) )
                A01->transpose( A10 );
            if ( transpose && !blockA(0,1) && blockA(1,0) )
                A10->transpose( A01 );

            auto u=Xh0->elementPtr();
            auto p=Xh1->elementPtr();

            BlocksBaseSparseMatrix<value_type> block_lhs(2,2);
            block_lhs(0,0)=A00;
            block_lhs(0,1)=A01;
            block_lhs(1,0)=A10;
            block_lhs(1,1)=A11;
            auto blockLHS = backend()->newBlockMatrix( _block=block_lhs, _copy_values=true );

            BlocksBaseVector<value_type> block_rhs(2);
            block_rhs(0)=F0;
            block_rhs(1)=F1;
            auto blockRHS = backend()->newBlockVector( _block=block_rhs, _copy_values=true );

            BlocksBaseVector<value_type> block_sol(2);
            block_sol(0)=u;
            block_sol(1)=p;
            auto blockSOL = backend()->newBlockVector( _block=block_sol, _copy_values=false );

            backend( _rebuild=true )->solve( _matrix=blockLHS, _rhs=blockRHS, _solution=blockSOL );

            M_U=*blockSOL;
            return M_U;
        }


    virtual void initDerived()
        {
            M_U = this->Xh->element();

            A00 = backend()->newMatrix(_test=Xh0, _trial=Xh0 );
            A01 = backend()->newMatrix(_test=Xh0, _trial=Xh1, _buildGraphWithTranspose=true );
            A10 = backend()->newMatrix(_test=Xh1, _trial=Xh0, _buildGraphWithTranspose=true );
            A11 = backend()->newMatrix(_test=Xh1, _trial=Xh1 );

            F0 = backend()->newVector( Xh0 );
            F1 = backend()->newVector( Xh1 );
        }

    virtual void setFunctionSpaces( functionspace_ptrtype Vh )
        {
            this->Xh=Vh;
            Xh0 = this->Xh->template functionSpace<0>();
            Xh1 = this->Xh->template functionSpace<1>();

            XN0 = rbspace_type<0>::New( _space=Xh0 );
            XN1 = rbspace_type<1>::New( _space=Xh1 );
            if (Environment::isMasterRank() )
                std::cout << "Number of dof : " << this->Xh->nDof() << std::endl
                          << "Number of local dof : " << this->Xh->nLocalDof() << std::endl
                          << "Number of dof Xh0,Xh1 = "<< Xh0->nDof() << "," << Xh1->nDof() << std::endl;
        }

    template<int T>
    rbspace_ptrtype<T> rBFunctionSpace()
        {
            fusion::vector<rbspace_ptrtype<0>,rbspace_ptrtype<1> > v( XN0, XN1 );
            return fusion::at_c<T>( v );
        }

protected :
    template<typename BlockType>
    void mergeBlock( sparse_matrix_ptrtype mat, BlockType b, parameter_type const& mu )
        {
            mat->zero();
            if ( b )
            {
                auto beta = b->computeBetaQm( mu );
                auto Aqm = b->compute();
                for( int q=0; q<b->Q(); q++ )
                    for( int m=0; m<b->mMax(q); m++ )
                        mat->addMatrix( beta[q][m], Aqm[q][m] );
            }
        }
    template<typename BlockType>
    void mergeBlock( vector_ptrtype vec, BlockType b, parameter_type const& mu )
        {
            vec->zero();
            if ( b )
            {
                auto beta = b->computeBetaQm( mu );
                auto Aqm = b->compute();
                for( int q=0; q<b->Q(); q++ )
                    for( int m=0; m<b->mMax(q); m++ )
                        vec->add( beta[q][m], Aqm[q][m] );
            }
        }

    sparse_matrix_ptrtype A00, A01, A10, A11;
    element_type M_U;
    vector_ptrtype F0, F1;

    subspace_ptrtype<0> Xh0;
    subspace_ptrtype<1> Xh1;

    rbspace_ptrtype<0> XN0;
    rbspace_ptrtype<1> XN1;
}; // class CRBModelSaddlepoint

} // namespace Feel

#endif // __CRBModelSaddlePoint_H
