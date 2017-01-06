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
    typedef boost::shared_ptr<ModelType> model_ptrtype;

    typedef typename model_type::value_type value_type;
    typedef typename model_type::parameter_type parameter_type;

    typedef typename ModelType::mesh_type mesh_type;
    typedef typename ModelType::mesh_ptrtype mesh_ptrtype;

    typedef typename model_type::space_type space_type;
    template <int T>
    using subspace_type = typename space_type::template sub_functionspace<T>::type;
    template <int T>
    using subspace_ptrtype = boost::shared_ptr<subspace_type<T>>;
    template<int T>
    using subelement_type = typename subspace_type<T>::element_type;

    CRBModelSaddlePoint( CRBModelMode mode = CRBModelMode::PFEM, int level=0, bool doInit = true ) :
        super ( mode, level, doInit )
    {
        this->initModelSaddlePoint();
    }

    CRBModelSaddlePoint( model_ptrtype const& model , CRBModelMode mode = CRBModelMode::PFEM, bool doInit = true ) :
        super ( model, mode, doInit )
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

        M_Xh0_indices.resize( Xh0->nLocalDofWithGhost() );
        M_Xh1_indices.resize( Xh1->nLocalDofWithGhost() );
        std::iota( M_Xh0_indices.begin(), M_Xh0_indices.end(), 0 );
        std::iota( M_Xh1_indices.begin(), M_Xh1_indices.end(), Xh0->nLocalDofWithGhost() );

        this->M_inner_product_matrix_vec.resize(2);
        auto u = Xh0->element();
        this->M_inner_product_matrix_vec[0] = backend()->newMatrix( _test=Xh0, _trial=Xh0 );
        form2( _trial=Xh0, _test=Xh0, _matrix=this->M_inner_product_matrix_vec[0] )
            = integrate( elements(Xh0->mesh()),
                         inner( idt(u), id(u) )
                         +inner( gradt(u), grad(u) ) );
        this->M_inner_product_matrix_vec[0]->close();

        auto p = Xh1->element();
        this->M_inner_product_matrix_vec[1] = backend()->newMatrix( _test=Xh1, _trial=Xh1 );
        form2( _test=Xh1, _trial=Xh1, _matrix=this->M_inner_product_matrix_vec[1] )
            =integrate( elements(Xh1->mesh()),
                        inner( idt(p), id(p) ) );
        this->M_inner_product_matrix_vec[1]->close();
    }

    subelement_type<0> supremizer( parameter_type const& mu, vector_ptrtype const& vec )
    {
        auto Xh0 = this->functionSpace()->template functionSpace<0>();
        auto us = Xh0->element();

        vector_ptrtype rhs = this->M_backend_l2_vec[0]->newVector( Xh0 );
        rhs->zero();

        sparse_matrix_ptrtype A01 = offlineMergeForSupremizer(mu);
        rhs->addVector( vec, A01 );
        this->M_backend_l2_vec[0]->solve( _matrix=this->M_inner_product_matrix_vec[0],
                                          _rhs=rhs, _solution=us );
        return us;
    }

    sparse_matrix_ptrtype offlineMergeForSupremizer( parameter_type const& mu )
    {
        sparse_matrix_ptrtype A;
        boost::tie( boost::tuples::ignore, A, boost::tuples::ignore ) = this->update(mu);

        return A->createSubMatrix( M_Xh0_indices, M_Xh1_indices );
    }


protected :
    std::vector<sparse_matrix_ptrtype> M_inner_product_matrix_vec;
    std::vector< backend_ptrtype > M_backend_l2_vec;
    std::vector<size_type> M_Xh0_indices;
    std::vector<size_type> M_Xh1_indices;

}; // class CRBModelSaddlepoint


} // namespace Feel

#endif // __CRBModelSaddlePoint_H
