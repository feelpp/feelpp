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

#include <feel/feelmor/crbmodelblock.hpp>

namespace Feel
{

template<typename ModelType>
class CRBModelSaddlePoint :
        public CRBModelBlock<ModelType>
{
    typedef CRBModelBlock<ModelType> super;
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

    CRBModelSaddlePoint( std::string modelName, crb::stage stage, int level = 0 ) :
        super ( modelName, stage, level ),
        M_addSupremizer(boption(_prefix=this->M_prefix,_name="crb.saddlepoint.add-supremizer"))
    {}

    CRBModelSaddlePoint( std::string modelName, model_ptrtype const& model , crb::stage stage, int level = 0 ) :
        super ( modelName, model, stage, level ),
        M_addSupremizer(boption(_prefix=this->M_prefix,_name="crb.saddlepoint.add-supremizer"))
    {}

    bool addSupremizerInSpace( int const& n_space ) const override
    {
        if ( n_space==0 )
            return M_addSupremizer;
        return false;
    }

    int QTri() { return 0; }

    element_type supremizer( parameter_type const& mu, element_type const& U, int n_space ) override
    {
        auto Us = this->functionSpace()->element();

        if ( n_space==0 )
        {
            auto vec = U.template elementPtr<1>();
            auto Xh0 = this->functionSpace()->template functionSpace<0>();
            auto Xh1 = this->functionSpace()->template functionSpace<1>();
            auto us = Us.template element<0>();
            //auto us = Xh0->element();

            vector_ptrtype rhs = this->M_backend_l2_vec[0]->newVector( Xh0 );
            rhs->zero();

            sparse_matrix_ptrtype A01 = this->M_backend_l2_vec[0]->newMatrix( _test=Xh0, _trial=Xh1 );
            A01 = offlineMergeForSupremizer(mu,U);

            rhs->addVector( vec, A01 );
            this->M_backend_l2_vec[0]->solve( _matrix=this->M_inner_product_matrix_vec[0],
                                              _rhs=rhs, _solution=us );
            //u = us.container();
        }

        return Us;
    }

    sparse_matrix_ptrtype offlineMergeForSupremizer( parameter_type const& mu,
                                                     element_type const& U )
    {
        sparse_matrix_ptrtype A;
        if ( this->isLinear() && !this->isTrilinear() )
             boost::tie( boost::tuples::ignore, A, boost::tuples::ignore ) = this->update(mu);
        else
            boost::tie( boost::tuples::ignore, A, boost::tuples::ignore ) = this->update(mu,U);

        auto const& Xh0_indices = A->mapRow().dofIdToContainerId( 0 );
        auto const& Xh1_indices = A->mapRow().dofIdToContainerId( 1 );
        auto A01 = A->createSubMatrix( Xh0_indices, Xh1_indices );
        A01->close();
        return A01;
    }
protected:
    bool M_addSupremizer;


}; // class CRBModelSaddlepoint


} // namespace Feel

#endif // __CRBModelSaddlePoint_H
