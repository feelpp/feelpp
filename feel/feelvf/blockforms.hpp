/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 24 Jul 2016

 Copyright (C) 2016 Feel++ Consortium

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
#ifndef FEELPP_BLOCKFORMS_HPP
#define FEELPP_BLOCKFORMS_HPP 1

#include <feel/feelvf/form.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/product.hpp>


namespace Feel {

/**
 * Handles bilinear form over a product of spaces
 */
template<typename PS>
class BlockBilinearForm
{
public :
    using product_space_t = PS;
    using value_type = double;
    BlockBilinearForm(product_space_t&& ps )
        :
        M_ps(ps),
        M_matrix(backend()->newBlockMatrix(_block=csrGraphBlocks(M_ps)))
        {}
    BlockBilinearForm(product_space_t&& ps, sparse_matrix_ptrtype&& m)
        :
        M_ps(ps),
        M_matrix(m)
        {}

    template<typename N1,typename N2>
    decltype(auto) operator()( N1 n1, N2 n2 )
        {
            cout << "filling out matrix block (" << n1 << "," << n2 << ")\n";
            return form2(_test=M_ps[n1],_trial=M_ps[n2], _matrix=M_matrix, _rowstart=int(n1), _colstart=int(n2) );
        }
    void close()
        {
            M_matrix->close();
        }
    sparse_matrix_ptrtype matrixPtr() const { return M_matrix; }
    sparse_matrix_ptrtype matrixPtr() { return M_matrix; }
    using pre_solve_type = typename Backend<value_type>::pre_solve_type;
    using post_solve_type = typename Backend<value_type>::post_solve_type;
    BOOST_PARAMETER_MEMBER_FUNCTION( ( typename Backend<double>::solve_return_type ),
                                     solve,
                                     tag,
                                     ( required
                                       ( in_out( solution ),* )
                                       ( rhs, * ) )
                                     ( optional
                                       ( name,           ( std::string ), "" )
                                       ( kind,           ( std::string ), soption(_prefix=name,_name="backend") )
                                       ( rebuild,        ( bool ), boption(_prefix=name,_name="backend.rebuild") )
                                       ( pre, (pre_solve_type), pre_solve_type() )
                                       ( post, (post_solve_type), post_solve_type() )
                                         ) )
        {
            auto U = backend()->newBlockVector(_block=solution, _copy_values=false);
            auto r = backend( _name=name, _kind=kind, _rebuild=rebuild,
                              _worldcomm=Environment::worldComm() )->solve( _matrix=this->matrixPtr(),
                                                                   _rhs=rhs.vectorPtr(),
                                                                   _solution=U,
                                                                   _pre=pre,
                                                                   _post=post
                                                                   );
            U->printMatlab("U.m");
            solution.localize(U);
            solution(0,0)->printMatlab("u.m");
            solution(1,0)->printMatlab("v.m");
            return r;
        }
    product_space_t M_ps;
    sparse_matrix_ptrtype M_matrix;
};

template<typename PS>
BlockBilinearForm<PS>
blockform2( PS&& ps )
{
    return BlockBilinearForm<PS>( std::forward<PS>(ps) );
}
template<typename PS>
BlockBilinearForm<PS>
blockform2( PS&& ps, sparse_matrix_ptrtype m )
{
    return BlockBilinearForm<PS>( std::forward<PS>(ps), std::forward<sparse_matrix_ptrtype>(m) );
}
/**
 * Handles linear form over a product of spaces
 */
template<typename PS>
class BlockLinearForm
{
public :
    using product_space_t = PS;

    BlockLinearForm() = default;
    BlockLinearForm(product_space_t&& ps )
        :
        M_ps(ps),
        M_vector(backend()->newBlockVector(_block=blockVector(M_ps), _copy_values=false))
        {}
    BlockLinearForm(product_space_t&& ps, vector_ptrtype&& v )
        :
        M_ps(ps),
        M_vector(v)
        {}
    BlockLinearForm( BlockLinearForm&& ) = default;
    BlockLinearForm( BlockLinearForm& ) = default;
    BlockLinearForm& operator=( BlockLinearForm && lf ) = default;
    BlockLinearForm& operator=( BlockLinearForm const& lf ) = default;

    template<typename N1>
    decltype(auto) operator()( N1 n1 )
        {
            cout << "filling out vector block (" << n1 << ")\n";
            return form1(_test=M_ps[n1],_vector=M_vector, _rowstart=int(n1) );
        }
    void close()
        {
            M_vector->close();
        }
    product_space_t functionSpace() const { return M_ps; }

    vector_ptrtype vectorPtr() const { return M_vector; }
    vector_ptrtype vectorPtr() { return M_vector; }

    /**
     * set linear form to 0
     */
    void zero() { M_vector->zero(); }

    BlockLinearForm& operator+=( BlockLinearForm const& l )
        {
            if ( this == &l )
            {
                M_vector->scale( 2. );
                return *this;
            }

            *M_vector += *l.M_vector;

            return *this;
        }
    product_space_t M_ps;
    vector_ptrtype M_vector;
};

template<typename PS>
BlockLinearForm<PS>
blockform1( PS&& ps )
{
    return BlockLinearForm<PS>( std::forward<PS>(ps) );
}


}
#endif
