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
   \file preassembleobject.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2012-04-23
 */
#ifndef __PreAssembleObject_H
#define __PreAssembleObject_H 1

#include <feel/feel.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel
{

/**
 * \class PreAssembleObject
 * \brief Pre Assemble Object class
 *
 * @author Christophe Prud'homme
 * @author Stephane Veys
 * @see
 */

class PreAssembleObjectBase
{
public :

    PreAssembleObjectBase()
    :
        M_backend ( backend_type::build( BACKEND_PETSC ) )
    {}

    typedef double value_type;

    /* backend */
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    /* vectors */
    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    backend_ptrtype backend()
    {
        return M_backend;
    }

    //assemble a matrix
    virtual sparse_matrix_ptrtype assembleMatrix()
    {
        return M_matrix;
    }

    //assemble a vector
    virtual vector_ptrtype assembleVector()
    {
        return M_vector;
    }

private :
    vector_ptrtype M_vector;
    sparse_matrix_ptrtype M_matrix;
    backend_ptrtype M_backend;
};//PreAssembleObjectBase

template< typename RangeType , typename ExprType , typename TestFunctionSpaceType , typename TrialFunctionSpaceType >
class PreAssembleMatrixObject
    : public PreAssembleObjectBase
{

public :
    typedef PreAssembleObjectBase super;
    typedef typename super::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef RangeType range_type;
    typedef boost::shared_ptr<range_type> range_ptrtype;

    typedef ExprType expr_type;
    typedef boost::shared_ptr<expr_type> expr_ptrtype;

    typedef TestFunctionSpaceType test_functionspace_type;
    typedef boost::shared_ptr<test_functionspace_type> test_functionspace_ptrtype;

    typedef TrialFunctionSpaceType trial_functionspace_type;
    typedef boost::shared_ptr<trial_functionspace_type> trial_functionspace_ptrtype;

    typedef PreAssembleMatrixObject< RangeType , ExprType , TestFunctionSpaceType , TrialFunctionSpaceType > preassemble_type;
    typedef boost::shared_ptr<preassemble_type> preassemble_ptrtype;

    //constructor
    PreAssembleMatrixObject( range_type range ,
                             expr_type expr   ,
                             test_functionspace_type test ,
                             trial_functionspace_type trial )
        :
        super(),
        M_range( range ),
        M_expr( expr ),
        M_test( test ),
        M_trial( trial )
    {}

    range_type range()
    {
        return M_range;
    }

    expr_type expr()
    {
        return M_expr;
    }

    test_functionspace_type test()
    {
        return M_test;
    }

    test_functionspace_type trial()
    {
        return M_trial;
    }

    virtual sparse_matrix_ptrtype assembleMatrix()
    {
        auto backend = super::backend();
        auto matrix = backend->newMatrix( _test=M_test , _trial=M_trial );
        form2( _test=M_test, _trial=M_trial , _matrix=matrix ) = integrate( _range=M_range, _expr=M_expr );
        return matrix;
    }

    static preassemble_ptrtype New ( range_type const& range , expr_type const& expr , test_functionspace_type const& test, trial_functionspace_type const& trial )
    {
        return preassemble_ptrtype( new PreAssembleMatrixObject( range , expr , test , trial ) );
    }


private :

    range_type M_range;
    expr_type M_expr;
    test_functionspace_type M_test;
    trial_functionspace_type M_trial;
};//PreAssembleMatrixObject class


template< typename RangeType , typename ExprType , typename TestFunctionSpaceType >
class PreAssembleVectorObject
    : public PreAssembleObjectBase
{

public :
    typedef PreAssembleObjectBase super;
    typedef typename super::vector_ptrtype vector_ptrtype;

    typedef RangeType range_type;
    typedef boost::shared_ptr<range_type> range_ptrtype;

    typedef ExprType expr_type;
    typedef boost::shared_ptr<expr_type> expr_ptrtype;

    typedef TestFunctionSpaceType test_functionspace_type;
    typedef boost::shared_ptr<test_functionspace_type> test_functionspace_ptrtype;

    typedef PreAssembleVectorObject< RangeType , ExprType , TestFunctionSpaceType > preassemble_type;
    typedef boost::shared_ptr<preassemble_type> preassemble_ptrtype;

    //constructor
    PreAssembleVectorObject()
    {}

    PreAssembleVectorObject( range_type range ,
                             expr_type expr   ,
                             test_functionspace_type test)
        :
        M_range( range ),
        M_expr( expr ),
        M_test( test )
    {}

    range_type range()
    {
        return M_range;
    }

    expr_type expr()
    {
        return M_expr;
    }

    test_functionspace_type test()
    {
        return M_test;
    }

    virtual vector_ptrtype assembleVector()
    {
        auto backend = super::backend();
        auto vector = backend->newVector( M_test );
        form1( _test=M_test,_vector=vector ) = integrate( _range=M_range, _expr=M_expr );
        return vector;
    }


    static preassemble_ptrtype New ( range_type const& range , expr_type const& expr , test_functionspace_type const& test )
    {
        return preassemble_ptrtype( new PreAssembleVectorObject( range , expr , test ) );
    }


private :

    range_type M_range;
    expr_type M_expr;
    test_functionspace_type M_test;
};//PreAssembleVectorObject class

template< typename RangeType, typename ExprType, typename TestFunctionSpaceType, typename TrialFunctionSpaceType >
inline
boost::shared_ptr< PreAssembleMatrixObject<RangeType, ExprType, TestFunctionSpaceType, TrialFunctionSpaceType >  >
PreAssembleMatrix( RangeType const& range ,
                   ExprType const& expr ,
                   TestFunctionSpaceType const& test,
                   TrialFunctionSpaceType const& trial )
{
    return PreAssembleMatrixObject< RangeType , ExprType , TestFunctionSpaceType , TrialFunctionSpaceType >::New( range , expr , test , trial );
}

template< typename RangeType, typename ExprType, typename TestFunctionSpaceType >
inline
boost::shared_ptr< PreAssembleVectorObject<RangeType, ExprType, TestFunctionSpaceType >  >
PreAssembleVector( RangeType const& range ,
                   ExprType const& expr ,
                   TestFunctionSpaceType const& test )
{
    return PreAssembleVectorObject< RangeType , ExprType , TestFunctionSpaceType >::New( range , expr , test );
}

}//namespace Feel
#endif


