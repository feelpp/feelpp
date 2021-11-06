/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-03-20

  Copyright (C) 2013-2016 Feel++ Consortium

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
   \file evaluator.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-03-20
 */
#ifndef FEELPP_VF_EVALUATORCONTEXT_H
#define FEELPP_VF_EVALUATORCONTEXT_H

#include <feel/feelcore/parameter.hpp>
#include <feel/feelcore/commobject.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/projectors.hpp>

namespace Feel
{
struct GeometricSpaceBase;

namespace vf
{
namespace details
{
/**
 * \class EvaluatorContext
 * \brief work class to evaluate expressions at sets of points
 *
 * @author Christophe Prud'homme
 */
template<typename CTX, typename ExprT, typename CTX2=CTX>
class EvaluatorContext : public CommObject
{
public:


    /** @name Typedefs
     */
    //@{
    using super = CommObject;
    
    //we don't use it, we use this delcared in function space Context
    //static const size_type context = ExprT::context|vm::POINT ;

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;

    typedef CTX context_type;
    typedef CTX2 context2_type;
    typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> element_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    EvaluatorContext( context_type const& ctx,
                      context2_type const& ctx2,
                      expression_type const& __expr,
                      int max_points_used,
                      GeomapStrategyType geomap_strategy,
                      bool mpi_communications,
                      bool projection,
                      worldcomm_ptr_t const& worldComm )
        :
        super( worldComm ),
        M_ctx( ctx ),
        M_ctx2( ctx2 ),
        M_expr( __expr ),
        M_max_points_used( max_points_used ),
        M_geomap_strategy( geomap_strategy ),
        M_mpi_communications( mpi_communications ),
        M_projection( projection )
    {
        DVLOG(2) << "EvaluatorContext constructor from expression\n";
    }


    EvaluatorContext( EvaluatorContext const& __vfi )
        :
        super( __vfi ),
        M_ctx( __vfi.M_ctx ),
        M_ctx2( __vfi.M_ctx2 ),
        M_expr( __vfi.M_expr ),
        M_max_points_used( __vfi.M_max_points_used ),
        M_geomap_strategy( __vfi.M_geomap_strategy ),
        M_mpi_communications( __vfi.M_mpi_communications ),
        M_projection( __vfi.M_projection )
    {
        DVLOG(2) << "EvaluatorContext copy constructor\n";
    }

    ~EvaluatorContext() override {}

    //@}

    /** @name Operator overloads
     */
    //@{

    element_type operator()() const;

    /**
     * instead of evaluate an expression, evaluate the projection
     * of the expression of the function space.
     * But warning, here the projection is made only on elements
     * that contains nodes given by the context and not on all
     * elements of the mesh.
     * As we project the expression on the function space linked
     * the the context, we can't project gradient of more
     * generally shape::N > 1
     */
    element_type evaluateProjection( ) const;
    element_type evaluateProjection( mpl::bool_< true > ) const;
    element_type evaluateProjection( mpl::bool_< false > ) const;
    //@}

    /** @name Accessors
     */
    //@{

    /**
     * get the variational expression
     *
     *
     * @return the variational expression
     */
    expression_type const& expression() const
    {
        return M_expr;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    //! polynomial order
    uint16_type polynomialOrder() const { return 1; }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    //@}

private:

    context_type M_ctx;
    context2_type M_ctx2;
    expression_type const&  M_expr;
    int M_max_points_used;
    GeomapStrategyType M_geomap_strategy;
    bool M_mpi_communications;
    bool M_projection;
};

template<typename CTX, typename ExprT, typename CTX2>
typename EvaluatorContext<CTX, ExprT, CTX2>::element_type
EvaluatorContext<CTX, ExprT, CTX2>::operator()() const
{
    if( M_projection )
        return evaluateProjection();

    //rank of the current processor
    int proc_number = this->worldComm().globalRank();

    //total number of processors
    int nprocs = this->worldComm().globalSize();

    auto it = M_ctx.begin();
    auto en = M_ctx.end();

    typedef typename CTX::mapped_type::element_type::geometric_mapping_context_ptrtype gm_context_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    typedef expression_type the_expression_type;
    typedef typename boost::remove_reference<typename boost::remove_const<the_expression_type>::type >::type iso_expression_type;
    typedef typename iso_expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename t_expr_type::value_type value_type;
    typedef typename t_expr_type::shape shape;

    //in case of scalar unknown, shape::M should be 1
    //CHECK( shape::M == 1 ) << "Invalid expression shape " << shape::M << " should be 1";

    int local_max_size = 0;
    int global_max_size = 0;
    int npoints = M_ctx.nPoints();
    if( M_max_points_used > 0 )
        local_max_size = M_max_points_used;
    else
        local_max_size = npoints;
    mpi::all_reduce( this->worldComm(), local_max_size, global_max_size, mpi::maximum<int>() );

    element_type __globalv( global_max_size*shape::M*shape::N );
    __globalv.setZero();

    //local version of __v on each proc
    element_type __localv( global_max_size*shape::M*shape::N );
    __localv.setZero();

    if ( !M_ctx.empty() )
    {
        /**
         * be careful there is no guarantee that the set of contexts will
         * have the reference points. We should probably have a flag set by
         * the programmer so that we don't have to re-create the expression
         * context if the reference points are the same
         */
        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >(it->second->gmContext() ) );

        t_expr_type tensor_expr( M_expr, mapgmc );

        auto Xh = M_ctx.ptrFunctionSpace();
        auto Xh2 = M_ctx2.ptrFunctionSpace();
        //loop on local points
        for ( int p = 0; it!=en ; ++it, ++p )
        {
            auto const& ctx = *it;

            //if( CTX::is_rb_context )
            //    LOG( INFO ) << "we have a RB context ";
            //else
            //    LOG( INFO ) << "we have a FEM context";
            int global_p = it->first;

            if( global_p < global_max_size )
            {
                if ( !M_ctx2.empty() && M_ctx2.find( global_p ) != M_ctx2.end() )
                {
                    auto const& ctx2 = *M_ctx2.find( global_p );
                    tensor_expr.updateContext( Xh->contextBasis( ctx, M_ctx ), Xh2->contextBasis( ctx2, M_ctx2 ) );
                }
                else
                    tensor_expr.updateContext( Xh->contextBasis( ctx, M_ctx ) );

                //LOG( INFO ) << "Xh->contextBasis returns a context of type \n"<< typeid( decltype( Xh->contextBasis( ctx, M_ctx ) )  ).name();

                for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
                {
                    for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                    {
                        //__localv(shape::M*p+c1) = tensor_expr.evalq( c1, 0, 0 );
                        __localv(global_p*shape::M*shape::N+c1+c2*shape::M) = tensor_expr.evalq( c1, c2, 0 );
                        //LOG( INFO ) << "__localv("<<shape::M*p+c1<<") = "<<tensor_expr.evalq( c1, 0, 0 )<<" and global p = "<<global_p;
                    }
                }
            }//only if globalp < max_size

        }//loop over local points
    }

    if( ! M_mpi_communications || nprocs == 1 )
    {
        //in this case, we don't call mpi::all_reduce
        //to fill __globalv
        return __localv;
    }

    //bring back each proc contribution in __globalv
    //mpi::all_reduce( worldComm , __localv, __globalv, std::plus< element_type >() );
    mpi::all_reduce( this->worldComm(), __localv, __globalv, [](element_type x, element_type y)
                     {
                         return x + y;
                     } );

    //LOG( INFO ) << "__globalv : "<<__globalv;

    return __globalv;
}

template<typename CTX, typename ExprT, typename CTX2>
typename EvaluatorContext<CTX, ExprT, CTX2>::element_type
EvaluatorContext<CTX, ExprT, CTX2>::evaluateProjection(  ) const
{
    typedef typename CTX::mapped_type::element_type::geometric_mapping_context_ptrtype gm_context_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    typedef expression_type the_expression_type;
    typedef typename boost::remove_reference<typename boost::remove_const<the_expression_type>::type >::type iso_expression_type;
    typedef typename iso_expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename t_expr_type::value_type value_type;
    typedef typename t_expr_type::shape shape;
    static const bool shapeN = (shape::N==1);

    typedef typename boost::remove_reference<typename boost::remove_const< decltype(*M_ctx.ptrFunctionSpace()) >::type >::type ctxspace_type;
    static const bool ctxspace_is_geometricspace = boost::is_base_of<GeometricSpaceBase,ctxspace_type>::type::value;
    return evaluateProjection( mpl::bool_<shapeN && !ctxspace_is_geometricspace>() );
}
template<typename CTX, typename ExprT, typename CTX2>
typename EvaluatorContext<CTX, ExprT, CTX2>::element_type
EvaluatorContext<CTX, ExprT, CTX2>::evaluateProjection( mpl::bool_<true> ) const
{
    auto Xh = M_ctx.ptrFunctionSpace();
    auto const& worldComm = Xh->worldComm();
    int proc_number = worldComm.globalRank();
    int nprocs = worldComm.globalSize();

    auto it = M_ctx.begin();
    auto en = M_ctx.end();

    typedef typename CTX::mapped_type::element_type::geometric_mapping_context_ptrtype gm_context_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    typedef expression_type the_expression_type;
    typedef typename boost::remove_reference<typename boost::remove_const<the_expression_type>::type >::type iso_expression_type;
    typedef typename iso_expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename t_expr_type::value_type value_type;
    typedef typename t_expr_type::shape shape;

    int max_size = 0;
    int npoints = M_ctx.nPoints();
    if( M_max_points_used > 0 )
        max_size = M_max_points_used;
    else
        max_size = npoints;

    element_type __globalv( max_size*shape::M );
    __globalv.setZero();

    //local version of __v on each proc
    element_type __localv( max_size*shape::M );
    __localv.setZero();

    if ( !M_ctx.empty() )
    {
        /**
         * be careful there is no guarantee that the set of contexts will
         * have the reference points. We should probably have a flag set by
         * the programmer so that we don't have to re-create the expression
         * context if the reference points are the same
         */
        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >(it->second->gmContext() ) );

        t_expr_type tensor_expr( M_expr, mapgmc );

        //loop on local points
        for ( int p = 0; it!=en ; ++it, ++p )
        {
            auto const& ctx = *it;

            int global_p = it->first;

            if( global_p < max_size )
            {
                auto const& e = ctx.second->gmContext()->element();
                //Xh is a pointer, not a shared ptr
                //functionspace is a shared ptr
                auto functionspace = M_ctx.functionSpace();
                auto projected_expression = vf::project( _space=functionspace, _expr=M_expr , _range=idedelements( Xh->mesh(), e.id() ) );
                auto myctx=functionspace->context();
                myctx.addCtx(  it->second , proc_number );
                bool do_communications=false;//we don't want that each proc have the result now ( but latter )
                auto val = projected_expression.evaluate( myctx , do_communications );
                //loop over components
                for(int comp=0; comp<shape::M; comp++)
                {
                    __localv( global_p*shape::M+comp ) = val( comp );
                }
            }//only if globalp < max_size

        }//loop over local points
    }


    if( ! M_mpi_communications || nprocs == 1 )
    {
        return __localv;
    }
    mpi::all_reduce( worldComm , __localv, __globalv, std::plus< element_type >() );
    return __globalv;

}
template<typename CTX, typename ExprT, typename CTX2>
typename EvaluatorContext<CTX, ExprT, CTX2>::element_type
EvaluatorContext<CTX, ExprT, CTX2>::evaluateProjection( mpl::bool_<false> ) const
{
    //here, the expression is a gradient
    //so wa can't project it on the function space
    bool go = false;
    CHECK( go ) << "evaluateFromContext( _expr=..., _projection=true ) can't be used if _expr is a gradient ! Or more generally if shape::N > 1 ! \n";
    element_type __globalv;
    return __globalv;
}

}
/// \endcond

#if 0
/// \cond DETAIL
namespace detail
{
template<typename Args>
struct evaluate_context
{
    typedef clean_type<Args,tag::expr> _expr_type;
    typedef clean_type<Args,tag::context> _context_type;
    typedef clean2_type<Args, tag::context2,_context_type> _context2_type;
    typedef details::EvaluatorContext<_context_type, Expr<_expr_type>, _context2_type > eval_t;
    typedef typename eval_t::element_type element_type;
};
}
/// \endcond
#endif


template<typename Ctx, typename ExprT, typename Ctx2>
typename details::EvaluatorContext<Ctx, Expr<ExprT>, Ctx2 >::element_type
evaluatecontext_impl( Ctx const& ctx,
                      Ctx2 const& ctx2,
                      Expr<ExprT> const& __expr,
                      int max_points_used = -1,
                      GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_HO,
                      bool mpi_communications = true,
                      bool projection = false,
                      worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
{
    typedef details::EvaluatorContext<Ctx, Expr<ExprT>, Ctx2 > proj_t;
    proj_t p( ctx, ctx2, __expr, max_points_used, geomap , mpi_communications , projection, worldComm );
    return p();
}


/**
 *
 * \brief evaluate an expression over a range of element at a set of points defined in the reference element
 *
 * \arg range the range of mesh elements to apply the projection (the remaining parts are set to 0)
 * \arg pset set of points (e.g. quadrature points) in the reference elements to be transformed in the real elements
 * \arg expr the expression to project
 * \arg geomap the type of geomap to use (make sense only using high order meshes)
 * \arg mpi_communications a bool that indicates if all proc communicate or not
 * \arg projection a bool that indicates if we project the expression on function space or not (usefull for EIM)
 */

template <typename ... Ts>
auto evaluateFromContext( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && context = args.get(_context);
    auto && expr = args.get(_expr);
    auto && context2 = args.get_else(_context2,context.functionSpace()->context());
    int max_points_used = args.get_else(_max_points_used,-1);
    GeomapStrategyType geomap = args.get_else(_geomap,GeomapStrategyType::GEOMAP_OPT );
    bool mpi_communications = args.get_else(_mpi_communications,true);
    bool projection = args.get_else(_projection,false);
    worldcomm_ptr_t worldcomm = args.get_else_invocable(_worldcomm,[&mpi_communications,&context](){
         return (mpi_communications && !context.ctxHaveBeenMpiBroadcasted() )? context.functionSpace()->worldCommPtr() : context.functionSpace()->worldComm().subWorldCommSeqPtr(); } );

    bool doMpiComm = mpi_communications && !context.ctxHaveBeenMpiBroadcasted();
    //LOG(INFO) << "evaluate expression..." << std::endl;
    return evaluatecontext_impl( context, context2, expr, max_points_used, geomap , doMpiComm/*mpi_communications*/, projection, worldcomm );
    //LOG(INFO) << "evaluate expression done." << std::endl;
}


} // vf
} // feel


#endif /* __FEELPP_EVALUATORS_H */


