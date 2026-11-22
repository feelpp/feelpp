/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2008 Universite Joseph Fourier (Grenoble I)
 * SPDX-FileCopyrightText: 2008-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file integratordirac.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-10
 */
#ifndef __INTEGRATORDIRAC_HPP
#define __INTEGRATORDIRAC_HPP 1

#include <feel/feelalg/enums.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/ranges.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class IntegratorDirac
 * \brief Dirac Delta function integrator
 *
 *
 *
 * @author Christophe Prud'homme
 */
template<typename ElementRange, typename Pts, typename DiracExpr >
class IntegratorDirac
{
public:


    /** @name Typedefs
     */
    //@{
    static const size_type context = DiracExpr::context | vm::POINT;

    template<typename Func>
    struct HasTestFunction
    {
        static inline const bool result = DiracExpr::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = DiracExpr::template HasTrialFunction<Func>::result;
    };

    static const size_type iDim = ElementRange::iDim();

    typedef IntegratorDirac<ElementRange, Pts, DiracExpr> self_type;

    using mesh_type = typename ElementRange::mesh_t;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using element_iterator = typename ElementRange::iterator_t;

    typedef Pts points_type;
    typedef DiracExpr expression_type;
    typedef typename expression_type::value_type value_type;
    using test_basis = typename expression_type::test_basis;
    using trial_basis = typename expression_type::trial_basis;

    struct eval
    {
        //
        // some typedefs
        //
        using element_type = typename mesh_type::element_type;
        using gm_type = typename mesh_type::gm_type;
        typedef std::shared_ptr<gm_type> gm_ptrtype;
        typedef typename gm_type::template Context<element_type> gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gm_type::precompute_ptrtype gmcpc_ptrtype;
        //typedef typename eval_expr_type::value_type value_type;
        //typedef typename strongest_numeric_type<typename Pts::value_type, typename expression_type::value_type>::type value_type;
        typedef typename expression_type::value_type value_type;

        //
        // Precompute some data in the reference element for
        // geometric mapping and reference finite element
        //
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
        typedef typename eval_expr_type::shape shape;
        //typedef typename shape_type::storage<value_type> storage_type;
        /*
          typedef mpl::if_<mpl::bool_<shape_type::is_scalar>,
          mpl::identity<value_type>,
          mpl::identity<ublas::vector<value_type,storage_type> > >::type::type ret_type;*/
        typedef Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> ret_type;
    };

    using evaluate_type = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>;
    using ret_type = evaluate_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    IntegratorDirac( ElementRange const& __elts,
                     Pts const& pts,
                     expression_type const& __expr )
        :
        M_eltbegin( __elts.begin() ),
        M_eltend( __elts.end() ),
        M_mesh( __elts.mesh() ),
        M_pts( pts ),
        M_expr( __expr )
    {
    }
    IntegratorDirac( IntegratorDirac const& ioe )
        :
        M_eltbegin( ioe.M_eltbegin ),
        M_eltend( ioe.M_eltend ),
        M_mesh( ioe.M_mesh ),
        M_pts( ioe.M_pts ),
        M_expr( ioe.M_expr )
    {
    }

    ~IntegratorDirac() {}

    //@}

    /** @name Accessors
    */
    //@{


    /**
     * iterator that points at the beginning of the container that
     * holds the data that will apply the Dirichlet condition upon
     */
    element_iterator beginElement() const
    {
        return M_eltbegin;
    }

    /**
     * iterator that points at the end of the container that
     * holds the data that will apply the Dirichlet condition upon
     */
    element_iterator endElement() const
    {
        return M_eltend;
    }


    //@}
    /** @name  Methods
     */
    //@{

    //! polynomial order
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    expression_type const& expression() const { return M_expr; }

    /**
     * assembly routine for Dirichlet condition
     *
     */
    template<typename Elem1, typename Elem2, typename FormType>
    void assemble( std::shared_ptr<Elem1> const& __u,
                   std::shared_ptr<Elem2> const& __v,
                   FormType& __f ) const
    {
        typedef typename Elem1::functionspace_type functionspace_type;
        DVLOG(2) << "[IntegratorDirac::assemble()] is_same: "
                      << mpl::bool_<boost::is_same<functionspace_type,Elem1>::value>::value << "\n";
        assemble( __u, __v, __f, mpl::bool_<boost::is_same<functionspace_type,Elem1>::value>() );
    }

    //typename expression_type::template tensor<Geo_t>::value_type
    ret_type
    evaluate( bool /*parallel*/ = true ) const
    {
        return evaluate( mpl::int_<iDim>() );
    }

    //@}
private:

    template<typename Elem1, typename Elem2, typename FormType>
    void assemble( std::shared_ptr<Elem1> const& /*__u*/,
                   std::shared_ptr<Elem2> const& /*__v*/,
                   FormType& /*__f*/, mpl::bool_<false> ) const {}

    template<typename Elem1, typename Elem2, typename FormType>
    void assemble( std::shared_ptr<Elem1> const& __u,
                   std::shared_ptr<Elem2> const& __v,
                   FormType& __f, mpl::bool_<true> ) const;

    ret_type evaluate( mpl::int_<MESH_ELEMENTS> ) const;
    ret_type evaluate( mpl::int_<MESH_FACES> ) const;

private:

    element_iterator M_eltbegin;
    element_iterator M_eltend;
    mesh_type const* M_mesh;

    Pts const& M_pts;
    expression_type M_expr;
};

template<typename ElementRange, typename Pts, typename DiracExpr>
template<typename Elem1, typename Elem2, typename FormType>
void
IntegratorDirac<ElementRange, Pts,  DiracExpr>::assemble( std::shared_ptr<Elem1> const& /*__u*/,
        std::shared_ptr<Elem2> const& /*__v*/,
        FormType& __form,
        mpl::bool_<true> ) const
{
#if 0
    std::map<std::string,std::pair<boost::timer,value_type> > timer;
    timer["init intvrho"].first.restart();

    element_type v( M_Space, "v" );
    std::fill( v.begin(), v.end(), .0 );

    molecule_type::atoms_const_iterator_type atom( molecule.begin() );

    if ( atom == molecule.end() ) return;

    mesh_type::Inverse meshinv( M_mesh );
    std::vector<value_type> atomcharges( molecule.size() );

    /* initialisation of the mesh::inverse data structure */
    for ( size_type atomid = 0; atom != molecule.end(); ++atom, ++atomid )
    {
        meshinv.addPointWithId( element_prod( M_invStretch , atom->center() - M_translation ),
                                atomid );
        atomcharges[atomid] = atom->charge();
    }

    meshinv.distribute();

    std::vector<bool> dof_done( molecule.size() );
    std::fill( dof_done.begin(), dof_done.end(), false );
    std::vector<size_type> itab;
    matrix_node<value_type>::type pts( mesh_type::nDim, 1 );
    typedef mesh_type::element_type geoelement_type;
    typedef mesh_type::element_iterator mesh_element_iterator;

    mesh_element_iterator it = M_mesh->beginElementWithProcessId( M_mesh->comm().rank() );
    mesh_element_iterator en = M_mesh->endElementWithProcessId( M_mesh->comm().rank() );

    // geometric mapping context
    using gm_type = typename eval::gm_type;
    using gmc_ptrtype = typename eval::gmc_ptrtype;

    // basis
    typedef space_type::basis_type basis_type;
    basis_type const* __basis = M_Space->basis().get();
    gm_ptrtype __gm = M_Space->gm();
    typedef gm_type::precompute_ptrtype geopc_ptrtype;
    typedef gm_type::precompute_type geopc_type;
    geopc_ptrtype __geopc( new geopc_type( __gm, __basis->dual().points() ) );
    gmc_ptrtype __c( new gmc_type( __gm, *it, __geopc ) );


    /* note: meshmover.hpp has this two lines: */
    if ( !v.areGlobalValuesUpdated() )
        v.updateGlobalValues();

    /*shall we use this ?
    */

    timer["init intvrho"].second = timer["init intvrho"].first.elapsed();
    LOG(INFO) << "[timer] init intvrho(): " << timer["init intvrho"].second << "\n";

    timer["intvrho"].first.restart();

    // size_type first_dof = M_Space->dof()->firstDof();
    for ( ; it != en; ++ it )
    {
        __c->update( *it, __geopc );
        meshinv.pointsInConvex( it->id(), itab );

        if ( itab.size() == 0 )
            continue;

        for ( size_type i = 0; i < itab.size(); ++i )
        {
            // get dof id in target dof table
            size_type dof = itab[i];

            if ( !dof_done[dof] )
            {
                dof_done[dof]=true;
                ublas::column( pts, 0 ) = meshinv.referenceCoords()[dof];
                element_type::pc_type pc( v.functionSpace()->fe(), pts );
                __geopc->update( pts );
                __c->update( *it, __geopc );

                for ( uint16_type loc_ind = 0;
                      loc_ind < __basis->localDof(); ++loc_ind )
                {
                    for ( uint16_type comp = 0; comp < basis_type::nComponents; ++comp )
                    {
                        size_type globaldof = boost::get<0>( M_Space->dof()->localToGlobal( it->id(), loc_ind, comp ) );

                        // update only values on the processor
                        if ( globaldof >= v.firstLocalIndex() &&
                                globaldof < v.lastLocalIndex() )
                        {
                            v.set( globaldof, 1 );

                            element_type::id_type interpfunc( v.id( *__c, pc ) );
                            //std::cout << "interpfunc :  " << interpfunc << "\n";

                            rhs->add( globaldof, atomcharges[dof] * interpfunc( comp, 0, 0 ) );
                            // DVLOG(2) << "rhs( " << globaldof << ")=" << (*rhs)( globaldof )
                            //           << " (just added " << atomcharges[dof] * interpfunc( comp, 0, 0 ) << " )" << "\n";
                            v.set( globaldof, 0 );

                        }
                    }
                }

            }
        }
    } // element

    timer["intvrho"].second = timer["intvrho"].first.elapsed();
    LOG(INFO) << "[timer] intvrho(): " << timer["intvrho"].second << "\n";
#endif
}
template<typename Elements, typename Pts, typename DiracExpr>
typename IntegratorDirac<Elements, Pts, DiracExpr>::ret_type
IntegratorDirac<Elements, Pts, DiracExpr>::evaluate( mpl::int_<MESH_ELEMENTS> ) const
{
    typename eval::ret_type res( eval::shape::M, eval::shape::N );
    res.setZero();
    if ( M_eltbegin == M_eltend || M_pts.size() == 0 )
        return res;

    CHECK( M_mesh ) << "Dirac integration requires an owning mesh range";
    auto mesh = std::const_pointer_cast<mesh_type>( M_mesh->shared_from_this() );
    using mesh_inverse_type = MeshInverse<mesh_type>;
    using mesh_inverse_size_type = typename mesh_inverse_type::size_type;
    mesh_inverse_type meshinv( mesh );

    auto ptit = M_pts.begin();
    auto pten = M_pts.end();

    /* initialisation of the mesh::inverse data structure */
    for ( size_type ptid = 0; ptit != pten; ++ptit, ++ptid )
    {
        meshinv.addPointWithId( ptit->node(), ptid, 0 );
    }

    meshinv.distribute();

    std::vector<bool> dof_done( M_pts.size() );
    std::fill( dof_done.begin(), dof_done.end(), false );
    std::vector<boost::tuple<mesh_inverse_size_type, uint16_type>> itab;
    typename matrix_node<value_type>::type pts( mesh_type::nDim, 1 );
    element_iterator it = M_eltbegin;
    element_iterator en = M_eltend;

    // geometric mapping context
    using gm_type = typename eval::gm_type;
    using gmc_ptrtype = typename eval::gmc_ptrtype;

    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    auto const& initialElement = boost::unwrap_ref( *it );
    auto gm = initialElement.gm();
    CHECK( gm )
        << "Dirac integration requires a geometric mapping";
    auto __geopc = gm->preCompute( pts );

    // make sure that we have elements to iterate over (return 0
    // otherwise)
    if ( it == en )
        return typename eval::ret_type( eval::shape::M, eval::shape::N );;

    gmc_ptrtype __c = gm->template context<context>(
        initialElement, __geopc, expression().dynamicContext() );
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

    typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
    eval_expr_type expr( expression(), mapgmc );
    typedef typename eval_expr_type::shape shape;

    for ( ; it != en; ++ it )
    {
        auto const& element = boost::unwrap_ref( *it );
        meshinv.pointsInConvex( element.id(), itab );

        if ( itab.size() == 0 )
            continue;

        for ( size_type i = 0; i < itab.size(); ++i )
        {
            // get dof id in target dof table
            size_type const pointId = boost::get<0>( itab[i] );

            if ( !dof_done[pointId] )
            {
                dof_done[pointId] = true;
                auto const referencePoint =
                    meshinv.referenceCoords().find( pointId );
                CHECK( referencePoint != meshinv.referenceCoords().end() )
                    << "missing Dirac reference coordinate";
                ublas::column( pts, 0 ) = referencePoint->second;
                __geopc->update( pts );
                __c->template update<context>( element );
                expr.update( mapgmc );

                for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                {
                    for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                        res( c1, c2 ) += expr.evalq( c1, c2, 0 );
                }

            }
        }
    } // element

    return res;
}


/// \endcond

/**
 * integrate the Dirac delta function
 *
 * \param elem_range tuple of iterators over geometric entities to which the points belong
 * \param pt_range pair of iterators over points that define the delta functions
 * \param expr expression to integrate
 *
 */
template<typename ElementRange, typename PointRange, typename DiracExpr>
Expr<IntegratorDirac<ElementRange, PointRange, DiracExpr> >
dirac( ElementRange const& elem_range, PointRange const& pt_range, DiracExpr const& expr )
{
    typedef IntegratorDirac<ElementRange, PointRange, DiracExpr> expr_t;
    return Expr<expr_t>( expr_t( elem_range, pt_range, expr ) );

}

} // vf
} // feel
#endif
