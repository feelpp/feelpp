/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Jan 2016

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
#define BOOST_TEST_MODULE test_forms
#include <feel/feelcore/testsuite.hpp>

#include <array>
#include <memory>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelvf/dirichletconstraints.hpp>

/** use Feel namespace */
using namespace Feel;
namespace bdata = boost::unit_test::data;

namespace
{
struct DeferredDirichletStrategyCase
{
    char const* label;
    char const* type;
    double tolerance;
};

std::ostream& operator<<( std::ostream& os, DeferredDirichletStrategyCase const& strategy )
{
    os << strategy.label;
    return os;
}

auto const deferred_dirichlet_strategies = std::array{
    DeferredDirichletStrategyCase{ .label = "elimination", .type = "elimination", .tolerance = 1e-10 },
    DeferredDirichletStrategyCase{ .label = "elimination_symmetric", .type = "elimination_symmetric", .tolerance = 1e-10 },
    DeferredDirichletStrategyCase{ .label = "penalisation", .type = "penalisation", .tolerance = 1e-8 }
};

auto const elimination_strategy = deferred_dirichlet_strategies.front();
auto const form2_test_dims = std::array{ 2, 3 };

struct DeferredDirichletMaterializationCounters
{
    int baseCloseCalls = 0;
    int zeroRowsCalls = 0;

    void reset() noexcept
    {
        baseCloseCalls = 0;
        zeroRowsCalls = 0;
    }
};

template<typename T>
class CountingMatrixSparse : public MatrixSparse<T>
{
public:
    using super = MatrixSparse<T>;
    using value_type = typename super::value_type;
    using real_type = typename super::real_type;
    using size_type = typename super::size_type;
    using graph_ptrtype = typename super::graph_ptrtype;
    using clone_ptrtype = typename super::clone_ptrtype;
    using vector_type = typename super::vector_type;

    CountingMatrixSparse( clone_ptrtype inner,
                          std::shared_ptr<DeferredDirichletMaterializationCounters> counters,
                          bool countBaseClose = true )
        :
        super( inner->mapRowPtr(), inner->mapColPtr(), inner->worldCommPtr() ),
        M_inner( std::move( inner ) ),
        M_counters( std::move( counters ) ),
        M_countBaseClose( countBaseClose )
    {
        this->syncState();
    }

    clone_ptrtype clone() const override
    {
        return std::make_shared<CountingMatrixSparse>( M_inner->clone(), M_counters, false );
    }

    void init( size_type m, size_type n, size_type m_l, size_type n_l, size_type nnz = 30, size_type noz = 10 ) override
    {
        M_inner->init( m, n, m_l, n_l, nnz, noz );
        this->syncState();
    }

    void init( size_type m, size_type n, size_type m_l, size_type n_l, graph_ptrtype const& graph ) override
    {
        M_inner->init( m, n, m_l, n_l, graph );
        this->syncState();
    }

    size_type nnz() const override { return M_inner->nnz(); }

    void clear() override
    {
        M_inner->clear();
        this->syncState();
    }

    void zero() override
    {
        M_inner->zero();
        this->syncState();
    }

    void zero( size_type start1, size_type stop1, size_type start2, size_type stop2 ) override
    {
        M_inner->zero( start1, stop1, start2, stop2 );
        this->syncState();
    }

    void close() const override
    {
        if ( M_countBaseClose )
            ++M_counters->baseCloseCalls;
        M_inner->close();
        const_cast<CountingMatrixSparse*>( this )->syncState();
    }

    void closeIfNeeded() const override
    {
        if ( !M_inner->closed() && M_countBaseClose )
            ++M_counters->baseCloseCalls;
        M_inner->closeIfNeeded();
        const_cast<CountingMatrixSparse*>( this )->syncState();
    }

    bool closed() const override
    {
        return M_inner->closed();
    }

    size_type size1() const override { return M_inner->size1(); }
    size_type size2() const override { return M_inner->size2(); }
    size_type rowStart() const override { return M_inner->rowStart(); }
    size_type rowStop() const override { return M_inner->rowStop(); }

    void set( size_type i, size_type j, value_type const& value ) override
    {
        M_inner->set( i, j, value );
        this->syncState();
    }

    void add( size_type i, size_type j, value_type const& value ) override
    {
        M_inner->add( i, j, value );
        this->syncState();
    }

    void addMatrix( const ublas::matrix<value_type>& dm,
                    std::vector<size_type> const& rows,
                    std::vector<size_type> const& cols ) override
    {
        M_inner->addMatrix( dm, rows, cols );
        this->syncState();
    }

    void addMatrix( int* rows, int nrows, int* cols, int ncols,
                    value_type* data, size_type K, size_type K2 ) override
    {
        M_inner->addMatrix( rows, nrows, cols, ncols, data, K, K2 );
        this->syncState();
    }

    void addMatrix( const ublas::matrix<value_type>& dm,
                    std::vector<size_type> const& dof_indices ) override
    {
        M_inner->addMatrix( dm, dof_indices );
        this->syncState();
    }

    void addMatrix( T const alpha, MatrixSparse<T> const& X, Feel::MatrixStructure matStruc = Feel::SAME_NONZERO_PATTERN ) override
    {
        M_inner->addMatrix( alpha, X, matStruc );
        this->syncState();
    }

    void scale( T const alpha ) override
    {
        M_inner->scale( alpha );
        this->syncState();
    }

    void multVector( Vector<T> const& arg, Vector<T>& dest, bool transpose ) const override
    {
        M_inner->multVector( arg, dest, transpose );
    }

    value_type operator()( size_type i, size_type j ) const override
    {
        return ( *M_inner )( i, j );
    }

    MatrixSparse<T>& operator=( MatrixSparse<value_type> const& M ) override
    {
        *M_inner = M;
        this->syncState();
        return *this;
    }

    void diagonal( Vector<T>& dest ) const override
    {
        M_inner->diagonal( dest );
    }

    void transpose( MatrixSparse<value_type>& Mt, size_type options = MATRIX_TRANSPOSE_ASSEMBLED ) const override
    {
        M_inner->transpose( Mt, options );
    }

    real_type energy( vector_type const& v, vector_type const& u, bool transpose = false ) const override
    {
        return M_inner->energy( v, u, transpose );
    }

    real_type l1Norm() const override { return M_inner->l1Norm(); }
    real_type linftyNorm() const override { return M_inner->linftyNorm(); }

    void zeroRows( std::vector<int> const& rows,
                   Vector<value_type> const& values,
                   Vector<value_type>& rhs,
                   Context const& on_context,
                   value_type value_on_diagonal ) override
    {
        ++M_counters->zeroRowsCalls;
        M_inner->zeroRows( rows, values, rhs, on_context, value_on_diagonal );
        this->syncState();
    }

    void updateBlockMat( std::shared_ptr<MatrixSparse<T>> const& m,
                         std::vector<size_type> const& start_i,
                         std::vector<size_type> const& start_j ) override
    {
        M_inner->updateBlockMat( m, start_i, start_j );
        this->syncState();
    }

private:
    void syncState()
    {
        this->setMapRow( M_inner->mapRowPtr() );
        this->setMapCol( M_inner->mapColPtr() );
        this->setInitialized( M_inner->isInitialized() );
        this->setIsClosed( M_inner->closed() );
        this->setGraph( M_inner->graph() );
        this->setIndexSplit( M_inner->indexSplit() );
    }

    clone_ptrtype M_inner;
    std::shared_ptr<DeferredDirichletMaterializationCounters> M_counters;
    bool M_countBaseClose = true;
};

template<typename MeshPtrType, typename UElementType, typename VElementType, typename BilinearFormType, typename LinearFormType>
void assembleScalarForm2System( MeshPtrType const& mesh,
                                UElementType const& u,
                                VElementType const& v,
                                BilinearFormType& a,
                                LinearFormType& l )
{
    a += integrate( _range=elements( mesh ),
                    _expr=inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
    l += integrate( _range=elements( mesh ),
                    _expr=cst( 1.0 ) * id( v ) );
}

template<typename MeshPtrType, typename BilinearFormType, typename LinearFormType, typename ElementType>
void applyDirichletCondition( BilinearFormType& a,
                              LinearFormType& l,
                              MeshPtrType const& mesh,
                              ElementType const& u,
                              char const* marker,
                              double value,
                              char const* type )
{
    a += on( _range=markedfaces( mesh, marker ),
             _rhs=l,
             _element=u,
             _expr=cst( value ),
             _type=type );
}

template<typename MeshPtrType, typename BilinearFormType, typename LinearFormType, typename ElementType>
void applyBoundaryPair( BilinearFormType& a,
                        LinearFormType& l,
                        MeshPtrType const& mesh,
                        ElementType const& u,
                        DeferredDirichletStrategyCase const& strategy,
                        double eastValue )
{
    applyDirichletCondition( a, l, mesh, u, "WEST", 0.0, strategy.type );
    applyDirichletCondition( a, l, mesh, u, "EAST", eastValue, strategy.type );
}

template<typename MeshPtrType, typename SolutionElementType>
void checkBoundaryPair( MeshPtrType const& mesh,
                        SolutionElementType const& solution,
                        double eastValue,
                        double tolerance )
{
    double const westError = normL2( _range=markedfaces( mesh, "WEST" ),
                                     _expr=idv( solution ) );
    double const eastError = normL2( _range=markedfaces( mesh, "EAST" ),
                                     _expr=idv( solution ) - cst( eastValue ) );

    BOOST_CHECK_SMALL( westError, tolerance );
    BOOST_CHECK_SMALL( eastError, tolerance );
}

template<typename BilinearFormType>
void enableDeferredDirichlet( BilinearFormType& a )
{
    a.deferDirichlet();
    BOOST_CHECK( a.useDeferredDirichlet() );
}

template<int Dim, typename Fn>
void withScalarForm2System( Fn&& fn )
{
    auto mesh = unitHypercube<Dim>();
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    backend( _rebuild=true );

    auto a = form2( _test=Vh, _trial=Vh );
    auto l = form1( _test=Vh );
    assembleScalarForm2System( mesh, u, v, a, l );

    std::forward<Fn>( fn )( mesh, Vh, u, v, a, l );
}

template<typename Fn>
void withScalarForm2Dim( int dim, Fn&& fn )
{
    switch ( dim )
    {
    case 2:
        std::forward<Fn>( fn ).template operator()<2>();
        break;
    case 3:
        std::forward<Fn>( fn ).template operator()<3>();
        break;
    default:
        BOOST_FAIL( "unsupported dimension" );
    }
}

template<int Dim>
void runRepeatedDeferredDirichletOnForm2( DeferredDirichletStrategyCase const& strategy )
{
    withScalarForm2System<Dim>( [&]( auto const& mesh, auto const& Vh, auto& u, auto&, auto& a, auto& l )
    {
        enableDeferredDirichlet( a );
        applyBoundaryPair( a, l, mesh, u, strategy, 1.0 );
        applyBoundaryPair( a, l, mesh, u, strategy, 1.0 );

        BOOST_CHECK( a.hasPendingDirichletConstraints() );
        BOOST_CHECK( a.hasDirichletConstraints() );

        auto solution = Vh->element();
        a.solve( _solution=solution, _rhs=l );
        BOOST_CHECK( !a.hasPendingDirichletConstraints() );
        BOOST_CHECK( a.hasDirichletConstraints() );

        checkBoundaryPair( mesh, solution, 1.0, strategy.tolerance );
    } );
}

template<int Dim>
void runManualApplyDeferredDirichletOnForm2( DeferredDirichletStrategyCase const& strategy )
{
    withScalarForm2System<Dim>( [&]( auto const& mesh, auto const& Vh, auto& u, auto&, auto& a, auto& l )
    {
        enableDeferredDirichlet( a );
        applyBoundaryPair( a, l, mesh, u, strategy, 1.0 );

        BOOST_CHECK( a.hasPendingDirichletConstraints() );
        BOOST_CHECK( a.hasDirichletConstraints() );
        BOOST_CHECK( a.supportsConstrainedOperatorView() );
        BOOST_CHECK( !a.hasMaterializedConstrainedOperator() );
        auto const baseMatrix = a.baseMatrixPtr();
        a.materializeConstrainedOperator( l );
        BOOST_CHECK( !a.hasPendingDirichletConstraints() );
        BOOST_CHECK( a.hasMaterializedConstrainedOperator() );
        auto const constrainedMatrix = a.constrainedMatrixPtr( l );
        auto const constrainedVector = a.constrainedVectorPtr( l );
        BOOST_CHECK( constrainedMatrix != baseMatrix );
        BOOST_CHECK( constrainedVector );
        BOOST_CHECK( !a.hasPendingDirichletConstraints() );
        BOOST_CHECK( a.hasDirichletConstraints() );
        BOOST_CHECK( a.activeMatrixPtr( l ) == constrainedMatrix );
        BOOST_CHECK( a.activeVectorPtr( l ) == constrainedVector );
        a.applyDeferredDirichlet( l );
        BOOST_CHECK( !a.hasPendingDirichletConstraints() );
        BOOST_CHECK( a.hasDirichletConstraints() );
        BOOST_CHECK( a.closed() );

        auto solution = Vh->element();
        a.solve( _solution=solution, _rhs=l );

        checkBoundaryPair( mesh, solution, 1.0, strategy.tolerance );
    } );
}

template<int Dim>
void runDeferredDirichletRhsReuseRefreshOnForm2( DeferredDirichletStrategyCase const& strategy )
{
    withScalarForm2System<Dim>( [&]( auto const& mesh, auto const& Vh, auto& u, auto&, auto& a, auto& l )
    {
        enableDeferredDirichlet( a );
        applyBoundaryPair( a, l, mesh, u, strategy, 1.0 );

        auto constrainedVector = a.constrainedVectorPtr( l );
        auto constrainedSnapshot = vf::cloneVectorWithValues( constrainedVector );
        auto rhsVector = l.vectorPtr();
        auto const initialRevision = rhsVector->revision();
        rhsVector->add( 3.0 );
        BOOST_CHECK_GT( rhsVector->revision(), initialRevision );

        auto refreshedConstrainedVector = a.constrainedVectorPtr( l );
        BOOST_CHECK( refreshedConstrainedVector != constrainedVector );

        auto constrainedDelta = vf::cloneVectorWithValues( refreshedConstrainedVector );
        constrainedDelta->add( -1.0, *constrainedSnapshot );
        if ( !constrainedDelta->closed() )
            constrainedDelta->close();
        BOOST_CHECK_GT( constrainedDelta->linftyNorm(), 1e-12 );

        a.applyDeferredDirichlet( l );
        BOOST_CHECK( a.activeVectorPtr( l ) == refreshedConstrainedVector );

        auto solution = Vh->element();
        a.solve( _solution=solution, _rhs=l );
        checkBoundaryPair( mesh, solution, 1.0, strategy.tolerance );
    } );
}

template<int Dim>
void runClearDeferredDirichletOnForm2()
{
    withScalarForm2System<Dim>( [&]( auto const& mesh, auto const& Vh, auto& u, auto&, auto& a, auto& l )
    {
        enableDeferredDirichlet( a );
        applyBoundaryPair( a, l, mesh, u, elimination_strategy, 1.0 );

        BOOST_CHECK( a.hasPendingDirichletConstraints() );
        BOOST_CHECK( a.hasDirichletConstraints() );
        a.clearDeferredDirichlet();
        BOOST_CHECK( !a.hasPendingDirichletConstraints() );
        BOOST_CHECK( !a.hasDirichletConstraints() );

        applyBoundaryPair( a, l, mesh, u, elimination_strategy, 2.0 );

        BOOST_CHECK( a.hasPendingDirichletConstraints() );
        BOOST_CHECK( a.hasDirichletConstraints() );

        auto solution = Vh->element();
        a.solve( _solution=solution, _rhs=l );
        BOOST_CHECK( !a.hasPendingDirichletConstraints() );
        BOOST_CHECK( a.hasDirichletConstraints() );

        checkBoundaryPair( mesh, solution, 2.0, elimination_strategy.tolerance );
    } );
}

template<int Dim>
void runImmediateDirichletCompatibilityOnForm2( DeferredDirichletStrategyCase const& strategy )
{
    withScalarForm2System<Dim>( [&]( auto const& mesh, auto const& Vh, auto& u, auto&, auto& a, auto& l )
    {
        a.immediateDirichlet();
        BOOST_CHECK( !a.useDeferredDirichlet() );

        applyBoundaryPair( a, l, mesh, u, strategy, 1.0 );
        applyBoundaryPair( a, l, mesh, u, strategy, 1.0 );

        BOOST_CHECK( !a.hasPendingDirichletConstraints() );
        BOOST_CHECK( !a.hasDirichletConstraints() );

        auto solution = Vh->element();
        a.solve( _solution=solution, _rhs=l );

        BOOST_CHECK( !a.hasPendingDirichletConstraints() );
        BOOST_CHECK( !a.hasDirichletConstraints() );
        checkBoundaryPair( mesh, solution, 1.0, strategy.tolerance );
    } );
}

template<int Dim>
void runDeferredDirichletCallCountTargetOnForm2()
{
    auto mesh = unitHypercube<Dim>();
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    backend( _rebuild=true );

    auto counters = std::make_shared<DeferredDirichletMaterializationCounters>();
    auto countingMatrix = std::make_shared<CountingMatrixSparse<double>>( backend()->newMatrix( _test=Vh, _trial=Vh ),
                                                                          counters );
    auto a = form2( _test=Vh, _trial=Vh, _matrix=countingMatrix );
    auto l = form1( _test=Vh );

    assembleScalarForm2System( mesh, u, v, a, l );
    enableDeferredDirichlet( a );
    applyBoundaryPair( a, l, mesh, u, elimination_strategy, 1.0 );

    BOOST_CHECK( !a.hasMaterializedConstrainedOperator() );
    counters->reset();
    auto constrainedMatrix = a.activeMatrixPtr( l );
    auto constrainedVector = a.activeVectorPtr( l );
    auto constrainedMatrixAgain = a.activeMatrixPtr( l );
    auto constrainedVectorAgain = a.activeVectorPtr( l );

    BOOST_CHECK( constrainedMatrix );
    BOOST_CHECK( constrainedVector );
    BOOST_CHECK( a.hasMaterializedConstrainedOperator() );
    BOOST_CHECK( constrainedMatrix != a.baseMatrixPtr() );
    BOOST_CHECK( dynamic_cast<CountingMatrixSparse<double>*>( constrainedMatrix.get() ) != nullptr );
    BOOST_CHECK_EQUAL( counters->baseCloseCalls, 1 );
    BOOST_CHECK_EQUAL( constrainedMatrix.get(), constrainedMatrixAgain.get() );
    BOOST_CHECK_EQUAL( constrainedVector.get(), constrainedVectorAgain.get() );
}

template<int Dim>
void runDeferredDirichletSingleApplyZeroRowsCount()
{
    if ( Environment::worldComm().globalSize() != 1 )
        return;

    auto mesh = unitHypercube<Dim>();
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    auto probe = Vh->element();

    backend( _rebuild=true );

    auto counters = std::make_shared<DeferredDirichletMaterializationCounters>();
    auto countingMatrix = std::make_shared<CountingMatrixSparse<double>>( backend()->newMatrix( _test=Vh, _trial=Vh ),
                                                                          counters );
    auto a = form2( _test=Vh, _trial=Vh, _matrix=countingMatrix );
    auto l = form1( _test=Vh );

    assembleScalarForm2System( mesh, u, v, a, l );

    auto const& rowMap = countingMatrix->mapRow();
    BOOST_REQUIRE_GT( rowMap.nLocalDofWithoutGhost(), 0 );

    std::vector<int> localDofs;
    localDofs.push_back( static_cast<int>( rowMap.firstDofGlobalCluster() ) );
    std::vector<double> values( localDofs.size(), 0.0 );
    vf::DeferredDirichletSet<double> constraints;
    constraints.append( localDofs, values, Feel::Context( ContextOn::ELIMINATION ), 1.0 );

    countingMatrix->close();
    l.vectorPtr()->close();
    counters->reset();

    auto const merged = constraints.mergedEntries();
    BOOST_REQUIRE_EQUAL( merged.size(), 1 );
    vf::applyDeferredDirichletEntries( merged, countingMatrix, l.vectorPtr() );

    BOOST_CHECK_EQUAL( counters->zeroRowsCalls, 1 );
}

template<int Dim>
void runDeferredDirichletApplyZeroRowsCountWithEmptyLocalRows()
{
    if ( Environment::worldComm().globalSize() < 2 )
        return;

    auto mesh = unitHypercube<Dim>();
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    backend( _rebuild=true );

    auto counters = std::make_shared<DeferredDirichletMaterializationCounters>();
    auto countingMatrix = std::make_shared<CountingMatrixSparse<double>>( backend()->newMatrix( _test=Vh, _trial=Vh ),
                                                                          counters );
    auto a = form2( _test=Vh, _trial=Vh, _matrix=countingMatrix );
    auto l = form1( _test=Vh );

    assembleScalarForm2System( mesh, u, v, a, l );

    auto const& rowMap = countingMatrix->mapRow();
    BOOST_REQUIRE_GT( rowMap.nLocalDofWithoutGhost(), 0 );

    std::vector<int> localDofs;
    if ( Environment::worldComm().globalRank() == 0 )
        localDofs.push_back( 0 );
    std::vector<double> values( localDofs.size(), 0.0 );

    vf::DeferredDirichletSet<double> constraints;
    constraints.append( localDofs, values, Feel::Context( ContextOn::ELIMINATION ), 1.0 );

    countingMatrix->close();
    l.vectorPtr()->close();
    counters->reset();

    auto const merged = constraints.mergedEntries();
    BOOST_REQUIRE_EQUAL( merged.size(), 1 );
    vf::applyDeferredDirichletEntries( merged, countingMatrix, l.vectorPtr() );

    BOOST_CHECK_EQUAL( counters->zeroRowsCalls, 1 );
}
} // namespace

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Forms  Options" );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_forms" ,
                     "test_forms" ,
                     "0.2",
                     "nD(n=2,3) test forms",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );

    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}



FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( forms_suite )

using dim_t = boost::mpl::list<boost::mpl::int_<2>, boost::mpl::int_<3> >;
//using dim_t = boost::mpl::list<boost::mpl::int_<2>>;
//using dim_t = boost::mpl::list<boost::mpl::int_<3>>;


BOOST_AUTO_TEST_CASE_TEMPLATE( test_form2_faces, T, dim_t )
{
    using Feel::cout;
    BOOST_MESSAGE( "test_form2_faces starts for dim=" << T::value);
    auto meshnd = unitHypercube<T::value>();
    auto mesh = createSubmesh( _mesh=meshnd, _range=faces(meshnd), _update=0 );
    size_type nFaceInParallelMeshnd = nelements(faces(meshnd),true);
    size_type nInternalFaceInPara1llelMeshnd = nelements(faces(meshnd),true);
    BOOST_CHECK_EQUAL( nelements(elements(mesh),true), nFaceInParallelMeshnd  );
    auto Vh=Pdhv<1>(meshnd);
    auto u = Vh->element();
    auto Wh=Pdh<1>(meshnd);
    auto p = Wh->element();
    p.on(_range=elements(meshnd),_expr=cst(1.),_close=true);
    auto Mh=Pdh<1>(mesh);
    auto l=Mh->element();
    l.on(_range=elements(mesh),_expr=cst(1.));

    // all the 4 integrals below should be the same
    // compute \int 1 over internalfaces
    auto I = integrate( _range=internalfaces(meshnd), _expr=cst(1.)).evaluate()(0,0)
        +integrate( _range=boundaryfaces(meshnd), _expr=cst(1.)).evaluate()(0,0);
    // compute \int 1 over d-1 mesh
    auto I1 = integrate( _range=elements(mesh), _expr=cst(1.)).evaluate()(0,0);
    // compute \int 1 over internalfaces using left element
    auto I2 = integrate( _range=internalfaces(meshnd), _expr=leftfacev(cst(1.))).evaluate()(0,0)
        +integrate( _range=boundaryfaces(meshnd), _expr=cst(1.)).evaluate()(0,0);
    // compute \int 1 over internalfaces using right element
    auto I3 = integrate( _range=internalfaces(meshnd), _expr=rightfacev(cst(1.))).evaluate()(0,0)
        +integrate( _range=boundaryfaces(meshnd), _expr=cst(1.)).evaluate()(0,0);
    cout << "I=" << I << " I1=" << I1 << " I2=" << I2 << " I3=" << I3 << std::endl;
    BOOST_CHECK_CLOSE( I, I1, 1e-10 );
    BOOST_CHECK_CLOSE( I, I2, 1e-10 );
    BOOST_CHECK_CLOSE( I, I3, 1e-10 );

    auto e=inner(P(),one());
    LOG(INFO) << "a start";
    auto a = form2(_test=Mh, _trial=Wh );
    a = integrate( _range=internalfaces(meshnd), _expr=(e*id(l))*(leftfacet(idt(p)/e)))
        + integrate( _range=boundaryfaces(meshnd), _expr=(e*id(l))*(idt(p)/e));
    a.close();
    cout << "a(1,1) = " << a(l,p) << std::endl;
    BOOST_CHECK_CLOSE( a(l,p), I1, 1e-10 );
    LOG(INFO) << "a done";

    LOG(INFO) << "ai start";
    a = integrate( _range=internalfaces(meshnd), _expr=(e*id(l))*(leftfacet(idt(p)/e)));
    a.close();
    cout << "ai(1,1) = " << a(l,p) << std::endl;
    LOG(INFO) << "ai done";

    LOG(INFO) << "ab start";
    a = integrate( _range=boundaryfaces(meshnd), _expr=(e*id(l))*(leftfacet(idt(p)/e)));
    a.close();
    cout << "ab(1,1) = " << a(l,p) << std::endl;
    LOG(INFO) << "ab done";

    LOG(INFO) << "a1 start";
    auto a1 = form2(_test=Mh, _trial=Mh );
    a1 = integrate( _range=internalfaces(meshnd), _expr=e*id(l)*idt(l)/(2*e))
        + integrate( _range=boundaryfaces(meshnd), _expr=e*id(l)*idt(l)/e);
    a1.close();
    auto a1en = a1(l,l);
    BOOST_CHECK_CLOSE( a1en, I1, 1e-10 );
    cout << "a1(1,1)=" << a1en << std::endl;
    LOG(INFO) << "a1 done";

    auto a11 = form2(_test=Mh, _trial=Mh );
    a11 = integrate( _range=elements(mesh), _expr=id(l)*idt(l));
    a11.close();
    auto a11en = a11(l,l);
    BOOST_CHECK_CLOSE( a11en, I1, 1e-10 );
    cout << "a11(1,1)=" << a11en << " int =" << I1 << std::endl;

    auto a111 = form2(_test=Mh, _trial=Mh );
    a111 = integrate( _range=internalfaces(meshnd), _expr=id(l)*idt(l)/2)
        + integrate( _range=boundaryfaces(meshnd), _expr=id(l)*idt(l));
    a111.close();
    auto a111en = a111(l,l);
    BOOST_CHECK_CLOSE( a111en, I1, 1e-10 );
    cout << "a111(1,1)=" << a111en << " int =" << I1 << std::endl;

    // - Tests A22 = <ph, w> with ph, w \in Wh
    auto a2 = form2(_test=Wh, _trial=Wh );
    a2 = integrate( _range=internalfaces(meshnd), _expr=leftface(e*id(p))*leftfacet(idt(p)/e))
        +integrate( _range=boundaryfaces(meshnd), _expr=(e*id(p))*(idt(p)/e));
    a2.close();
    auto a2en = a2(p,p);
    BOOST_CHECK_CLOSE( a2en, I1, 1e-10 );
    cout << "a2(1,1)=" << a2en << std::endl;

    // - Same as a2, but with rightface instead of leftface.
    auto a3 = form2(_test=Wh, _trial=Wh );
    a3 = integrate(_range=internalfaces(meshnd), _expr=rightface(e*id(p))*rightfacet(idt(p)/e))
        +integrate( _range=boundaryfaces(meshnd), _expr=(e*id(p))*(idt(p)/e));
    a3.close();
    auto a3en = a3(p,p);
    BOOST_CHECK_CLOSE( a3en, I1, 1e-10 );
    cout << "a3(1,1)=" << a3en << std::endl;

    auto a23 = form2(_test=Wh, _trial=Wh );
    a23 = integrate(_range=internalfaces(meshnd), _expr=rightface(e*id(p))*leftfacet(idt(p)/e));
    a23.close();
    auto a23en = a23(p,p);
    BOOST_CHECK_SMALL( a23en, 1e-10 );
    cout << "a23(1,1)=" << a23en << std::endl;

    auto a231 = form2(_test=Wh, _trial=Wh, _pattern=size_type(Pattern::EXTENDED) );
    a231 = integrate(_range=internalfaces(meshnd), _expr=rightface(e*id(p))*leftfacet(idt(p)/e))
        +integrate( _range=boundaryfaces(meshnd), _expr=(e*id(p))*(idt(p)/e));
    a231.close();
    auto a231en = a231(p,p);
    BOOST_CHECK_CLOSE( a231en, I1, 1e-10 );
    cout << "a231(1,1)=" << a231en << std::endl;

    // - Tests A23 = <phat, w>, with phat \in Mh, w \in Wh
    //
    auto a4 = form2(_test=Wh, _trial=Mh );
    //auto a4 = form2(_test=Wh, _trial=Mh );
    a4 = integrate( _range=internalfaces(meshnd), _expr=leftface(e*id(p))*idt(l)/e );
    a4.close();
    auto a4en = a4(p, l);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a4l(1,1)=" << a4en );
    a4 = integrate( _range=internalfaces(meshnd), _expr=rightface(e*id(p))*idt(l)/e );
    a4.close();
    a4en = a4(p, l);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a4r(1,1)=" << a4en );
    a4 = integrate( _range=internalfaces(meshnd), _expr=0.5*(leftface(e*id(p))+rightface(e*id(p)))*idt(l)/e );
    a4.close();
    a4en = a4(p, l);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a4lr(1,1)=" << a4en );

    // - Tests A32 = <ph, \mu>, with ph \in Wh, \mu \in Mh
    auto a5 = form2(_test=Mh, _trial=Wh );
    a5 = integrate(_range=internalfaces(meshnd), _expr=leftfacet(e*idt(p))*id(l)/e );
    a5.close();
    auto a5en = a5(l, p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a5l(1,1)=" << a5en );
    a5 = integrate(_range=internalfaces(meshnd), _expr=rightfacet(e*idt(p))*id(l)/e );
    a5.close();
    a5en = a5(l, p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a5r(1,1)=" << a5en );

    // - Tests (p, p)_Omega. Should give the measure of the domain
    auto a6 = form2(_test=Wh, _trial=Wh);
    a6 = integrate(_range=elements(meshnd), _expr=idt(p)*id(p));
    a6.close();
    auto a6en = a6(p, p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a6(1,1)=" << a6en );

    double Ib = integrate(_range=boundaryfaces(meshnd), _expr=cst(1.)).evaluate()(0,0);

    // The five following tests should all provide the (n-1)-Lebesgue
    // measure of the boundary of the domain.
    // The goal is to check whether we need to multiply by 0.5 or 0.25 or nothing in
    // boundary integrals. Also
    auto a7 = form2( _test=Wh, _trial=Wh);
    a7 = integrate(_range=boundaryfaces(meshnd), _expr=idt(p)*id(p));
    a7.close();
    auto a7eval = a7(p,p);
    BOOST_CHECK_CLOSE( a7eval, Ib, 1e-10 );
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a7(1,1)=" << a7eval );

    // a8
    auto a8 = form2( _test=Mh, _trial=Mh);
    a8 = integrate(_range=boundaryfaces(meshnd), _expr=idt(l)*id(l));
    a8.close();
    auto a8eval = a8(l,l);
    BOOST_CHECK_CLOSE( a8eval, Ib, 1e-10 );
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a8(1,1)=" << a8eval );

    auto a9 = form2(_test=Mh, _trial=Wh);
    a9 = integrate(_range=boundaryfaces(meshnd), _expr=idt(p)*id(l));
    a9.close();
    auto a9eval = a9(l,p);
    BOOST_CHECK_CLOSE( a9eval, Ib, 1e-10 );
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a9(1,1)=" << a9eval );

    auto a10 = form1(_test=Mh);
    a10 = integrate(_range=boundaryfaces(meshnd), _expr=id(l));
    a10.close();
    auto a10eval = a10(l);
    BOOST_CHECK_CLOSE( a10eval, Ib, 1e-10 );
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a10(1)=" << a10eval );

    auto a12 = form1(_test=Wh);
    a12 = integrate(_range=boundaryfaces(meshnd), _expr=id(p));
    a12.close();
    auto a12eval = a12(p);
    BOOST_CHECK_CLOSE( a12eval, Ib, 1e-10 );
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a12(1)=" << a12eval );



    // ------------------------------------------------------------------------------------

    // The following tests help to understand what we have to do when mixing
    // integrals that require the extended pattern with others that do not

    LOG(INFO) << "a7b starts";
    // Works fine
    auto a7b = form2( _test=Wh, _trial=Wh);
    a7b = integrate(_range=boundaryfaces(meshnd), _expr=idt(p)*id(p));
    a7b.close();
    auto a7beval = a7b(p,p);
    BOOST_CHECK_CLOSE( a7beval, Ib, 1e-10 );
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a7b(1,1)=" << a7beval );
    LOG(INFO) << "a7b ends";

    auto a8b = form2( _test=Mh, _trial=Mh);
    a8b = integrate(_range=boundaryfaces(meshnd), _expr=idt(l)*id(l));
    a8b.close();
    auto a8beval = a8b(l,l);
    BOOST_CHECK_CLOSE( a8beval, Ib, 1e-10 );
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a8b(1,1)=" << a8beval );

    LOG(INFO) << "a9b starts";
    auto a9b = form2(_test=Mh, _trial=Wh, _pattern=size_type(Pattern::EXTENDED));
    a9b = integrate(_range=boundaryfaces(meshnd), _expr=idt(p)*id(l));
    a9b.close();
    auto a9beval = a9b(l,p);
    BOOST_CHECK_CLOSE( a9beval, Ib, 1e-10 );
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a9b(1,1)=" << a9beval );
    LOG(INFO) << "a9b ends";

    auto a10b = form1(_test=Mh, _pattern=size_type(Pattern::EXTENDED));
    a10b = integrate(_range=boundaryfaces(meshnd), _expr=id(l));
    a10b.close();
    auto a10beval = a10b(l);
    BOOST_CHECK_CLOSE( a10beval, Ib, 1e-10 );
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a10b(1)=" << a10beval );

    LOG(INFO) << "a12b starts";
    auto a12b = form1(_test=Wh, _pattern=size_type(Pattern::EXTENDED));
    a12b = integrate(_range=boundaryfaces(meshnd), _expr=id(p));
    a12b.close();
    auto a12beval = a12b(p);
    BOOST_CHECK_CLOSE( a12beval, Ib, 1e-10 );
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a12b(1)=" << a12beval );
    LOG(INFO) << "a12b ends";
    BOOST_MESSAGE( "test_form2_faces ends for dim=" << T::value);
}

BOOST_AUTO_TEST_CASE( test_deferred_dirichlet_collective_zero_rows_with_empty_local_rows )
{
    runDeferredDirichletApplyZeroRowsCountWithEmptyLocalRows<2>();
}

BOOST_DATA_TEST_CASE( test_repeated_dirichlet_on_form2,
                      bdata::make( form2_test_dims ) * bdata::make( deferred_dirichlet_strategies ),
                      dim,
                      strategy )
{
    withScalarForm2Dim( dim,
                        [&]<int Dim>()
                        {
                            runRepeatedDeferredDirichletOnForm2<Dim>( strategy );
                        } );
}

BOOST_DATA_TEST_CASE( test_immediate_dirichlet_default_on_form2,
                      bdata::make( form2_test_dims ) * bdata::make( deferred_dirichlet_strategies ),
                      dim,
                      strategy )
{
    withScalarForm2Dim( dim,
                        [&]<int Dim>()
                        {
                            runImmediateDirichletCompatibilityOnForm2<Dim>( strategy );
                        } );
}

BOOST_DATA_TEST_CASE( test_manual_apply_deferred_dirichlet_on_form2,
                      bdata::make( form2_test_dims ) * bdata::make( deferred_dirichlet_strategies ),
                      dim,
                      strategy )
{
    withScalarForm2Dim( dim,
                        [&]<int Dim>()
                        {
                            runManualApplyDeferredDirichletOnForm2<Dim>( strategy );
                        } );
}

BOOST_DATA_TEST_CASE( test_deferred_dirichlet_rhs_reuse_refresh_on_form2,
                      bdata::make( form2_test_dims ) * bdata::make( deferred_dirichlet_strategies ),
                      dim,
                      strategy )
{
    withScalarForm2Dim( dim,
                        [&]<int Dim>()
                        {
                            runDeferredDirichletRhsReuseRefreshOnForm2<Dim>( strategy );
                        } );
}

BOOST_DATA_TEST_CASE( test_clear_deferred_dirichlet_on_form2,
                      bdata::make( form2_test_dims ),
                      dim )
{
    withScalarForm2Dim( dim,
                        [&]<int Dim>()
                        {
                            runClearDeferredDirichletOnForm2<Dim>();
                        } );
}

BOOST_DATA_TEST_CASE( test_deferred_dirichlet_call_count_target_on_form2,
                      bdata::make( form2_test_dims ),
                      dim )
{
    withScalarForm2Dim( dim,
                        [&]<int Dim>()
                        {
                            runDeferredDirichletCallCountTargetOnForm2<Dim>();
                        } );
}

BOOST_DATA_TEST_CASE( test_deferred_dirichlet_single_apply_zero_rows_count,
                      bdata::make( form2_test_dims ),
                      dim )
{
    withScalarForm2Dim( dim,
                        [&]<int Dim>()
                        {
                            runDeferredDirichletSingleApplyZeroRowsCount<Dim>();
                        } );
}

BOOST_AUTO_TEST_CASE( test_deferred_dirichlet_merge )
{
    vf::DeferredDirichletSet<double> constraints;
    constraints.append( { 5, 1, 5 }, { 2.0, 1.0, 2.0 }, Feel::Context( ContextOn::ELIMINATION ), 1.0 );
    constraints.append( { 3, 1 }, { 4.0, 1.0 }, Feel::Context( ContextOn::ELIMINATION ), 1.0 );

    auto const merged = constraints.mergedEntries();
    BOOST_REQUIRE_EQUAL( merged.size(), 1 );

    auto const expectedDofs = std::array<int, 3>{ 1, 3, 5 };
    auto const expectedValues = std::array<double, 3>{ 1.0, 4.0, 2.0 };
    BOOST_CHECK_EQUAL_COLLECTIONS( merged.front().dofs.begin(), merged.front().dofs.end(),
                                   expectedDofs.begin(), expectedDofs.end() );
    BOOST_CHECK_EQUAL_COLLECTIONS( merged.front().values.begin(), merged.front().values.end(),
                                   expectedValues.begin(), expectedValues.end() );
}

BOOST_AUTO_TEST_CASE( test_deferred_dirichlet_groups_by_application_context )
{
    vf::DeferredDirichletSet<double> constraints;
    constraints.append( { 1 }, { 0.0 }, Feel::Context( ContextOn::ELIMINATION ), 1.0 );
    constraints.append( { 2 }, { 0.0 }, Feel::Context( ContextOn::ELIMINATION ), 2.0 );
    constraints.append( { 3 }, { 0.0 }, Feel::Context( ContextOn::ELIMINATION | ContextOn::SYMMETRIC ), 1.0 );

    auto const merged = constraints.mergedEntries();
    BOOST_REQUIRE_EQUAL( merged.size(), 3 );
}

BOOST_AUTO_TEST_CASE( test_deferred_dirichlet_last_write_wins )
{
    vf::DeferredDirichletSet<double> constraints;
    constraints.append( { 1 }, { 0.0 }, Feel::Context( ContextOn::ELIMINATION ), 1.0 );
    constraints.append( { 1 }, { 2.0 }, Feel::Context( ContextOn::ELIMINATION ), 1.0 );

    BOOST_REQUIRE_NO_THROW( static_cast<void>( constraints.mergedEntries() ) );
    auto const merged = constraints.mergedEntries();
    BOOST_REQUIRE_EQUAL( merged.size(), 1 );
    BOOST_REQUIRE_EQUAL( merged.front().dofs.size(), 1 );
    BOOST_REQUIRE_EQUAL( merged.front().values.size(), 1 );
    BOOST_CHECK_EQUAL( merged.front().dofs.front(), 1 );
    BOOST_CHECK_EQUAL( merged.front().values.front(), 2.0 );
}

BOOST_AUTO_TEST_CASE( test_deferred_dirichlet_precedence_target )
{
    vf::DeferredDirichletSet<double> constraints;

    constraints.append( { 20, 50 }, { 32.0, 32.0 }, Feel::Context( ContextOn::ELIMINATION ), 1.0,
                        vf::deferredDirichletEntityPriority( vf::DeferredDirichletEntity::point ) );
    constraints.append( { 10, 20 }, { 2.0, 2.0 }, Feel::Context( ContextOn::ELIMINATION ), 1.0,
                        vf::deferredDirichletEntityPriority( vf::DeferredDirichletEntity::element ) );
    constraints.append( { 20, 40 }, { 22.0, 22.0 }, Feel::Context( ContextOn::ELIMINATION ), 1.0,
                        vf::deferredDirichletEntityPriority( vf::DeferredDirichletEntity::edge ) );
    constraints.append( { 20, 30 }, { 12.0, 12.0 }, Feel::Context( ContextOn::ELIMINATION ), 1.0,
                        vf::deferredDirichletEntityPriority( vf::DeferredDirichletEntity::face ) );

    BOOST_REQUIRE_NO_THROW( static_cast<void>( constraints.mergedEntries() ) );
    auto const merged = constraints.mergedEntries();
    BOOST_REQUIRE_EQUAL( merged.size(), 1 );

    auto const expectedDofs = std::array<int, 5>{ 10, 20, 30, 40, 50 };
    auto const expectedValues = std::array<double, 5>{ 2.0, 32.0, 12.0, 22.0, 32.0 };
    BOOST_CHECK_EQUAL_COLLECTIONS( merged.front().dofs.begin(), merged.front().dofs.end(),
                                   expectedDofs.begin(), expectedDofs.end() );
    BOOST_CHECK_EQUAL_COLLECTIONS( merged.front().values.begin(), merged.front().values.end(),
                                   expectedValues.begin(), expectedValues.end() );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( test_form1, T, dim_t )
{
    BOOST_MESSAGE( "test_form2_faces starts for dim=" << T::value);
    auto meshnd = unitHypercube<T::value>();
    auto mesh = createSubmesh( _mesh=meshnd, _range=faces(meshnd), _update=0 );
    auto Vh=Pdhv<1>(meshnd);
    auto u = Vh->element();
    auto Wh=Pdh<1>(meshnd);
    auto p = Wh->element();

    auto a = form1(_test=Wh);
    a = integrate( _range=elements(meshnd), _expr=id(p));
    a.close();
    auto b = a;
    b+=a;
    auto c = a+b;
    c-=3*a;
    BOOST_CHECK_SMALL( c.vectorPtr()->linftyNorm(), 1e-12 );
    b /= 2;
    BOOST_CHECK_SMALL( (b-a).vectorPtr()->linftyNorm(), 1e-12 );
}
BOOST_AUTO_TEST_SUITE_END()
