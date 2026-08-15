/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
#define BOOST_TEST_MODULE test_product
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelalg/backendpetsc.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <array>
#include <cstdlib>

#if defined( __unix__ )
#include <sys/wait.h>
#include <unistd.h>
#endif

/** use Feel namespace */
using namespace Feel;
namespace bdata = boost::unit_test::data;

namespace
{
struct DirichletStrategyCase
{
    std::string label;
    std::string type;
    double tolerance;
};

std::ostream& operator<<( std::ostream& os, DirichletStrategyCase const& strategy )
{
    os << strategy.label;
    return os;
}

auto const monolithicDirichletStrategies = std::array{
    DirichletStrategyCase{ .label = "elimination", .type = "elimination", .tolerance = 1e-10 },
    DirichletStrategyCase{ .label = "elimination_symmetric", .type = "elimination_symmetric", .tolerance = 1e-10 },
    DirichletStrategyCase{ .label = "elimination_keep_diagonal", .type = "elimination_keep_diagonal", .tolerance = 1e-10 },
    DirichletStrategyCase{ .label = "elimination_symmetric_keep_diagonal", .type = "elimination_symmetric_keep_diagonal", .tolerance = 1e-10 },
    DirichletStrategyCase{ .label = "penalisation", .type = "penalisation", .tolerance = 1e-8 }
};

auto const eliminationStrategy = monolithicDirichletStrategies.front();

#if defined( FEELPP_HAS_PETSC_H )
struct DeferredDirichletMaterializationCounters
{
    int baseCloseCalls = 0;
    int zeroRowsCalls = 0;
    int cloneCalls = 0;

    void reset() noexcept
    {
        baseCloseCalls = 0;
        zeroRowsCalls = 0;
        cloneCalls = 0;
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
        ++M_counters->cloneCalls;
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

class CountingBackendPetsc : public BackendPetsc<double>
{
public:
    using super = BackendPetsc<double>;
    using size_type = typename super::size_type;
    using sparse_matrix_ptrtype = typename super::sparse_matrix_ptrtype;
    using graph_ptrtype = typename super::graph_ptrtype;
    using datamap_ptrtype = typename super::datamap_ptrtype;

    explicit CountingBackendPetsc( std::shared_ptr<DeferredDirichletMaterializationCounters> counters,
                                   worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        :
        super( worldComm ),
        M_counters( std::move( counters ) )
    {}

    sparse_matrix_ptrtype newMatrix() override
    {
        return this->wrapMatrix( super::newMatrix() );
    }

    sparse_matrix_ptrtype newMatrix( size_type m,
                                     size_type n,
                                     size_type m_l,
                                     size_type n_l,
                                     size_type nnz = 30,
                                     size_type noz = 10,
                                     size_type matrix_properties = NON_HERMITIAN ) override
    {
        return this->wrapMatrix( super::newMatrix( m, n, m_l, n_l, nnz, noz, matrix_properties ) );
    }

    sparse_matrix_ptrtype newMatrix( datamap_ptrtype const& domainmap,
                                     datamap_ptrtype const& imagemap,
                                     size_type matrix_properties = NON_HERMITIAN,
                                     bool init = true ) override
    {
        return this->wrapMatrix( super::newMatrix( domainmap, imagemap, matrix_properties, init ) );
    }

    sparse_matrix_ptrtype newMatrix( size_type m,
                                     size_type n,
                                     size_type m_l,
                                     size_type n_l,
                                     graph_ptrtype const& graph,
                                     size_type matrix_properties = NON_HERMITIAN ) override
    {
        return this->wrapMatrix( super::newMatrix( m, n, m_l, n_l, graph, matrix_properties ) );
    }

    sparse_matrix_ptrtype newZeroMatrix( datamap_ptrtype const& domainmap,
                                         datamap_ptrtype const& imagemap ) override
    {
        return this->wrapMatrix( super::newZeroMatrix( domainmap, imagemap ) );
    }

    sparse_matrix_ptrtype newZeroMatrix( size_type m,
                                         size_type n,
                                         size_type m_l,
                                         size_type n_l ) override
    {
        return this->wrapMatrix( super::newZeroMatrix( m, n, m_l, n_l ) );
    }

private:
    sparse_matrix_ptrtype wrapMatrix( sparse_matrix_ptrtype inner ) const
    {
        return std::make_shared<CountingMatrixSparse<double>>( std::move( inner ), M_counters );
    }

    std::shared_ptr<DeferredDirichletMaterializationCounters> M_counters;
};
#endif

template<typename MeshPtrType, typename WElementType, typename TElementType, typename BlockFormType, typename BlockLinearFormType>
void assembleMonolithicBlockformSystem( MeshPtrType const& mesh,
                                        WElementType const& W,
                                        TElementType const& T,
                                        BlockFormType& a,
                                        BlockLinearFormType& l )
{
    auto u = W( 0_c );
    auto p = W( 1_c );
    auto v = T( 0_c );
    auto q = T( 1_c );

    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
    a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( p ), grad( q ) ) + idt( p ) * id( q ) );

    l( 0_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 1.0 ) * id( v ) );
    l( 1_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 2.0 ) * id( q ) );
}

template<typename MeshPtrType, typename BlockFormType, typename BlockLinearFormType, typename ElementType>
void applyBoundaryPairOnBlockform( BlockFormType& a,
                                   BlockLinearFormType& l,
                                   MeshPtrType const& mesh,
                                   ElementType const& u,
                                   DirichletStrategyCase const& strategy,
                                   double eastValue )
{
    a.row( 0_c ) += on( _range=markedfaces( mesh, "WEST" ),
                        _rhs=l( 0_c ),
                        _element=u,
                        _expr=cst( 0.0 ),
                        _type=strategy.type );
    a.row( 0_c ) += on( _range=markedfaces( mesh, "EAST" ),
                        _rhs=l( 0_c ),
                        _element=u,
                        _expr=cst( eastValue ),
                        _type=strategy.type );
}

template<typename MeshPtrType, typename WElementType, typename TElementType, typename BlockFormType, typename BlockLinearFormType>
void assembleStaticCondensationSystem( MeshPtrType const& mesh,
                                       WElementType const& W,
                                       TElementType const& T,
                                       BlockFormType& a,
                                       BlockLinearFormType& l )
{
    auto u = W( 0_c );
    auto alpha = W( 1_c );
    auto v = T( 0_c );
    auto beta = T( 1_c );

    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
    a( 0_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=idt( alpha ) * id( v ) );
    a( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=idt( u ) * id( beta ) );
    a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=cst( 2.0 ) * idt( alpha ) * id( beta ) );

    l( 0_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 1.0 ) * id( v ) );
}

template<typename MeshPtrType, typename BlockFormType, typename BlockLinearFormType, typename ElementType>
void applyBoundaryZeroOnStaticCondensationBlockform( BlockFormType& a,
                                                     BlockLinearFormType& l,
                                                     MeshPtrType const& mesh,
                                                     ElementType const& u,
                                                     std::string const& type = "elimination" )
{
    a.row( 0_c ) += on( _range=boundaryfaces( mesh ),
                        _rhs=l( 0_c ),
                        _element=u,
                        _expr=cst( 0.0 ),
                        _type=type );
}
} // namespace

inline
po::options_description makeOptions()
{
    po::options_description options( "Test space product Options" );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_product" ,
                     "test_product" ,
                     "0.2",
                     "test product of spaces",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );

    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

 FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
 BOOST_AUTO_TEST_SUITE( productspace_suite )


 BOOST_AUTO_TEST_CASE( test1 )
 {
     using namespace Feel;
     using Feel::cout;

     auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
     auto Xh=Pch<1>(mesh);
     auto Wh=Pchv<2>(mesh);
     auto sp = product(Xh,Wh);
     BOOST_CHECK_EQUAL( sp.numberOfSpaces(), 2 );
     cout << tc::red << "number of spaces " << tc::reset << sp.numberOfSpaces() << std::endl;
     BOOST_CHECK_EQUAL( sp.nDof(), Xh->nDof()+Wh->nDof() );
     cout << tc::red << "total number of dof " << tc::reset << sp.nDof() << std::endl;
     BOOST_CHECK_EQUAL( sp.nLocalDof(), Xh->nLocalDof()+Wh->nLocalDof() );
     cout << tc::red << "local number of dof " << tc::reset << sp.nLocalDof() << std::endl;
     BOOST_CHECK_EQUAL( sp[0_c], Xh );
     BOOST_CHECK_EQUAL( sp[1_c], Wh );

     auto cp = hana::cartesian_product(hana::make_tuple(sp.tupleSpaces(),sp.tupleSpaces()));

     BOOST_CHECK_EQUAL( cp[0_c][0_c], Xh );
     BOOST_CHECK_EQUAL( cp[0_c][1_c], Xh );
     BOOST_CHECK_EQUAL( cp[1_c][0_c], Xh );
     BOOST_CHECK_EQUAL( cp[1_c][1_c], Wh );
     BOOST_CHECK_EQUAL( cp[2_c][0_c], Wh );
     BOOST_CHECK_EQUAL( cp[2_c][1_c], Xh );
     BOOST_CHECK_EQUAL( cp[3_c][0_c], Wh );
     BOOST_CHECK_EQUAL( cp[3_c][1_c], Wh );

     auto ps = product(Xh,Wh);
     auto u = Xh->element();
     auto v = Wh->element();
     auto bbf = blockform2( ps );
     bbf( 0_c, 0_c ) = integrate( _range=elements(mesh), _expr=id(u)*idt(u) );
     bbf( 0_c, 1_c ) += integrate( _range=elements(mesh), _expr=id(u)*(trans(idt(v))*one() ) ) ;
     bbf( 1_c, 0_c ) += integrate( _range=elements(mesh), _expr=1./2*(trans(id(v))*one())*idt(u) );
     bbf( 1_c, 1_c ) += integrate( _range=elements(mesh), _expr=trans(id(v))*idt(v) );
     bbf.close();

     auto blf = blockform1( ps );
     blf( 0_c) = integrate( _range=elements(mesh), _expr=id(u) );
     blf( 1_c) += integrate( _range=elements(mesh), _expr=trans(id(v))*one() );
     blf.close();

     auto U=ps.element();
     bbf.solve( _solution=U, _rhs=blf );
     auto ex = exporter(_mesh=mesh);
     ex->add("u",U(0_c));
     ex->add("v",U(1_c));
     ex->save();
 }

BOOST_AUTO_TEST_CASE( test_productspace_dof_counts )
{
    using namespace Feel;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Vh0 = Pdh<0>( mesh );
    auto Uh1 = Pch<1>( mesh );

    auto staticProduct = product( Vh0, Uh1 );
    BOOST_CHECK_EQUAL( staticProduct.nDof(), Vh0->nDof() + Uh1->nDof() );
    BOOST_CHECK_EQUAL( staticProduct.nLocalDof(), Vh0->nLocalDof() + Uh1->nLocalDof() );

    std::vector<decltype( Uh1 )> dynamicEntries = { Uh1, Uh1 };
    ProductSpace<decltype( Uh1 ), false> dynamicProduct( dynamicEntries );
    BOOST_CHECK_EQUAL( dynamicProduct.numberOfSpaces(), 2 );
    BOOST_CHECK_EQUAL( dynamicProduct.nDof(), 2*Uh1->nDof() );
    BOOST_CHECK_EQUAL( dynamicProduct.nLocalDof(), 2*Uh1->nLocalDof() );

    auto repeatedUh1 = std::make_shared<ProductSpace<decltype( Uh1 ), true>>( 2, mesh );
    auto mixedProduct = product2( repeatedUh1, Vh0 );
    BOOST_CHECK_EQUAL( mixedProduct.numberOfSpaces(), 3 );
    BOOST_CHECK_EQUAL( mixedProduct.nDof(), Vh0->nDof() + 2*Uh1->nDof() );
    BOOST_CHECK_EQUAL( mixedProduct.nLocalDof(), Vh0->nLocalDof() + 2*Uh1->nLocalDof() );
}

BOOST_AUTO_TEST_CASE( test_blockform2_square_api_compatibility )
{
    using namespace Feel;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto Yh = Pdh<0>( mesh );
    auto ps = product( Xh, Yh );

    backend( _rebuild=true );

    auto checkBlockMap = [&ps]( auto const& a )
    {
        BOOST_CHECK_EQUAL( a.matrixPtr()->mapRow().nDof(), ps.nDof() );
        BOOST_CHECK_EQUAL( a.matrixPtr()->mapCol().nDof(), ps.nDof() );
    };

    auto aDefault = blockform2( ps );
    checkBlockMap( aDefault );

    auto aBackend = blockform2( ps, backend() );
    checkBlockMap( aBackend );

    auto aMonolithic = blockform2( ps, solve::strategy::monolithic, backend(), Pattern::COUPLED );
    checkBlockMap( aMonolithic );
    BOOST_CHECK( aMonolithic.matrix().monolithic() );

    std::vector<size_type> patterns = {
        Pattern::COUPLED, Pattern::ZERO,
        Pattern::ZERO, Pattern::COUPLED
    };
    auto aPatterns = blockform2( ps, solve::strategy::monolithic, backend(), patterns );
    checkBlockMap( aPatterns );
    BOOST_CHECK( aPatterns.matrix().monolithic() );

    auto const& cps = ps;
    auto aConst = blockform2( cps, solve::strategy::monolithic, backend(), Pattern::COUPLED );
    checkBlockMap( aConst );
}

BOOST_AUTO_TEST_CASE( test_csrgraphblocks_rectangular_productspaces )
{
    using namespace Feel;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Vh0 = Pdh<0>( mesh );
    auto Vh1 = Pdh<1>( mesh );
    auto Uh1 = Pch<1>( mesh );
    auto Uh2 = Pch<2>( mesh );

    backend( _rebuild=true );

    auto test = product( Vh0, Vh1 );
    auto trial = product( Uh1, Uh2 );

    std::vector<size_type> patterns = {
        Pattern::COUPLED, Pattern::ZERO,
        Pattern::ZERO, Pattern::COUPLED
    };
    auto graph = csrGraphBlocks( test, trial, patterns );
    MatrixCondensed<double> matrix( solve::strategy::monolithic, graph, backend(), false );

    BOOST_CHECK_EQUAL( matrix.mapRow().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( matrix.mapCol().nDof(), trial.nDof() );

    auto repeatedTrial = std::make_shared<ProductSpace<decltype( Uh1 ), true>>( 2, mesh );
    auto nestedTrial = product2( repeatedTrial, Vh0 );
    auto nestedGraph = csrGraphBlocks( test, nestedTrial, Pattern::COUPLED );
    MatrixCondensed<double> nestedMatrix( solve::strategy::monolithic, nestedGraph, backend(), false );

    BOOST_CHECK_EQUAL( nestedMatrix.mapRow().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( nestedMatrix.mapCol().nDof(), Vh0->nDof() + 2*Uh1->nDof() );
}

BOOST_AUTO_TEST_CASE( test_blockform2_rectangular_core_allocation )
{
    using namespace Feel;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Vh0 = Pdh<0>( mesh );
    auto Vh1 = Pdh<1>( mesh );
    auto Uh1 = Pch<1>( mesh );
    auto Uh2 = Pch<2>( mesh );

    backend( _rebuild=true );

    auto test = product( Vh0, Vh1 );
    auto trial = product( Uh1, Uh2 );

    BlockBilinearForm<decltype( test ), decltype( trial )> a( test, trial, solve::strategy::monolithic, backend(), Pattern::COUPLED );

    BOOST_CHECK( a.isRectangular() );
    BOOST_CHECK_EQUAL( a.testFunctionSpace().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( a.trialFunctionSpace().nDof(), trial.nDof() );
    BOOST_CHECK_EQUAL( a.matrixPtr()->mapRow().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( a.matrixPtr()->mapCol().nDof(), trial.nDof() );
    BOOST_CHECK( a.matrix().monolithic() );

    a.allocateMatrix( solve::strategy::monolithic, backend() );
    BOOST_CHECK_EQUAL( a.matrixPtr()->mapRow().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( a.matrixPtr()->mapCol().nDof(), trial.nDof() );
    BOOST_CHECK( a.matrix().monolithic() );
}

BOOST_AUTO_TEST_CASE( test_blockform2_named_argument_api )
{
    using namespace Feel;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto Yh = Pdh<0>( mesh );
    auto Vh0 = Pdh<0>( mesh );
    auto Vh1 = Pdh<1>( mesh );
    auto Uh1 = Pch<1>( mesh );
    auto Uh2 = Pch<2>( mesh );

    backend( _rebuild=true );

    auto squarePs = product( Xh, Yh );
    auto squareNamed = blockform2( _test=squarePs );

    BOOST_CHECK( !squareNamed.isRectangular() );
    BOOST_CHECK_EQUAL( squareNamed.testFunctionSpace().nDof(), squarePs.nDof() );
    BOOST_CHECK_EQUAL( squareNamed.trialFunctionSpace().nDof(), squarePs.nDof() );
    BOOST_CHECK_EQUAL( squareNamed.matrixPtr()->mapRow().nDof(), squarePs.nDof() );
    BOOST_CHECK_EQUAL( squareNamed.matrixPtr()->mapCol().nDof(), squarePs.nDof() );

    auto squareShapeTrial = product( Yh, Xh );
    auto squareShapeNamed = blockform2( _test=squarePs,
                                        _trial=squareShapeTrial,
                                        _strategy=solve::strategy::monolithic,
                                        _backend=backend(),
                                        _pattern=Pattern::COUPLED );

    BOOST_CHECK( !squareShapeNamed.isRectangular() );
    BOOST_CHECK_EQUAL( squareShapeNamed.testFunctionSpace().nDof(), squarePs.nDof() );
    BOOST_CHECK_EQUAL( squareShapeNamed.trialFunctionSpace().nDof(), squareShapeTrial.nDof() );
    BOOST_CHECK_EQUAL( squareShapeNamed.matrixPtr()->mapRow().nDof(), squarePs.nDof() );
    BOOST_CHECK_EQUAL( squareShapeNamed.matrixPtr()->mapCol().nDof(), squareShapeTrial.nDof() );

    auto test = product( Vh0, Vh1 );
    auto trial = product( Uh1, Uh2 );

    auto rectNamed = blockform2( _test=test,
                                 _trial=trial,
                                 _backend=backend(),
                                 _pattern=Pattern::COUPLED );

    BOOST_CHECK( rectNamed.isRectangular() );
    BOOST_CHECK_EQUAL( rectNamed.testFunctionSpace().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( rectNamed.trialFunctionSpace().nDof(), trial.nDof() );
    BOOST_CHECK_EQUAL( rectNamed.matrixPtr()->mapRow().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( rectNamed.matrixPtr()->mapCol().nDof(), trial.nDof() );
    BOOST_CHECK( rectNamed.matrix().monolithic() );

    auto rectStrategy = blockform2( _test=test,
                                    _trial=trial,
                                    _strategy=solve::strategy::monolithic,
                                    _backend=backend(),
                                    _pattern=Pattern::COUPLED );

    BOOST_CHECK_EQUAL( rectStrategy.matrixPtr()->mapRow().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( rectStrategy.matrixPtr()->mapCol().nDof(), trial.nDof() );
    BOOST_CHECK( rectStrategy.matrix().monolithic() );

    std::vector<size_type> patterns = {
        Pattern::COUPLED, Pattern::ZERO,
        Pattern::ZERO, Pattern::COUPLED
    };
    auto rectPatterns = blockform2( _test=test,
                                    _trial=trial,
                                    _strategy=solve::strategy::monolithic,
                                    _backend=backend(),
                                    _pattern=patterns );

    BOOST_CHECK_EQUAL( rectPatterns.matrixPtr()->mapRow().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( rectPatterns.matrixPtr()->mapCol().nDof(), trial.nDof() );
    BOOST_CHECK( rectPatterns.matrix().monolithic() );
}

BOOST_AUTO_TEST_CASE( test_blockform2_rectangular_monolithic_assembly )
{
    using namespace Feel;
    using namespace boost::hana::literals;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Vh0 = Pdh<0>( mesh );
    auto Vh1 = Pdh<1>( mesh );
    auto Uh1 = Pch<1>( mesh );
    auto Uh2 = Pch<2>( mesh );

    backend( _rebuild=true );

    auto test = product( Vh0, Vh1 );
    auto trial = product( Uh1, Uh2 );
    auto A = blockform2( _test=test,
                         _trial=trial,
                         _strategy=solve::strategy::monolithic,
                         _backend=backend(),
                         _pattern=Pattern::COUPLED );

    auto v0 = Vh0->element();
    auto v1 = Vh1->element();
    auto u1 = Uh1->element();
    auto u2 = Uh2->element();

    A( 0_c, 0_c ) += integrate( _range=elements( mesh ), _expr=id( v0 )*idt( u1 ) );
    A( 0_c, 1_c ) += integrate( _range=elements( mesh ), _expr=id( v0 )*idt( u2 ) );
    A( 1_c, 0_c ) += integrate( _range=elements( mesh ), _expr=id( v1 )*idt( u1 ) );
    A( 1_c, 1_c ) += integrate( _range=elements( mesh ), _expr=id( v1 )*idt( u2 ) );
    A.close();

    BOOST_CHECK_EQUAL( A.matrixPtr()->mapRow().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( A.matrixPtr()->mapCol().nDof(), trial.nDof() );
    BOOST_CHECK( A.matrix().monolithic() );
    BOOST_CHECK_GT( A.nnz(), 0 );
}

BOOST_AUTO_TEST_CASE( test_blockform2_box_scheme_rt_cr_trial_p0_test_quarter_turn_3d )
{
    using namespace Feel;
    using namespace boost::hana::literals;

    using mesh_type = Mesh<Simplex<3>>;
    using cr_space_type = FunctionSpace<mesh_type, bases<CrouzeixRaviart<1>>>;

    auto mesh = loadMesh( _mesh=new mesh_type,
                          _filename=Environment::expand( "${top_srcdir}/feelpp/quickstart/laplacian/cases/quarter-turn/quarter-turn3D.geo" ) );
    auto XhRT = Dh<0>( mesh );
    auto XhCR = cr_space_type::New( _mesh=mesh,
                                    _worldscomm=makeWorldsComm( 1, mesh->worldComm() ) );
    auto VhP0v = Pdhv<0>( mesh );
    auto VhP0s = Pdh<0>( mesh );

    backend( _rebuild=true );

    auto trial = product( XhRT, XhCR );
    auto test = product( VhP0v, VhP0s );
    auto A = blockform2( _test=test,
                         _trial=trial,
                         _strategy=solve::strategy::monolithic,
                         _backend=backend(),
                         _pattern=Pattern::COUPLED );

    auto sigma = XhRT->element( "sigma" );
    auto phi = XhCR->element( "phi" );
    auto v = VhP0v->element( "v" );
    auto q = VhP0s->element( "q" );

    A( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( idt( sigma ), id( v ) ) );
    A( 0_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=gradt( phi )*id( v ) );
    A( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=divt( sigma )*id( q ) );
    A( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=idt( phi )*id( q ) );
    A.close();

    BOOST_CHECK( A.isRectangular() );
    BOOST_CHECK_EQUAL( A.testFunctionSpace().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( A.trialFunctionSpace().nDof(), trial.nDof() );
    BOOST_CHECK_EQUAL( A.matrixPtr()->mapRow().nDof(), test.nDof() );
    BOOST_CHECK_EQUAL( A.matrixPtr()->mapCol().nDof(), trial.nDof() );
    BOOST_CHECK_GT( A.nnz(), 0 );
    BOOST_CHECK_GT( A.l1Norm(), 0.0 );
    BOOST_CHECK_GT( A.linftyNorm(), 0.0 );
}

BOOST_AUTO_TEST_CASE( test_blockform2_row_uses_flattened_test_space_index )
{
    using namespace Feel;
    using namespace boost::hana::literals;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto repeatedXh = std::make_shared<ProductSpace<decltype( Xh ), true>>( 2, mesh );
    auto ps = product2( repeatedXh, Xh );

    backend( _rebuild=true );

    auto A = blockform2( ps, solve::strategy::monolithic, backend(), Pattern::COUPLED );

    auto staticRow = A.row( 0_c );
    auto const& staticActual = staticRow.dofIdToContainerIdTest();
    auto const& staticExpected = A.matrix().mapRowPtr()->dofIdToContainerId( 0 );
    BOOST_REQUIRE_EQUAL( staticActual.size(), staticExpected.size() );
    BOOST_CHECK_EQUAL_COLLECTIONS( staticActual.begin(), staticActual.end(),
                                   staticExpected.begin(), staticExpected.end() );

    auto nestedFirstRow = A.row( 1_c );
    auto const& nestedFirstActual = nestedFirstRow.dofIdToContainerIdTest();
    auto const& nestedFirstExpected = A.matrix().mapRowPtr()->dofIdToContainerId( 1 );
    BOOST_REQUIRE_EQUAL( nestedFirstActual.size(), nestedFirstExpected.size() );
    BOOST_CHECK_EQUAL_COLLECTIONS( nestedFirstActual.begin(), nestedFirstActual.end(),
                                   nestedFirstExpected.begin(), nestedFirstExpected.end() );

    auto nestedRow = A.row( 1_c, 1 );
    auto const& nestedActual = nestedRow.dofIdToContainerIdTest();
    auto const& nestedExpected = A.matrix().mapRowPtr()->dofIdToContainerId( 2 );
    BOOST_REQUIRE_EQUAL( nestedActual.size(), nestedExpected.size() );
    BOOST_CHECK_EQUAL_COLLECTIONS( nestedActual.begin(), nestedActual.end(),
                                   nestedExpected.begin(), nestedExpected.end() );
}

BOOST_AUTO_TEST_CASE( test_blockform2_rectangular_row_dirichlet_worker )
{
    if ( !std::getenv( "FEELPP_BLOCKFORM2_RECTANGULAR_ROW_DIRICHLET_WORKER" ) )
        return;

    using namespace Feel;
    using namespace boost::hana::literals;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Vh0 = Pdh<0>( mesh );
    auto Vh1 = Pdh<1>( mesh );
    auto Uh1 = Pch<1>( mesh );
    auto Uh2 = Pch<2>( mesh );

    backend( _rebuild=true );

    auto test = product( Vh0, Vh1 );
    auto trial = product( Uh1, Uh2 );
    auto A = blockform2( _test=test,
                         _trial=trial,
                         _strategy=solve::strategy::monolithic,
                         _backend=backend(),
                         _pattern=Pattern::COUPLED );

    (void)A.row( 0_c );
}

BOOST_AUTO_TEST_CASE( test_blockform2_rectangular_row_dirichlet_rejected )
{
#if defined( __unix__ )
    if ( Environment::worldComm().globalSize() != 1 )
    {
        BOOST_TEST_MESSAGE( "rectangular row Dirichlet rejection death test is exercised only in sequential runs" );
        BOOST_CHECK( true );
        return;
    }

    std::string const executable = boost::unit_test::framework::master_test_suite().argv[0];
    pid_t pid = fork();
    BOOST_REQUIRE_NE( pid, -1 );

    if ( pid == 0 )
    {
        setenv( "FEELPP_BLOCKFORM2_RECTANGULAR_ROW_DIRICHLET_WORKER", "1", 1 );
        execl( executable.c_str(), executable.c_str(),
               "--run_test=productspace_suite/test_blockform2_rectangular_row_dirichlet_worker",
               "--log_level=test_suite",
               "--",
               "--directory=testsuite/test_productspaces_rectangular_row_dirichlet_worker",
               static_cast<char*>( nullptr ) );
        _exit( 127 );
    }

    int status = 0;
    BOOST_REQUIRE_EQUAL( waitpid( pid, &status, 0 ), pid );
    BOOST_CHECK_MESSAGE( WIFSIGNALED( status ) || ( WIFEXITED( status ) && WEXITSTATUS( status ) != 0 ),
                         "rectangular row Dirichlet worker exited successfully; expected same-space CHECK failure" );
#else
    BOOST_TEST_MESSAGE( "rectangular row Dirichlet rejection death test requires fork/exec support" );
#endif
}

BOOST_AUTO_TEST_CASE( test_blockform2_unsupported_solve_policy_worker )
{
    char const* mode = std::getenv( "FEELPP_BLOCKFORM2_UNSUPPORTED_SOLVE_POLICY_WORKER" );
    if ( !mode )
        return;

    BOOST_CHECK( true );

    using namespace Feel;
    using namespace boost::hana::literals;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Vh0 = Pdh<0>( mesh );
    auto Vh1 = Pdh<1>( mesh );
    auto Uh1 = Pch<1>( mesh );
    auto Uh2 = Pch<2>( mesh );

    backend( _rebuild=true );

    std::string const modeValue = mode;
    if ( modeValue == "solve_rectangular" )
    {
        auto test = product( Vh0, Vh1 );
        auto trial = product( Uh1, Uh2 );
        auto A = blockform2( _test=test,
                             _trial=trial,
                             _strategy=solve::strategy::monolithic,
                             _backend=backend(),
                             _pattern=Pattern::COUPLED );
        auto rhs = blockform1( test, solve::strategy::monolithic, backend() );
        auto solution = trial.element();
        A.solve( _solution=solution, _rhs=rhs );
        return;
    }

    if ( modeValue == "solve_condense_distinct_square" ||
         modeValue == "solve_local_distinct_square" )
    {
        auto Xh = Pch<1>( mesh );
        auto Yh = Pdh<0>( mesh );
        auto test = product( Xh, Yh );
        auto trial = product( Yh, Xh );
        auto A = blockform2( _test=test,
                             _trial=trial,
                             _strategy=solve::strategy::monolithic,
                             _backend=backend(),
                             _pattern=Pattern::COUPLED );
        auto rhs = blockform1( test, solve::strategy::monolithic, backend() );
        auto solution = trial.element();
        if ( modeValue == "solve_condense_distinct_square" )
            A.solve( _solution=solution, _rhs=rhs, _condense=true );
        else
            A.solve( _solution=solution, _rhs=rhs, _local=true );
        return;
    }

    if ( modeValue == "static_condensation_strategy_rectangular" )
    {
        auto test = product( Vh0, Vh1 );
        auto trial = product( Uh1, Uh2 );
        auto A = blockform2( _test=test,
                             _trial=trial,
                             _strategy=solve::strategy::static_condensation,
                             _backend=backend(),
                             _pattern=Pattern::COUPLED );
        (void)A;
        return;
    }

    BOOST_FAIL( "unknown unsupported solve policy worker mode: " << modeValue );
}

BOOST_AUTO_TEST_CASE( test_blockform2_unsupported_solve_paths_rejected )
{
#if defined( __unix__ )
    if ( Environment::worldComm().globalSize() != 1 )
    {
        BOOST_TEST_MESSAGE( "unsupported solve policy death tests are exercised only in sequential runs" );
        BOOST_CHECK( true );
        return;
    }

    std::array<std::string, 4> const modes = {
        "solve_rectangular",
        "solve_condense_distinct_square",
        "solve_local_distinct_square",
        "static_condensation_strategy_rectangular"
    };

    std::string const executable = boost::unit_test::framework::master_test_suite().argv[0];
    for ( auto const& mode : modes )
    {
        pid_t pid = fork();
        BOOST_REQUIRE_NE( pid, -1 );

        if ( pid == 0 )
        {
            setenv( "FEELPP_BLOCKFORM2_UNSUPPORTED_SOLVE_POLICY_WORKER", mode.c_str(), 1 );
            execl( executable.c_str(), executable.c_str(),
                   "--run_test=productspace_suite/test_blockform2_unsupported_solve_policy_worker",
                   "--log_level=test_suite",
                   "--",
                   "--directory=testsuite/test_productspaces_unsupported_solve_policy_worker",
                   static_cast<char*>( nullptr ) );
            _exit( 127 );
        }

        int status = 0;
        BOOST_REQUIRE_EQUAL( waitpid( pid, &status, 0 ), pid );
        BOOST_CHECK_MESSAGE( WIFSIGNALED( status ) || ( WIFEXITED( status ) && WEXITSTATUS( status ) != 0 ),
                             "unsupported solve policy worker mode " << mode << " exited successfully; expected CHECK failure" );
    }
#else
    BOOST_TEST_MESSAGE( "unsupported solve policy death tests require fork/exec support" );
#endif
}

BOOST_AUTO_TEST_CASE( test3 )
{
    using namespace Feel;
    using Feel::cout;
    using namespace boost::hana::literals;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);

    backend(_rebuild=true);
    int n = int(doption("parameters.n"));
    if ( n <= 0 || n >= 10 ) return;
    ProductSpace<Pch_ptrtype<Mesh<Simplex<2>>,1>, true> ps( n, mesh );
    BOOST_CHECK_EQUAL( ps.numberOfSpaces(), n );
    cout << tc::red << "number of spaces " << tc::reset << ps.numberOfSpaces() << std::endl;


    auto U = ps.element();
    auto u = U[0];
    auto b = blockform2( ps );
    auto l = blockform1( ps );
    std::vector<std::string> alphabet { "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega" };

    for( int i = 0; i < n; ++i )
    {
        b(i,i) += integrate( _range=elements(mesh), _expr=idt(u)*id(u));
        l(i) += integrate( _range=elements(mesh), _expr=expr(soption("functions."+alphabet[i]))*id(u));
    }
    b.solve( _rhs=l, _solution=U );
    auto ex = exporter(_mesh=mesh);
    for(int i = 0; i < n; ++i )
    {
        ex->add(alphabet[i],U[i]);
        U[i].printMatlab(alphabet[i]+".m");
    }
    ex->save();
}

BOOST_AUTO_TEST_CASE( test4 )
{
    using namespace Feel;
    using Feel::cout;
    using namespace boost::hana::literals;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);

    backend(_rebuild=true);
    int n = int(doption("parameters.n"));
    auto Xh = Pch<1>(mesh);
    auto Zh = Pch<1>(mesh);
    auto Yh = Pchv<3>(mesh);
    auto ps = std::make_shared<ProductSpace<decltype(Pch<2>(mesh)), true>>( n, mesh );
    auto p = product2( ps, Xh, Yh, Zh );

    BOOST_CHECK_EQUAL( p.numberOfSpaces(), n+3 );
    cout << tc::red << "number of spaces " << tc::reset << p.numberOfSpaces() << std::endl;

    auto U = p.element();

    std::vector<std::string> alphabet { "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega" };

    auto u = U(0_c);
    auto w = U(1_c);
    auto z = U(2_c);
    auto l = blockform1( p );
    auto b = blockform2( p );


    l(0_c) = integrate( _range=elements(mesh), _expr=expr(soption("functions."+alphabet[0]))*id(u));
    b(0_c,0_c) += integrate( _range=elements(mesh), _expr=idt(u)*id(u));
    l(1_c) = integrate( _range=elements(mesh), _expr=expr(soption("functions."+alphabet[1]))*trans(id(w))*one());
    b(1_c,1_c) += integrate( _range=elements(mesh), _expr=trans(idt(w))*id(w));
    l(2_c) = integrate( _range=elements(mesh), _expr=expr(soption("functions."+alphabet[0]))*id(z));
    b(2_c,2_c) += integrate( _range=elements(mesh), _expr=idt(z)*id(z));
    for( int i = 0; i < n; ++i )
    {
        auto v = U(3_c,i);

        b(3_c,3_c,i,i) += integrate( _range=elements(mesh), _expr=idt(v)*id(v));
        l(3_c,i) += integrate( _range=elements(mesh), _expr=expr(soption("functions."+alphabet[2+i]))*id(v));
    }
    b.solve( _rhs=l, _solution=U );
    auto ex = exporter(_mesh=mesh);
    ex->add(alphabet[0],U(0_c));
    ex->add(alphabet[1],U(1_c));
    ex->add(alphabet[2],U(2_c));
    for(int i = 0; i < ps->numberOfSpaces(); ++i )
    {
        ex->add(alphabet[i+2],U(3_c,i));
        //U[i].printMatlab(alphabet[i]+".m");
    }
    ex->save();

}

BOOST_AUTO_TEST_CASE( test_row_dirichlet_static_condensation )
{
    using namespace Feel;
    using namespace vf;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Uh = Pch<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto ps = product( Uh, Ah );

    auto W = ps.element();
    auto T = ps.element();
    auto u = W( 0_c );
    auto alpha = W( 1_c );
    auto v = T( 0_c );
    auto beta = T( 1_c );

    auto assembleMixedSystem = [&]( auto& a, auto& l )
    {
        a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _expr=inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
        a( 0_c, 1_c ) += integrate( _range=elements( mesh ),
                                    _expr=idt( alpha ) * id( v ) );
        a( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _expr=idt( u ) * id( beta ) );
        a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                    _expr=cst( 2.0 ) * idt( alpha ) * id( beta ) );

        l( 0_c ) += integrate( _range=elements( mesh ),
                               _expr=cst( 1.0 ) * id( v ) );

        l.close();
        a.close();

        a.row( 0_c ) += on( _range=boundaryfaces( mesh ),
                            _rhs=l( 0_c ),
                            _element=u,
                            _expr=cst( 0.0 ),
                            _type="elimination" );
    };

    backend( _rebuild=true );

    auto aMonolithic = blockform2( ps, solve::strategy::monolithic, backend() );
    auto lMonolithic = blockform1( ps, solve::strategy::monolithic, backend() );
    assembleMixedSystem( aMonolithic, lMonolithic );

    auto aCondensed = blockform2( ps, solve::strategy::static_condensation, backend() );
    auto lCondensed = blockform1( ps, solve::strategy::static_condensation, backend() );
    assembleMixedSystem( aCondensed, lCondensed );

    auto UMonolithic = ps.element();
    auto UCondensed = ps.element();

    aMonolithic.solve( _solution=UMonolithic, _rhs=lMonolithic );
    aCondensed.solve( _solution=UCondensed, _rhs=lCondensed,
                      _condense=true, _condenser=condenser_sb9() );

    double const uError = normL2( _range=elements( mesh ),
                                  _expr=idv( UMonolithic( 0_c ) ) - idv( UCondensed( 0_c ) ) );
    double const alphaError = normL2( _range=elements( mesh ),
                                      _expr=idv( UMonolithic( 1_c ) ) - idv( UCondensed( 1_c ) ) );
    double const boundaryError = normL2( _range=boundaryfaces( mesh ),
                                         _expr=idv( UCondensed( 0_c ) ) );

    // With the MPI gasm/umfpack path used by ctest, monolithic and condensed
    // solves agree to about 1e-10 rather than bitwise-identical precision.
    BOOST_CHECK_SMALL( uError, 2e-10 );
    BOOST_CHECK_SMALL( alphaError, 2e-10 );
    BOOST_CHECK_SMALL( boundaryError, 1e-12 );
}

BOOST_DATA_TEST_CASE( test_row_dirichlet_repeated_blockform,
                      bdata::make( monolithicDirichletStrategies ),
                      strategy )
{
    using namespace Feel;
    using namespace vf;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto Yh = Pch<1>( mesh );
    auto ps = product( Xh, Yh );

    auto W = ps.element();
    auto T = ps.element();
    auto u = W( 0_c );
    auto p = W( 1_c );
    auto v = T( 0_c );
    auto q = T( 1_c );

    backend( _rebuild=true );

    auto a = blockform2( ps, solve::strategy::monolithic, backend() );
    auto l = blockform1( ps, solve::strategy::monolithic, backend() );
    a.deferDirichlet();

    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
    a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( p ), grad( q ) ) + idt( p ) * id( q ) );

    l( 0_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 1.0 ) * id( v ) );
    l( 1_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 2.0 ) * id( q ) );

    l.close();
    a.close();

    auto applyDirichlet = [&]( std::string const& marker, double value )
    {
        a.row( 0_c ) += on( _range=markedfaces( mesh, marker ),
                            _rhs=l( 0_c ),
                            _element=u,
                            _expr=cst( value ),
                            _type=strategy.type );
    };

    applyDirichlet( "WEST", 0.0 );
    applyDirichlet( "EAST", 1.0 );
    applyDirichlet( "WEST", 0.0 );
    applyDirichlet( "EAST", 1.0 );

    BOOST_CHECK( a.hasPendingDirichletConstraints() );
    BOOST_CHECK( a.hasDirichletConstraints() );

    auto solution = ps.element();
    a.solve( _solution=solution, _rhs=l );
    BOOST_CHECK( !a.hasPendingDirichletConstraints() );
    BOOST_CHECK( a.hasDirichletConstraints() );

    double const westError = normL2( _range=markedfaces( mesh, "WEST" ),
                                     _expr=idv( solution( 0_c ) ) );
    double const eastError = normL2( _range=markedfaces( mesh, "EAST" ),
                                     _expr=idv( solution( 0_c ) ) - cst( 1.0 ) );
    double const secondFieldNorm = normL2( _range=elements( mesh ),
                                           _expr=idv( solution( 1_c ) ) );

    BOOST_CHECK_SMALL( westError, strategy.tolerance );
    BOOST_CHECK_SMALL( eastError, strategy.tolerance );
    BOOST_CHECK( std::isfinite( secondFieldNorm ) );
}

BOOST_DATA_TEST_CASE( test_row_dirichlet_manual_apply_blockform,
                      bdata::make( monolithicDirichletStrategies ),
                      strategy )
{
    using namespace Feel;
    using namespace vf;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto Yh = Pch<1>( mesh );
    auto ps = product( Xh, Yh );

    auto W = ps.element();
    auto T = ps.element();
    auto u = W( 0_c );
    auto p = W( 1_c );
    auto v = T( 0_c );
    auto q = T( 1_c );

    backend( _rebuild=true );

    auto a = blockform2( ps, solve::strategy::monolithic, backend() );
    auto l = blockform1( ps, solve::strategy::monolithic, backend() );
    a.deferDirichlet();
    a.setDirichletInPlace( false );

    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
    a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( p ), grad( q ) ) + idt( p ) * id( q ) );

    l( 0_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 1.0 ) * id( v ) );
    l( 1_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 2.0 ) * id( q ) );

    l.close();
    a.close();

    a.row( 0_c ) += on( _range=markedfaces( mesh, "WEST" ),
                        _rhs=l( 0_c ),
                        _element=u,
                        _expr=cst( 0.0 ),
                        _type=strategy.type );
    a.row( 0_c ) += on( _range=markedfaces( mesh, "EAST" ),
                        _rhs=l( 0_c ),
                        _element=u,
                        _expr=cst( 1.0 ),
                        _type=strategy.type );

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

    auto solution = ps.element();
    a.solve( _solution=solution, _rhs=l );

    double const westError = normL2( _range=markedfaces( mesh, "WEST" ),
                                     _expr=idv( solution( 0_c ) ) );
    double const eastError = normL2( _range=markedfaces( mesh, "EAST" ),
                                     _expr=idv( solution( 0_c ) ) - cst( 1.0 ) );

    BOOST_CHECK_SMALL( westError, strategy.tolerance );
    BOOST_CHECK_SMALL( eastError, strategy.tolerance );
}

BOOST_DATA_TEST_CASE( test_row_dirichlet_rhs_reuse_refresh_blockform,
                      bdata::make( monolithicDirichletStrategies ),
                      strategy )
{
    using namespace Feel;
    using namespace vf;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto Yh = Pch<1>( mesh );
    auto ps = product( Xh, Yh );

    auto W = ps.element();
    auto T = ps.element();
    auto u = W( 0_c );
    auto p = W( 1_c );
    auto v = T( 0_c );
    auto q = T( 1_c );

    backend( _rebuild=true );

    auto a = blockform2( ps, solve::strategy::monolithic, backend() );
    auto l = blockform1( ps, solve::strategy::monolithic, backend() );
    a.deferDirichlet();

    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
    a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( p ), grad( q ) ) + idt( p ) * id( q ) );

    l( 0_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 1.0 ) * id( v ) );
    l( 1_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 2.0 ) * id( q ) );

    l.close();
    a.close();

    a.row( 0_c ) += on( _range=markedfaces( mesh, "WEST" ),
                        _rhs=l( 0_c ),
                        _element=u,
                        _expr=cst( 0.0 ),
                        _type=strategy.type );
    a.row( 0_c ) += on( _range=markedfaces( mesh, "EAST" ),
                        _rhs=l( 0_c ),
                        _element=u,
                        _expr=cst( 1.0 ),
                        _type=strategy.type );

    auto constrainedVector = a.constrainedVectorPtr( l );
    auto constrainedSnapshot = vf::cloneVectorWithValues( constrainedVector );
    auto rhsVector = l.vectorPtr()->getVector();
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

    auto solution = ps.element();
    a.solve( _solution=solution, _rhs=l );

    double const westError = normL2( _range=markedfaces( mesh, "WEST" ),
                                     _expr=idv( solution( 0_c ) ) );
    double const eastError = normL2( _range=markedfaces( mesh, "EAST" ),
                                     _expr=idv( solution( 0_c ) ) - cst( 1.0 ) );

    BOOST_CHECK_SMALL( westError, strategy.tolerance );
    BOOST_CHECK_SMALL( eastError, strategy.tolerance );
}

#if defined( FEELPP_HAS_PETSC_H )
BOOST_AUTO_TEST_CASE( test_row_dirichlet_call_count_target_blockform_monolithic )
{
    using namespace Feel;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto Yh = Pch<1>( mesh );
    auto ps = product( Xh, Yh );

    auto W = ps.element();
    auto T = ps.element();
    auto u = W( 0_c );

    backend( _rebuild=true );

    auto counters = std::make_shared<DeferredDirichletMaterializationCounters>();
    auto countingBackend = std::make_shared<CountingBackendPetsc>( counters, Environment::worldCommPtr() );
    auto a = blockform2( ps, solve::strategy::monolithic, countingBackend );
    auto l = blockform1( ps, solve::strategy::monolithic, countingBackend );

    assembleMonolithicBlockformSystem( mesh, W, T, a, l );
    l.close();
    a.close();

    BOOST_CHECK( a.useDeferredDirichlet() );
    applyBoundaryPairOnBlockform( a, l, mesh, u, eliminationStrategy, 1.0 );
    BOOST_CHECK( a.hasPendingDirichletConstraints() );
    BOOST_CHECK( a.hasDirichletConstraints() );

    auto const baseMatrix = a.baseMatrixPtr();
    BOOST_CHECK( dynamic_cast<CountingMatrixSparse<double>*>( baseMatrix.get() ) != nullptr );
    BOOST_CHECK( !a.hasMaterializedConstrainedOperator() );

    counters->reset();
    auto rhsSnapshot = vf::cloneVectorWithValues( l.vectorPtr()->getVector() );
    auto constrainedMatrix = a.activeMatrixPtr( l );
    auto constrainedVector = a.activeVectorPtr( l );
    auto constrainedMatrixAgain = a.activeMatrixPtr( l );
    auto constrainedVectorAgain = a.activeVectorPtr( l );

    BOOST_CHECK( constrainedMatrix );
    BOOST_CHECK( constrainedVector );
    BOOST_CHECK( a.hasMaterializedConstrainedOperator() );
    BOOST_CHECK( constrainedMatrix == baseMatrix );
    BOOST_CHECK( !a.hasUnconstrainedMatrix() );
    BOOST_CHECK( dynamic_cast<CountingMatrixSparse<double>*>( constrainedMatrix.get() ) != nullptr );
    BOOST_CHECK_EQUAL( counters->baseCloseCalls, 2 );
    BOOST_CHECK_EQUAL( counters->cloneCalls, 0 );
    BOOST_CHECK_EQUAL( counters->zeroRowsCalls, 1 );
    BOOST_CHECK_EQUAL( constrainedMatrix.get(), constrainedMatrixAgain.get() );
    BOOST_CHECK_EQUAL( constrainedVector.get(), constrainedVectorAgain.get() );
    auto rhsDelta = vf::cloneVectorWithValues( l.vectorPtr()->getVector() );
    rhsDelta->add( -1.0, *rhsSnapshot );
    if ( !rhsDelta->closed() )
        rhsDelta->close();
    BOOST_CHECK_SMALL( rhsDelta->linftyNorm(), 1e-14 );
}

BOOST_AUTO_TEST_CASE( test_row_dirichlet_single_apply_zero_rows_count_blockform_monolithic )
{
    if ( Environment::worldComm().globalSize() != 1 )
        return;

    using namespace Feel;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto Yh = Pch<1>( mesh );
    auto ps = product( Xh, Yh );

    auto W = ps.element();
    auto T = ps.element();

    backend( _rebuild=true );

    auto counters = std::make_shared<DeferredDirichletMaterializationCounters>();
    auto countingBackend = std::make_shared<CountingBackendPetsc>( counters, Environment::worldCommPtr() );
    auto a = blockform2( ps, solve::strategy::monolithic, countingBackend );
    auto l = blockform1( ps, solve::strategy::monolithic, countingBackend );

    assembleMonolithicBlockformSystem( mesh, W, T, a, l );
    auto countingMatrix = a.baseMatrixPtr();
    auto rhsVector = l.vectorPtr()->getVector();

    auto const& rowMap = countingMatrix->mapRow();
    BOOST_REQUIRE_GT( rowMap.nLocalDofWithoutGhost(), 0 );

    std::vector<int> localDofs;
    localDofs.push_back( static_cast<int>( rowMap.firstDofGlobalCluster() ) );
    if ( rowMap.nLocalDofWithoutGhost() > 1 )
        localDofs.push_back( static_cast<int>( rowMap.firstDofGlobalCluster() + 1 ) );
    std::vector<double> values( localDofs.size(), 0.0 );

    vf::DeferredDirichletSet<double> constraints;
    constraints.append( localDofs, values, Feel::Context( ContextOn::ELIMINATION ), 1.0 );

    countingMatrix->close();
    rhsVector->close();
    counters->reset();
    auto const merged = constraints.mergedEntries();

    BOOST_REQUIRE_EQUAL( merged.size(), 1 );
    vf::applyDeferredDirichletEntries( merged, countingMatrix, rhsVector );
    BOOST_CHECK_EQUAL( counters->zeroRowsCalls, 1 );
}

BOOST_AUTO_TEST_CASE( test_row_dirichlet_static_condensation_close_count_blockform )
{
    using namespace Feel;
    using namespace vf;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Uh = Pch<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto ps = product( Uh, Ah );

    auto W = ps.element();
    auto T = ps.element();
    auto u = W( 0_c );

    backend( _rebuild=true );

    auto a = blockform2( ps, solve::strategy::static_condensation, backend() );
    auto l = blockform1( ps, solve::strategy::static_condensation, backend() );

    assembleStaticCondensationSystem( mesh, W, T, a, l );
    l.close();
    a.close();

    applyBoundaryZeroOnStaticCondensationBlockform( a, l, mesh, u );

    BOOST_CHECK( a.hasPendingDirichletConstraints() );
    BOOST_CHECK( a.hasDirichletConstraints() );
    BOOST_CHECK( !a.supportsConstrainedOperatorView() );

    auto solution = ps.element();
    auto sc = a.matrix().sc();
    auto psS = product( solution( 0_c ).functionSpace() );

    auto counters = std::make_shared<DeferredDirichletMaterializationCounters>();
    auto countingBackend = std::make_shared<CountingBackendPetsc>( counters, Environment::worldCommPtr() );
    auto S = blockform2( psS, solve::strategy::monolithic, countingBackend );
    auto V = blockform1( psS, solve::strategy::monolithic, countingBackend );

    a.syncLocalMatrix();
    sc->condense( l.vectorPtr()->sc(), solution, S, V );

    auto const baseMatrix = S.baseMatrixPtr();
    BOOST_CHECK( dynamic_cast<CountingMatrixSparse<double>*>( baseMatrix.get() ) != nullptr );

    counters->reset();
    a.applyDeferredDirichlet( S, V );

    BOOST_CHECK_EQUAL( counters->baseCloseCalls, 2 );
    BOOST_CHECK( S.baseOperatorClosed() );
    BOOST_CHECK( V.vectorPtr()->getVector()->closed() );
    BOOST_CHECK( !a.hasPendingDirichletConstraints() );
    BOOST_CHECK( a.hasDirichletConstraints() );
}

BOOST_AUTO_TEST_CASE( test_row_dirichlet_static_condensation_single_apply_zero_rows_count_blockform )
{
    if ( Environment::worldComm().globalSize() != 1 )
        return;

    using namespace Feel;
    using namespace vf;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Uh = Pch<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto ps = product( Uh, Ah );

    auto W = ps.element();
    auto T = ps.element();
    auto u = W( 0_c );

    backend( _rebuild=true );

    auto a = blockform2( ps, solve::strategy::static_condensation, backend() );
    auto l = blockform1( ps, solve::strategy::static_condensation, backend() );

    assembleStaticCondensationSystem( mesh, W, T, a, l );
    l.close();
    a.close();

    applyBoundaryZeroOnStaticCondensationBlockform( a, l, mesh, u );

    auto solution = ps.element();
    auto sc = a.matrix().sc();
    auto psS = product( solution( 0_c ).functionSpace() );

    auto counters = std::make_shared<DeferredDirichletMaterializationCounters>();
    auto countingBackend = std::make_shared<CountingBackendPetsc>( counters, Environment::worldCommPtr() );
    auto S = blockform2( psS, solve::strategy::monolithic, countingBackend );
    auto V = blockform1( psS, solve::strategy::monolithic, countingBackend );

    a.syncLocalMatrix();
    sc->condense( l.vectorPtr()->sc(), solution, S, V );

    counters->reset();
    a.applyDeferredDirichlet( S, V );

    BOOST_CHECK_EQUAL( counters->baseCloseCalls, 2 );
    BOOST_CHECK_EQUAL( counters->zeroRowsCalls, 1 );
    BOOST_CHECK( S.baseOperatorClosed() );
    BOOST_CHECK( V.vectorPtr()->getVector()->closed() );
}
#endif

BOOST_AUTO_TEST_CASE( test_row_dirichlet_static_condensation_manual_apply )
{
    using namespace Feel;
    using namespace vf;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Uh = Pch<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto ps = product( Uh, Ah );

    auto W = ps.element();
    auto T = ps.element();
    auto u = W( 0_c );
    auto alpha = W( 1_c );
    auto v = T( 0_c );
    auto beta = T( 1_c );

    backend( _rebuild=true );

    auto a = blockform2( ps, solve::strategy::static_condensation, backend() );
    auto l = blockform1( ps, solve::strategy::static_condensation, backend() );

    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
    a( 0_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=idt( alpha ) * id( v ) );
    a( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=idt( u ) * id( beta ) );
    a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=cst( 2.0 ) * idt( alpha ) * id( beta ) );

    l( 0_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 1.0 ) * id( v ) );

    l.close();
    a.close();

    a.row( 0_c ) += on( _range=boundaryfaces( mesh ),
                        _rhs=l( 0_c ),
                        _element=u,
                        _expr=cst( 0.0 ),
                        _type="elimination" );

    BOOST_CHECK( a.hasPendingDirichletConstraints() );
    BOOST_CHECK( a.hasDirichletConstraints() );
    BOOST_CHECK( !a.supportsConstrainedOperatorView() );

    auto solution = ps.element();
    auto sc = a.matrix().sc();
    auto psS = product( solution( 0_c ).functionSpace() );
    auto S = blockform2( psS, solve::strategy::monolithic, backend() );
    auto V = blockform1( psS, solve::strategy::monolithic, backend() );

    a.syncLocalMatrix();
    sc->condense( l.vectorPtr()->sc(), solution, S, V );
    a.applyDeferredDirichlet( S, V );
    BOOST_CHECK( !a.hasPendingDirichletConstraints() );
    BOOST_CHECK( a.hasDirichletConstraints() );
    BOOST_CHECK( S.baseOperatorClosed() );
    BOOST_CHECK( V.vectorPtr()->getVector()->closed() );

    auto U = psS.element();
    S.solve( _solution=U, _rhs=V );
    solution( 0_c ) = U( 0_c );
    sc->localSolve( l.vectorPtr()->sc(), solution );

    double const boundaryError = normL2( _range=boundaryfaces( mesh ),
                                         _expr=idv( solution( 0_c ) ) );

    BOOST_CHECK_SMALL( boundaryError, 1e-12 );
    BOOST_CHECK( a.hasDirichletConstraints() );
}

BOOST_AUTO_TEST_SUITE_END()
