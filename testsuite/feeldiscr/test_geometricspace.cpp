#define BOOST_TEST_MODULE test_geometricspace
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/geometricspace.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( geometricspace )

BOOST_AUTO_TEST_CASE( context )
{
    using _mesh_type = Mesh<Simplex<2>>;
    auto m = loadMesh( _mesh = new _mesh_type );

    using geometricspace_type = GeometricSpace<_mesh_type>;
    auto geospace = std::make_shared<geometricspace_type>(m);
    auto geospacectx = std::make_shared<typename geometricspace_type::Context>( geospace );
    std::vector<double> checkResults;
    node_type ptCoord(_mesh_type::nRealDim);

    for ( auto const& [x,y] : { std::make_pair(0.5,0.5),std::make_pair(0.45,0.65), std::make_pair(0.50001,0.50001), std::make_pair(0.1,0.9), std::make_pair(0.95,0.05) } )
    {
        ptCoord[0]=x;
        ptCoord[1]=y;
        geospacectx->add( ptCoord,false );
        checkResults.push_back( x*y );
    }
    // update geo ctx for use (necessary after adding a point in geoctx)
    geospacectx->updateForUse();


    auto Vh = Pch<2>( m );
    auto u = Vh->element();
    u.on(_range=elements(m),_expr=Px()*Py());

    auto exprUsed = id(u);
    geospacectx->updateGmcContext<std::decay_t<decltype(exprUsed)>::context>();

    std::vector<std::shared_ptr<Vector<double>>> lfsVec;
    for ( int k=0;k<geospacectx->nPoints();++k )
        lfsVec.push_back( backend()->newVector(Vh) );

    using space_type = std::decay_t<decltype(unwrap_ptr(Vh))>;
    using gmc_type = typename geometricspace_type::gmc_type;
    using gmc_ptrtype = std::shared_ptr<gmc_type>;
    using gm_type = typename gmc_type::gm_type;
    using geoelement_type = typename gmc_type::element_type;
    using fe_type = typename space_type::fe_type;

    using expr_type = std::decay_t<decltype(exprUsed)>;
    using expr_basis_t = typename expr_type::test_basis;
    static constexpr size_type expr_context = expr_type::context;

    // fe context
    using fecontext_type = typename space_type::fe_type::template Context< expr_context, fe_type, gm_type, geoelement_type>;
    using fecontext_ptrtype = std::shared_ptr<fecontext_type>;

    for ( auto& [geoCtxId,geoCtx] : *geospacectx )
    {
        auto gmc = std::get<0>(geoCtx)->gmContext();
        auto const& curCtxIdToPointIds = std::get<1>(geoCtx);
        auto const& elt = gmc->element();
        auto fepc = Vh->fe()->preCompute( Vh->fe(), gmc->xRefs() );
        auto fec = std::make_shared<fecontext_type>( Vh->fe(), gmc, fepc );
        auto mapgmc = Feel::vf::mapgmc( gmc );
        auto mapfec = Feel::vf::mapfec( fec );
        auto tExpr = exprUsed.evaluator(mapgmc, mapfec);

        // TODO check if next lines are really usefull?
        fec->update( gmc );
        tExpr.update( mapgmc, mapfec );

        //Eigen::MatrixXd M_IhLoc = Vh->fe()->localInterpolants(1);
        Eigen::MatrixXd IhLoc = Eigen::MatrixXd::Zero(expr_basis_t::nLocalDof,gmc->nPoints());
        Vh->fe()->interpolateBasisFunction( tExpr, IhLoc );


        for( auto const& ldof : Vh->dof()->localDof( elt.id() ) )
        {
            index_type index = ldof.second.index();

            for ( uint16_type q=0;q<curCtxIdToPointIds.size();++q )
            {
                lfsVec.at( curCtxIdToPointIds[q] )->set( index, IhLoc( ldof.first.localDof(), q ) ); // TODO use datamap mapping
            }
        }
    }



    for ( int k=0;k<lfsVec.size();++k )
    {
        auto & lfVec = lfsVec[k];
        lfVec->close();
        double res = inner_product( *lfVec, u );
        BOOST_CHECK_CLOSE( checkResults[k], res, 1e-8 );
    }

}

BOOST_AUTO_TEST_SUITE_END()
