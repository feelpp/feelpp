#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_add_matrix
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( add_matrix_suite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    typedef Backend<double> backend_type;
    typedef typename std::shared_ptr<backend_type> backend_ptrtype ;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;

    typedef Mesh<Simplex<2,1> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef Mesh<Simplex<1,1,2> > face_mesh_type;
    typedef std::shared_ptr<face_mesh_type> face_mesh_ptrtype;

    auto mesh = loadMesh( _mesh=new mesh_type);
    auto face_mesh = createSubmesh( _mesh=mesh, _range=faces(mesh), _update=0 );

    auto Vh = Pdh<1>( mesh, true );
    auto Mh = Pdh<1>( face_mesh, true );

    backend_ptrtype b = backend( _rebuild=true);

    BlocksBaseGraphCSR hdg_graph(2,2);
    hdg_graph(0,0) = stencil( _test=Vh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,0) = stencil( _test=Mh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(0,1) = stencil( _test=Vh,_trial=Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,1) = stencil( _test=Mh,_trial=Mh, _diag_is_nonzero=false, _close=false)->graph();

    sparse_matrix_ptrtype A = b->newBlockMatrix(_block=hdg_graph);
    sparse_matrix_ptrtype A2 = b->newBlockMatrix(_block=hdg_graph);

    auto v = Vh->element( "v" );
    auto phat = Mh->element( "phat" );
    auto a13 = form2( _trial=Mh, _test=Vh,_matrix=A,
                      _rowstart=0, _colstart=1);
    a13 += integrate( _range=internalfaces(mesh),
                      _expr=idt(phat)*id(v) );

    A2->addMatrix(1.,A);
}

BOOST_AUTO_TEST_SUITE_END()
