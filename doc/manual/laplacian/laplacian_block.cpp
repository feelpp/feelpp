// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

#include <algorithm>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <string>

#include <feel/feelalg/graphcsr.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelmesh/meshbase.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

using laplacian_block_mesh_type = Mesh<Simplex<2>>;
using laplacian_block_mesh_ptrtype = std::shared_ptr<laplacian_block_mesh_type>;
using laplacian_block_primal_space_type = Pch_type<laplacian_block_mesh_type,2>;
using laplacian_block_primal_space_ptrtype = Pch_ptrtype<laplacian_block_mesh_type,2>;
using laplacian_block_primal_element_type = typename laplacian_block_primal_space_type::element_type;

struct BlockSolveResult
{
    std::string name;
    laplacian_block_primal_space_ptrtype space;
    laplacian_block_primal_element_type u1;
    double normL1_A = 0.0;
    double normL2_F = 0.0;
    double normL2_U = 0.0;
};

po::options_description
makeLaplacianBlockOptions()
{
    auto opts = feel_options();
    opts.add_options()
        ( "variant", po::value<std::string>()->default_value( "all" ), "run one of: v1, v2, v3, all" );
    return opts;
}

std::string
normalizeVariant( std::string variant )
{
    std::transform( variant.begin(), variant.end(), variant.begin(),
                    []( unsigned char c ) { return static_cast<char>( std::tolower( c ) ); } );
    return variant;
}

void
printBlockSolveResult( BlockSolveResult const& result )
{
    if ( Environment::worldComm().globalRank() == 0 )
        std::cout << result.name
                  << ": normL1_A " << std::setprecision( 9 ) << result.normL1_A
                  << " normL2_F " << std::setprecision( 9 ) << result.normL2_F
                  << " normL2_U " << std::setprecision( 9 ) << result.normL2_U << std::endl;
}

void
printBlockSolveComparison( laplacian_block_mesh_ptrtype const& mesh,
                           BlockSolveResult const& r1,
                           BlockSolveResult const& r2,
                           BlockSolveResult const& r3 )
{
    auto printPair = [&]( BlockSolveResult const& lhs, BlockSolveResult const& rhs )
    {
        double const diff = normL2( _range=elements( mesh ),
                                    _expr=idv( lhs.u1 ) - idv( rhs.u1 ) );
        double const ref = normL2( _range=elements( mesh ),
                                   _expr=idv( lhs.u1 ) );
        double const rel = ( ref > 1e-16 ) ? diff / ref : diff;
        std::cout << "||u1_" << lhs.name << " - u1_" << rhs.name << "||_L2 = " << diff
                  << ", relative = " << rel << std::endl;
    };

    if ( Environment::worldComm().globalRank() == 0 )
    {
        std::cout << "Comparison:" << std::endl;
        printPair( r1, r2 );
        printPair( r1, r3 );
        printPair( r2, r3 );
    }
}

std::shared_ptr<Mesh<Simplex<2> > >
createMeshLaplacianBlock()
{
    double const meshSize = doption( _name="gmsh.hsize" );
    GeoTool::Rectangle R( meshSize, "Rectangle",
                          GeoTool::Node(-1,-1),
                          GeoTool::Node( 1, 1) );
    R.setMarker(_type="line",_name="BoundaryA",_marker1=true);
    R.setMarker(_type="line",_name="BoundaryB",_marker2=true,_marker3=true,_marker4=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh(_mesh=new Mesh<Simplex<2> >,
                             _name="meshrect",
                             _hmax=meshSize );
    return mesh;
}


std::shared_ptr<Mesh<Simplex<2> > >
createSubMeshLaplacianBlock(std::shared_ptr<Mesh<Simplex<2> > > mesh)
{
    auto P0d = Pdh<0>(mesh);
    double r=0.2;
    auto proj = vf::project(_space=P0d,_range=elements(mesh),
                            _expr=vf::chi( (Px()*Px()+Py()*Py()) < r*r ) );
    mesh->updateMarker3( proj );
    auto submesh = createSubmesh( _mesh=mesh, _range=marked3elements(mesh,1) );
    //saveGMSHMesh(_mesh=submesh,_filename="mysubmesh.msh");
    return submesh;
}

}


Feel::BlockSolveResult runLaplacianBlockV1( Feel::laplacian_block_mesh_ptrtype const& mesh,
                                            Feel::laplacian_block_mesh_ptrtype const& submesh )
{
    using namespace Feel;

    auto Vh1 = Pch<2>( mesh );
    auto Vh2 = Pch<1>( submesh );
    if (Environment::worldComm().isMasterRank())
    {
        std::cout << "mesh->numGlobalElements() "<< mesh->numGlobalElements() << std::endl;
        std::cout << "submesh->numGlobalElements() "<< submesh->numGlobalElements() << std::endl;
        std::cout << "Vh1->nDof() "<<Vh1->nDof() << std::endl;
        std::cout << "Vh2->nDof() "<<Vh2->nDof() << std::endl;
    }

    CHECK( Vh2->nDof() > 0 ) << "not take into account\n";

    auto u1 = Vh1->elementPtr();
    auto u2 = Vh2->elementPtr();

    BlocksBaseGraphCSR myblockGraph(2,2);
    myblockGraph(0,0) = stencil(_test=Vh1,_trial=Vh1, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(0,1) = stencil(_test=Vh1,_trial=Vh2, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(1,0) = stencil(_test=Vh2,_trial=Vh1, _diag_is_nonzero=false, _close=false)->graph();
    auto A = backend()->newBlockMatrix(_block=myblockGraph);

    BlocksBaseVector<double> myblockVec(2);
    myblockVec(0,0) = backend()->newVector( Vh1 );
    myblockVec(1,0) = backend()->newVector( Vh2 );
    auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

    BlocksBaseVector<double> myblockVecSol(2);
    myblockVecSol(0,0) = u1;
    myblockVecSol(1,0) = u2;
    auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);

    form2( _trial=Vh1, _test=Vh1 ,_matrix=A )
        += integrate(_range=elements(mesh),
                    _expr=gradt(u1)*trans(grad(u1)) );

    form2( _trial=Vh2, _test=Vh1 ,_matrix=A,
           _rowstart=0, _colstart=1 )
        += integrate( _range=elements(submesh),
                      _expr=idt(u2)*id(u1) );

    form2( _trial=Vh1, _test=Vh2 ,_matrix=A,
           _rowstart=1, _colstart=0 )
        += integrate( _range=elements(submesh),
                      _expr=idt(u1)*id(u2) );

    form1( _test=Vh1, _vector=F )
        = integrate(_range=elements(mesh),
                    _expr=id(u1));

    form2( _trial=Vh1, _test=Vh1 ,_matrix=A )
        +=on(_range=boundaryfaces(mesh), _rhs=F, _element=*u1,
        _expr=constant(0.) );

    backend(_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=U );

    double normL1_A = A->l1Norm();
    double normL2_F = F->l2Norm();
    double normL2_U = U->l2Norm();

    myblockVecSol.localize(U);

    BlockSolveResult result{ "V1", Vh1, Vh1->element( "u1V1" ), normL1_A, normL2_F, normL2_U };
    result.u1 = *u1;
    printBlockSolveResult( result );

    auto e = exporter( _mesh=mesh, _name="exportV1" );
    e->add( "u1", result.u1 );
    e->save();

    return result;
}

Feel::BlockSolveResult runLaplacianBlockV2( Feel::laplacian_block_mesh_ptrtype const& mesh,
                                            Feel::laplacian_block_mesh_ptrtype const& submesh )
{
    using namespace Feel;

    auto Vh1 = Pch<2>( mesh );
    auto Vh2 = Pch<1>( submesh );

    auto u1 = Vh1->elementPtr();
    auto u2 = Vh2->elementPtr();

    auto A11 = backend()->newMatrix(_test=Vh1,_trial=Vh1);
    auto A12 = backend()->newMatrix(_test=Vh1,_trial=Vh2 );
    auto A21 = backend()->newMatrix(_test=Vh2,_trial=Vh1);

    auto F1 = backend()->newVector( Vh1 );
    auto F2 = backend()->newVector( Vh2 );

    form2( _trial=Vh1, _test=Vh1 ,_matrix=A11 )
        += integrate(_range=elements(mesh),
                    _expr=gradt(u1)*trans(grad(u1)) );

    form2( _trial=Vh2, _test=Vh1 ,_matrix=A12 )
        += integrate( _range=elements(submesh),
                      _expr=idt(u2)*id(u1) );

    form2( _trial=Vh1, _test=Vh2 ,_matrix=A21 )
        += integrate( _range=elements(submesh),
                      _expr=idt(u1)*id(u2) );

    form1( _test=Vh1, _vector=F1 )
        = integrate(_range=elements(mesh),
                    _expr=id(u1));

    BlocksBaseVector<double> myblockVec(2);
    myblockVec(0,0) = F1;
    myblockVec(1,0) = F2;
    auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=true);

    BlocksBaseVector<double> myblockVecSol(2);
    myblockVecSol(0,0) = u1;
    myblockVecSol(1,0) = u2;
    auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);

    BlocksBaseSparseMatrix<double> myblockMat(2,2);
    myblockMat(0,0) = A11;
    myblockMat(0,1) = A12;
    myblockMat(1,0) = A21;
    auto A = backend()->newBlockMatrix(_block=myblockMat, _copy_values=true);

    form2( _trial=Vh1, _test=Vh1 ,_matrix=A )
        +=on(_range=boundaryfaces(mesh), _rhs=F, _element=*u1,
        _expr=constant(0.) );

    backend(_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=U );

    double normL1_A = A->l1Norm();
    double normL2_F = F->l2Norm();
    double normL2_U = U->l2Norm();

    myblockVecSol.localize(U);

    BlockSolveResult result{ "V2", Vh1, Vh1->element( "u1V2" ), normL1_A, normL2_F, normL2_U };
    result.u1 = *u1;
    printBlockSolveResult( result );

    auto e = exporter( _mesh=mesh, _name="exportV2" );
    e->add( "u1", result.u1 );
    e->save();

    return result;
}

Feel::BlockSolveResult runLaplacianBlockV3( Feel::laplacian_block_mesh_ptrtype const& mesh,
                                            Feel::laplacian_block_mesh_ptrtype const& submesh )
{
    using namespace Feel;

    auto Vh1 = Pch<2>( mesh );
    auto Vh2 = Pch<1>( submesh );

    auto Xh = product( Vh1, Vh2 );
    auto a = blockform2( Xh, solve::strategy::monolithic, backend() );
    auto l = blockform1( Xh, solve::strategy::monolithic, backend() );

    auto u1 = Vh1->element( "u1" );
    auto v1 = Vh1->element( "v1" );
    auto u2 = Vh2->element( "u2" );
    auto v2 = Vh2->element( "v2" );

    a( 0_c, 0_c ) = integrate( _range=elements( mesh ),
                               _expr=inner( gradt( u1 ), grad( v1 ) ) );
    a( 0_c, 1_c ) += integrate( _range=elements( submesh ),
                                _expr=idt( u2 )*id( v1 ) );
    a( 1_c, 0_c ) += integrate( _range=elements( submesh ),
                                _expr=idt( u1 )*id( v2 ) );

    l( 0_c ) = integrate( _range=elements( mesh ),
                          _expr=id( v1 ) );

    auto U = Xh.element();
    auto u1h = U( 0_c );

    a.row( 0_c ) += on( _range=boundaryfaces( mesh ),
                        _rhs=l( 0_c ),
                        _element=u1h,
                        _expr=cst( 0.0 ),
                        _type="elimination" );

    a.solve( _rhs=l, _solution=U );

    double normL1_A = a.l1Norm();
    double normL2_F = l.l2Norm();
    double normL2_U = U.l2Norm();

    BlockSolveResult result{ "V3", Vh1, Vh1->element( "u1V3" ), normL1_A, normL2_F, normL2_U };
    result.u1 = U( 0_c );
    printBlockSolveResult( result );

    auto e = exporter( _mesh=mesh, _name="exportV3" );
    e->add( "u1", result.u1 );
    e->save();

    return result;
}


int main(int argc, char**argv )
{
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=makeLaplacianBlockOptions(),
                     _about=about(_name="laplacian_block",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto variant = normalizeVariant( soption( _name="variant" ) );
    auto mesh = createMeshLaplacianBlock();
    auto submesh = createSubMeshLaplacianBlock( mesh );

    if ( variant == "v1" )
        runLaplacianBlockV1( mesh, submesh );
    else if ( variant == "v2" )
        runLaplacianBlockV2( mesh, submesh );
    else if ( variant == "v3" )
        runLaplacianBlockV3( mesh, submesh );
    else if ( variant == "all" )
    {
        auto r1 = runLaplacianBlockV1( mesh, submesh );
        auto r2 = runLaplacianBlockV2( mesh, submesh );
        auto r3 = runLaplacianBlockV3( mesh, submesh );
        printBlockSolveComparison( mesh, r1, r2, r3 );
    }
    else
    {
        if ( Environment::worldComm().globalRank() == 0 )
            std::cerr << "Invalid variant '" << variant << "'. Expected v1, v2, v3 or all.\n";
        return 1;
    }

    return 0;
}
