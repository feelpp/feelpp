/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2018-05-21

  Copyright (C) 2018 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <feel/feelmodels/modelmesh/winslow.hpp>

#include <feel/feelvf/vf.hpp>

//#define FSI_WINSLOW_USE_OPT_EXPR 0

namespace Feel
{
namespace FeelModels
{


template< typename MeshType, int Order >
Winslow<MeshType,Order>::Winslow( mesh_ptrtype mesh, std::string const& prefix,
                                  ModelBaseRepository const& modelRep )
    :
    super_type( prefix, mesh->worldCommPtr(),"", modelRep ),
    ModelBase( prefix, mesh->worldCommPtr(),"", modelRep ),
    //M_backend( backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr() ) ),
    M_solverType( soption(_prefix=this->prefix(),_name="solver") ),
    M_mesh(mesh),
    M_flagSet(/*flagSet*/),
    M_Xh( space_type::New(_mesh=mesh ) )
{
    this->initAlgebraicBackend();
}

//----------------------------------------------------------------------------//

template< typename MeshType, int Order >
Winslow<MeshType,Order>::Winslow( space_ptrtype const& space, std::string const& prefix,
                                  ModelBaseRepository const& modelRep )
    :
    super_type( prefix, space->worldCommPtr(),"",modelRep ),
    ModelBase( prefix, space->worldCommPtr(),"",modelRep ),
    M_solverType( soption(_prefix=this->prefix(),_name="solver") ),
    M_mesh(space->mesh()),
    M_flagSet(/*flagSet*/),
    M_Xh(space)
{
    this->initAlgebraicBackend();
}

//----------------------------------------------------------------------------//

template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::init()
{
    this->log("Winslow","init", "start" );


    M_displacement.reset( new element_type(M_Xh) );
    M_displacementOld.reset( new element_type(M_Xh) );
    M_dispImposedOnBoundary.reset( new element_type(M_Xh) );
    M_identity.reset( new element_type(M_Xh) );
    M_identity->on( _range=elements(M_mesh), _expr=vf::P() );

    //M_XhScalP0Disc = space_p0_type::New(_mesh=M_mesh);
#if 0
    M_XhScalM1 = space_scal_m1_type::New(_mesh=M_mesh);
    M_l2projector = opProjection(_domainSpace=this->functionSpaceScalM1(),
                                 _imageSpace=this->functionSpaceScalM1(),
                                 _backend=backend_type::build( soption( _name="backend" ), prefixvm(this->prefix(),"l2proj"), this->worldCommPtr() ),
                                 _type=Feel::L2 ) );
#endif


    // update dofsWithValueImposed
    std::set<flag_type> flagsMovingAndFixed;
    if ( this->hasFlagSet("moving" ) )
        flagsMovingAndFixed.insert( this->flagSet("moving").begin(), this->flagSet("moving").end() );
    if ( this->hasFlagSet("fixed" ) )
        flagsMovingAndFixed.insert( this->flagSet("fixed").begin(), this->flagSet("fixed").end() );
    M_dofsWithValueImposed.clear();
    for ( auto const& faceWrap : markedfaces(M_Xh->mesh(), flagsMovingAndFixed ) )
    {
        auto const& face = unwrap_ref( faceWrap );
        auto facedof = M_Xh->dof()->faceLocalDof( face.id() );
        for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
        {
            M_dofsWithValueImposed.insert( it->index() );
        }
    }


    M_fieldMetricDerivative.reset( new element_type( M_Xh ) );
    M_backendMetricDerivative = backend_type::build( soption( _name="backend" ), prefixvm(this->prefix(),"metric-derivative"), this->worldCommPtr() );
    M_matrixMetricDerivative = M_backendMetricDerivative->newMatrix(_test=M_Xh,_trial=M_Xh);
    M_vectorMetricDerivative = M_backendMetricDerivative->newVector(M_Xh);

    auto Xh = this->functionSpace();
    auto mesh = this->functionSpace()->mesh();
    auto rangeElt = elements(support(Xh));
    auto const& u = *M_displacement;
    form2( _trial=Xh, _test=Xh,_matrix=M_matrixMetricDerivative)
    += integrate(_range=rangeElt,
                     _expr=inner(idt(u),id(u)) );
    M_matrixMetricDerivative->close();

//M_useMeshAdapation = boption(_name="mesh-adaptation",_prefix=this->prefix() );
//  M_useMeshAdapationScalar = boption(_name="mesh-adaptation.scalar-weight",_prefix=this->prefix() );

    M_weightFunctionScalar = M_Xh->compSpace()->elementPtr();
    M_weightFunctionScalar->setConstant( 1.0 );
#if 0
if ( Environment::vm().count( prefixvm(this->prefix(),"mesh-adaptation.scalar-weight.expr").c_str() ) )
    {
        auto eScal = expr( soption(_name="mesh-adaptation.scalar-weight.expr", _prefix=this->prefix() ) );
        M_weightFunctionScalar->on(_range=elements(mesh),_expr=eScal);
    }
    if ( M_useMeshAdapation && !M_useMeshAdapationScalar )
    {
        M_XhTensor2P0Disc = space_p0_tensor2_type::New(_mesh=mesh);
        M_weightFunctionTensor2.reset( new element_p0_tensor2_type( M_XhTensor2P0Disc ) );
        M_weightFunctionTensor2->on(_range=elements(mesh),_expr=Id<mesh_type::nDim>());
    }
#endif

    this->setRebuildCstPartInLinearSystem(false);
    this->setUseLinearJacobianInResidual(false);
    this->setRebuildLinearPartInJacobian(false);
    this->setUseCstMatrix(false);
    this->setUseCstVector(false);

    if ( true/*false*//*true*/ )
    {
        int nBlock=2;

        M_blockVectorSolution.resize( nBlock );
        int cptBlock = 0;
        M_blockVectorSolution(cptBlock++) = M_displacement;//this->backend()->newVector( this->functionSpace() );
        *M_blockVectorSolution(0) = *M_identity;
        M_blockVectorSolution(cptBlock++) = M_fieldMetricDerivative;//this->backend()->newVector( this->functionSpace() );
        M_blockVectorSolution.buildVector( this->backend() );

        BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
        int indexBlock=0;
        auto graph = stencil( _test=M_Xh, _trial=M_Xh )->graph();
        myblockGraph(0,0) = graph;
        myblockGraph(0,1) = graph;
        myblockGraph(1,0) = graph;
        myblockGraph(1,1) = graph;
        myblockGraph.close();

        graph_ptrtype fullGraph( new graph_type( myblockGraph ) );
        M_algebraicFactoryMixedFormulation.reset( new model_algebraic_factory_type( this->shared_from_this(),Feel::backend(_rebuild=true,_name=this->backend()->prefix())/*this->backend()*/,
                                                                                    fullGraph, fullGraph->mapRow().indexSplit() ) );
    }
    if ( true )//else
    {
        auto graph = stencil( _test=M_Xh, _trial=M_Xh )->graph();
        M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),Feel::backend(_rebuild=true,_name=this->backend()->prefix())/*this->backend()*/,
                                                                    graph, graph->mapRow().indexSplit() ) );
        // mesh not move so not rebuild cst part
        M_vectorSolution = this->backend()->newVector( this->functionSpace() );
        *M_vectorSolution = *M_identity;
    }

    this->log("Winslow","init", "finish" );
}

template< typename MeshType, int Order >
std::shared_ptr<std::ostringstream>
Winslow<MeshType,Order>::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );

    *_ostr << "\n   Winslow "
           << "\n     -- solver type : " << M_solverType;
        //<<  "\n    -- mesh adaptation : " << M_useMeshAdapation;
        //if ( M_useMeshAdapation )
        //*_ostr <<  "\n       + type : " << ( ( M_useMeshAdapationScalar)? "Scalar" : "Tensor2" );

    return _ostr;
}

//----------------------------------------------------------------------------//

template< typename MeshType,int Order >
void
Winslow<MeshType,Order>::solve()
{
    this->log("Winslow","solve", "start" );

    if ( M_solverType=="Picard" || M_solverType == "FixPoint" )
    {
        //ModelAlgebraic::DataNewtonInitialGuess dataInitialGuess( M_vectorSolution );
        //this->updateNewtonInitialGuess( dataInitialGuess );
#if 1
        M_algebraicFactory->solvePicard( M_vectorSolution );
        *M_displacement = *M_vectorSolution;
#else
        M_blockVectorSolution.updateVectorFromSubVectors();
        M_algebraicFactoryMixedFormulation->solvePicard( M_blockVectorSolution.vectorMonolithic() );
        M_blockVectorSolution.localize();
#endif
    }
    else if ( M_solverType=="Newton" )
    {
#if 0
        M_algebraicFactory->solveNewton( M_vectorSolution );
#else
        M_blockVectorSolution.updateVectorFromSubVectors();
        M_algebraicFactoryMixedFormulation->solveNewton( M_blockVectorSolution.vectorMonolithic() );
        M_blockVectorSolution.localize();
#endif
    }
    else if (M_solverType=="Picard-Newton")
    {
        int saveMaxIt = M_algebraicFactory->backend()->maxIterationsSNES();
        double rtol = M_algebraicFactory->backend()->rToleranceSNES();
        double atol = M_algebraicFactory->backend()->aToleranceSNES();
        double stol = M_algebraicFactory->backend()->sToleranceSNES();
        int picardMaxIt = ioption(_prefix=this->prefix(),_name="Picard-Newton.maxit-Picard");
        M_algebraicFactory->backend()->setTolerancesSNES(_rtolerance=rtol,_maxit=picardMaxIt,_atolerance=atol,_stolerance=stol);
        M_algebraicFactory->solvePicard( M_vectorSolution );
        M_algebraicFactory->backend()->setTolerancesSNES(_rtolerance=rtol,_maxit=saveMaxIt,_atolerance=atol,_stolerance=stol);
        //M_algebraicFactory->solveNewton( M_vectorSolution );
        *M_displacement = *M_vectorSolution;
        M_blockVectorSolution.updateVectorFromSubVectors();
        M_algebraicFactoryMixedFormulation->solveNewton( M_blockVectorSolution.vectorMonolithic() );
        M_blockVectorSolution.localize();
    }

    //update displacement
#if 0
    *M_displacement = *M_vectorSolution;
#endif
    *M_displacement -= *M_identity;

    this->log("Winslow","solve", "finish" );
}


//----------------------------------------------------------------------------//

template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateLinearPDE( DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto Xh = this->functionSpace();
    auto mesh = this->functionSpace()->mesh();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;

    auto rangeElt = elements(support(Xh));
    auto rangeBoundaryFaces = boundaryfaces(support(Xh));
#if 0
    const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
    auto u = Xh->element( vecCurrentPicardSolution, 0 );
    auto alpha = Xh->element( vecCurrentPicardSolution, 1 );

    auto G = trans(gradv(u))*gradv(u);
    auto invG = inv(G);

    form2( _test=Xh, _trial=Xh, _matrix=A,
           _rowstart=0,_colstart=0 ) +=
        integrate( _range=elements(mesh),
                   _expr= inner(invG*trans(gradt(u)),trans(grad(u))) );

    form2( _test=Xh, _trial=Xh, _matrix=A,
           _rowstart=0,_colstart=0 ) +=
        integrate( _range=elements(mesh),
                   _expr= -inner( gradt(u)*idv(alpha), id(u) ) );


    form2( _test=Xh, _trial=Xh, _matrix=A,
           _rowstart=1,_colstart=1 ) +=
        integrate( _range=elements(mesh),
                   _expr= inner( idt(alpha), id(alpha) ) );

    form1( _test=Xh, _vector=F,
           _rowstart=1 ) +=
        integrate( _range=elements(mesh),
                   _expr= inner( invG, grad(alpha) ) );

    form1( _test=Xh, _vector=F,
           _rowstart=1 ) +=
        integrate( _range=boundaryfaces(mesh),
                   _expr= -inner( invG*N(), id(alpha) ) );

#endif

#if 1
    const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
    auto u = Xh->element( vecCurrentPicardSolution );
    auto G = trans(gradv(u))*gradv(u);
    auto invG = inv(G);

    M_vectorMetricDerivative->zero();
    auto lf = form1( _test=Xh,_vector=M_vectorMetricDerivative );
    lf += integrate(_range=rangeElt,
                    _expr=inner(invG, grad(u)) );
    lf += integrate(_range=rangeBoundaryFaces,//boundaryfaces(mesh),
                    _expr=-inner(invG*N(), id(u)) );
    //M_vectorMetricDerivative->close();
    M_backendMetricDerivative->solve(_matrix=M_matrixMetricDerivative,_solution=*M_fieldMetricDerivative,_rhs=M_vectorMetricDerivative);


    if ( false /*M_useMeshAdapation && !M_useMeshAdapationScalar*/ )
    {
#if 0
        auto tau = idv(M_weightFunctionTensor2);
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            integrate( _range=rangeElt,
                       _expr= inner(invG*trans(tau*gradt(u)),trans(grad(u))) );
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            integrate( _range=rangeElt,
                       _expr=-inner(tau*gradt(u)*idv(M_fieldMetricDerivative),id(u)) );
#else
        CHECK( false ) << "not work";
#endif
    }
    else
    {
        auto tau = idv( M_weightFunctionScalar );
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            integrate( _range=rangeElt,
                       _expr= (tau)*inner(invG*trans(gradt(u)),trans(grad(u))) );
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            integrate( _range=rangeElt,
                       _expr=-(tau)*inner(gradt(u)*idv(M_fieldMetricDerivative),id(u)) );
    }
#endif
}


template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto Xh = this->functionSpace();
    auto mesh = this->functionSpace()->mesh();
    auto const& u = *M_displacement;

    if ( this->hasFlagSet("moving" ) )
    {
        auto aux_pol = idv( this->dispImposedOnBoundary() ) + idv(this->identity());
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            on( _range=markedfaces( mesh, this->flagSet("moving") ),
                _element=u, _rhs=F,
                _expr=aux_pol/*, ON_ELIMINATION*/ );
    }

    if ( this->hasFlagSet("fixed" ) )
    {
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            on( _range=markedfaces(mesh, this->flagSet("fixed") ),
                _element=u, _rhs=F,
                _expr=idv(this->identity())/*, ON_ELIMINATION*/ );

    }
}

//----------------------------------------------------------------------------//

template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    this->log("Winslow","updateNewtonInitialGuess", "start" );

    vector_ptrtype& U = data.initialGuess();
    auto Xh = this->functionSpace();
    auto mesh = Xh->mesh();
    auto u = Xh->element( U );

    u.on(_range=markedfaces(mesh, this->flagSet("moving") ),
         _expr=idv(M_dispImposedOnBoundary)+idv(M_identity) );
    u.on(_range=markedfaces(mesh, this->flagSet("fixed") ),
         _expr=idv(M_identity) );

    sync( u, "=", M_dofsWithValueImposed );

    this->log("Winslow","updateNewtonInitialGuess", "finish" );
}


template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool BuildNonCstPart = !buildCstPart;
    if ( buildCstPart )
        return;

    this->log("Winslow","updateResidual", "start" );

    auto Xh = this->functionSpace();
    auto mesh = Xh->mesh();
    auto u = Xh->element( XVec, 0 );
    auto alpha = Xh->element( XVec, 1 );
    auto G = trans(gradv(u))*gradv(u);
    auto invG = inv(G);

    auto linearFormMapping = form1( _test=Xh, _vector=R,_rowstart=0 );
    auto linearFormLagMult = form1( _test=Xh, _vector=R,_rowstart=1 );

    auto rangeElt = elements(support(Xh));
    auto rangeBoundaryFaces = boundaryfaces(support(Xh));

    if ( false/*M_useMeshAdapation && !M_useMeshAdapationScalar*/ )
    {
#if 0
        auto tau = idv(M_weightFunctionTensor2);
        linearFormMapping +=
            integrate( _range=rangeElt,
                       _expr= inner(invG*trans(tau*gradv(u)),trans(grad(u))) );
        linearFormMapping +=
            integrate( _range=rangeElt,
                       _expr= -inner( tau*gradv(u)*idv(alpha), id(u) ) );
#else
        CHECK( false ) << "not work";
#endif
    }
    else
    {
        auto tau = idv( M_weightFunctionScalar );
        linearFormMapping +=
            integrate( _range=rangeElt,
                       _expr= tau*inner(invG*trans(gradv(u)),trans(grad(u))) );
        linearFormMapping +=
            integrate( _range=rangeElt,
                       _expr= -tau*inner( gradv(u)*idv(alpha), id(u) ) );
    }
    linearFormLagMult +=
        integrate( _range=rangeElt,
                   _expr= inner( idv(alpha), id(alpha) ) );
    linearFormLagMult +=
        integrate( _range=rangeElt,
                   _expr= -inner( invG, grad(alpha) ) );
    linearFormLagMult +=
        integrate( _range=rangeBoundaryFaces,//boundaryfaces(mesh),
                   _expr= inner( invG*N(), id(alpha) ) );

    this->log("Winslow","updateResidual", "finish" );
}


template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    this->log("Winslow","updateResidualDofElimination", "start" );

    vector_ptrtype& R = data.residual();
    auto Xh = this->functionSpace();
    auto mesh = Xh->mesh();
    auto residualView = Xh->element( R );

    for ( size_type thedof : M_dofsWithValueImposed )
        residualView.set( thedof,0. );
    sync( residualView, "=", M_dofsWithValueImposed );

    this->log("Winslow","updateResidualDofElimination", "finish" );
}


template<typename ElementMappingType>
auto
MyInvGt( ElementMappingType const& u,  mpl::int_<2> /**/ )
{
    auto G = trans(gradv(u))*gradv(u);

    auto Gt00 = 2*gradt(u)(0,0)*gradv(u)(0,0) +  2*gradt(u)(1,0)*gradv(u)(1,0);
    auto Gt11 = 2*gradt(u)(0,1)*gradv(u)(0,1) +  2*gradt(u)(1,1)*gradv(u)(1,1);
    auto Gt01 = gradt(u)(0,1)*gradv(u)(0,0) + gradv(u)(0,1)*gradt(u)(0,0) + gradt(u)(1,1)*gradv(u)(1,0) + gradv(u)(1,1)*gradt(u)(1,0);
    auto Gt10 = Gt01;
    auto detGt = Gt00*G(1,1)+G(0,0)*Gt11 - Gt01*G(1,0) - G(0,1)*Gt10;
    auto detG = G(0,0)*G(1,1) - G(0,1)*G(1,0);
    auto invGt = (cst(1.)/detG)*mat<2,2>( Gt11, -Gt01, -Gt10, Gt00 ) - detGt*(cst(1.)/pow(detG,2))*mat<2,2>( G(1,1), -G(0,1), -G(1,0), G(0,0) );
    return invGt;
}
template<typename ElementMappingType>
auto
MyInvGt( ElementMappingType const& u,  mpl::int_<3> /**/ )
{
    CHECK( false ) << "TODO";
    return vf::zero<3,3>();
}

template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    bool buildCstPart = data.buildCstPart();

    if ( buildCstPart )
        return;

    this->log("Winslow","updateJacobian", "start" );
    auto Xh = this->functionSpace();
    auto mesh = Xh->mesh();
    auto rangeElt = elements(support(Xh));
    auto rangeBoundaryFaces = boundaryfaces(support(Xh));

    auto u = Xh->element( XVec, 0 );
    auto alpha = Xh->element( XVec, 1 );

    auto G = trans(gradv(u))*gradv(u);
    auto invG = inv(G);
    auto invGt = MyInvGt( u, mpl::int_<mesh_type::nDim>() );

    if ( false/*M_useMeshAdapation && !M_useMeshAdapationScalar*/ )
    {
#if 0
        auto tau = idv(M_weightFunctionTensor2);
        form2( _test=Xh, _trial=Xh, _matrix=J,
               _rowstart=0,_colstart=0 ) +=
            integrate( _range=rangeElt,
                       _expr= inner(invG*trans(tau*gradt(u)),trans(grad(u))) );
        form2( _test=Xh, _trial=Xh, _matrix=J,
               _rowstart=0,_colstart=0 ) +=
            integrate( _range=rangeElt,
                       _expr= inner(invGt*trans(tau*gradv(u)),trans(grad(u))) );

        form2( _test=Xh, _trial=Xh, _matrix=J,
               _rowstart=0,_colstart=0 ) +=
            integrate( _range=rangeElt,
                       _expr= -inner( tau*gradt(u)*idv(alpha), id(u) ) );
        form2( _test=Xh, _trial=Xh, _matrix=J,
               _rowstart=0,_colstart=1 ) +=
            integrate( _range=rangeElt,
                       _expr= -inner( tau*gradv(u)*idt(alpha), id(u) ) );
#else
        CHECK( false ) << "not work";
#endif
    }
    else
    {
        auto tau = idv( M_weightFunctionScalar );
        form2( _test=Xh, _trial=Xh, _matrix=J,
               _rowstart=0,_colstart=0 ) +=
            integrate( _range=rangeElt,
                       _expr= tau*inner(invG*trans(gradt(u)),trans(grad(u))) );
        form2( _test=Xh, _trial=Xh, _matrix=J,
               _rowstart=0,_colstart=0 ) +=
            integrate( _range=rangeElt,
                       _expr= tau*inner(invGt*trans(gradv(u)),trans(grad(u))) );

        form2( _test=Xh, _trial=Xh, _matrix=J,
               _rowstart=0,_colstart=0 ) +=
            integrate( _range=rangeElt,
                       _expr= -tau*inner( gradt(u)*idv(alpha), id(u) ) );
        form2( _test=Xh, _trial=Xh, _matrix=J,
               _rowstart=0,_colstart=1 ) +=
            integrate( _range=rangeElt,
                       _expr= -tau*inner( gradv(u)*idt(alpha), id(u) ) );
    }

    form2( _test=Xh, _trial=Xh, _matrix=J,
           _rowstart=1,_colstart=1 ) +=
        integrate( _range=rangeElt,
                   _expr= inner( idt(alpha), id(alpha) ) );

    form2( _test=Xh, _trial=Xh, _matrix=J,
               _rowstart=1,_colstart=0 ) +=
        integrate( _range=rangeElt,
                   _expr= -inner( invGt*eye<mesh_type::nDim,mesh_type::nDim>(), grad(alpha) ) );
    form2( _test=Xh, _trial=Xh, _matrix=J,
               _rowstart=1,_colstart=0 ) +=
        integrate( _range=rangeBoundaryFaces,
                   _expr= inner( invGt*N(), id(alpha) ) );


    this->log("Winslow","updateJacobian", "finish" );
}

template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    this->log("Winslow","updateJacobianDofElimination", "start" );

    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    auto Xh = this->functionSpace();
    auto mesh = Xh->mesh();
    auto const& u = *M_displacement;

    if ( this->hasFlagSet("moving" ) )
    {
        form2( _test=Xh, _trial=Xh, _matrix=J ) +=
            on( _range=markedfaces( mesh, this->flagSet("moving") ),
                _element=u, _rhs=RBis, _expr=vf::zero<mesh_type::nDim,1>() );
    }
    if ( this->hasFlagSet("fixed" ) )
    {
        form2( _test=Xh, _trial=Xh, _matrix=J ) +=
            on( _range=markedfaces(mesh, this->flagSet("fixed") ),
                _element=u, _rhs=RBis, _expr=vf::zero<mesh_type::nDim,1>() );
    }

    this->log("Winslow","updateJacobianDofElimination", "finish" );
}











#if 0


template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateJacobian( DataUpdateJacobian & data, mpl::int_<2> /**/ ) const
{

#if 0
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    bool buildCstPart = data.buildCstPart();

    if ( buildCstPart )
        return;

    this->log("Winslow","updateJacobian<2>", "start" );
    auto Xh = this->functionSpace();
    auto mesh = Xh->mesh();
    auto u = Xh->element( XVec );

    auto du1dxv = gradv(u)(0,0);
    auto du1dyv = gradv(u)(0,1);
    auto du2dxv = gradv(u)(1,0);
    auto du2dyv = gradv(u)(1,1);

    auto du1dxt = gradt(u)(0,0);
    auto du1dyt = gradt(u)(0,1);
    auto du2dxt = gradt(u)(1,0);
    auto du2dyt = gradt(u)(1,1);

    auto alpha1 = du1dyv*du1dyv + du2dyv*du2dyv;
    auto beta1 = du1dxv*du1dyv + du2dxv*du2dyv;
    auto gamma1 = du1dxv*du1dxv + du2dxv*du2dxv;

    /*auto l2p = projector(WinslowFactory.functionSpaceScalM1(),
     WinslowFactory.functionSpaceScalM1(),
     winslow_type::backend_type::build( WinslowFactory.vm() ),
     Feel::L2);*/

    boost::mpi::timer thetimer;


    auto proj_alpha = this->l2projector()->operator()(alpha1);
    auto proj_beta = this->l2projector()->operator()(beta1);
    auto proj_gamma = this->l2projector()->operator()(gamma1);

    auto proj_du1dx = this->l2projector()->operator()(du1dxv);
    auto proj_du1dy = this->l2projector()->operator()(du1dyv);
    auto proj_du2dx = this->l2projector()->operator()(du2dxv);
    auto proj_du2dy = this->l2projector()->operator()(du2dyv);

    double tProj = thetimer.elapsed();
    thetimer.restart();

    //auto Metric1 = vf::mat<2,2>(alpha1,-beta1,-beta1,gamma1);
#if 0
    auto alpha2 = 2*(du1dyv*du1dyt + du2dyv*du2dyt);
    auto beta2 = du1dxv*du1dyt + du1dxt*du1dyv + du2dxv*du2dyt + du2dxt*du2dyv;
    auto gamma2 = 2*(du1dxv*du1dxt + du2dxv*du2dxt);
#else
    auto alpha2 = 2*(idv(proj_du1dy)*du1dyt + idv(proj_du2dy)*du2dyt);
    auto beta2 = idv(proj_du1dx)*du1dyt + du1dxt*idv(proj_du1dy) + idv(proj_du2dx)*du2dyt + du2dxt*idv(proj_du2dy);
    auto gamma2 = 2*(idv(proj_du1dx)*du1dxt + idv(proj_du2dx)*du2dxt);
#endif


#if !FSI_WINSLOW_USE_OPT_EXPR
    form2( _test=Xh, _trial=Xh, _matrix=J ) +=
        integrate( _range=elements(mesh),
                   _expr=
                   //comp 1
                   - gradt(u)(0,0)*(dxv(proj_alpha)*id(u)(0,0)+ idv(proj_alpha)*grad(u)(0,0) )
                   + gradt(u)(0,0)*(dyv(proj_beta)*id(u)(0,0) + idv(proj_beta)*grad(u)(0,1) )
                   + gradt(u)(0,1)*(dxv(proj_beta)*id(u)(0,0) + idv(proj_beta)*grad(u)(0,0) )
                   - gradt(u)(0,1)*(dyv(proj_gamma)*id(u)(0,0) + idv(proj_gamma)*grad(u)(0,1) )
                   //comp 2
                   - gradt(u)(1,0)*(dxv(proj_alpha)*id(u)(1,0)+ idv(proj_alpha)*grad(u)(1,0) )
                   + gradt(u)(1,0)*(dyv(proj_beta)*id(u)(1,0) + idv(proj_beta)*grad(u)(1,1) )
                   + gradt(u)(1,1)*(dxv(proj_beta)*id(u)(1,0) + idv(proj_beta)*grad(u)(1,0) )
                   - gradt(u)(1,1)*(dyv(proj_gamma)*id(u)(1,0) + idv(proj_gamma)*grad(u)(1,1) ) );

    form2( _test=Xh, _trial=Xh, _matrix=J ) +=
        integrate( _range=elements(mesh),
                   _expr= (alpha2*dxv(proj_du1dx)-/*2**/beta2*dxv(proj_du1dy) - beta2*dyv(proj_du1dx)   + gamma2*dyv(proj_du1dy) )*id(u)(0,0) );

    form2( _test=Xh, _trial=Xh, _matrix=J ) +=
        integrate( _range=elements(mesh),
                   _expr= (alpha2*dxv(proj_du2dx)-/*2**/beta2*dxv(proj_du2dy) - beta2*dyv(proj_du2dx)   +gamma2*dyv(proj_du2dy) )*id(u)(1,0) );

#else
    std::map<std::string, element_scal_m1_type> myProj;
    myProj["alpha"]=proj_alpha;
    myProj["beta"]=proj_beta;
    myProj["gamma"]=proj_gamma;
    myProj["du1dxv"]=proj_du1dx;
    myProj["du1dyv"]=proj_du1dy;
    myProj["du2dxv"]=proj_du2dx;
    myProj["du2dyv"]=proj_du2dy;
    const uint16_type nQuadOrderElement = 2*winslow_type::space_type::basis_type::nOrder+winslow_type::space_scal_m1_type::basis_type::nOrder-2;
    auto const winslowJacobianExpr = Feel::vf::FSI::winslowJacobian<nQuadOrderElement>(u,myProj);
    form2( _test=Xh, _trial=Xh, _matrix=J ) +=
        integrate( _range=elements(mesh),
                   _expr=winslowJacobianExpr );
#endif

#if 1
    // update markedfaces free
    if ( this->hasFlagSet("free" ) )
    {
        form2( _test=Xh, _trial=Xh, _matrix=J ) +=
            integrate( _range=markedfaces( mesh, this->flagSet("free") ),
                       //_range=boundaryfaces(mesh),
                       _expr=
                       // comp 1
                       + gradt(u)(0,0)*vf::Nx()*idv(proj_alpha)*id(u)(0,0)
                       //- gradt(u)(0,0)*vf::Nx()*idv(proj_gamma)*id(u)(0)
                       + gradt(u)(0,1)*vf::Ny()*idv(proj_gamma)*id(u)(0,0)
                       - gradt(u)(0,1)*vf::Nx()*idv(proj_beta)*id(u)(0,0)
                       - gradt(u)(0,0)*vf::Ny()*idv(proj_beta)*id(u)(0,0)
                       // comp 2
                       + gradt(u)(1,0)*vf::Nx()*idv(proj_alpha)*id(u)(1,0)
                       //- gradt(u)(1,0)*vf::Nx()*idv(proj_gamma)*id(u)(1)
                       + gradt(u)(1,1)*vf::Ny()*idv(proj_gamma)*id(u)(1,0)
                       - gradt(u)(1,1)*vf::Nx()*idv(proj_beta)*id(u)(1,0)
                       - gradt(u)(1,0)*vf::Ny()*idv(proj_beta)*id(u)(1,0)
                       );
    }
#endif

    double tAssembly = thetimer.elapsed();

    this->log("Winslow","updateJacobian<2>", "finish" );


#endif
}
template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateJacobian( DataUpdateJacobian & data, mpl::int_<3> /**/ ) const
{
    CHECK( false ) << "TODO";
}

template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateResidual( DataUpdateResidual & data, mpl::int_<2> /**/ ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool BuildNonCstPart = !BuildCstPart;
    if ( BuildCstPart )
        return;

    this->log("Winslow","updateResidual<2>", "start" );

    auto Xh = this->functionSpace();
    auto mesh = Xh->mesh();
    auto u = Xh->element( XVec );

    auto du1dxv = gradv(u)(0,0);
    auto du1dyv = gradv(u)(0,1);
    auto du2dxv = gradv(u)(1,0);
    auto du2dyv = gradv(u)(1,1);

    auto alpha = pow(du1dyv,2) + pow(du2dyv,2);
    auto beta = du1dxv*du1dyv + du2dxv*du2dyv;
    auto gamma = pow(du1dxv,2) + pow(du2dxv,2);
#if 1
    auto proj_alpha = this->l2projector()->operator()(alpha);
    auto proj_beta = this->l2projector()->operator()(beta);
    auto proj_gamma = this->l2projector()->operator()(gamma);

#if !FSI_WINSLOW_USE_OPT_EXPR
    form1( _test=Xh, _vector=R ) +=
        integrate( _range=elements(mesh),
                   _expr=
                   //comp 1
                   - gradv(u)(0,0)*(dxv(proj_alpha)*id(u)(0,0)+ idv(proj_alpha)*grad(u)(0,0) )
                   + gradv(u)(0,0)*(dyv(proj_beta)*id(u)(0,0) + idv(proj_beta)*grad(u)(0,1) )
                   + gradv(u)(0,1)*(dxv(proj_beta)*id(u)(0,0) + idv(proj_beta)*grad(u)(0,0) )
                   - gradv(u)(0,1)*(dyv(proj_gamma)*id(u)(0,0) + idv(proj_gamma)*grad(u)(0,1) )
                   //comp 2
                   - gradv(u)(1,0)*(dxv(proj_alpha)*id(u)(1,0)+ idv(proj_alpha)*grad(u)(1,0) )
                   + gradv(u)(1,0)*(dyv(proj_beta)*id(u)(1,0) + idv(proj_beta)*grad(u)(1,1) )
                   + gradv(u)(1,1)*(dxv(proj_beta)*id(u)(1,0) + idv(proj_beta)*grad(u)(1,0) )
                   - gradv(u)(1,1)*(dyv(proj_gamma)*id(u)(1,0) + idv(proj_gamma)*grad(u)(1,1) ) );
#else
    std::map<std::string, element_scal_m1_type> myProj;
    myProj["alpha"]=proj_alpha;
    myProj["beta"]=proj_beta;
    myProj["gamma"]=proj_gamma;
    const uint16_type nQuadOrderElement = 2*winslow_type::space_type::basis_type::nOrder+winslow_type::space_scal_m1_type::basis_type::nOrder-2;
    //auto const winslowResidualExpr = Feel::vf::FSI::winslowResidual<nQuadOrderElement>(u,proj_alpha,proj_beta,proj_gamma);
    auto const winslowResidualExpr = Feel::vf::FSI::winslowResidual<nQuadOrderElement>(u,myProj);
    form1( _test=Xh, _vector=R ) +=
        integrate( _range=elements(mesh),
                   _expr=winslowResidualExpr );
#endif
#else
    auto ux = u[Component::X];
    auto uy = u[Component::Y];
    auto uxhess = hessv(ux);
    auto uyhess = hessv(uy);
    auto g11 = alpha;
    auto g12 = beta;
    auto g22 = gamma;
    auto dxvg11 = 2*(uxhess(0,1)*dyv(ux)+uyhess(0,1)*dyv(uy));
    auto dyvg22 = 2*(uxhess(1,0)*dxv(ux)+uyhess(1,0)*dxv(uy));
    auto dxvg12 = uxhess(0,0)*dyv(ux)+dxv(ux)*uxhess(0,1) + uyhess(0,0)*dyv(uy)+dxv(uy)*uyhess(0,1);
    auto dyvg12 = uxhess(0,1)*dyv(ux)+dxv(ux)*uxhess(1,1) + uyhess(0,1)*dyv(uy)+dxv(uy)*uyhess(1,1);
    form1( _test=Xh, _vector=R ) +=
        integrate( _range=elements(mesh),
                   _expr=
                   //comp 1
                   - gradv(u)(0,0)*(dxvg11*id(u)(0,0)+ g11*grad(u)(0,0) )
                   + gradv(u)(0,0)*(dyvg12*id(u)(0,0) + g12*grad(u)(0,1) )
                   + gradv(u)(0,1)*(dxvg12*id(u)(0,0) + g12*grad(u)(0,0) )
                   - gradv(u)(0,1)*(dyvg22*id(u)(0,0) + g22*grad(u)(0,1) )
                   //comp 2
                   - gradv(u)(1,0)*(dxvg11*id(u)(1,0) + g11*grad(u)(1,0) )
                   + gradv(u)(1,0)*(dyvg12*id(u)(1,0) + g12*grad(u)(1,1) )
                   + gradv(u)(1,1)*(dxvg12*id(u)(1,0) + g12*grad(u)(1,0) )
                   - gradv(u)(1,1)*(dyvg22*id(u)(1,0) + g22*grad(u)(1,1) ) );
#endif
#if 0
    // update markedfaces free
    if ( this->hasFlagSet("free" ) )
    {
        form1( _test=Xh, _vector=R ) +=
            integrate( _range=markedfaces( mesh, this->flagSet("free") ),
                       //_range=boundaryfaces(mesh),
                       _expr=
                       // comp 1
                       + gradv(u)(0,0)*vf::Nx()*idv(proj_alpha)*id(u)(0,0)
                       //- (gradv(u)(0,0)*vf::Nx()/*+vf::Nx()+vf::Ny()*/)*idv(proj_gamma)*id(u)(0)
                       + gradv(u)(0,1)*vf::Ny()*idv(proj_gamma)*id(u)(0,0)
                       - gradv(u)(0,1)*vf::Nx()*idv(proj_beta)*id(u)(0,0)
                       - gradv(u)(0,0)*vf::Ny()*idv(proj_beta)*id(u)(0,0)
                       // comp 2
                       + gradv(u)(1,0)*vf::Nx()*idv(proj_alpha)*id(u)(1,0)
                       //- (gradv(u)(1,0)*vf::Nx()/*+vf::Nx()+vf::Ny()*/ )*idv(proj_gamma)*id(u)(1)
                       + gradv(u)(1,1)*vf::Ny()*idv(proj_gamma)*id(u)(1,0)
                       - gradv(u)(1,1)*vf::Nx()*idv(proj_beta)*id(u)(1,0)
                       - gradv(u)(1,0)*vf::Ny()*idv(proj_beta)*id(u)(1,0)
                       );
    }
#if 0
    form1( _test=Xh, _vector=R ) +=
        integrate( _range=boundaryfaces(mesh),
                   _expr=
                   // comp 1
                   + vf::Nx()*idv(proj_alpha)*id(u)(0)
                   // comp 2
                   + vf::Ny()*idv(proj_gamma)*id(u)(1)
                   );
#endif
#endif

    this->log("Winslow","updateResidual<2>", "finish" );
}
template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateResidual( DataUpdateResidual & data, mpl::int_<3> /**/ ) const
{
    CHECK( false ) << "TODO";
}


template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateLinearPDE( DataUpdateLinear & data, mpl::int_<2> /**/ ) const
{
#if 1
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto Xh = this->functionSpace();
    auto mesh = this->functionSpace()->mesh();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;

    const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
    auto u = Xh->element( vecCurrentPicardSolution );
    auto ux = u[Component::X];
    auto uy = u[Component::Y];

    auto XhScalM1 = this->functionSpaceScalM1();

    auto G = trans(gradv(u))*gradv(u);
    auto invG = inv(G);
#if 0
    auto g11 = invG(0,0);
    auto g12 = invG(0,1);
    auto g22 = invG(1,1);
#else
    auto g11 = G(1,1);
    auto g12 = G(0,1);
    auto g22 = G(0,0);
#endif

#if 1
    //auto l2p = projector(XhScalM1,XhScalM1,winslow_type::backend_type::build( WinslowFactory.vm() ),Feel::L2);
    auto l2p = this->l2projector();
    auto proj_alpha/*proj_g11*/ = (*l2p)(g11);
    auto proj_beta/*proj_g12*/ = (*l2p)(g12);
    auto proj_gamma/* proj_g22*/ = (*l2p)(g22);
    form2( _test=Xh, _trial=Xh, _matrix=A ) +=
        integrate( _range=elements(mesh),
                   _expr=
                   //comp 1
                   - gradt(u)(0,0)*(dxv(proj_alpha)*id(u)(0,0) + idv(proj_alpha)*grad(u)(0,0) )
                   + gradt(u)(0,0)*(dyv(proj_beta)*id(u)(0,0) + idv(proj_beta)*grad(u)(0,1) )
                   + gradt(u)(0,1)*(dxv(proj_beta)*id(u)(0,0) + idv(proj_beta)*grad(u)(0,0) )
                   - gradt(u)(0,1)*(dyv(proj_gamma)*id(u)(0,0) + idv(proj_gamma)*grad(u)(0,1) )
                   //comp 2
                   - gradt(u)(1,0)*(dxv(proj_alpha)*id(u)(1,0) + idv(proj_alpha)*grad(u)(1,0) )
                   + gradt(u)(1,0)*(dyv(proj_beta)*id(u)(1,0) + idv(proj_beta)*grad(u)(1,1) )
                   + gradt(u)(1,1)*(dxv(proj_beta)*id(u)(1,0) + idv(proj_beta)*grad(u)(1,0) )
                   - gradt(u)(1,1)*(dyv(proj_gamma)*id(u)(1,0) + idv(proj_gamma)*grad(u)(1,1) ) );
#else
    tic();
    auto uxhess = hessv(ux);
    auto uyhess = hessv(uy);
    auto dxvg11 = 2*(uxhess(0,1)*dyv(ux)+uyhess(0,1)*dyv(uy));
    auto dyvg22 = 2*(uxhess(1,0)*dxv(ux)+uyhess(1,0)*dxv(uy));
    auto dxvg12 = uxhess(0,0)*dyv(ux)+dxv(ux)*uxhess(0,1) + uyhess(0,0)*dyv(uy)+dxv(uy)*uyhess(0,1);
    auto dyvg12 = uxhess(0,1)*dyv(ux)+dxv(ux)*uxhess(1,1) + uyhess(0,1)*dyv(uy)+dxv(uy)*uyhess(1,1);
    form2( _test=Xh, _trial=Xh, _matrix=A ) +=
        integrate( _range=elements(mesh),
                   _expr=
                   //comp 1
                   - gradt(u)(0,0)*(dxvg11*id(u)(0,0) + g11*grad(u)(0,0) )
                   + gradt(u)(0,0)*(dyvg12*id(u)(0,0) + g12*grad(u)(0,1) )
                   + gradt(u)(0,1)*(dxvg12*id(u)(0,0) + g12*grad(u)(0,0) )
                   - gradt(u)(0,1)*(dyvg22*id(u)(0,0) + g22*grad(u)(0,1) )
                   //comp 2
                   - gradt(u)(1,0)*(dxvg11*id(u)(1,0) + g11*grad(u)(1,0) )
                   + gradt(u)(1,0)*(dyvg12*id(u)(1,0) + g12*grad(u)(1,1) )
                   + gradt(u)(1,1)*(dxvg12*id(u)(1,0) + g12*grad(u)(1,0) )
                   - gradt(u)(1,1)*(dyvg22*id(u)(1,0) + g22*grad(u)(1,1) ) );
    double telapsed = toc("HOLA");
    std::cout << "telapsed="<<telapsed<<"\n";
#endif


#if 0
    // update markedfaces free
    if ( this->hasFlagSet("free" ) )
    {
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            integrate( _range=markedfaces( mesh, this->flagSet("free") ),
                       //_range=boundaryfaces(mesh),
                       _expr=
                       // comp 1
                       + gradt(u)(0,0)*vf::Nx()*idv(proj_alpha)*id(u)(0,0)
                       //- gradt(u)(0,0)*vf::Nx()*idv(proj_gamma)*id(u)(0)
                       + gradt(u)(0,1)*vf::Ny()*idv(proj_gamma)*id(u)(0,0)
                       - gradt(u)(0,1)*vf::Nx()*idv(proj_beta)*id(u)(0,0)
                       - gradt(u)(0,0)*vf::Ny()*idv(proj_beta)*id(u)(0,0)
                       // comp 2
                       + gradt(u)(1,0)*vf::Nx()*idv(proj_alpha)*id(u)(1,0)
                       //- gradt(u)(1,0)*vf::Nx()*idv(proj_gamma)*id(u)(1)
                       + gradt(u)(1,1)*vf::Ny()*idv(proj_gamma)*id(u)(1,0)
                       - gradt(u)(1,1)*vf::Nx()*idv(proj_beta)*id(u)(1,0)
                       - gradt(u)(1,0)*vf::Ny()*idv(proj_beta)*id(u)(1,0)
                       );
    }
#endif

#endif
}
template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateLinearPDE( DataUpdateLinear & data, mpl::int_<3> /**/ ) const
{
    CHECK( false ) << "TODO";
}

#endif

} // namespace FeelModels
} // namespace Feel
