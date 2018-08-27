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

//#include <feel/feelalg/aitken.hpp>

#include <feel/feelvf/vf.hpp>


#define FSI_WINSLOW_USE_OPT_EXPR 0

namespace Feel
{
namespace FeelModels
{


template< typename MeshType, int Order >
Winslow<MeshType,Order>::Winslow( mesh_ptrtype mesh, std::string const& prefix, worldcomm_ptr_t const& worldcomm,
                                  bool useGhostEltFromExtendedStencil, ModelBaseRepository const& modelRep )
    :
    super_type( prefix, worldcomm,"", modelRep ),
    M_backend( backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr() ) ),
    M_solverType( soption(_prefix=this->prefix(),_name="solver") ),
    M_mesh(mesh),
    M_flagSet(/*flagSet*/),
    M_Xh(space_type::New(_mesh=mesh,_worldscomm=makeWorldsComm(1,this->worldCommPtr()),
                         _extended_doftable=std::vector<bool>(1,useGhostEltFromExtendedStencil) )),
    M_displacement( new element_type(M_Xh) ),
    M_displacementOld( new element_type(M_Xh) ),
    M_dispImposedOnBoundary( new element_type(M_Xh) ),
    M_identity( new element_type(M_Xh) ),
    M_XhScalM1(space_scal_m1_type::New(_mesh=mesh,_worldscomm=makeWorldsComm(1,this->worldCommPtr()))),
    M_XhScalP0Disc(space_p0_type::New(_mesh=mesh,_worldscomm=makeWorldsComm(1,this->worldCommPtr()))),
    M_l2projector(opProjection(_domainSpace=this->functionSpaceScalM1(),
                               _imageSpace=this->functionSpaceScalM1(),
                               _backend=backend_type::build( soption( _name="backend" ), prefixvm(this->prefix(),"l2proj"), this->worldCommPtr() ),
                               _type=Feel::L2 ) )
{
    EntityProcessType entityProcess = (useGhostEltFromExtendedStencil)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    *M_identity = vf::project( _space=M_Xh, _range=elements(M_mesh,entityProcess), _expr=vf::P() );
}

//----------------------------------------------------------------------------//

template< typename MeshType, int Order >
Winslow<MeshType,Order>::Winslow( space_ptrtype const& space, std::string const& prefix,
                                  ModelBaseRepository const& modelRep )
    :
    super_type( prefix, space->worldCommPtr(),"",modelRep ),
    M_backend( backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr() ) ),
    M_solverType( soption(_prefix=this->prefix(),_name="solver") ),
    M_mesh(space->mesh()),
    M_flagSet(/*flagSet*/),
    M_Xh(space),
    M_displacement( new element_type(M_Xh) ),
    M_displacementOld( new element_type(M_Xh) ),
    M_dispImposedOnBoundary( new element_type(M_Xh) ),
    M_identity( new element_type(M_Xh) ),
    M_XhScalM1(space_scal_m1_type::New(_mesh=M_mesh,_worldscomm=makeWorldsComm(1,this->worldCommPtr()))),
    M_XhScalP0Disc(space_p0_type::New(_mesh=M_mesh,_worldscomm=makeWorldsComm(1,this->worldCommPtr()))),
    M_l2projector(opProjection(_domainSpace=this->functionSpaceScalM1(),
                               _imageSpace=this->functionSpaceScalM1(),
                               _backend=backend_type::build( soption( _name="backend" ), prefixvm(this->prefix(),"l2proj"), this->worldCommPtr() ),
                               _type=Feel::L2 ) )
{
    *M_identity = vf::project( _space=M_Xh, _range=elements(M_mesh), _expr=vf::P() );
}

//----------------------------------------------------------------------------//
template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::init()
{
    this->log("Winslow","init", "start" );

    this->setRebuildCstPartInLinearSystem(false);
    this->setUseLinearJacobianInResidual(false);
    this->setRebuildLinearPartInJacobian(false);
    this->setUseCstMatrix(false);
    this->setUseCstVector(false);

    auto graph = stencil( _test=M_Xh, _trial=M_Xh )->graph();
    M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),Feel::backend(_rebuild=true,_name=this->backend()->prefix())/*this->backend()*/,
                                                                graph, graph->mapRow().indexSplit() ) );
    // mesh not move so not rebuild cst part
    M_vectorSolution = this->backend()->newVector( this->functionSpace() );
    *M_vectorSolution = *M_identity;

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
    auto const& u = *M_displacement;
    form2( _trial=Xh, _test=Xh,_matrix=M_matrixMetricDerivative)
        += integrate(_range=elements(mesh),
                     _expr=inner(idt(u),id(u)) );
    M_matrixMetricDerivative->close();

    M_useMeshAdapation = boption(_name="mesh-adaptation",_prefix=this->prefix() );
    M_useMeshAdapationScalar = boption(_name="mesh-adaptation.scalar-weight",_prefix=this->prefix() );

    M_weightFunctionScalar.reset( new element_p0_type( M_XhScalP0Disc ) );
    M_weightFunctionScalar->on(_range=elements(mesh),_expr=cst(1.));
    if ( M_useMeshAdapation && !M_useMeshAdapationScalar )
    {
        M_XhTensor2P0Disc = space_p0_tensor2_type::New(_mesh=mesh);
        M_weightFunctionTensor2.reset( new element_p0_tensor2_type( M_XhTensor2P0Disc ) );
        M_weightFunctionTensor2->on(_range=elements(mesh),_expr=Id<mesh_type::nDim>());
    }

#if 0
    M_hMinRadius.reset( new element_p0_type(M_XhScalP0Disc));
    for ( auto const& eltWrap : elements(mesh) )
    {
        double dmin = 1e10;
        auto const& elt = unwrap_ref( eltWrap );
        auto bary = elt.barycenter();
        em_node_type<double> bary2( bary.data().begin(), bary.size() );
        for( uint16_type p = 0; p < elt.numVertices; ++p ) {
            auto const& pt = elt.point( p ).node();
            em_node_type<double> pt2( const_cast<double *>(pt.data().begin()), pt.size() );
            dmin = std::min(dmin,(bary2-pt2).norm());
        }
        for ( auto const& ldof : M_XhScalP0Disc->dof()->localDof( elt.id() ) )
        {
            M_hMinRadius->set( ldof.second.index(), dmin );
        }
    }
#endif
    this->log("Winslow","init", "finish" );
}

template< typename MeshType, int Order >
std::shared_ptr<std::ostringstream>
Winslow<MeshType,Order>::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );

    *_ostr << "\n   Winslow "
           << "\n     -- solver type : " << M_solverType
           <<  "\n    -- mesh adaptation : " << M_useMeshAdapation;
    if ( M_useMeshAdapation )
        *_ostr <<  "\n       + type : " << ( ( M_useMeshAdapationScalar)? "Scalar" : "Tensor2" );

    return _ostr;
}

//----------------------------------------------------------------------------//

template< typename MeshType,int Order >
void
Winslow<MeshType,Order>::solve()
{
    this->log("Winslow","solve", "start" );

    if ( M_useMeshAdapation )
        this->updateMeshAdaptation();

    if ( M_solverType=="Picard" || M_solverType == "FixPoint" )
    {
        //this->updateNewtonInitialGuess( M_vectorSolution );
        M_algebraicFactory->solvePicard( M_vectorSolution );
    }
    else if ( M_solverType=="Newton" )
    {
        M_algebraicFactory->solveNewton( M_vectorSolution );
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
        M_algebraicFactory->solveNewton( M_vectorSolution );
    }

    //update displacement
    *M_displacement = *M_vectorSolution;
    *M_displacement -= *M_identity;

    this->log("Winslow","solve", "finish" );
}


//----------------------------------------------------------------------------//

template< typename MeshType,int Order >
void
Winslow<MeshType,Order>::updateMeshAdaptation()
{
#if 0
    auto curMapBB = M_Xh->element(M_vectorSolution);
    auto curDispBB = M_Xh->element( idv(curMapBB)-idv(M_identity) );
    auto curMinRadius = M_XhScalP0Disc->element();
    auto newBary = M_XhTensor2P0Disc->element();
    auto newBaryx = newBary.comp(Component::X,Component::X);
    auto newBaryy = newBary.comp(Component::Y,Component::Y);
    newBaryx.on(_range=elements(M_Xh->mesh()),_expr=Cx()+idv(curDispBB)(0,0));
    newBaryy.on(_range=elements(M_Xh->mesh()),_expr=Cy()+idv(curDispBB)(1,0));
    for ( auto const& eltWrap : elements(M_Xh->mesh()) )
    {
        double dmin = 1e10;
        auto const& elt = unwrap_ref( eltWrap );
        size_type index = M_XhScalP0Disc->dof()->localToGlobal( elt.id(),0,0 ).index();
        double baryX = newBaryx(index);double baryY = newBaryy(index);
        for( uint16_type p = 0; p < elt.numVertices; ++p ) {
            double ptX = curMapBB( M_Xh->dof()->localToGlobal( elt.id(),p,0 ).index() );
            double ptY = curMapBB( M_Xh->dof()->localToGlobal( elt.id(),p,1 ).index() );
            dmin = std::min(dmin,std::sqrt(std::pow(ptX-baryX,2)+std::pow(ptY-baryY,2)));
        }
        curMinRadius.set( index, dmin );
    }
#endif


#if 0
    auto wfxx = M_weightFunction->comp(Component::X,Component::X);
    auto wfyy = M_weightFunction->comp(Component::Y,Component::Y);
    auto exprg = expr(soption(_name="functions.g") );
    auto exprh = expr(soption(_name="functions.h") );
    std::map<std::string,double> mp;mp["k"]=M_weigthFunctionScaling;
    exprg.setParameterValues( mp );
    exprh.setParameterValues( mp );

    wfxx.on(_range=elements(M_Xh->mesh()),_expr=exprg );
    wfyy.on(_range=elements(M_Xh->mesh()),_expr=exprh );
#elif 1
    auto Xh = this->functionSpace();
    auto mesh = this->functionSpace()->mesh();
    auto const& u = *M_displacement;
    auto aa = form2( _trial=Xh, _test=Xh);
    aa = integrate(_range=elements(mesh),
                  _expr=inner(gradt(u),grad(u)) );
    auto ll = form1( _test=Xh );
    ll = integrate(_range=elements(mesh),
                   _expr=0*inner(one(),id(u)));
    if ( this->hasFlagSet("moving" ) )
    {
        aa +=
            on( _range=markedfaces( mesh, this->flagSet("moving") ),
                _element=u, _rhs=ll,
                _expr=idv( this->dispImposedOnBoundary() ) );
    }
    if ( this->hasFlagSet("fixed" ) )
    {
        aa +=
            on( _range=markedfaces(mesh, this->flagSet("fixed") ),
                _element=u, _rhs=ll,
                _expr=0*idv(this->identity())/*, ON_ELIMINATION*/ );

    }
    auto uHarmonic = Xh->element();
    aa.solve(_rhs=ll,_solution=uHarmonic);

    auto uHarmonicMagnitude = Xh->compSpace()->element( sqrt(inner(idv(uHarmonic)) ));
#if 0
    auto wfx = (*M_weightFunction)[Component::X];
    auto Volume = integrate( elements(this->functionSpace()->mesh()), cst(1.) ).broken( XhP0 );
    //auto Volume = integrate( elements(this->functionSpace()->mesh()), det( gradv(u) ) ).broken( XhP0 );
    double Vmin = Volume.min();
    double Vmax = Volume.max();
#endif

#if 0
    double gf=doption(_name="parameters.a");//1.4;
    double nLayers=doption(_name="parameters.b");//2;
    //auto dispMag = sqrt(inner(idv(uHarmonic)));
    double uMax = uHarmonicMagnitude.max();
    double uMin = uMax-uMax/2.;
    //auto myweight = max(cst(1.),pow(gf,nLayers*idv(uHarmonicMagnitude)/uMax) );
    auto myweight = max(cst(1.),pow(cst(gf),nLayers*(idv(uHarmonicMagnitude)-uMin)/(uMax-uMin) ) );
    M_weightFunctionScalar->on(_range=elements(M_Xh->mesh()),_expr=myweight);
#endif


#if 1
#if 1
    //auto curDisp = uHarmonic;//uHarmonicMagnitude;
    auto curDisp = Xh->element();
    curDisp.on(_range=elements(mesh),_expr=idv(uHarmonic)*exp(inner(idv(uHarmonic))));
    //double newArea = integrate(_range=elements(mesh),_expr=det(Id<mesh_type::nDim>()+gradv(uHarmonic)) ).evaluate()(0,0);
#else
    auto curMap = M_Xh->element(M_vectorSolution);
    //auto uHarmonicMagnitude = M_Xh->compSpace()->element( sqrt(inner(idv(curMap)-idv(M_identity)) ));
    auto curDisp = M_Xh->element( idv(curMap)-idv(M_identity) );
#endif

    double alpha = integrate(_range=elements(M_Xh->mesh()),_expr=pow(inner(gradv(curDisp)),cst(1./2.) )).evaluate()(0,0);
    alpha /= M_Xh->mesh()->measure();
    //alpha /= newArea;
    alpha = std::pow( alpha,2);
    //alpha *= doption(_name="parameters.a");
    //std::cout << "\n\nalpha="<<alpha<<"\n\n";
    auto sqrtofg = sqrt(cst(1.)+ (1./alpha)*inner(gradv(curDisp)));

    if ( M_useMeshAdapationScalar )
    {
        if ( alpha > 1e-8 )
        {
            //auto newScaling = idv(curMinRadius)/idv(M_hMinRadius);
            auto scalarMatrixCoeff = sqrt(cst(1.)+(1./alpha)*inner(gradv(curDisp)));
            M_weightFunctionScalar->on(_range=elements(M_Xh->mesh()),_expr=scalarMatrixCoeff);
        }
        else
            M_weightFunctionScalar->on(_range=elements(M_Xh->mesh()),_expr=cst(1.));
    }
    else
    {
        if ( alpha > 1e-8 )
        {
            auto blabla = Id<mesh_type::nDim>()+(1./alpha)*(gradv(curDisp)*trans(gradv(curDisp)));
            auto blablaRootSquare=(1./sqrt(trace(blabla)+2*sqrt(det(blabla))))*(blabla+sqrt(det(blabla))*Id<mesh_type::nDim>());
            M_weightFunctionTensor2->on(_range=elements(M_Xh->mesh()),_expr=(1./sqrt(det(blablaRootSquare)))*blablaRootSquare);
        }
        else
            M_weightFunctionTensor2->on(_range=elements(M_Xh->mesh()),_expr=Id<mesh_type::nDim>());
    }
#endif
#endif
}

//----------------------------------------------------------------------------//

template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateNewtonInitialGuess(vector_ptrtype& U) const
{
    this->log("Winslow","updateNewtonInitialGuess", "start" );

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
Winslow<MeshType,Order>::updateJacobian( DataUpdateJacobian & data ) const
{
    this->updateJacobian( data, mpl::int_<mesh_type::nDim>() );
}

template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateJacobian( DataUpdateJacobian & data, mpl::int_<2> /**/ ) const
{
#if 1
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

template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateResidual( DataUpdateResidual & data ) const
{
    this->updateResidual( data, mpl::int_<mesh_type::nDim>() );
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
#if 0
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

template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateLinearPDE( DataUpdateLinear & data ) const
{
    if ( false )
    {
        this->updateLinearPDE( data, mpl::int_<mesh_type::nDim>() );
        return;
    }


    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto Xh = this->functionSpace();
    auto mesh = this->functionSpace()->mesh();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;

    const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
    auto u = Xh->element( vecCurrentPicardSolution );
    auto G = trans(gradv(u))*gradv(u);
    auto invG = inv(G);

    M_vectorMetricDerivative->zero();
    auto lf = form1( _test=Xh,_vector=M_vectorMetricDerivative );
    lf += integrate(_range=elements(mesh),
                    _expr=inner(invG, grad(u)) );
    lf += integrate(_range=boundaryfaces(mesh),
                    _expr=-inner(invG*N(), id(u)) );
    //M_vectorMetricDerivative->close();
    M_backendMetricDerivative->solve(_matrix=M_matrixMetricDerivative,_solution=*M_fieldMetricDerivative,_rhs=M_vectorMetricDerivative);


    if ( M_useMeshAdapation && !M_useMeshAdapationScalar )
    {
        auto tau = idv(M_weightFunctionTensor2);
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            integrate( _range=elements(mesh),
                       _expr= inner(invG*trans(tau*gradt(u)),trans(grad(u))) );
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            integrate( _range=elements(mesh),
                   _expr=-inner(tau*gradt(u)*idv(M_fieldMetricDerivative),id(u)) );
    }
    else
    {
        auto tau = idv( M_weightFunctionScalar );
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            integrate( _range=elements(mesh),
                       _expr= (tau)*inner(invG*trans(gradt(u)),trans(grad(u))) );
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            integrate( _range=elements(mesh),
                       _expr=-(tau)*inner(gradt(u)*idv(M_fieldMetricDerivative),id(u)) );
    }
}

template< typename MeshType, int Order >
void
Winslow<MeshType,Order>::updateLinearPDE( DataUpdateLinear & data, mpl::int_<2> /**/ ) const
{
#if 0
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

} // namespace FeelModels
} // namespace Feel
