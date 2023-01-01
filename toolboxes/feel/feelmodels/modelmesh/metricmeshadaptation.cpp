#include <feel/feelmodels/modelmesh/metricmeshadaptation.hpp>

#include <feel/feeldiscr/projector.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelmodels/modelvf/eigendecomposition.hpp>

namespace Feel
{
namespace FeelModels
{

template<typename ConvexType>
MetricMeshAdaptation<ConvexType>::MetricMeshAdaptation( space_vectorial_ptrtype space, std::string const& prefix )
    :
    M_vectorialSpace( space )
{}

template<typename ConvexType>
void
MetricMeshAdaptation<ConvexType>::init()
{
    M_scalarMetric = M_vectorialSpace->compSpace()->elementPtr();
}






namespace detail
{

template<typename ExprType>
auto
myHessianFromSymbolicExpr( ExprType const& e,  mpl::int_<2> /**/ )
{
    //std::cout << "myHessianFromSymbolicExpr\n";
#if 1
    auto hessian_xx = diff( e, "x", 2 );
    auto hessian_yy = diff( e, "y", 2 );
    auto hessian_xy = diff( diff( e, "x", 1 ), "y", 1);
    auto hessian_yx = hessian_xy;
    return mat<2,2>( hessian_xx,hessian_xy,
                     hessian_yx,hessian_yy );
#else
    return eye<2,2>()*eye<2,2>();// not compile if only one eye
#endif
}
template<typename ExprType>
auto
myHessianFromSymbolicExpr( ExprType const& e,  mpl::int_<3> /**/ )
{
    auto hessian_xx = diff( e, "x", 2 );
    auto hessian_xy = diff( diff( e, "x", 1 ), "y", 1);
    auto hessian_xz = diff( diff( e, "x", 1 ), "z", 1);
    auto hessian_yx = hessian_xy;
    auto hessian_yy = diff( e, "y", 2 );
    auto hessian_yz = diff( diff( e, "y", 1 ), "z", 1);
    auto hessian_zx = hessian_xz;
    auto hessian_zy = hessian_yz;
    auto hessian_zz = diff( e, "z", 2 );
    return mat<3,3>( hessian_xx,hessian_xy,hessian_xz,
                     hessian_yx,hessian_yy,hessian_yz,
                     hessian_zx,hessian_zy,hessian_zz );
}


template<typename FieldType, typename VectorialSpaceType, typename VectorialElementType>
auto
myHessianFromField( FieldType const& u, VectorialSpaceType const& Vh, std::vector<VectorialElementType> & tmpField, mpl::int_<2> /**/ )
{
    //std::cout << "myHessianFromField\n";
    auto & grad_dxu_P0toP1 = *tmpField[0];
    auto & grad_dyu_P0toP1 = *tmpField[1];

    if ( true ) // fast
    {
        auto sumMeasElt = sum( Vh, one()*meas() );
        auto gradu_P0toP1 = div( sum( Vh, trans(gradv(u))*meas() ), sumMeasElt );
        grad_dxu_P0toP1 = div( sum( Vh, trans(gradv( gradu_P0toP1.comp(Component::X) ))*meas() ), sumMeasElt );
        grad_dyu_P0toP1 = div( sum( Vh, trans(gradv( gradu_P0toP1.comp(Component::Y) ))*meas() ), sumMeasElt );
    }
    else
    {
         auto opP = opProjection(_domainSpace=Vh,
                                 _imageSpace=Vh,
                                 //_backend=backend_type::build( soption( _name="backend" ), prefixvm(this->prefix(),"l2proj"), this->worldCommPtr() ),
                                 _type=Feel::L2 );
         auto gradu_proj = opP->operator()(trans(gradv(u)));
         grad_dxu_P0toP1 = opP->operator()(trans(gradv(gradu_proj.comp(Component::X))));
         grad_dyu_P0toP1 = opP->operator()(trans(gradv(gradu_proj.comp(Component::Y))));
    }

#if 0
    auto hessian_xx_field = grad_dxu_P0toP1.comp(Component::X);
    auto hessian_xy_field = grad_dxu_P0toP1.comp(Component::Y);
    auto hessian_yy_field = grad_dyu_P0toP1.comp(Component::Y);
    auto hessian_xx = idv(hessian_xx_field);
    auto hessian_xy = idv(hessian_xy_field);
    auto hessian_yx = hessian_xy;
    auto hessian_yy = idv(hessian_yy_field);
#else
    auto hessian_xx = idv(grad_dxu_P0toP1)(0);
    auto hessian_xy = idv(grad_dxu_P0toP1)(1);
    auto hessian_yx = hessian_xy;
    auto hessian_yy = idv(grad_dyu_P0toP1)(1);
#endif
    return mat<2,2>( hessian_xx,hessian_xy,
                     hessian_yx,hessian_yy );
}
template<typename FieldType, typename VectorialSpaceType, typename VectorialElementType>
auto
myHessianFromField( FieldType const& u, VectorialSpaceType const& Vh, std::vector<VectorialElementType> & tmpField, mpl::int_<3> /**/ )
{
    auto sumMeasElt = sum( Vh, one()*meas() );
    auto gradu_P0toP1 = div( sum( Vh, trans(gradv(u))*meas() ), sumMeasElt );
    auto & grad_dxu_P0toP1 = *tmpField[0];
    auto & grad_dyu_P0toP1 = *tmpField[1];
    auto & grad_dzu_P0toP1 = *tmpField[2];
    grad_dxu_P0toP1 = div( sum( Vh, trans(gradv( gradu_P0toP1.comp(Component::X) ))*meas() ), sumMeasElt );
    grad_dyu_P0toP1 = div( sum( Vh, trans(gradv( gradu_P0toP1.comp(Component::Y) ))*meas() ), sumMeasElt );
    grad_dzu_P0toP1 = div( sum( Vh, trans(gradv( gradu_P0toP1.comp(Component::Z) ))*meas() ), sumMeasElt );
    auto hessian_xx = idv(grad_dxu_P0toP1)(0);
    auto hessian_xy = idv(grad_dxu_P0toP1)(1);
    auto hessian_xz = idv(grad_dxu_P0toP1)(2);
    auto hessian_yx = hessian_xy;
    auto hessian_yy = idv(grad_dyu_P0toP1)(1);
    auto hessian_yz = idv(grad_dyu_P0toP1)(2);
    auto hessian_zx = hessian_xz;
    auto hessian_zy = hessian_yz;
    auto hessian_zz = idv(grad_dzu_P0toP1)(2);
    return mat<3,3>( hessian_xx,hessian_xy,hessian_xz,
                     hessian_yx,hessian_yy,hessian_yz,
                     hessian_zx,hessian_zy,hessian_zz );
}


template <typename MetricMeshAdaptationType, typename ExprHessianType>
void
updateFromHessian( MetricMeshAdaptationType & mma, ExprHessianType const& eHessian )
{
    if ( true/*M_useMeshAdapationScalar*/ )
    {
        int k=1;// k=1 ->hessian, k=0->grad
        int m=0; // m=0-> L2, m=1->semi H1
        int n=2;//dimension

        double exposant1 = 1.0*n/(2*(k+1-m)+n); // warning 1.0 and not 2.0 because we compute square of norm
        double exposant2 = (2.0*(k+1-m)+n)/n;

        auto mesh = mma.mesh();
        double alpha = integrate(_range=elements(mesh),_expr=pow(inner(eHessian),cst(exposant1) )).evaluate()(0,0);
        alpha /= mesh->measure();
        alpha = std::pow( alpha,exposant2);

        if ( alpha > 1e-8 )
        {
            double exposant3 = 2.0/(2*(k+1-m)+n);
            auto exprMetricScalar = pow( cst(1.)+(1./alpha)*inner(eHessian), cst(exposant3));
            mma.updateScalarMetric( exprMetricScalar );
        }
        else
            mma.updateScalarMetric( cst(1.) );
    }
    else
    {
#if 0
        auto eHessionModifed = eigenDecomposition( eHessian );
        auto Id = eye<mesh_type::nDim,mesh_type::nDim>();
        double alpha = 1;
        M_weightFunctionTensor2->on(_range=elements(M_Xh->mesh()),
                                    //_expr=eHessionModifed );
                                    //_expr= /*inv*/(  pow(det(eHessionModifed),cst(-1./(mesh_type::nDim+4)) )*eHessionModifed ) );
                                    _expr=pow( det( Id + (1./alpha)*eHessionModifed ),cst(-1./6) )*(Id + (1./alpha)*eHessionModifed ) );
#endif
    }
}

}

template<typename ConvexType>
void
MetricMeshAdaptation<ConvexType>::update( Expr<GinacExVF<2>> const& uScal )
{
    auto eHessian =  Feel::FeelModels::detail::myHessianFromSymbolicExpr( uScal, mpl::int_<mesh_type::nDim>() );
    Feel::FeelModels::detail::updateFromHessian( *this, eHessian );
}
template<typename ConvexType>
void
MetricMeshAdaptation<ConvexType>::update( element_scalar_type const& uScal )
{
    std::vector<element_vectorial_ptrtype> tmpField( mesh_type::nDim ); // need to store this field else the expression crash
    for ( int k=0;k<tmpField.size();++k )
        tmpField[k] = M_vectorialSpace->elementPtr();
    auto eHessian =  Feel::FeelModels::detail::myHessianFromField( uScal, M_vectorialSpace, tmpField, mpl::int_<mesh_type::nDim>() );
    Feel::FeelModels::detail::updateFromHessian( *this, eHessian );
}

#if 0
template<typename ConvexType>
void
MetricMeshAdaptation<ConvexType>::update( element_scalar_type const& uScal )
{
    if ( true/*M_useMeshAdapationScalar*/ )
    {
        int k=1;// k=1 ->hessian, k=0->grad
        int m=0; // m=0-> L2, m=1->semi H1
        int n=2;//dimension

        double exposant1 = 1.0*n/(2*(k+1-m)+n); // warning 1.0 and not 2.0 because we compute square of norm
        double exposant2 = (2.0*(k+1-m)+n)/n;

        //auto eScal = expr( soption(_name="mesh-adaptation.scalar-weight.expr", _prefix=this->prefix() ) );
#if 0
        auto eHessian = myHessianFromSymbolicExpr( eScal, mpl::int_<mesh_type::nDim>() );
#else
        //auto uScal = M_Xh->compSpace()->elementPtr();
        //uScal->on(_range=elements(M_Xh->mesh()),_expr=eScal );
        std::vector<element_vectorial_ptrtype> tmpField( mesh_type::nDim ); // need to store this field else the expression crash
        for ( int k=0;k<tmpField.size();++k )
            tmpField[k] = M_vectorialSpace->elementPtr();
        auto eHessian = myHessianFromField( uScal, M_vectorialSpace, tmpField, mpl::int_<mesh_type::nDim>() );
#endif
        auto mesh = M_vectorialSpace->mesh();
        double alpha = integrate(_range=elements(mesh),_expr=pow(inner(eHessian),cst(exposant1) )).evaluate()(0,0);
        alpha /= mesh->measure();
        alpha = std::pow( alpha,exposant2);

        if ( alpha > 1e-8 )
        {
            double exposant3 = 2.0/(2*(k+1-m)+n);
            auto exprMetricScalar = pow( cst(1.)+(1./alpha)*inner(eHessian), cst(exposant3));
            M_scalarMetric->on(_range=elements(mesh),_expr=exprMetricScalar);
        }
        else
            M_scalarMetric->on(_range=elements(mesh),_expr=cst(1.));
    }
    else
    {
#if 0
        auto eScal = expr( soption(_name="mesh-adaptation.scalar-weight.expr", _prefix=this->prefix() ) );
        auto eHessian = myHessianFromSymbolicExpr( eScal, mpl::int_<mesh_type::nDim>() );
        auto eHessionModifed = eigenDecomposition( eHessian );
        auto Id = eye<mesh_type::nDim,mesh_type::nDim>();
        double alpha = 1;
        M_weightFunctionTensor2->on(_range=elements(M_Xh->mesh()),
                                    //_expr=eHessionModifed );
                                    //_expr= /*inv*/(  pow(det(eHessionModifed),cst(-1./(mesh_type::nDim+4)) )*eHessionModifed ) );
                                    _expr=pow( det( Id + (1./alpha)*eHessionModifed ),cst(-1./6) )*(Id + (1./alpha)*eHessionModifed ) );
#endif
    }

}
#endif





#if 0
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
#endif



} // namespace FeelModels
} // namespace Feel


