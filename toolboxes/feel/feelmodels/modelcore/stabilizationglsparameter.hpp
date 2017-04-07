
#ifndef FEELPP_MODELS_STABILIZATIONGLSPARAMETER_HPP
#define FEELPP_MODELS_STABILIZATIONGLSPARAMETER_HPP 1

#include <feel/feelalg/matrixeigendense.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>


namespace Feel {

namespace FeelModels
{


template<typename MeshType, uint16_type OrderPoly>
class StabilizationGLSParameter
{
public :
    static const uint16_type nOrder = OrderPoly;
    using mesh_t = MeshType;
    using mesh_ptr_t = boost::shared_ptr<mesh_t>;

    using space_t = Pdh_type<mesh_t,nOrder>;
    using space_ptr_t = boost::shared_ptr<space_t>;
    using element_t = typename space_t::element_type;
    using space_P0_t = Pdh_type<mesh_t,0>;
    using space_P0_ptr_t = boost::shared_ptr<space_P0_t>;
    using element_P0_t = typename space_P0_t::element_type;

    StabilizationGLSParameter( mesh_ptr_t const& mesh, std::string const& prefix )
        :
        M_mesh( mesh ),
        M_spaceP0( Pdh<0>(mesh) ),
        M_Xh( Pdh<nOrder>(mesh) ),
        M_hSize( M_spaceP0, "hSize" ),
        M_lambdaK( M_spaceP0, "lambdak" ),
        M_sqrtLambdaK( M_spaceP0, "sqrtLambdaK"),
        M_mK( M_spaceP0, "mk"),
        M_method( soption(_name="method",_prefix=prefix ) ),
        M_hSizeMethod( soption(_name="hsize.method",_prefix=prefix ) ),
        M_penalLambdaK( doption( _name="eigenvalue.penal-lambdaK",_prefix=prefix ) )
        {
            CHECK( M_method == "eigenvalue" || M_method == "doubly-asymptotic-approximation" )  << "invalid method";
        }

    void init();

    template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
    auto localPeRe( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<0> /**/ ) const;
    template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
    auto localPeRe( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<1> /**/ ) const;
    template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
    auto localPeRe( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<2> /**/ ) const;

    template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
    auto tau( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion ) const
        {
            return this->tau( exprConvection, exprCoeffDiffusion, mpl::int_<0>() );
        }
    template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
    auto tau( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<0> /**/ ) const;
    template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
    auto tau( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<1> /**/ ) const;
    template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
    auto tau( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<2> /**/ ) const;

    template<int ParameterType, typename ExprConvectiontype, typename ExprCoeffDiffusionType>
    auto delta( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion ) const;

    element_P0_t const& hSize() const { return M_hSize; }
    element_P0_t const& sqrtLambdaK() const { return M_sqrtLambdaK; }
    element_P0_t const& lambdaK() const { return M_lambdaK; }
    element_P0_t const& mK() const { return M_mK; }

    std::string const& method() const { return M_method; }
    double penalLambdaK() const { return M_penalLambdaK; }
    void setPenalLambdaK( double val ) { M_penalLambdaK = val; }

private :
    mesh_ptr_t M_mesh;
    space_P0_ptr_t M_spaceP0;
    space_ptr_t M_Xh;
    element_P0_t M_hSize, M_lambdaK, M_sqrtLambdaK, M_mK;
    std::string M_method;
    std::string M_hSizeMethod;
    double M_penalLambdaK;
}; // class StabilizationGLSParameter


template<typename MeshType, uint16_type OrderPoly>
void
StabilizationGLSParameter<MeshType,OrderPoly>::init()
{
    auto mesh = M_mesh;

    if ( M_hSizeMethod == "hmin" )
        M_hSize.on(_range=elements(mesh), _expr=hMin() );
    else if ( M_hSizeMethod == "h" )
        M_hSize.on(_range=elements(mesh), _expr=h() );
    else if ( M_hSizeMethod == "meas" )
        M_hSize.on(_range=elements(mesh), _expr=pow(meas(),1./mesh_t::nDim) );

    // Stabilization : computation of lambdaK
    if ( M_method == "eigenvalue" || M_method == "eigenvalue-simplified" )
    {
        auto Ph = M_Xh;
        auto w = Ph->element();

        auto matA = backend()->newMatrix(_test=Ph,_trial=Ph );
        auto matB = backend()->newMatrix(_test=Ph,_trial=Ph );
        auto fA = form2( _trial=Ph, _test=Ph, _matrix=matA );
        fA = integrate( _range=elements(mesh),
                        _expr=inner(laplacian(w),laplaciant(w)) );
        matA->close();
        auto fB = form2( _trial=Ph, _test=Ph, _matrix=matB );
        fB = integrate( _range=elements(mesh),
                        _expr=inner(grad(w),gradt(w)) );
        if ( math::abs( M_penalLambdaK ) > 1e-12 )
            fB += integrate( elements(mesh),
                             M_penalLambdaK*inner(id(w),idt(w)) );
        matB->close();

        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> localEigenMatA(0,0);
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> localEigenMatB(0,0);

        for ( auto const& eltWrap : elements(M_mesh) )
        {
            auto const& elt = unwrap_ref( eltWrap );
            auto extractIndices = Ph->dof()->getIndices( elt.id() );
            int nLocalIndices = extractIndices.size();
            if ( localEigenMatA.rows() != nLocalIndices)
            {
                localEigenMatA.resize( nLocalIndices,nLocalIndices );
                localEigenMatB.resize( nLocalIndices,nLocalIndices );
            }
            for ( int i=0;i<nLocalIndices;++i )
            {
                for ( int j=0;j<nLocalIndices;++j )
                {
                    localEigenMatA(i,j) = matA->operator()(extractIndices[i],extractIndices[j]);
                    localEigenMatB(i,j) = matB->operator()(extractIndices[i],extractIndices[j]);
                }
            }
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es( localEigenMatA,localEigenMatB );
            Eigen::VectorXd eigen_vls= es.eigenvalues();
            double lambda = eigen_vls.maxCoeff();
            M_lambdaK.on( _range=idedelements( mesh,elt.id() ), _expr=cst(lambda) );
        }

        M_sqrtLambdaK.on(_range=elements(mesh), _expr=sqrt(idv(M_lambdaK)) );
        auto cK = cst(1.)/(idv(M_lambdaK)*pow(idv(M_hSize),2));
        M_mK.on(_range=elements(mesh),_expr=min( cst(1./3.),2*cK ) );
    }
}



template<typename MeshType, uint16_type OrderPoly>
template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
auto
StabilizationGLSParameter<MeshType,OrderPoly>::localPeRe( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<0> /**/ ) const
{
    auto norm_u = norm2( exprConvection );
    auto Re = idv(M_mK)*norm_u*idv(M_hSize)/(2*exprCoeffDiffusion);
    auto xi_Re = min( cst(1.), Re );
    return xi_Re;
}

template<typename MeshType, uint16_type OrderPoly>
template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
auto
StabilizationGLSParameter<MeshType,OrderPoly>::localPeRe( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<1> /**/ ) const
{
    auto norm_u = norm2( exprConvection );
    auto Re = norm_u/(4*idv(M_sqrtLambdaK)*exprCoeffDiffusion);
    auto xi_Re = min( cst(1.), Re );
    return xi_Re;
}

template<typename MeshType, uint16_type OrderPoly>
template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
auto
StabilizationGLSParameter<MeshType,OrderPoly>::localPeRe( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<2> /**/ ) const
{
    auto psi=min(cst(1.),(1./3.)*idv(M_hSize)*norm2(exprConvection)/(2*exprCoeffDiffusion));
    return psi;
}


template<typename MeshType, uint16_type OrderPoly>
template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
auto
StabilizationGLSParameter<MeshType,OrderPoly>::tau( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<0> /**/ ) const
{
    auto norm_u = norm2( exprConvection );
    auto xi_Re = this->localPeRe( exprConvection, exprCoeffDiffusion, mpl::int_<0>() );
    auto tau = idv(M_hSize)*xi_Re/( max( 1e-10, 2*norm_u) );
    return tau;
}

template<typename MeshType, uint16_type OrderPoly>
template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
auto
StabilizationGLSParameter<MeshType,OrderPoly>::tau( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<1> /**/ ) const
{
    auto norm_u = norm2( exprConvection );
    auto xi_Re = this->localPeRe( exprConvection, exprCoeffDiffusion, mpl::int_<1>() );
    auto tau = 2*xi_Re/( max( 1e-10, idv(M_sqrtLambdaK)*norm_u) );
    return tau;
}

template<typename MeshType, uint16_type OrderPoly>
template<typename ExprConvectiontype, typename ExprCoeffDiffusionType>
auto
StabilizationGLSParameter<MeshType,OrderPoly>::tau( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion, mpl::int_<2> /**/ ) const
{
    auto psi = this->localPeRe( exprConvection, exprCoeffDiffusion, mpl::int_<2>() );
#if 0
    auto chiConv = chi(norm2(exprConvection)>cst(1e-5));
    auto tau = (psi*h()/(2*norm2(exprConvection) + (1-chiConv) ) )*chiConv;
#else
    auto tau = psi*idv(M_hSize)/( max( 2*norm2(exprConvection), 1e-10) );
#endif
    return tau;
}

template<typename MeshType, uint16_type OrderPoly>
template<int ParameterType, typename ExprConvectiontype, typename ExprCoeffDiffusionType>
auto
StabilizationGLSParameter<MeshType,OrderPoly>::delta( ExprConvectiontype const& exprConvection, ExprCoeffDiffusionType const& exprCoeffDiffusion ) const
{
    auto norm_u = norm2( exprConvection );
    auto xi_Re = localPeRe( exprConvection, exprCoeffDiffusion, mpl::int_<ParameterType>() );
    auto myDelta = cst(1.)*norm_u*idv(M_hSize)*xi_Re;
    return myDelta;
}

}  // namespace FeelModels

} // namespace Feel

#endif
