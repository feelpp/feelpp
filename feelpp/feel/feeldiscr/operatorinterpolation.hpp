/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-01

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2016 Feel++ Consortium

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
/**
   \file operatorinterpolation.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Chabannes Vincent <vincent.chabannes@imag.fr>
   \date 2008-02-01
 */
#ifndef FEELPP_DISCR_OPERATORINTERPOLATION_H
#define FEELPP_DISCR_OPERATORINTERPOLATION_H 1

#include <feel/feeldiscr/operatorlinear.hpp>

#include <feel/feeldiscr/stencil.hpp>

namespace Feel
{
using conforming_t = mpl::bool_<true>;
using nonconforming_t = mpl::bool_<false>;
//extern constexpr conforming_t conforming;
//extern constexpr nonconforming_t nonconforming;
template<typename C>
constexpr bool is_conforming = C::value;

enum class interpolation_operand_type { ID = 0, GRADIENT, CURL, DIV, EXPR };

inline std::ostream&
operator<<( std::ostream& os, interpolation_operand_type const& o )
{
    if ( o == interpolation_operand_type::ID )
        os << "ID";
    if ( o == interpolation_operand_type::GRADIENT )
        os << "GRADIENT";
    if ( o == interpolation_operand_type::CURL )
        os << "CURL";
    if ( o == interpolation_operand_type::DIV )
        os << "DIV";
    if ( o == interpolation_operand_type::EXPR )
        os << "EXPR";
    return os;
}


template<bool C, interpolation_operand_type O>
class InterpolationTypeBaseBase
{
public :
    static constexpr bool is_conforming = C;
    static constexpr interpolation_operand_type interpolation_operand = O;
    using operand_t = mpl::int_<static_cast<int>( interpolation_operand )>;

    constexpr InterpolationTypeBaseBase()
        :
        M_searchWithCommunication( true ),
        M_componentsAreSamePoint( true ),
        M_onlyLocalizeOnBoundary( false ),
        M_nbNearNeighborInKdTree( 15)
        {}
    constexpr InterpolationTypeBaseBase( bool useComm,
                                         bool compAreSamePt=true,
                                         bool onlyLocalizeOnBoundary=false,
                                         int nbNearNeighborInKdTree=15 )
    :
        M_searchWithCommunication( useComm ),
        M_componentsAreSamePoint( compAreSamePt ),
        M_onlyLocalizeOnBoundary( onlyLocalizeOnBoundary ),
        M_nbNearNeighborInKdTree( nbNearNeighborInKdTree )
    {}

    InterpolationTypeBaseBase( InterpolationTypeBaseBase const& a) = default;
    InterpolationTypeBaseBase( InterpolationTypeBaseBase && a) = default;

    InterpolationTypeBaseBase& operator=( InterpolationTypeBaseBase const& a) = default;
    InterpolationTypeBaseBase& operator=( InterpolationTypeBaseBase && a) = default;

    constexpr bool searchWithCommunication() const noexcept { return M_searchWithCommunication; }
    constexpr bool componentsAreSamePoint() const noexcept { return M_componentsAreSamePoint; }
    constexpr bool onlyLocalizeOnBoundary() const noexcept { return M_onlyLocalizeOnBoundary; }
    constexpr int nbNearNeighborInKdTree() const noexcept { return M_nbNearNeighborInKdTree; }

    static constexpr bool isConforming() noexcept { return is_conforming; }
    static constexpr interpolation_operand_type interpolationOperand() noexcept { return interpolation_operand; }

private :

    bool M_searchWithCommunication;
    bool M_componentsAreSamePoint;
    bool M_onlyLocalizeOnBoundary;
    int M_nbNearNeighborInKdTree;
};


template<bool C, interpolation_operand_type O >
class InterpolationTypeBase : public InterpolationTypeBaseBase<C,O>
{
    using super_type = InterpolationTypeBaseBase<C,O>;
public :
    static constexpr bool is_conforming = super_type::is_conforming;
    static constexpr interpolation_operand_type interpolation_operand = super_type::interpolation_operand;
    using operand_t = typename super_type::operand_t; //mpl::int_<static_cast<int>( interpolation_operand )>;
    using id_t = mpl::int_<static_cast<int>( interpolation_operand_type::ID )>;
    using gradient_t = mpl::int_<static_cast<int>( interpolation_operand_type::GRADIENT )>;
    using curl_t = mpl::int_<static_cast<int>( interpolation_operand_type::CURL )>;
    using div_t = mpl::int_<static_cast<int>( interpolation_operand_type::DIV )>;

    constexpr InterpolationTypeBase() = default;
    constexpr InterpolationTypeBase( bool useComm,
                                     bool compAreSamePt=true,
                                     bool onlyLocalizeOnBoundary=false,
                                     int nbNearNeighborInKdTree=15 )
        :
        super_type( useComm,compAreSamePt,onlyLocalizeOnBoundary,nbNearNeighborInKdTree )
    {}

    InterpolationTypeBase( InterpolationTypeBase const& a) = default;
    InterpolationTypeBase( InterpolationTypeBase && a) = default;

    InterpolationTypeBase& operator=( InterpolationTypeBase const& a) = default;
    InterpolationTypeBase& operator=( InterpolationTypeBase && a) = default;

    static constexpr bool requiresDerivatives() noexcept
        {
            return interpolation_operand == interpolation_operand_type::GRADIENT ||
                interpolation_operand == interpolation_operand_type::CURL ||
                interpolation_operand == interpolation_operand_type::DIV;
        }

    template<typename Elt>
    static constexpr decltype(auto) operand( Elt&& e )
        {
            return operand( std::forward<Elt>(e), operand_t() );
        }
private:
    template<typename Elt>
    static constexpr decltype(auto) operand( Elt&& e, id_t )
        {
            return vf::id( std::forward<Elt>(e) );
        }
    template<typename Elt>
    static constexpr decltype(auto) operand( Elt&& e, gradient_t )
        {
            return trans(vf::grad( std::forward<Elt>(e) ));
        }
    template<typename Elt>
    static constexpr decltype(auto) operand( Elt&& e, curl_t )
        {
            return vf::curl( std::forward<Elt>(e) );
        }
    template<typename Elt>
    static constexpr decltype(auto) operand( Elt&& e, div_t )
        {
            return vf::div( std::forward<Elt>(e) );
        }
    /*private :

    bool M_searchWithCommunication;
    bool M_componentsAreSamePoint;
    bool M_onlyLocalizeOnBoundary;
     int M_nbNearNeighborInKdTree;*/
};

struct InterpolationNonConforming : public InterpolationTypeBase<false,interpolation_operand_type::ID>
{
    using super = InterpolationTypeBase<false,interpolation_operand_type::ID>;
    static const uint16_type value=super::isConforming();

    constexpr InterpolationNonConforming() = default;

    constexpr InterpolationNonConforming( bool useComm,
                                          bool compAreSamePt=true,
                                          bool onlyLocalizeOnBoundary=false,
                                          int nbNearNeighborInKdTree=15 )
        :
        super(useComm,compAreSamePt, onlyLocalizeOnBoundary,nbNearNeighborInKdTree)
    {}
    constexpr InterpolationNonConforming( nonconforming_t nc,
                                          bool useComm=true,
                                          bool compAreSamePt=true,
                                          bool onlyLocalizeOnBoundary=false,
                                          int nbNearNeighborInKdTree=15 )
        :
        super(useComm,compAreSamePt, onlyLocalizeOnBoundary,nbNearNeighborInKdTree)
        {}
    InterpolationNonConforming( InterpolationNonConforming const& ) = default;
    InterpolationNonConforming( InterpolationNonConforming && ) = default;
    InterpolationNonConforming& operator=( InterpolationNonConforming const& ) = default;
};
using InterpolationNonConforme = InterpolationNonConforming;
struct InterpolationConforming : public InterpolationTypeBase<true,interpolation_operand_type::ID>
{
    using super = InterpolationTypeBase<true,interpolation_operand_type::ID>;
    static const uint16_type value=super::isConforming();

    constexpr InterpolationConforming() = default;

    constexpr InterpolationConforming(bool useComm,
                                      bool compAreSamePt=true,
                                      bool onlyLocalizeOnBoundary=false,
                                      int nbNearNeighborInKdTree=15 )
        :
        super(useComm,compAreSamePt,onlyLocalizeOnBoundary,nbNearNeighborInKdTree)
    {}
    constexpr InterpolationConforming(conforming_t c,
                                      bool useComm=true,
                                      bool compAreSamePt=true,
                                      bool onlyLocalizeOnBoundary=false,
                                      int nbNearNeighborInKdTree=15 )
        :
        super(useComm,compAreSamePt,onlyLocalizeOnBoundary,nbNearNeighborInKdTree)
        {}
    InterpolationConforming( InterpolationConforming const& ) = default;
    InterpolationConforming( InterpolationConforming && ) = default;
    InterpolationConforming& operator=( InterpolationConforming const& ) = default;

};
using InterpolationConforme = InterpolationConforming;

template<typename C>
using InterpolationID = typename mpl::if_<C,
                                          mpl::identity<InterpolationConforming>,
                                          mpl::identity<InterpolationNonConforming>>::type::type;

template<typename C>
struct InterpolationGradient : public InterpolationTypeBase<is_conforming<C>,interpolation_operand_type::GRADIENT>
{
    using super = InterpolationTypeBase<is_conforming<C>,interpolation_operand_type::GRADIENT>;

    constexpr InterpolationGradient() = default;
    constexpr InterpolationGradient( C c,
                                     bool useComm = true,
                                     bool compAreSamePt=true,
                                     bool onlyLocalizeOnBoundary=false,
                                     int nbNearNeighborInKdTree=15 )
        :
        super(useComm,compAreSamePt,onlyLocalizeOnBoundary,nbNearNeighborInKdTree)
        {}
    InterpolationGradient( InterpolationGradient const& ) = default;
    InterpolationGradient( InterpolationGradient && ) = default;
    InterpolationGradient& operator=( InterpolationGradient const& ) = default;
};

template<typename C>
struct InterpolationCurl : public InterpolationTypeBase<is_conforming<C>,interpolation_operand_type::CURL>
{
    using super = InterpolationTypeBase<is_conforming<C>,interpolation_operand_type::CURL>;
    constexpr InterpolationCurl() = default;
    constexpr InterpolationCurl( C c,
                                     bool useComm = true,
                                     bool compAreSamePt=true,
                                     bool onlyLocalizeOnBoundary=false,
                                     int nbNearNeighborInKdTree=15 )
        :
        super(useComm,compAreSamePt,onlyLocalizeOnBoundary,nbNearNeighborInKdTree)
        {}
    InterpolationCurl( InterpolationCurl const& ) = default;
    InterpolationCurl( InterpolationCurl && ) = default;
    InterpolationCurl& operator=( InterpolationCurl const& ) = default;
};

template<typename C>
struct InterpolationDiv : public InterpolationTypeBase<is_conforming<C>,interpolation_operand_type::DIV>
{
    using super = InterpolationTypeBase<is_conforming<C>,interpolation_operand_type::DIV>;
    constexpr InterpolationDiv() = default;
    constexpr InterpolationDiv( C c,
                                bool useComm = true,
                                bool compAreSamePt=true,
                                bool onlyLocalizeOnBoundary=false,
                                int nbNearNeighborInKdTree=15 )
        :
        super(useComm,compAreSamePt,onlyLocalizeOnBoundary,nbNearNeighborInKdTree)
        {}
    InterpolationDiv( InterpolationDiv const& ) = default;
    InterpolationDiv( InterpolationDiv && ) = default;
    InterpolationDiv& operator=( InterpolationDiv const& ) = default;
};

template<typename ExprType,typename C>
struct InterpolationExpr : public InterpolationTypeBaseBase<is_conforming<C>,interpolation_operand_type::EXPR>
{
    using super = InterpolationTypeBaseBase<is_conforming<C>,interpolation_operand_type::EXPR>;
    using expr_type = ExprType;
    constexpr InterpolationExpr() = default;
    constexpr InterpolationExpr(expr_type const& expr,
                                C c,
                                bool useComm = true,
                                bool compAreSamePt=true,
                                bool onlyLocalizeOnBoundary=false,
                                int nbNearNeighborInKdTree=15 )
        :
        super(useComm,compAreSamePt,onlyLocalizeOnBoundary,nbNearNeighborInKdTree),
        M_expr( expr )
        {}
    InterpolationExpr( InterpolationExpr const& ) = default;
    InterpolationExpr( InterpolationExpr && ) = default;
    //InterpolationExpr& operator=( InterpolationExpr const& ) = default;

    template<typename Elt>
    constexpr decltype(auto) operand( Elt&& e )
        {
            return M_expr;
        }

private :
    expr_type M_expr;
};

template<typename C, typename ...Args>
InterpolationID<C>
makeInterpolation( C&& c, Args&&... args )
{
    return InterpolationID<C>( std::forward<C>(c), std::forward<Args>(args)... );
}

template<typename C, typename ...Args>
InterpolationGradient<C>
makeGradientInterpolation( C&& c, Args&&... args )
{
    return InterpolationGradient<C>( std::forward<C>(c), std::forward<Args>(args)... );
}

template<typename C, typename ...Args>
InterpolationCurl<C>
makeCurlInterpolation( C&& c, Args&&... args )
{
    return InterpolationCurl<C>( std::forward<C>(c), std::forward<Args>(args)... );
}

template<typename C, typename ...Args>
InterpolationDiv<C>
makeDivInterpolation( C&& c, Args&&... args )
{
    return InterpolationDiv<C>( std::forward<C>(c), std::forward<Args>(args)... );
}

template<typename ExprType, typename C, typename ...Args>
InterpolationExpr<ExprType,C>
makeExprInterpolation( ExprType const& expr, C&& c, Args&&... args )
{
    return InterpolationExpr<ExprType,C>( expr, std::forward<C>(c), std::forward<Args>(args)... );
}


template <typename T=double>
struct OperatorInterpolationMatrixSetup
{
    using sparse_matrix_type = MatrixSparse<T>;
    using sparse_matrix_ptrtype = std::shared_ptr<sparse_matrix_type>;

    OperatorInterpolationMatrixSetup( sparse_matrix_ptrtype m = sparse_matrix_ptrtype{}, Feel::MatrixStructure c = Feel::DIFFERENT_NONZERO_PATTERN,
                                      size_type indexBlockRow = 0, size_type indexBlockCol = 0 )
        :
        M_matrix( m ),
        M_indexBlockSpaceRow( indexBlockRow ),
        M_indexBlockSpaceCol( indexBlockCol ),
        M_context( c )
        {}

    OperatorInterpolationMatrixSetup( OperatorInterpolationMatrixSetup const& ) = default;
    OperatorInterpolationMatrixSetup( OperatorInterpolationMatrixSetup && ) = default;

    OperatorInterpolationMatrixSetup& operator=( OperatorInterpolationMatrixSetup const& ) = default;
    OperatorInterpolationMatrixSetup& operator=( OperatorInterpolationMatrixSetup && ) = default;

    sparse_matrix_ptrtype matrix() const { return M_matrix; }
    size_type indexBlockSpaceRow() const { return M_indexBlockSpaceRow; }
    size_type indexBlockSpaceCol() const { return M_indexBlockSpaceCol; }

    bool hasContext( Feel::MatrixStructure c ) const { return (M_context == c); }

    void setMatrix( sparse_matrix_ptrtype m ) { M_matrix = m;; }
private :
    sparse_matrix_ptrtype M_matrix;
    size_type M_indexBlockSpaceRow, M_indexBlockSpaceCol;
    Feel::MatrixStructure M_context;
};
//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
namespace detail
{



template <typename DomainSpaceType, typename ImageSpaceType, typename ExprType,int DomainEvalOnSubEntityCoDim/*=0*/>
struct PrecomputeDomainBasisFunction
{
    typedef std::shared_ptr<DomainSpaceType> domain_space_ptrtype;
    typedef std::shared_ptr<ImageSpaceType> image_space_ptrtype;
    //typedef GeoElementType geoelement_type;
    using domain_mesh_type = typename DomainSpaceType::mesh_type;
    using image_mesh_type = typename ImageSpaceType::mesh_type;
    using domain_geoelement_type = typename domain_mesh_type::element_type;
    using image_geoelement_type = typename image_mesh_type::element_type;
    typedef typename domain_geoelement_type::permutation_type domain_permutation_type;
    typedef typename DomainSpaceType::basis_type fe_type;
    typedef typename ImageSpaceType::basis_type image_fe_type;
    typedef typename image_fe_type::template ChangeDim<fe_type::nDim>::type image_fe_changedim_type;

    typedef ExprType expression_type;

    // geomap context
    typedef typename domain_geoelement_type::gm_type domain_gm_type;
    typedef typename domain_geoelement_type::gm_ptrtype domain_gm_ptrtype;
    typedef typename image_geoelement_type::gm_type image_gm_type;
    typedef typename image_geoelement_type::gm_ptrtype image_gm_ptrtype;
    typedef typename domain_gm_type::precompute_type domain_geopc_type;
    typedef typename domain_gm_type::precompute_ptrtype domain_geopc_ptrtype;
    using image_geopc_type = typename image_gm_type::precompute_type;

    static const size_type context2 = (is_hdiv_conforming_v<image_fe_type> || is_hcurl_conforming_v<image_fe_type> )?
        expression_type::context|vm::JACOBIAN|vm::KB :
        expression_type::context;
    // static const size_type context = ( DomainSpaceType::nDim == ImageSpaceType::nDim )? context2 : context2|vm::POINT;
    static const size_type context = context2|vm::POINT;//|vm::KB|vm::GRAD; // TOD REVERT
    typedef typename domain_gm_type::template Context<domain_geoelement_type,DomainEvalOnSubEntityCoDim> domain_gmc_type;
    typedef std::shared_ptr<domain_gmc_type> domain_gmc_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, domain_gmc_ptrtype> > map_gmc_type;


    using image_gmc_type = typename image_gm_type::template Context<image_geoelement_type,0>;


    // fe context
    typedef typename fe_type::template Context< context, fe_type, domain_gm_type, domain_geoelement_type,0, domain_gmc_type::subEntityCoDim> fecontext_type;
    typedef std::shared_ptr<fecontext_type> fecontext_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, fecontext_ptrtype> > map_fec_type;

    // tensor expression
    typedef typename expression_type::template tensor<map_gmc_type,map_fec_type> t_expr_type;


    PrecomputeDomainBasisFunction( domain_space_ptrtype const& XhDomain, image_space_ptrtype const& XhImage, ExprType const& expr )
        :
        M_XhDomain( XhDomain ),
        M_XhImage( XhImage ),
        M_expr( expr ),
        M_isSameMesh( M_XhDomain->mesh()->isSameMesh( M_XhImage->mesh() ) )
        {}

    //! return current interpolant
    Eigen::MatrixXd const& interpolant() const { return M_IhLoc; }

    domain_gmc_ptrtype & gmc() { return M_gmc; }

    template <typename InterpOnType,
              std::enable_if_t<(InterpOnType::onEntity() == ElementsType::MESH_ELEMENTS),bool> = true >
    void update( InterpOnType const& interpOnElts )
        {
            auto const& interpOnElt = interpOnElts.elements().front();
            auto const& elt = interpOnElt.element();
            this->init( interpOnElt );

            // update geomap context, finite element context and expr
            M_gmc->template update<context>( elt );
            M_fec->update( M_gmc );
            M_tensorExpr->update( mapgmc( M_gmc), mapfec( M_fec ) );

            // update local interpolation matrix
            if constexpr ( DomainSpaceType::nDim == ImageSpaceType::nDim )
            {
                M_XhImage->fe()->interpolateBasisFunction( *M_tensorExpr, M_IhLoc );
            }
            else
            {
                auto const& imageElt = interpOnElts.imageElement();
                if ( !M_gmcImageElt )
                    M_gmcImageElt = M_XhImage->gm()->template context<vm::POINT>( imageElt, M_geopcImageElt );
                M_gmcImageElt->template update<vm::POINT>( imageElt );

                std::vector<uint16_type> mapImageChangeDimFaceDofPointToImageDofPoint( M_gmc->nPoints(), invalid_v<uint16_type> );
                upMatchingDofPoints( M_gmc, M_gmcImageElt, mapImageChangeDimFaceDofPointToImageDofPoint );
                M_XhImage->fe()->interpolateBasisFunction( *M_tensorExpr, M_IhLoc, mapImageChangeDimFaceDofPointToImageDofPoint );
            }
        }


    template <typename InterpOnType,
              std::enable_if_t<(InterpOnType::onEntity() == ElementsType::MESH_FACES) && (DomainSpaceType::nDim == ImageSpaceType::nDim) ,bool> = true >
    void update( InterpOnType const& interpOnFaceFull )
        {
            auto const& interpOnFace = interpOnFaceFull.elements().front();
            auto const& imageElt = interpOnFaceFull.imageElement();

            this->init( interpOnFace );

            uint16_type __face_id = interpOnFace.faceIdInElement();
            auto const& domainElt = interpOnFace.element();
            //M_gmc->template update
            M_gmc->template update<context>( domainElt, __face_id );
#if 1
            auto image_fe = M_XhImage->fe();
            domain_gm_ptrtype domain_gm = M_XhDomain->gm();
            //auto fepc = std::make_shared<geopc_type>( domain_gm, image_fe->points( __face_id ) );
            auto fepc = M_XhDomain->fe()->preCompute( image_fe->points( __face_id ) );

            //std::cout << "elt.id()="<< elt.id() << "__face_id=" << __face_id << " pts: " << image_fe->points( __face_id ) << " xReal()= " << M_gmc->xReal() << std::endl;
            M_fec.reset( new fecontext_type( M_XhDomain->fe(), M_gmc, fepc /*geopc*/ ) );
            M_tensorExpr = std::make_shared<t_expr_type>( M_expr, mapgmc( M_gmc), mapfec( M_fec ) );
#else
            M_fec->update( M_gmc );
#endif
            M_tensorExpr->update( mapgmc( M_gmc), mapfec( M_fec ) );

            // WARNING : here this case is turn to off in waiting the fix in Lagrange::interpolateBasisFunction for example
            if ( false && M_isSameMesh && ( imageElt.id() == domainElt.id() ) )
            {
                M_XhImage->fe()->interpolateBasisFunction( *M_tensorExpr, M_IhLoc );
            }
            else
            {
                if ( !M_gmcImageElt )
                    M_gmcImageElt = M_XhImage->gm()->template context<vm::POINT>( imageElt, M_geopcImageElt );
                M_gmcImageElt->template update<vm::POINT>( imageElt );

                std::vector<uint16_type> mapImageChangeDimFaceDofPointToImageDofPoint( M_gmc->nPoints(), invalid_v<uint16_type> );
                upMatchingDofPoints( M_gmc, M_gmcImageElt, mapImageChangeDimFaceDofPointToImageDofPoint );
                M_XhImage->fe()->interpolateBasisFunction( *M_tensorExpr, M_IhLoc, mapImageChangeDimFaceDofPointToImageDofPoint );
            }

        }


    template <typename Gmc1Type,typename Gmc2Type>
    static void upMatchingDofPoints( std::shared_ptr<Gmc1Type> gmc1, std::shared_ptr<Gmc2Type> gmc2, std::vector<uint16_type>& mapGmc1ToGmc2 )
        {
            uint16_type nPoint1 = gmc1->nPoints();
            uint16_type nPoint2 = gmc2->nPoints();
            double dofPtCompareTol = std::max( 1e-15, std::min( gmc1->element().hMin(), gmc2->element().hMin() ) * 1e-5 );
            for ( uint16_type jloc_d = 0; jloc_d < nPoint1; ++jloc_d )
            {
                auto const& domainGlobDofPt = gmc1->xReal( jloc_d );
                bool find = false;
                size_type thelocDofToFind = invalid_v<size_type>;
                for ( uint16_type jloc_i = 0; jloc_i < nPoint2; ++jloc_i )
                {
                    auto const& imageGlobDofPt = gmc2->xReal( jloc_i );
                    bool find2 = true;
                    for ( uint16_type d = 0; d < DomainSpaceType::nRealDim; ++d )
                    {
                        find2 = find2 && ( std::abs( imageGlobDofPt[d] - domainGlobDofPt[d] ) < dofPtCompareTol );
                    }
                    if ( find2 )
                    {
                        mapGmc1ToGmc2[ jloc_d ] = jloc_i;
                        find = true;
                        break;
                    }
                }
                CHECK( find ) << "[OperatorInterpolation::update] Compatible dof not found";
            }
        }

    template <typename InterpOnType,//typename ImageEltType,
              std::enable_if_t<(InterpOnType::onEntity() == ElementsType::MESH_FACES) && (DomainSpaceType::nDim == (ImageSpaceType::nDim+1)) ,bool> = true >
    void update( InterpOnType const& interpOnFaceFull )//, ImageEltType const& imageElt, uint16_type imageLocDof, size_type imageGlobDof, uint16_type comp )
        {
            auto const& interpOnFace = interpOnFaceFull.elements().front();
            auto const& imageElt = interpOnFaceFull.imageElement(); // TODO CHECK face id in elt maybe?

            this->init( interpOnFace );

            if ( !M_gmcImageElt )
                M_gmcImageElt = M_XhImage->gm()->template context<vm::POINT>( imageElt, M_geopcImageElt );
            M_gmcImageElt->template update<vm::POINT>( imageElt );

            uint16_type __face_id = interpOnFace.faceIdInElement();
            auto const& elt = interpOnFace.element();
            M_gmc->template update<context>( elt, __face_id );

            std::vector<uint16_type> mapImageChangeDimFaceDofPointToImageDofPoint( M_gmc->nPoints()/*image_fe_changedim_type::nLocalDof*/, invalid_v<uint16_type> );
            upMatchingDofPoints( M_gmc, M_gmcImageElt, mapImageChangeDimFaceDofPointToImageDofPoint );



#if 1
            //auto image_fe = M_XhImage->fe();
            domain_gm_ptrtype domain_gm = M_XhDomain->gm();
            //auto fepc = std::make_shared<geopc_type>( domain_gm, image_fe->points( __face_id ) );
            auto fepc = M_XhDomain->fe()->preCompute( M_imageFeChangeDim->points( __face_id ) );

            //std::cout << "elt.id()="<< elt.id() << "__face_id=" << __face_id << " pts: " << image_fe->points( __face_id ) << " xReal()= " << M_gmc->xReal() << std::endl;
            M_fec.reset( new fecontext_type( M_XhDomain->fe(), M_gmc, fepc /*geopc*/ ) );

            M_tensorExpr = std::make_shared<t_expr_type>( M_expr, mapgmc( M_gmc), mapfec( M_fec ) );

#else
            M_fec->update( M_gmc );
#endif
            M_tensorExpr->update( mapgmc( M_gmc), mapfec( M_fec ) );


            M_XhImage->fe()->interpolateBasisFunction( *M_tensorExpr, M_IhLoc, mapImageChangeDimFaceDofPointToImageDofPoint );


            // CHECK(false) << "TODO";
        }

private :
    //template <typename MeshEntityType,std::enable_if_t< std::is_same_v<MeshEntityType,typename domain_mesh_type::element_type>, bool> = true >
    template <typename InterpOnType,
              std::enable_if_t<InterpOnType::onEntity() == ElementsType::MESH_ELEMENTS,bool> = true >
    void init( InterpOnType const& interpFromElt )//MeshEntityType const& elt )
        {
            if ( M_tensorExpr )
                return;

            if ( !M_geopcImageElt )
            {
                auto image_gm = M_XhImage->gm();
                auto image_fe = M_XhImage->fe();
                auto refPts = image_fe->dual().points();
                M_geopcImageElt = image_gm->preCompute( image_gm, refPts );
            }

            if constexpr ( DomainSpaceType::nDim != ImageSpaceType::nDim )
            {
                if ( !M_imageFeChangeDim )
                    M_imageFeChangeDim = std::make_shared<image_fe_changedim_type>();
                this->initCommon( interpFromElt, M_imageFeChangeDim );
            }
            else
            {
                this->initCommon( interpFromElt, M_XhImage->fe() );
            }
        }
    template <typename InterpOnType,typename TheImageFeType,
              std::enable_if_t<InterpOnType::onEntity() == ElementsType::MESH_ELEMENTS,bool> = true >
    void initCommon( InterpOnType const& interpFromElt, std::shared_ptr<TheImageFeType> image_fe )
        {
            auto const& elt = interpFromElt.element();
            domain_gm_ptrtype gm = M_XhDomain->gm();

            auto refPts = image_fe->dual().points();
            auto geopc = gm->preCompute( gm, refPts );
            auto fepc = M_XhDomain->fe()->preCompute( M_XhDomain->fe(), refPts /*gmc->xRefs()*/ );

            domain_gmc_ptrtype gmc = gm->template context<context>( elt, geopc );
            fecontext_ptrtype fec( new fecontext_type( M_XhDomain->fe(), gmc, fepc /*geopc*/ ) );

            M_gmc = gm->template context<context>( elt, geopc );
            M_fec.reset( new fecontext_type( M_XhDomain->fe(), gmc, fepc /*geopc*/ ) );

            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0>>( M_gmc ) );
            map_fec_type mapfec( fusion::make_pair<vf::detail::gmc<0>>( M_fec ) );
            M_tensorExpr = std::make_shared<t_expr_type>( M_expr, mapgmc, mapfec );

            using shape = typename t_expr_type::shape;
            M_IhLoc = Eigen::MatrixXd::Zero( fe_type::is_product ? fe_type::nComponents * fe_type::nLocalDof : fe_type::nLocalDof,
                                             image_fe_type::is_product ? image_fe_type::nComponents * image_fe_type::nLocalDof : image_fe_type::nLocalDof );
#if 0
            image_fe->interpolateBasisFunction( *M_tensorExpr, M_IhLoc );
#endif
        }

    template <typename InterpOnType,
              std::enable_if_t<InterpOnType::onEntity() == ElementsType::MESH_FACES,bool> = true >
    void init( InterpOnType const& interpOnFace )
        {
            auto image_fe = M_XhImage->fe();

            if ( !M_geopcImageElt )
            {
                auto image_gm = M_XhImage->gm();
                auto refPts = image_fe->dual().points();
                M_geopcImageElt = image_gm->preCompute( image_gm, refPts );
            }

            if constexpr ( DomainSpaceType::nDim != ImageSpaceType::nDim )
            {
                if ( !M_imageFeChangeDim )
                    M_imageFeChangeDim = std::make_shared<image_fe_changedim_type>();
                this->initCommon( interpOnFace, M_imageFeChangeDim );
            }
            else
            {
                this->initCommon( interpOnFace, M_XhImage->fe() );
            }
        }

    template <typename InterpOnType,typename TheImageFeType,
              std::enable_if_t<InterpOnType::onEntity() == ElementsType::MESH_FACES,bool> = true >
    void initCommon( InterpOnType const& interpOnFace, std::shared_ptr<TheImageFeType> image_fe )
        {
            domain_gm_ptrtype domain_gm = M_XhDomain->gm();

            //geopc_ptrtype geopc( new geopc_type( gm, imageSpace->fe()->dual().points() ) );
            //auto refPts = M_XhImage->fe()->dual().points();
            //auto geopc = gm->preCompute( gm, refPts );



            std::vector<std::map<domain_permutation_type, domain_geopc_ptrtype> > __geopc( domain_geoelement_type::numTopologicalFaces );
            //std::vector<std::map<permutation_type, geopc1_ptrtype> > __geopc1( geoelement_type::numTopologicalFaces );

            for ( uint16_type __f = 0; __f < domain_geoelement_type::numTopologicalFaces; ++__f )
            {
                for ( domain_permutation_type __p( domain_permutation_type::IDENTITY );
                      __p < domain_permutation_type( domain_permutation_type::N_PERMUTATIONS ); ++__p )
                {
                    __geopc[__f][__p] = domain_geopc_ptrtype(  new domain_geopc_type( domain_gm, image_fe->points( __f ) ) );
                    //__geopc1[__f][__p] = geopc1_ptrtype(  new geopc1_type( __gm1, __fe->points( __f ) ) );
                    //DVLOG(2) << "[geopc] FACE_ID = " << __f << " ref pts=" << __fe->dual().points( __f ) << "\n";
                    //FEELPP_ASSERT( __geopc[__f][__p]->nPoints() ).error( "invalid number of points for geopc" );
                    //FEELPP_ASSERT( __geopc1[__f][__p]->nPoints() ).error( "invalid number of points for geopc1" );
                }
            }

            // auto const& element() const { return unwrap_ref( M_element ); }
            // uint16_type faceIdInElement() const { return M_faceIdInElement; }
            // auto const& face() const { return this->element().face( faceIdInElement() ); }

            uint16_type __face_id = interpOnFace.faceIdInElement();
            auto const& elt = interpOnFace.element();

            M_gmc = domain_gm->template context<context>( elt, __geopc, __face_id, M_expr.dynamicContext() );

            //  gmc_ptrtype gmc = domain_gm->template context<context>( face.element(0), __geopc, __face_id, M_expr.dynamicContext() );
            //auto __c1 = __gm1->template context<context>( firstFace.element( 0 ), __geopc1, __face_id, M_expr.dynamicContext() );


            //if constexpr ( gmc_type::subEntityCoDim > 0 )
            //M_pc->update(fusion::at_key<key_type>( geom )->pc()->nodes() );

            //auto fepc = std::make_shared<geopc_type>( domain_gm, image_fe->points( __face_id ) );

            auto fepc = M_XhDomain->fe()->preCompute( image_fe->points( __face_id ) );
            //basis_0_type::PreCompute
            //pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );

            M_fec.reset( new fecontext_type( M_XhDomain->fe(), M_gmc, fepc /*geopc*/ ) );


            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0>>( M_gmc ) );
            map_fec_type mapfec( fusion::make_pair<vf::detail::gmc<0>>( M_fec ) );
            //t_expr_type texpr( M_expr, mapgmc, mapfec );
            M_tensorExpr = std::make_shared<t_expr_type>( M_expr, mapgmc, mapfec );

            using shape = typename t_expr_type::shape;
            M_IhLoc = Eigen::MatrixXd::Zero( fe_type::is_product ? fe_type::nComponents * fe_type::nLocalDof : fe_type::nLocalDof,
                                             image_fe_type::is_product ? image_fe_type::nComponents * image_fe_type::nLocalDof : image_fe_type::nLocalDof );
#if 0
            M_XhImage->fe()->interpolateBasisFunction( *M_tensorExpr, M_IhLoc );
#endif


        }


private:
    domain_space_ptrtype M_XhDomain;
    image_space_ptrtype M_XhImage;
    expression_type const& M_expr;
    bool M_isSameMesh;
    std::shared_ptr<t_expr_type> M_tensorExpr;
    std::shared_ptr<image_fe_changedim_type> M_imageFeChangeDim;
    domain_gmc_ptrtype M_gmc;
    std::shared_ptr<image_geopc_type> M_geopcImageElt;
    std::shared_ptr<image_gmc_type> M_gmcImageElt;
    fecontext_ptrtype M_fec;
    Eigen::MatrixXd M_IhLoc;
};

template <int DomainEvalOnSubEntityCoDim,typename DomainSpaceType, typename ImageSpaceType, typename ExprType>
std::shared_ptr<PrecomputeDomainBasisFunction<DomainSpaceType,ImageSpaceType,ExprType,DomainEvalOnSubEntityCoDim>>
precomputeDomainBasisFunction( std::shared_ptr<DomainSpaceType> const& domainSpace , std::shared_ptr<ImageSpaceType> const& imageSpace,
                               ExprType const& expr )
{
    typedef PrecomputeDomainBasisFunction<DomainSpaceType,ImageSpaceType,ExprType,DomainEvalOnSubEntityCoDim> res_type;
    return std::make_shared<res_type>( domainSpace, imageSpace, expr );
}







//--------------------------------------------------------------------------------------------------//
template <typename DomainSpaceType, typename ImageSpaceType>
struct DomainEltIdFromImageEltId
{
    using domain_space_type = DomainSpaceType;
    using domain_mesh_type = typename domain_space_type::mesh_type;
    using image_space_type = ImageSpaceType;
    using image_mesh_type = typename image_space_type::mesh_type;

    static constexpr uint16_type nDimDomain = domain_mesh_type::nDim;
    static constexpr uint16_type nDimImage = image_mesh_type::nDim;
    static constexpr int nDimDiffBetweenDomainImage = nDimDomain-nDimImage;


    using domain_vector_mesh_element_type = typename domain_mesh_type::super_elements::elements_reference_wrapper_type;
    using image_vector_mesh_element_type = typename image_mesh_type::super_elements::elements_reference_wrapper_type;

    using size_type = typename domain_mesh_type::size_type;
    using index_type = typename domain_mesh_type::index_type;

#if 0
    struct ReturnFromElement
    {
        static constexpr ElementsType onEntity() { return ElementsType::MESH_ELEMENTS; }

        ReturnFromElement( typename domain_vector_mesh_element_type::value_type const& elt )
            :
            M_element( elt )
            {}
        auto const& element() const { return unwrap_ref( M_element ); }

    private :
        typename domain_vector_mesh_element_type::value_type M_element;
    };

    struct ReturnFromFace
    {
        static constexpr ElementsType onEntity() { return ElementsType::MESH_FACES; }

        ReturnFromFace( typename domain_vector_mesh_element_type::value_type const& elt, uint16_type faceIdInElement )
            :
            M_element( elt ),
            M_faceIdInElement( faceIdInElement )
            {}
        ReturnFromFace( ReturnFromFace && ) = default;
        ReturnFromFace( ReturnFromFace const& ) = default;

        auto const& element() const { return unwrap_ref( M_element ); }
        uint16_type faceIdInElement() const { return M_faceIdInElement; }
        auto const& face() const { return this->element().face( faceIdInElement() ); }

    private :
        typename domain_vector_mesh_element_type::value_type M_element;
        uint16_type M_faceIdInElement;

    };


    struct ReturnFromElements
    {
        static constexpr ElementsType onEntity() { return ReturnFromElement::onEntity(); }

        ReturnFromElements( typename image_vector_mesh_element_type::value_type imageElement,
                            uint16_type faceIdInImageElement = invalid_v<uint16_type> )
            :
            M_imageElement( imageElement ),
            M_faceIdInImageElement( faceIdInImageElement )
            {}

        auto const& imageElement() const { return unwrap_ref( M_imageElement ); }
        uint16_type faceIdInImageElement() const { return M_faceIdInImageElement; }

        std::vector<ReturnFromElement> const& elements() const { return M_elements; }

        bool empty() const { return M_elements.empty(); }

        bool addElement( std::shared_ptr<MeshSupport<domain_mesh_type>> meshSupport, typename domain_vector_mesh_element_type::value_type const& elt )
            {
                if ( meshSupport->isPartialSupport() )
                {
                    index_type eltId = unwrap_ref( elt ).id();
                    if ( meshSupport->hasElement( eltId ) && !meshSupport->hasGhostElement( eltId ) )
                    {
                        M_elements.push_back( ReturnFromElement{elt} );
                        return true;
                    }
                }
                else
                {
                    M_elements.push_back( ReturnFromElement{elt} );
                    return true;
                }
                return false;
            }

        void keepOnly( index_type domainId )
            {
                auto itFindId = std::find_if(M_elements.begin(), M_elements.end(), [&domainId]( auto const& rfe ){ return rfe.element().id() == domainId; });
                if ( itFindId != M_elements.end() )
                {
                    auto rfeFound = std::move(*itFindId);
                    M_elements.clear();
                    M_elements.push_back( std::move( rfeFound ) );
                }
                else
                    CHECK( false ) << "elt id not found";
            }

    private :
        std::vector<ReturnFromElement> M_elements;
        typename image_vector_mesh_element_type::value_type M_imageElement;
        uint16_type M_faceIdInImageElement = invalid_v<uint16_type> ;
    };

    struct ReturnFromFaces
    {
        static constexpr ElementsType onEntity() { return ReturnFromFace::onEntity(); }

        ReturnFromFaces( typename image_vector_mesh_element_type::value_type imageElement,
                         uint16_type faceIdInImageElement = invalid_v<uint16_type> )
            :
            M_imageElement( imageElement ),
            M_faceIdInImageElement( faceIdInImageElement )
            {}


        auto const& imageElement() const { return unwrap_ref( M_imageElement ); }
        uint16_type faceIdInImageElement() const { return M_faceIdInImageElement; }

        std::vector<ReturnFromFace> const& elements() const { return M_elements; }

        bool empty() const { return M_elements.empty(); }

        bool addElement( std::shared_ptr<MeshSupport<domain_mesh_type>> meshSupport, typename domain_vector_mesh_element_type::value_type const& elt, uint16_type faceIdInElt )
            {
                if ( meshSupport->isPartialSupport() )
                {
                    index_type eltId = unwrap_ref( elt ).id();
                    if ( meshSupport->hasElement( eltId ) && !meshSupport->hasGhostElement( eltId ) )
                    {
                        M_elements.push_back( ReturnFromFace{elt, faceIdInElt} );
                        return true;
                    }
                }
                else
                {
                    M_elements.push_back( ReturnFromFace{elt, faceIdInElt} );
                    return true;
                }
                return false;
            }
        void keepOnly( index_type domainId )
            {
                auto itFindId = std::find_if(M_elements.begin(), M_elements.end(), [&domainId]( auto const& rfe ){ return rfe.element().id() == domainId; });
                if ( itFindId != M_elements.end() )
                {
                    auto rfeFound = std::move(*itFindId);
                    M_elements.clear();
                    M_elements.push_back( std::move( rfeFound ) );
                }
                else
                    CHECK( false ) << "elt id not found";
            }

    private:
        std::vector<ReturnFromFace> M_elements;
        typename image_vector_mesh_element_type::value_type M_imageElement;
        uint16_type M_faceIdInImageElement = invalid_v<uint16_type> ;
    };
#else

    template <int TheCoDim>
    struct ReturnDomainEntities
    {
        static constexpr ElementsType onEntity() { return TheCoDim==0?ElementsType::MESH_ELEMENTS:ElementsType::MESH_FACES; }

        ReturnDomainEntities( typename domain_vector_mesh_element_type::value_type const& elt, uint16_type faceIdInElement = invalid_v<uint16_type> )
            :
            M_element( elt ),
            M_faceIdInElement( faceIdInElement )
            {}
        ReturnDomainEntities( ReturnDomainEntities && ) = default;
        ReturnDomainEntities( ReturnDomainEntities const& ) = default;

        auto const& element() const { return unwrap_ref( M_element ); }
        uint16_type faceIdInElement() const { return M_faceIdInElement; }
        auto const& face() const { return this->element().face( faceIdInElement() ); }

    private :
        typename domain_vector_mesh_element_type::value_type M_element;
        uint16_type M_faceIdInElement;

    };

    using ReturnFromElement = ReturnDomainEntities<0>;
    using ReturnFromFace = ReturnDomainEntities<1>;



    template <int TheCoDim>
    struct ReturnType
    {
        using domain_entities_type = ReturnDomainEntities<TheCoDim>;
        static constexpr ElementsType onEntity() { return domain_entities_type::onEntity(); }

        ReturnType( typename image_vector_mesh_element_type::value_type imageElement,
                    uint16_type faceIdInImageElement = invalid_v<uint16_type> )
            :
            M_imageElement( imageElement ),
            M_faceIdInImageElement( faceIdInImageElement )
            {}

        auto const& imageElement() const { return unwrap_ref( M_imageElement ); }
        uint16_type faceIdInImageElement() const { return M_faceIdInImageElement; }

        std::vector<domain_entities_type> const& elements() const { return M_elements; }

        bool empty() const { return M_elements.empty(); }

        bool addElement( std::shared_ptr<MeshSupport<domain_mesh_type>> meshSupport, typename domain_vector_mesh_element_type::value_type const& elt, uint16_type faceIdInElement = invalid_v<uint16_type> )
            {
                if ( meshSupport->isPartialSupport() )
                {
                    index_type eltId = unwrap_ref( elt ).id();
                    if ( meshSupport->hasElement( eltId ) && !meshSupport->hasGhostElement( eltId ) )
                    {
                        M_elements.push_back( domain_entities_type{elt,faceIdInElement} );
                        return true;
                    }
                }
                else
                {
                    M_elements.push_back( domain_entities_type{elt,faceIdInElement} );
                    return true;
                }
                return false;
            }

        void keepOnly( index_type domainId )
            {
                auto itFindId = std::find_if(M_elements.begin(), M_elements.end(), [&domainId]( auto const& rfe ){ return rfe.element().id() == domainId; });
                if ( itFindId != M_elements.end() )
                {
                    auto rfeFound = std::move(*itFindId);
                    M_elements.clear();
                    M_elements.push_back( std::move( rfeFound ) );
                }
                else
                    CHECK( false ) << "elt id not found";
            }

    private :
        std::vector<domain_entities_type> M_elements;
        typename image_vector_mesh_element_type::value_type M_imageElement;
        uint16_type M_faceIdInImageElement = invalid_v<uint16_type> ;
    };

    using ReturnFromElements = ReturnType<0>;
    using ReturnFromFaces = ReturnType<1>;


#endif

    DomainEltIdFromImageEltId( std::shared_ptr<domain_space_type> XhDomain, std::shared_ptr<image_space_type> XhImage )
        :
        M_XhDomain( XhDomain ),
        M_XhImage( XhImage ),
        M_meshDomain( XhDomain->mesh() ),
        M_meshImage( XhImage->mesh() ),
        M_isSameMesh( M_meshImage->isSameMesh( M_meshDomain ) ),
        M_image_related_to_domain( M_meshImage->isSubMeshFrom( M_meshDomain ) ),
        M_domain_related_to_image( M_meshDomain->isSubMeshFrom( M_meshImage ) ),
        M_domain_sibling_of_image( M_meshDomain->isSiblingOf( M_meshImage ) ),
        M_meshSupportDomain( M_XhDomain->template meshSupport<0>() ),
        M_meshSupportImage( M_XhImage->template meshSupport<0>() ),
        M_hasMeshSupportPartialDomain( M_meshSupportDomain->isPartialSupport() ),
        M_hasMeshSupportPartialImage( M_meshSupportImage->isPartialSupport() )
        {}

    auto
    apply( typename image_mesh_type::element_type const& imageElt ) const
        {
            image_vector_mesh_element_type imageElts = {  boost::cref( imageElt ) };
            return this->apply<nDimDiffBetweenDomainImage>( imageElt, imageElts );
        }

    auto
    apply( typename image_mesh_type::face_type const& imageFace ) const
        {
            image_vector_mesh_element_type imageElts;
#if 0
            if ( imageFace.isConnectedTo0() && !imageFace.element0().isGhostCell() )
                if ( !M_hasMeshSupportPartialImage || M_meshSupportImage->hasElement( imageFace.element0().id() ) )
                    imageElts.push_back( boost::cref( imageFace.element0() ) );
            if ( imageFace.isConnectedTo1() && !imageFace.element1().isGhostCell() )
                if ( !M_hasMeshSupportPartialImage || M_meshSupportImage->hasElement( imageFace.element1().id() ) )
                    imageElts.push_back( boost::cref( imageFace.element1() ) );
#else
            if ( imageFace.isConnectedTo0() )
            {
                if ( M_hasMeshSupportPartialImage )
                {
                    //if ( !M_meshSupportImage->hasGhostElement( imageFace.element0().id() ) )
                    if ( M_meshSupportImage->hasElement( imageFace.element0().id() ) && !M_meshSupportImage->hasGhostElement( imageFace.element0().id() ) )
                        imageElts.push_back( boost::cref( imageFace.element0() ) );
                }
                else if ( !imageFace.element0().isGhostCell() )
                {
                    imageElts.push_back( boost::cref( imageFace.element0() ) );
                }
            }

            if ( imageFace.isConnectedTo1() )
            {
                if ( M_hasMeshSupportPartialImage )
                {
                    //if ( !M_meshSupportImage->hasGhostElement( imageFace.element1().id() ) )
                    if ( M_meshSupportImage->hasElement( imageFace.element1().id() ) && !M_meshSupportImage->hasGhostElement( imageFace.element1().id() ) )
                        imageElts.push_back( boost::cref( imageFace.element1() ) );
                }
                else if ( !imageFace.element1().isGhostCell() )
                {
                    imageElts.push_back( boost::cref( imageFace.element1() ) );
                }
            }

#endif

            return this->apply<nDimDiffBetweenDomainImage>( imageFace, imageElts );
        }
private :


    /*
     * - domain and image meshes haves same dimension
     * - interpolation on mesh element (belong to image mesh)
     */
    template <int ApplyType, typename MeshEntityType,
              std::enable_if_t< ApplyType == 0 && std::is_same_v<MeshEntityType,typename image_mesh_type::element_type>, bool> = true >
    auto apply( MeshEntityType const& imageElt, image_vector_mesh_element_type const& imageElts ) const
        {
            std::vector<ReturnFromElements> res;

            if ( imageElts.empty() )
                return res;

            auto imageEltWrap = boost::cref( imageElt );
            ReturnFromElements rfe{ imageEltWrap };

            size_type imageEltId = unwrap_ref(imageElt).id();
            std::optional<typename domain_vector_mesh_element_type::value_type> domainElt;
            if ( M_isSameMesh )
            {
                if constexpr ( std::is_same_v<typename domain_mesh_type::element_type,typename image_mesh_type::element_type> )
                    domainElt = imageEltWrap;
                else
                    CHECK( false ) << "should not go here";
            }
            else if ( M_image_related_to_domain )
            {
                index_type domainRelationEltId = M_meshImage->subMeshToMesh( imageEltId );
                if ( domainRelationEltId == invalid_v<index_type> )
                    return res;
                VLOG(2) << "[image_related_to_domain] image element id: "  << imageEltId << " domain element id : " << domainRelationEltId << "\n";
                domainElt = boost::cref( M_meshDomain->element( domainRelationEltId ) );

            }
            else if ( M_domain_related_to_image )
            {
                index_type domainRelationEltId = M_meshDomain->meshToSubMesh( imageEltId );
                if ( domainRelationEltId == invalid_v<index_type> )
                    return res;
                VLOG(2) << "[domain_related_to_image] image element id: "  << imageEltId << " domain element id : " << domainRelationEltId << "\n";
                domainElt = boost::cref( M_meshDomain->element( domainRelationEltId ) );
            }
            else if( M_domain_sibling_of_image )
            {
                index_type domainRelationEltId = M_meshDomain->meshToSubMesh( M_meshImage, imageEltId );
                if ( domainRelationEltId == invalid_v<index_type> )
                    return res;
                DVLOG(1) << "[domain_sibling_of_image] image element id: "  << imageEltId << " domain element id : " << domainRelationEltId << "\n";
                domainElt = boost::cref( M_meshDomain->element( domainRelationEltId ) );
            }

            if ( !domainElt )
                return res;

            rfe.addElement( M_meshSupportDomain, *domainElt );

            if ( !rfe.elements().empty() )
                res.push_back( std::move( rfe ) );

            return res;
        }

    /*
     * - domain and image meshes haves same dimension
     * - interpolation on mesh face (belong to image mesh)
     */
    template <int ApplyType, typename MeshEntityType,
              std::enable_if_t< ApplyType == 0 && std::is_same_v<MeshEntityType,typename image_mesh_type::face_type>, bool> = true >
    auto apply( MeshEntityType const& imageFace, image_vector_mesh_element_type const& imageElts ) const
        {
            std::vector<ReturnFromFaces> res;

            std::map<index_type,index_type> mapImageToDomainEltId;

            for ( auto const& imageElt : imageElts )
            {
                ReturnFromFaces rff{ imageElt };//imageElts.front() };

                size_type imageEltId = unwrap_ref(imageElt).id();

                uint16_type faceIdInElt = invalid_v<uint16_type>;
                if ( imageFace.isConnectedTo0() && imageFace.idElement0() == imageEltId )
                    faceIdInElt = imageFace.pos_first();
                else if ( imageFace.isConnectedTo1() && imageFace.idElement1() == imageEltId )
                    faceIdInElt = imageFace.pos_second();

                for ( uint16_type connectId = 0; connectId < 2 ;++connectId )
                {
                    if ( connectId == 0 && !imageFace.isConnectedTo0() )
                        continue;
                    else if ( connectId == 1 && !imageFace.isConnectedTo1() )
                        continue;

                    auto const& imageEltConnected = imageFace.element( connectId );
                    index_type imageEltConnectedId = imageEltConnected.id();

                    uint16_type faceIdInEltConnected = connectId == 0? imageFace.idInElement0() : imageFace.idInElement1();

                    std::optional<typename domain_vector_mesh_element_type::value_type> domainElt;
                    uint16_type domainFaceIdInElt = faceIdInEltConnected;//faceIdInElt;

                    if ( M_isSameMesh )
                    {
                        domainElt = boost::cref( imageEltConnected );
                        domainFaceIdInElt = faceIdInEltConnected;
                    }
                    else
                    {
                        // TODO CHECK domainFaceIdInElt is equal to faceIdInElt else find the good one
                        if ( M_image_related_to_domain )
                        {
                            index_type domainRelationEltId = M_meshImage->subMeshToMesh( imageEltConnectedId );
                            if ( domainRelationEltId == invalid_v<index_type> )
                                continue;
                            VLOG(2) << "[image_related_to_domain] image element id: "  << imageEltConnectedId << " domain element id : " << domainRelationEltId << "\n";
                            domainElt = boost::cref( M_meshDomain->element( domainRelationEltId ) );
                            domainFaceIdInElt = faceIdInEltConnected;
                        }
                        else if ( M_domain_related_to_image )
                        {
                            index_type domainRelationEltId = M_meshDomain->meshToSubMesh( imageEltConnectedId );
                            if ( domainRelationEltId == invalid_v<index_type> )
                                continue;
                            VLOG(2) << "[domain_related_to_image] image element id: "  << imageEltConnectedId << " domain element id : " << domainRelationEltId << "\n";
                            domainElt = boost::cref( M_meshDomain->element( domainRelationEltId ) );
                            domainFaceIdInElt = faceIdInEltConnected;
                        }
                        else if( M_domain_sibling_of_image )
                        {
                            index_type domainRelationEltId = M_meshDomain->meshToSubMesh( M_meshImage, imageEltConnectedId );
                            if ( domainRelationEltId == invalid_v<index_type> )
                                continue;
                            DVLOG(1) << "[domain_sibling_of_image] image element id: "  << imageEltConnectedId << " domain element id : " << domainRelationEltId << "\n";
                            domainElt = boost::cref( M_meshDomain->element( domainRelationEltId ) );
                            domainFaceIdInElt = faceIdInEltConnected;
                        }
                    }

                    if ( !domainElt )
                        continue;

                    bool isAdded = rff.addElement( M_meshSupportDomain, *domainElt, domainFaceIdInElt );
                    if ( isAdded )
                        mapImageToDomainEltId[imageEltConnectedId] = unwrap_ref( *domainElt ).id();
                }


                // if continue space and size > 1 : take only one (and if possible the same as image)
                if ( domain_space_type::is_continuous )
                {
                    if ( rff.elements().size() > 1 )
                    {
                        auto itFind = mapImageToDomainEltId.find( imageEltId );
                        index_type idKeep = itFind != mapImageToDomainEltId.end() ?
                            itFind->second : rff.elements().front().element().id();
                        rff.keepOnly( idKeep );
                    }
                    CHECK( rff.elements().size() <= 1 ) << "to many domain elements";
                }

                if ( !rff.elements().empty() )
                    res.push_back( std::move( rff ) );
            }

            return res;
        }


    /*
     * - domain dim = image dim + 1
     * - interpolation on mesh element (belong to image mesh)
     */
    template <int ApplyType, typename MeshEntityType,
              std::enable_if_t< ApplyType == 1 && std::is_same_v<MeshEntityType,typename image_mesh_type::element_type>, bool> = true >
    auto apply( MeshEntityType const& imageElt, image_vector_mesh_element_type const& imageElts ) const
        {
            std::vector<ReturnFromFaces> res;

            ReturnFromFaces rff{ boost::cref( imageElt ) };
            //ReturnFromFaces res;
            index_type imageEltId = imageElt.id();
            std::optional<boost::reference_wrapper<typename domain_mesh_type::face_type const>> domainFace;

            if ( M_image_related_to_domain )
            {
                index_type domainRelationFaceId = M_meshImage->subMeshToMesh( imageEltId );
                domainFace = boost::cref( M_meshDomain->face( domainRelationFaceId ) );
                VLOG(2) << "[image_related_to_domain] image element id: "  << imageEltId << " domain face id : " << domainRelationFaceId << "\n";
            }
            else if ( M_domain_related_to_image )
            {}

            if ( !domainFace )
                return res;

            auto const& thedomainFace = unwrap_ref( *domainFace );
            if ( thedomainFace.isConnectedTo0() )
                rff.addElement( M_meshSupportDomain, boost::cref(thedomainFace.element0()), thedomainFace.idInElement0() );
            if ( thedomainFace.isConnectedTo1() )
                rff.addElement( M_meshSupportDomain, boost::cref(thedomainFace.element1()), thedomainFace.idInElement1() );


            // if continue space and size > 1 : take only one (and if possible the same as image)
            if ( domain_space_type::is_continuous )
            {
                if ( rff.elements().size() > 1 )
                {
                    index_type idKeep = rff.elements().front().element().id();
                    rff.keepOnly( idKeep );
                }
                CHECK( rff.elements().size() <= 1 ) << "to many domain elements";
            }

            if ( !rff.elements().empty() )
                res.push_back( std::move( rff ) );


            return res;
        }

    /*
     * - domain dim = image dim - 1
     * - interpolation on mesh element (belong to image mesh)
     */
    template <int ApplyType, typename MeshEntityType,
              std::enable_if_t< ApplyType == -1 && std::is_same_v<MeshEntityType,typename image_mesh_type::face_type>, bool> = true >
    auto apply( MeshEntityType const& imageFace, image_vector_mesh_element_type const& imageElts ) const
        {
            std::vector<ReturnFromElements> res;

            for ( auto const& imageElt : imageElts )
            {
                ReturnFromElements rfe{ imageElt };

                std::optional<typename domain_vector_mesh_element_type::value_type> domainElt;

                if ( M_image_related_to_domain )
                {
                }
                else if ( M_domain_related_to_image )
                {
                    index_type domainElementId = M_meshDomain->meshToSubMesh( imageFace.id() );
                    CHECK( domainElementId != invalid_v<index_type> ) << "invalid relation";
                    domainElt = boost::cref( M_meshDomain->element( domainElementId ) );
                }
                if ( !domainElt )
                    return res;

                rfe.addElement( M_meshSupportDomain, *domainElt );

                if ( !rfe.elements().empty() )
                    res.push_back( std::move( rfe ) );
            }
            return res;

        }


    /*
     * - domain dim = image dim - 1
     * - interpolation on mesh element (belong to image mesh)
     */
    template <int ApplyType, typename MeshEntityType,
              std::enable_if_t< ApplyType == -1 && std::is_same_v<MeshEntityType,typename image_mesh_type::element_type>, bool> = true >
    auto apply( MeshEntityType const& imageElt, image_vector_mesh_element_type const& imageElts ) const
        {
            std::vector<ReturnFromElements> res;
            CHECK(false) << "not implemented\n";
            return res;
        }


    /*
     * - domain dim = image dim - 2
     * - interpolation on mesh element (belong to image mesh)
     */
    template <int ApplyType, typename MeshEntityType,
              std::enable_if_t< ApplyType == -2 /*&& std::is_same_v<MeshEntityType,typename image_mesh_type::face_type>*/, bool> = true >
    auto apply( MeshEntityType const& imageFace, image_vector_mesh_element_type const& imageElts ) const
        {
            std::vector<ReturnFromElements> res;
            CHECK(false) << "not implemented\n";
            return res;
        }
#if 0
       template <int ApplyType, typename EntityType>
    std::set<size_type> //domain_vector_mesh_element_type
       apply( EntityType const& imageEntity, image_vector_mesh_element_type const& imageElts,
           typename std::enable_if_t< ApplyType == 1 >* = nullptr ) const
        {
            std::set<size_type> idsFind;

            for ( auto const& imageElt : imageElts )
            {
                size_type imageEltId = unwrap_ref(imageElt).id();
                if ( M_image_related_to_domain )
                {
                    auto const& theface = M_meshDomain->face( M_meshImage->subMeshToMesh( imageEltId ) );
                    size_type domainEltId = invalid_v<size_type>;
                    if ( M_hasMeshSupportPartialDomain )
                    {
                        if ( M_meshSupportDomain->hasElement( theface.element0().id() ) && !M_meshSupportDomain->hasGhostElement( theface.element0().id() ) )
                            domainEltId = theface.element0().id();
                        else if ( theface.isConnectedTo1() && M_meshSupportDomain->hasElement( theface.element1().id() ) && !M_meshSupportDomain->hasGhostElement( theface.element1().id() ) )
                            domainEltId = theface.element1().id();
                        else
                            CHECK(false) << " error : maybe the faces is not on partition or invalid connection\n";
                    }
                    else
                    {
                        if ( !theface.element0().isGhostCell() )
                            domainEltId = theface.element0().id();
                        else if ( theface.isConnectedTo1() && !theface.element1().isGhostCell() )
                            domainEltId = theface.element1().id();
                        else
                            CHECK(false) << " error : maybe the faces is not on partition or invalid connection\n";
                    }
                    VLOG(2) << "[image_related_to_domain] image element id: "  << imageEltId << " domain element id : " << domainEltId << "\n";
                    if ( domainEltId != invalid_v<size_type> ) idsFind.insert( domainEltId );
                }
                else if( M_domain_related_to_image )
                {
#if 0
                    auto const& eltImage = M_meshImage->element(imageEltId);
                    for (uint16_type f=0;f< M_meshImage->numLocalFaces();++f)
                    {
                        const size_type idFind = M_meshDomain->meshToSubMesh( eltImage.face(f).id() );
                        if ( idFind != invalid_v<size_type> ) idsFind.insert( idFind );
                    }
                    DVLOG(2) << "[trial_related_to_test<1>] test element id: "  << imageEltId << " idsFind.size() "<< idsFind.size() << "\n";
#else
                    // imageEntity is face here
                    const size_type idFind = M_meshDomain->meshToSubMesh( imageEntity.id() );
                    if ( idFind == invalid_v<size_type> )
                        continue;
                    // warning this case can appears
                    if ( M_hasMeshSupportPartialImage && M_meshSupportDomain->hasGhostElement( idFind ) )
                        continue;
                    idsFind.insert( idFind );
#endif
                }
                else if( M_domain_sibling_of_image )
                {
                    if constexpr ( domain_mesh_type::nDim > image_mesh_type::nDim )
                                 {
                                     size_type domainEltId1 = invalid_v<size_type>;
                                     size_type domainEltId2 = invalid_v<size_type>;
                                     auto const& theface = dynamic_cast<domain_mesh_type const*>(M_meshImage->parentMesh().get())->face( M_meshImage->subMeshToMesh( imageEltId ) );
                                     if ( !theface.element0().isGhostCell() )
                                         domainEltId1 = theface.element0().id();
                                     if ( theface.isConnectedTo1() && !theface.element1().isGhostCell() )
                                         domainEltId2 = theface.element1().id();

                                     DVLOG(3) << "[image_related_to_domain] image element id: "  << imageEltId  << " : " << theface.element0().id() << " , " << theface.element1().id();
                                     // now recover the element id in domain mesh
                                     domainEltId1 = M_meshDomain->meshToSubMesh( domainEltId1 );
                                     domainEltId2 = M_meshDomain->meshToSubMesh( domainEltId2 );
                                     DVLOG(3) << "[image_related_to_domain] image element id: "  << imageEltId << " domain element id 1 : " << domainEltId1 << " domain element id 2 : " << domainEltId2 << "\n";
                                     if ( domainEltId1 != invalid_v<size_type> ) idsFind.insert( domainEltId1 );
                                     else if ( domainEltId2 != invalid_v<size_type> ) idsFind.insert( domainEltId2 );
                                 }
                    else
                    {
                        auto const& eltImage = M_meshImage->element(imageEltId);
                        for (uint16_type f=0;f< M_meshImage->numLocalFaces();++f)
                        {
                            const size_type id_in_parent_face = dynamic_cast<image_mesh_type const*>(M_meshImage->parentMesh().get())->subMeshToMesh( eltImage.face(f).id() );
                            // get now the id of the face in the domain mesh
                            const size_type idFind = M_meshDomain->meshToSubMesh( id_in_parent_face );
                            if ( idFind != invalid_v<size_type> ) idsFind.insert( idFind );
                        }
                        DVLOG(3) << "[trial_related_to_test<1>] test element id: "  << imageEltId << " idsFind.size() "<< idsFind.size() << "\n";
                    }
                }
                else // same mesh
                {
                    idsFind.insert( imageEltId );
                }
            }
            return idsFind;
        }

       template <int ApplyType, typename EntityType>
    std::set<size_type> //domain_vector_mesh_element_type
    apply( EntityType const& imageEntity, image_vector_mesh_element_type const& imageElts,
           typename std::enable_if_t< ApplyType == 2 || ApplyType == -2 >* = nullptr ) const
        {
            std::set<size_type> idsFind;
            CHECK(false) << "not implemented\n";
            return idsFind;
        }
#endif
private :
    std::shared_ptr<domain_space_type> M_XhDomain;
    std::shared_ptr<image_space_type> M_XhImage;
    std::shared_ptr<domain_mesh_type> M_meshDomain;
    std::shared_ptr<image_mesh_type> M_meshImage;
    bool M_isSameMesh;
    bool M_image_related_to_domain;
    bool M_domain_related_to_image;
    bool M_domain_sibling_of_image;
    std::shared_ptr<MeshSupport<domain_mesh_type>> M_meshSupportDomain;
    std::shared_ptr<MeshSupport<image_mesh_type>> M_meshSupportImage;
    bool M_hasMeshSupportPartialDomain, M_hasMeshSupportPartialImage;
};

} // namespace detail




/**
 * \class OperatorInterpolation
 * \brief Global interpolation operator
 *
 * @author Christophe Prud'homme
 * @see
 */

template<typename DomainSpaceType,
         typename ImageSpaceType,
         typename IteratorRange = elements_pid_t<typename ImageSpaceType::mesh_type>,
         typename InterpType = decltype(makeInterpolation<nonconforming_t>(nonconforming_t())) >
class OperatorInterpolation : public OperatorLinear<DomainSpaceType, ImageSpaceType >
{
    typedef OperatorLinear<DomainSpaceType, ImageSpaceType> super;
public:


    /** @name Typedefs
     */
    //@{

    /*
     * domain
     */
    typedef typename super::domain_space_type domain_space_type;
    typedef typename super::domain_space_ptrtype domain_space_ptrtype;
    typedef typename domain_space_type::mesh_type domain_mesh_type;
    typedef typename domain_mesh_type::element_type domain_geoelement_type;
    typedef typename domain_mesh_type::element_iterator domain_mesh_element_iterator;


    typedef typename super::backend_ptrtype backend_ptrtype;


    using size_type = typename domain_mesh_type::size_type;
    using datamap_type = DataMap<size_type>;
    using datamap_ptrtype = std::shared_ptr<datamap_type>;
    /*
     * image
     */
    typedef typename super::dual_image_space_type dual_image_space_type;
    typedef typename super::dual_image_space_ptrtype dual_image_space_ptrtype;
    typedef typename dual_image_space_type::value_type value_type;
    typedef typename dual_image_space_type::mesh_type image_mesh_type;
    typedef typename image_mesh_type::element_type image_geoelement_type;
    typedef typename image_mesh_type::element_iterator image_mesh_element_iterator;

    // geometric mapping context
    // typedef typename image_mesh_type::gm_type image_gm_type;
    // typedef typename image_mesh_type::gm_ptrtype image_gm_ptrtype;
    // using image_gmc_type = typename image_gm_type::template Context<typename image_mesh_type::element_type>;
    // using image_gmc_ptrtype = std::shared_ptr<image_gmc_type>;
    // typedef typename image_mesh_type::template gmc<vm::POINT>::type image_gmc_type;
    // typedef typename image_mesh_type::template gmc<vm::POINT>::ptrtype image_gmc_ptrtype;

    // dof
    typedef typename dual_image_space_type::dof_type dof_type;

    // basis
    typedef typename dual_image_space_type::basis_type image_basis_type;
    typedef typename domain_space_type::basis_type domain_basis_type;

    typedef typename boost::tuples::template element<0, IteratorRange>::type idim_type;
    typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_type;
    typedef IteratorRange range_iterator;


    static const uint16_type nLocalDofInDualImageElt = mpl::if_< mpl::equal_to< idim_type ,mpl::size_t<MESH_ELEMENTS> >,
                             mpl::int_< image_basis_type::nLocalDof > ,
                             mpl::int_< image_mesh_type::face_type::numVertices*dual_image_space_type::fe_type::nDofPerVertex +
                             image_mesh_type::face_type::numEdges*dual_image_space_type::fe_type::nDofPerEdge +
                             image_mesh_type::face_type::numFaces*dual_image_space_type::fe_type::nDofPerFace > >::type::value;

    // type conforme or non conforme
    using interpolation_type = std::remove_reference_t<std::remove_const_t<InterpType>>;

    // matrix graph
    typedef GraphCSR graph_type;
    typedef std::shared_ptr<graph_type> graph_ptrtype;

    // node type
    typedef typename matrix_node<typename image_mesh_type::value_type>::type matrix_node_type;

    using matrix_setup_type = OperatorInterpolationMatrixSetup<typename super::matrix_type::value_type>;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     */
    OperatorInterpolation() = default;

    /**
     * Construction the global interpolation operator from \p
     * domainspace to \p imagespace and represent it in matrix form
     * using the backend \p backend
     */
    OperatorInterpolation( domain_space_ptrtype const& domainspace,
                           dual_image_space_ptrtype const& imagespace,
                           backend_ptrtype const& backend,
                           InterpType const& interptype,
                           bool ddmethod=false,
                           matrix_setup_type const& matSetup = matrix_setup_type{} );

    OperatorInterpolation( domain_space_ptrtype const& domainspace,
                           dual_image_space_ptrtype const& imagespace,
                           IteratorRange const& r,
                           backend_ptrtype const& backend,
                           InterpType const& interptype,
                           bool ddmethod=false,
                           matrix_setup_type const& matSetup = matrix_setup_type{} );

    OperatorInterpolation( domain_space_ptrtype const& domainspace,
                           dual_image_space_ptrtype const& imagespace,
                           std::list<IteratorRange> const& r,
                           backend_ptrtype const& backend,
                           InterpType const& interptype,
                           bool ddmethod=false,
                           matrix_setup_type const& matSetup = matrix_setup_type{} );


    /**
     * copy constructor
     */
    OperatorInterpolation( OperatorInterpolation const& oi ) = default;
    OperatorInterpolation( OperatorInterpolation && oi ) = default;
    ~OperatorInterpolation() override = default;

    //@}

    /** @name Operator overloads
     */
    //@{

    OperatorInterpolation& operator=( OperatorInterpolation const& ) = default;
    OperatorInterpolation& operator=( OperatorInterpolation && ) = default;

    //@}

    /** @name Accessors
     */
    //@{

    WorldComm const& worldCommFusion() const { return *M_WorldCommFusion; }

    InterpType const& interpolationType() const { return M_interptype; }

    bool isDomainMeshRelatedToImageMesh() const { return this->domainSpace()->mesh()->isSubMeshFrom( this->dualImageSpace()->mesh() ); }

    bool isImageMeshRelatedToDomainMesh() const { return this->dualImageSpace()->mesh()->isSubMeshFrom( this->domainSpace()->mesh() ); }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * @return interpolation type
     */
    constexpr interpolation_type const& type() const noexcept { return M_interptype; }

    //@}



protected:

    virtual void update();

private:

    // determine if ghost dofs must used if the active dof in not present from the range
    std::set<size_type> defineGhostDofUsedToInterpolate( Feel::detail::DomainEltIdFromImageEltId<DomainSpaceType,ImageSpaceType> const& deifiei );

    void updateSameMesh();
    void updateNoRelationMesh();
    void updateNoRelationMeshMPI();
    void updateNoRelationMeshMPI_run(bool buildNonZeroMatrix=true);

    typedef std::vector<std::list<boost::tuple<int,
                                               typename image_mesh_type::node_type,
                                               typename image_mesh_type::node_type,
                                               std::vector<std::pair<size_type,size_type> >, // col gdof and gdofInGlobalCluster
                                               uint16_type // comp
                                               > > > extrapolation_memory_type;

    // point distrubition
    boost::tuple<std::vector< std::vector<size_type> >,
                 std::vector< std::vector<uint16_type> >,
                 std::vector<std::vector<typename image_mesh_type::node_type> >,
                 std::vector<std::vector< std::vector<typename image_mesh_type::node_type > > > >
    updateNoRelationMeshMPI_pointDistribution(const std::vector< std::list<boost::tuple<int,size_type,double> > > & memory_valueInMatrix,
                                              std::vector<std::set<size_type> > const& dof_searchWithProc,
                                              std::set<size_type> const& ghostDofUsedToInterpolate );

    // search in my world
    std::list<boost::tuple<size_type,uint16_type> >
    updateNoRelationMeshMPI_upWithMyWorld(const std::vector< std::vector<size_type> > & memmapGdof,
                                          const std::vector< std::vector<uint16_type> > & memmapComp,
                                          const std::vector<std::vector<typename image_mesh_type::node_type> > & pointsSearched,
                                          const std::vector<std::vector< std::vector<typename image_mesh_type::node_type > > > & memmap_vertices,
                                          graph_ptrtype & sparsity_graph,
                                          std::vector< std::list<boost::tuple<int,size_type,double> > > & memory_valueInMatrix,
                                          std::vector<std::map<size_type,size_type> > & memory_col_globalProcessToGlobalCluster,
                                          std::vector<std::set<size_type> > & dof_searchWithProc,
                                          bool extrapolation_mode,
                                          extrapolation_memory_type & dof_extrapolationData);

    // search in other world (MPI communication)
    std::list<boost::tuple<size_type,uint16_type> >
    updateNoRelationMeshMPI_upWithOtherWorld( boost::tuple<std::vector<rank_type>,std::vector<rank_type>,std::vector<boost::tuple<int,int> > > const& worldcommFusionProperties,
                                              std::vector< std::vector<size_type> > const& memmapGdof,
                                              std::vector< std::vector<uint16_type> > const& memmapComp,
                                              std::vector<std::vector<typename image_mesh_type::node_type> > const& pointsSearched,
                                              std::vector<std::vector< std::vector<typename image_mesh_type::node_type > > > const& memmap_vertices,
                                              graph_ptrtype & sparsity_graph,
                                              std::vector< std::list<boost::tuple<int,size_type,double> > > & memory_valueInMatrix,
                                              std::vector<std::map<size_type,size_type> > & memory_col_globalProcessToGlobalCluster,
                                              std::vector<std::set<size_type> > & dof_searchWithProc,
                                              bool extrapolation_mode,
                                              extrapolation_memory_type & dof_extrapolationData);

    std::list<boost::tuple<size_type,uint16_type> >
    updateNoRelationMeshMPI_upWithOtherWorld2( boost::tuple<std::vector<rank_type>,std::vector<rank_type>,std::vector<boost::tuple<int,int> > > const& worldcommFusionProperties,
                                               std::vector< std::vector<size_type> > const& memmapGdof,
                                               std::vector< std::vector<uint16_type> > const& memmapComp,
                                               std::vector<std::vector<typename image_mesh_type::node_type> > const& pointsSearched,
                                               std::vector<std::vector< std::vector<typename image_mesh_type::node_type > > > const& memmap_vertices,
                                               graph_ptrtype & sparsity_graph,
                                               std::vector< std::list<boost::tuple<int,size_type,double> > > & memory_valueInMatrix,
                                              std::vector<std::map<size_type,size_type> > & memory_col_globalProcessToGlobalCluster,
                                               std::vector<std::set<size_type> > & dof_searchWithProc,
                                               bool extrapolation_mode,
                                               extrapolation_memory_type & dof_extrapolationData);

    std::list<range_iterator> M_listRange;
    worldcomm_ptr_t M_WorldCommFusion;
    InterpType M_interptype;
    matrix_setup_type M_matrixSetup;
};

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::
OperatorInterpolation( domain_space_ptrtype const& domainspace,
                       dual_image_space_ptrtype const& imagespace,
                       backend_ptrtype const& backend,
                       InterpType const& interptype,
                       bool ddmethod,
                       matrix_setup_type const& matSetup )
    :
    super( domainspace, imagespace, backend, false ),
    M_listRange(),
    M_WorldCommFusion( (ddmethod || ( this->domainSpace()->worldCommPtr() == this->dualImageSpace()->worldCommPtr() ) ) ?
                       this->domainSpace()->worldCommPtr() :
                       this->domainSpace()->worldComm()+this->dualImageSpace()->worldComm() ),
    M_interptype(interptype),
    M_matrixSetup( matSetup )
{
    M_listRange.push_back( elements( imagespace->mesh() ) );
    update();
}


template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
OperatorInterpolation<DomainSpaceType,ImageSpaceType,IteratorRange,InterpType>::
OperatorInterpolation( domain_space_ptrtype const& domainspace,
                       dual_image_space_ptrtype const& imagespace,
                       IteratorRange const& r,
                       backend_ptrtype const& backend,
                       InterpType const& interptype,
                       bool ddmethod,
                       matrix_setup_type const& matSetup )
    :
    super( domainspace, imagespace, backend, false ),
    M_listRange(),
    M_WorldCommFusion( (ddmethod || ( this->domainSpace()->worldCommPtr() == this->dualImageSpace()->worldCommPtr() ) ) ?
                       this->domainSpace()->worldCommPtr() :
                       this->domainSpace()->worldComm()+this->dualImageSpace()->worldComm() ),
    M_interptype(interptype),
    M_matrixSetup( matSetup )
{
    M_listRange.push_back( r );
    update();
}

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::
OperatorInterpolation( domain_space_ptrtype const& domainspace,
                       dual_image_space_ptrtype const& imagespace,
                       std::list<IteratorRange> const& r,
                       backend_ptrtype const& backend,
                       InterpType const& interptype,
                       bool ddmethod,
                       matrix_setup_type const& matSetup )
    :
    super( domainspace, imagespace, backend, false ),
    M_listRange( r ),
    M_WorldCommFusion( (ddmethod || ( this->domainSpace()->worldCommPtr() == this->dualImageSpace()->worldCommPtr() ) ) ?
                       this->domainSpace()->worldCommPtr() :
                       this->domainSpace()->worldCommPtr()+this->dualImageSpace()->worldCommPtr() ),
    M_interptype(interptype),
    M_matrixSetup( matSetup )
{
    update();
}


template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::update()
{
    if ( this->dualImageSpace()->mesh()->numGlobalElements() == 0 || this->domainSpace()->mesh()->numGlobalElements() == 0 )
    {
        //std::cout << "OperatorInterpolation : update nothing!" << std::endl;
        //this->matPtr() = this->backend()->newZeroMatrix( this->domainSpace()->dofOnOff(),
        //                                                this->dualImageSpace()->dofOn() );
        return;
    }

    // if same mesh but not same function space (e.g. different polynomial
    // order, different basis) or if the image of domain mesh are related to
    // each other through an extraction (one of them is the sub mesh of the
    // other)
    VLOG(1) << "OperatorInterpolation: is image related to domain : " << this->dualImageSpace()->mesh()->isRelatedTo( this->domainSpace()->mesh() ) ;
    if ( this->dualImageSpace()->mesh()->isRelatedTo( this->domainSpace()->mesh() ) /*&&
                                                                                      boost::is_same<domain_mesh_type,image_mesh_type>::type::value*/ )
    {
        VLOG(1) << "OperatorInterpolation: use same mesh\n";
        VLOG(1) << "isDomainMeshRelatedToImageMesh: "  << isDomainMeshRelatedToImageMesh() << "\n";
        VLOG(1) << "isImageMeshRelatedToDomainMesh: "  << isImageMeshRelatedToDomainMesh() << "\n";
        this->updateSameMesh();
    }
    else // no relation between meshes
    {
        if ( this->dualImageSpace()->worldCommPtr()->localSize() > 1 ||
             this->domainSpace()->worldCommPtr()->localSize() > 1 )
            this->updateNoRelationMeshMPI();

        else
            this->updateNoRelationMesh();
    }

    // close matrix after build
    if ( this->matPtr() )
        this->mat().close();
}

//--------------------------------------------------------------------------------------------------//

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
std::set<typename OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::size_type>
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::defineGhostDofUsedToInterpolate( Feel::detail::DomainEltIdFromImageEltId<DomainSpaceType,ImageSpaceType> const& deifiei )
{
    std::set<size_type> ghostDofUsedToInterpolate;
    bool meshAreRelated = this->dualImageSpace()->mesh()->isRelatedTo( this->domainSpace()->mesh() );

    auto const& imagedof = this->dualImageSpace()->dof();

    std::map< rank_type, std::vector< size_type > > dataToSend, dataToRecv;
    // // init container used in send/recv
    // for ( rank_type p : this->dualImageSpace()->dof()->neighborSubdomains() )
    // {
    //     dataToSend[p].clear();
    //     dataToRecv[p].clear();
    // }

    std::vector<std::set<uint16_type> > dof_done( this->dualImageSpace()->nLocalDof(), std::set<uint16_type>() );
    std::set<size_type> activeDofSharedPresentInRange;
    auto itListRange = M_listRange.begin();
    auto const enListRange = M_listRange.end();
    for ( ; itListRange!=enListRange ; ++itListRange)
    {
        for( auto const& theImageEltWrap : *itListRange )
        {
            auto const& theImageElt = boost::unwrap_ref(theImageEltWrap);

            // if related -> check connection between domain and image element
            if ( meshAreRelated )
            {
                auto domains_eid_set = deifiei.apply( theImageElt );
                if ( domains_eid_set.empty() )
                    continue;
            }

            for ( uint16_type iloc = 0; iloc < nLocalDofInDualImageElt; ++iloc )
            {
                for ( uint16_type comp = 0; comp < image_basis_type::nComponents; ++comp )
                {
                    uint16_type compDofTableImage = (image_basis_type::is_product)? comp : 0;
                    auto thedofImage = imagedof->localToGlobal( theImageElt, iloc, compDofTableImage );
                    size_type igp = thedofImage.index();

                    if ( ( image_basis_type::is_product && dof_done[igp].empty() ) ||
                         ( !image_basis_type::is_product && dof_done[igp].find( comp ) == dof_done[igp].end() ) )
                    {
                        if ( imagedof->dofGlobalProcessIsGhost( igp ) )
                        {
                            const size_type igc = imagedof->mapGlobalProcessToGlobalCluster()[igp];
                            const rank_type theproc = imagedof->procOnGlobalCluster( igc );
                            dataToSend[theproc].push_back( igc );

                            dof_done[igp].insert( comp );
                        }
                        else if ( imagedof->activeDofSharedOnCluster().find( igp ) != imagedof->activeDofSharedOnCluster().end() )
                        {
                            activeDofSharedPresentInRange.insert(igp);
                        }
                    } // dof_done
                } // comp
            }
        }
    }

    int nbRequest=2*this->dualImageSpace()->dof()->neighborSubdomains().size();
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;

    // get size of data to transfer (first phase)
    std::map<rank_type,size_type> sizeRecv;
    for ( rank_type neighborRank : this->dualImageSpace()->dof()->neighborSubdomains() )
    {
        reqs[cptRequest++] = this->dualImageSpace()->worldCommPtr()->localComm().isend( neighborRank , 0, (size_type)dataToSend[neighborRank].size() );
        reqs[cptRequest++] = this->dualImageSpace()->worldCommPtr()->localComm().irecv( neighborRank , 0, sizeRecv[neighborRank] );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);

    // send/recv first phase
    cptRequest=0;
    for ( rank_type neighborRank : this->dualImageSpace()->dof()->neighborSubdomains() )
    {
        int nSendData = dataToSend[neighborRank].size();
        if ( nSendData > 0 )
            reqs[cptRequest++] = this->dualImageSpace()->worldCommPtr()->localComm().isend( neighborRank , 0, &(dataToSend[neighborRank][0]), nSendData );

        int nRecvData = sizeRecv[neighborRank];
        dataToRecv[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[cptRequest++] = this->dualImageSpace()->worldCommPtr()->localComm().irecv( neighborRank , 0, &(dataToRecv[neighborRank][0]), nRecvData );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);


    std::map< size_type, std::set<rank_type> > dataToTreat;
    for ( auto const& dataR : dataToRecv )
    {
        rank_type theproc = dataR.first;
        for ( auto const& dofInProc : dataR.second )
            dataToTreat[dofInProc].insert( theproc );
    }

    std::map< rank_type, std::vector< size_type > > dataToReSend, dataToReRecv;
    for ( auto const& dofAsked : dataToTreat )
    {
        size_type thedofGC = dofAsked.first;
        CHECK( imagedof->dofGlobalClusterIsOnProc( thedofGC ) ) << " thedofGC "<< thedofGC << "is not on proc\n";
        size_type thedofGP = thedofGC - imagedof->firstDofGlobalCluster();

        CHECK ( imagedof->mapGlobalProcessToGlobalCluster( thedofGP ) == thedofGC ) << "error " << imagedof->mapGlobalProcessToGlobalCluster( thedofGP ) << "must be equal to " << thedofGC;

        if ( activeDofSharedPresentInRange.find(thedofGP) == activeDofSharedPresentInRange.end() )
        {
            // define one rank to interpolate the dof
            dataToReSend[ *(dofAsked.second.begin()) ].push_back( thedofGC );
        }

    }

    // get size of data to transfer (second phase)
    cptRequest=0;
    for ( rank_type neighborRank : this->dualImageSpace()->dof()->neighborSubdomains() )
    {
        reqs[cptRequest++] = this->dualImageSpace()->worldCommPtr()->localComm().isend( neighborRank , 0, (size_type)dataToReSend[neighborRank].size() );
        reqs[cptRequest++] = this->dualImageSpace()->worldCommPtr()->localComm().irecv( neighborRank , 0, sizeRecv[neighborRank] );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);

    // send/recv second phase
    cptRequest=0;
    for ( rank_type neighborRank : this->dualImageSpace()->dof()->neighborSubdomains() )
    {
        int nSendData = dataToReSend[neighborRank].size();
        if ( nSendData > 0 )
            reqs[cptRequest++] = this->dualImageSpace()->worldCommPtr()->localComm().isend( neighborRank , 0, &(dataToReSend[neighborRank][0]), nSendData );
        int nRecvData = sizeRecv[neighborRank];
        dataToReRecv[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[cptRequest++] = this->dualImageSpace()->worldCommPtr()->localComm().irecv( neighborRank , 0, &(dataToReRecv[neighborRank][0]), nRecvData );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);
    delete [] reqs;

    for ( rank_type p : this->dualImageSpace()->dof()->neighborSubdomains() )
    {
        for ( size_type k : dataToReRecv[p] )
            ghostDofUsedToInterpolate.insert( k );
    }

    return ghostDofUsedToInterpolate;
}

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::updateSameMesh()
{

#if 0
    std::cout << "Interpolation operator (Basis info) \n"
              << "  domain_basis_type::nDof         : " << domain_basis_type::nDof << "\n"
              << "  domain_basis_type::nLocalDof    : " << domain_basis_type::nLocalDof << "\n"
              << "  domain_basis_type::nComponents  : " << domain_basis_type::nComponents << "\n"
              << "  domain_basis_type::nComponents1 : " << domain_basis_type::nComponents1 << "\n"
              << "  domain_basis_type::is_product   : " << domain_basis_type::is_product << "\n"
              << "  image_basis_type::nDof          : " << image_basis_type::nDof << "\n"
              << "  image_basis_type::nLocalDof     : " << image_basis_type::nLocalDof << "\n"
              << "  image_basis_type::nComponents   : " << image_basis_type::nComponents << "\n"
              << "  image_basis_type::nComponents1  : " << image_basis_type::nComponents1 << "\n"
              << "  image_basis_type::is_product    : " << image_basis_type::is_product << "\n"
              << "\n";
#endif
    auto const& imagedof = this->dualImageSpace()->dof();
    auto const& domaindof = this->domainSpace()->dof();

    bool needToUpdateGraph = true;
    bool needToCopyMatrix = false;
    std::shared_ptr<graph_type> sparsity_graph;
    if ( M_matrixSetup.matrix() && M_matrixSetup.matrix()->hasGraph() &&
         ( M_matrixSetup.hasContext( Feel::SAME_NONZERO_PATTERN ) || M_matrixSetup.hasContext( Feel::SUBSET_NONZERO_PATTERN ) ) )
    {
        sparsity_graph = M_matrixSetup.matrix()->graph();
        if ( M_matrixSetup.hasContext( Feel::SAME_NONZERO_PATTERN ) )
            needToUpdateGraph = false;
        else
        {
            sparsity_graph->unlock();
            needToCopyMatrix = true;
        }
    }
    else
    {
        auto dmRow = M_matrixSetup.matrix()?  M_matrixSetup.matrix()->mapRowPtr() : this->dualImageSpace()->mapPtr();
        auto dmCol = M_matrixSetup.matrix()?  M_matrixSetup.matrix()->mapColPtr() : this->domainSpace()->mapPtr();
        sparsity_graph = std::make_shared<graph_type>( dmRow, dmCol );
    }

    auto const& dofIdToContainerIdRow = sparsity_graph->mapRowPtr()->dofIdToContainerId( M_matrixSetup.indexBlockSpaceRow() );
    auto const& dofIdToContainerIdCol = sparsity_graph->mapColPtr()->dofIdToContainerId( M_matrixSetup.indexBlockSpaceCol() );

    // Local assembly: compute matrix by evaluating
    // the domain space basis function at the dual image space
    // dof points (nodal basis)
    // for Lagrange we have only computation
    // in the ref elements and the basis and dof points in ref
    // element are the same, this compute are done only one times
    // For other fe like Nedelec,Raviart-Thomas this assembly must
    // done for each element


    Feel::detail::DomainEltIdFromImageEltId<DomainSpaceType,ImageSpaceType> deifiei( this->domainSpace(), this->dualImageSpace() );

    //using range_entity_type =  std::decay_t<decltype( unwrap_ref( std::declval<typename iterator_type::value_type>() ) )>; //filter_entity_t<decltype(M_listRange.begin())>;
    using range_entity_type = typename boost::unwrap_reference<std::decay_t<typename iterator_type::value_type>>::type; //filter_entity_t<decltype(M_listRange.begin())>;
    static constexpr ElementsType deifiei_apply_on_entity = std::decay_t<decltype( deifiei.apply( std::declval<range_entity_type>() ).front() )>::onEntity();
    static constexpr int codimUsed = (deifiei_apply_on_entity == ElementsType::MESH_ELEMENTS)? 0 : 1;
    auto uDomain = this->domainSpace()->element();
    //auto expr = interpolation_type::operand(uDomain);
    auto expr = M_interptype.operand(uDomain);
    auto MlocEvalBasisNEW = Feel::detail::precomputeDomainBasisFunction<codimUsed>( this->domainSpace(), this->dualImageSpace(), expr );

    Eigen::MatrixXd IhLoc;


    // determine if ghost dofs must used if the active dof in not present from the range
    std::set<size_type> ghostDofUsedToInterpolate;
    if ( this->dualImageSpace()->worldCommPtr()->localSize() > 0 )
    {
        ghostDofUsedToInterpolate = this->defineGhostDofUsedToInterpolate( deifiei );
    }


    // we perfom 2 pass : first build matrix graph, second assembly matrix
    enum OpToApplyEnum { BUILD_GRAPH, ASSEMBLY_MATRIX };
    std::vector<OpToApplyEnum> opToApplySet;
    if ( needToUpdateGraph )
        opToApplySet.push_back( OpToApplyEnum::BUILD_GRAPH );
    opToApplySet.push_back( OpToApplyEnum::ASSEMBLY_MATRIX );
    for ( OpToApplyEnum opToApply : opToApplySet )
    {
        std::vector<std::set<uint16_type> > dof_done( this->dualImageSpace()->nLocalDof(), std::set<uint16_type>() );

        if ( opToApply == OpToApplyEnum::ASSEMBLY_MATRIX )
        {
            // create matrix
            if ( needToUpdateGraph )
            {
                VLOG(1) << "Building interpolation matrix ( " << this->domainSpace()->dofOnOff()->nDof() << "," << this->domainSpace()->dofOnOff()->nLocalDof()
                        << "," << this->dualImageSpace()->dofOn()->nDof() << ", " << this->dualImageSpace()->dofOn()->nLocalDof() << ")";
                google::FlushLogFiles(google::INFO);
                CHECK( !needToCopyMatrix ) << "TODO : copy matrix from setup";
                M_matrixSetup.setMatrix( this->backend()->newMatrix( sparsity_graph->mapColPtr(), sparsity_graph->mapRowPtr(),
                                                                     sparsity_graph ) );
                 // this->matPtr() = this->backend()->newMatrix( this->domainSpace()->dofOnOff(),
                 //                                              this->dualImageSpace()->dofOn(),
                 //                                              sparsity_graph  );
            }
            this->setMatrix( M_matrixSetup.matrix() );
        }

#if 0
        for ( auto& itListRange : M_listRange )
        {
            for( auto const& theImageEltWrap : itListRange )
            {
                auto const& theImageElt = boost::unwrap_ref(theImageEltWrap);
                auto domains_eid_set = deifiei.apply( theImageElt );
                if ( domains_eid_set.size() == 0 )
                    continue;

                for ( uint16_type iloc = 0; iloc < nLocalDofInDualImageElt; ++iloc )
                {
                    int nc = (image_basis_type::is_product)? image_basis_type::nComponents : 1;
                    for ( uint16_type comp = 0; comp < nc; ++comp )
                    {
                        uint16_type compDofTableImage = (image_basis_type::is_product)? comp : 0;
                        auto const& thedofImage = imagedof->localToGlobal( theImageElt, iloc, compDofTableImage );
                        size_type i = thedofImage.index();

                        if ( ( image_basis_type::is_product && dof_done[i].empty() ) ||
                             ( !image_basis_type::is_product && dof_done[i].find( comp ) == dof_done[i].end() ) )
                        {
                            const auto ig1 = imagedof->mapGlobalProcessToGlobalCluster()[i];
                            const auto theproc = imagedof->procOnGlobalCluster( ig1 );

                            if ( imagedof->dofGlobalProcessIsGhost( i ) &&
                                 ( ghostDofUsedToInterpolate.find(ig1) == ghostDofUsedToInterpolate.end() ) ) continue;


                            if ( opToApply == OpToApplyEnum::BUILD_GRAPH )
                            {
                                // define row in graph
                                auto& row = sparsity_graph->row( ig1 );
                                row.template get<0>() = theproc;
                                const size_type il1 = ig1 - imagedof->firstDofGlobalCluster( theproc );
                                row.template get<1>() = il1;
                            }

                            for ( auto const& domain_eid : domains_eid_set )
                            {

                                auto const& s = domaindof->localToGlobalSigns( domain_eid );

                                uint16_type ilocprime = invalid_v<uint16_type>;

                                if ( opToApply == OpToApplyEnum::ASSEMBLY_MATRIX )
                                {
                                    ilocprime = MlocEvalBasisNEW->update( this->domainSpace()->mesh()->element( domain_eid ),  theImageElt, iloc, i,comp );
                                    IhLoc = MlocEvalBasisNEW->interpolant();
                                }

                                for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                {
                                    uint16_type nCompDomain = (domain_basis_type::is_product)?domain_basis_type::nComponents:1;
                                    for(int cdomain = 0;cdomain < nCompDomain;++cdomain )
                                    {
                                        if ( domain_basis_type::is_product &&
                                             image_basis_type::is_product &&
                                             ( M_interptype.interpolationOperand() == interpolation_operand_type::ID )
                                             && cdomain != comp) continue;

                                        // get column
                                        const size_type j = domaindof->localToGlobal( domain_eid, jloc, cdomain ).index();

                                        if ( opToApply == OpToApplyEnum::BUILD_GRAPH )
                                        {
                                            //up the pattern graph
                                            auto& row = sparsity_graph->row( ig1 );
                                            row.template get<2>().insert( domaindof->mapGlobalProcessToGlobalCluster()[j] );
                                        }
                                        else if ( opToApply == OpToApplyEnum::ASSEMBLY_MATRIX )
                                        {
                                            const value_type val = s(jloc)*IhLoc( cdomain*domain_basis_type::nLocalDof+jloc,
                                                                                  ilocprime );
#if 0
                                            this->matPtr()->set( i,j,val );
#else
                                            this->matPtr()->set( dofIdToContainerIdRow[i], dofIdToContainerIdCol[j], val );
#endif
                                        }
                                    }
                                }

                            }  // for ( ; it_trial!=en_trial ; ++it_trial )

                            dof_done[i].insert( comp );
                        } // if ( !dof_done[i] )
                    } // for ( uint16_type comp ... )
                } // for ( uint16_type iloc ... )

            } // for ( ; it != en; ++ it )
        } // for ( ; itListRange!=enListRange ; ++itListRange)
#else // NEWWWWWWWWWWWWWWWWW
        for ( auto& itListRange : M_listRange )
        {
            for( auto const& theImageEltWrap : itListRange )
            {
                auto const& theImageElt = boost::unwrap_ref(theImageEltWrap);


                auto meshEntitiesImageToDomain = deifiei.apply( theImageElt );
                for ( auto const& domains_eid_set : meshEntitiesImageToDomain )
                {
                    if ( domains_eid_set.empty() )
                        continue;

                    //auto const& theImageElt = domains_eid_set.imageElement();

                    if ( opToApply == OpToApplyEnum::ASSEMBLY_MATRIX )
                    {
                        MlocEvalBasisNEW->update( domains_eid_set );
                        IhLoc = MlocEvalBasisNEW->interpolant();
                    }


                    // TODO get local dof from elt and maybe with face/edge/point id
                    for( auto const& ldof : imagedof->localDof( theImageElt ) )
                    {
                        index_type i = invalid_v<index_type>;
                        uint16_type iloc = invalid_v<uint16_type>;
                        if constexpr ( idim_type::value == MESH_FACES )
                        {
                            i = ldof.index();
                            iloc = ldof.localDof(); // WRONG HERE IN VECTORIAL, SEE FIX BELOW
                            //comp = ldof.first.component( image_basis_type::nLocalDof ); // TODO

                            uint16_type ccc = ldof.localDofInFace()/nLocalDofInDualImageElt;//image_basis_type::nLocalEdgeDof;
                            iloc = ldof.localDof()+ccc*image_basis_type::nLocalDof;
                        }
                        else
                        {
                            i = ldof.second.index();
                            iloc = ldof.first.localDof();
                            //comp = ldof.first.component( image_basis_type::nLocalDof ); // TODO
                        }
                        uint16_type comp = iloc/image_basis_type::nLocalDof;

                        if ( dof_done[i].empty() )
                        {
                            const auto ig1 = imagedof->mapGlobalProcessToGlobalCluster()[i];
                            const auto theproc = imagedof->procOnGlobalCluster( ig1 );

                            if ( imagedof->dofGlobalProcessIsGhost( i ) &&
                                 ( ghostDofUsedToInterpolate.find(ig1) == ghostDofUsedToInterpolate.end() ) ) continue;


                            if ( opToApply == OpToApplyEnum::BUILD_GRAPH )
                            {
                                // define row in graph
                                auto& row = sparsity_graph->row( ig1 );
                                row.template get<0>() = theproc;
                                const size_type il1 = ig1 - imagedof->firstDofGlobalCluster( theproc );
                                row.template get<1>() = il1;
                            }

                            for ( auto const& domain_elt : domains_eid_set.elements() )
                            {
                                index_type domain_eid = domain_elt.element().id();
                                auto const& s = domaindof->localToGlobalSigns( domain_eid );
                                for( auto const& ldof_domain : domaindof->localDof( domain_eid ) )
                                {
                                    index_type j = ldof_domain.second.index();
                                    uint16_type jloc = ldof_domain.first.localDof();
                                    uint16_type cdomain = ldof_domain.first.component( domain_basis_type::nLocalDof );

                                    if ( domain_basis_type::is_product &&
                                         image_basis_type::is_product &&
                                         ( M_interptype.interpolationOperand() == interpolation_operand_type::ID )
                                         && cdomain != comp ) continue;

                                    if ( opToApply == OpToApplyEnum::BUILD_GRAPH )
                                    {
                                        //up the pattern graph
                                        auto& row = sparsity_graph->row( ig1 );
                                        row.template get<2>().insert( domaindof->mapGlobalProcessToGlobalCluster()[j] );
                                    }
                                    else if ( opToApply == OpToApplyEnum::ASSEMBLY_MATRIX )
                                    {
                                        const value_type val = s(jloc)*IhLoc( jloc,iloc );
                                        this->matPtr()->set( dofIdToContainerIdRow[i], dofIdToContainerIdCol[j], val );
                                    }
                                }
                                //break;/////WARING ONLY FOR CONTINUE!!!!
                            }  // for ( auto const& domain_elt : domains_eid_set.elements() )

                            dof_done[i].insert( comp );
                        } // if ( !dof_done[i] )
                    } // for( auto const& ldof : imagedof->localDof( theImageElt ) )
                } /// for ( auto const& domains_eid_set : meshEntitiesImageToDomain )

            } // for ( ; it != en; ++ it )
        } // for ( ; itListRange!=enListRange ; ++itListRange)

#endif
        if ( opToApply == OpToApplyEnum::BUILD_GRAPH )
        {
            // compute graph
            sparsity_graph->close();
        }
    } // opToApply

}

//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::updateNoRelationMesh()
{
    DVLOG(2) << "[interpolate] different meshes\n";
    //std::cout << "OperatorInterpolation::updateNoRelationMesh start " << std::endl;

    const size_type proc_id = this->dualImageSpace()->mesh()->worldCommPtr()->localRank();
    const size_type n1_dof_on_proc = this->dualImageSpace()->nLocalDof();
    const size_type firstrow_dof_on_proc = this->dualImageSpace()->dof()->firstDof( proc_id );
    const size_type lastrow_dof_on_proc = this->dualImageSpace()->dof()->lastDof( proc_id );
    const size_type firstcol_dof_on_proc = this->domainSpace()->dof()->firstDof( proc_id );
    const size_type lastcol_dof_on_proc = this->domainSpace()->dof()->lastDof( proc_id );
#if 0
    graph_ptrtype sparsity_graph( new graph_type( n1_dof_on_proc,
                                                  firstrow_dof_on_proc, lastrow_dof_on_proc,
                                                  firstcol_dof_on_proc, lastcol_dof_on_proc,
                                                  this->dualImageSpace()->mesh()->worldCommPtr()->subWorldCommSeq() ) );
#else
    graph_ptrtype sparsity_graph( new graph_type( this->dualImageSpace()->dof(), this->domainSpace()->dof() ) );
#endif

    auto const* imagedof = this->dualImageSpace()->dof().get();
    auto const* domaindof = this->domainSpace()->dof().get();
    auto const* domainbasis = this->domainSpace()->basis().get();

    //-----------------------------------------
    //init the localization tool
    //auto locTool = this->domainSpace()->mesh()->tool_localization();
    auto locTool = support(this->domainSpace())->tool_localization();
    if ( this->interpolationType().onlyLocalizeOnBoundary() )
        locTool->updateForUseBoundaryFaces();
    else
        locTool->updateForUse();
    // kdtree parameter
    locTool->kdtree()->nbNearNeighbor(this->interpolationType().nbNearNeighborInKdTree());
#if 0
    bool notUseOptLocTest = domain_mesh_type::nDim!=domain_mesh_type::nRealDim;
#else
    bool notUseOptLocTest = true;
#endif
    //if (notUseOptLocTest) locTool->kdtree()->nbNearNeighbor(domain_mesh_type::element_type::numPoints);

    //locTool->kdtree()->nbNearNeighbor(3);
    //locTool->kdtree()->nbNearNeighbor(this->domainSpace()->mesh()->numElements());
    //locTool->setExtrapolation(false);

    //-----------------------------------------
    // usefull data
    matrix_node_type ptsReal( image_mesh_type::nRealDim, 1 );
    matrix_node_type ptsRef( domain_mesh_type::nDim , 1 );
    typename Localization<domain_mesh_type>::container_search_iterator_type itanal,itanal_end;
    typename Localization<domain_mesh_type>::container_output_iterator_type itL,itL_end;
    matrix_node_type MlocEval( domain_basis_type::nLocalDof*domain_basis_type::nComponents1,1 );

    std::vector<bool> dof_done( this->dualImageSpace()->nLocalDof(), false );
    std::vector< std::list<std::pair<size_type,double> > > memory_valueInMatrix( this->dualImageSpace()->nLocalDof() );

    //-----------------------------------------
    size_type eltIdLocalised = invalid_v<size_type>;

    // for each element in range
    auto itListRange = M_listRange.begin();
    auto const enListRange = M_listRange.end();
    for ( ; itListRange!=enListRange ; ++itListRange)
    {
        //iterator_type it, en;
        //boost::tie( boost::tuples::ignore, it, en ) = *itListRange;
        //for ( ; it != en; ++ it )
        for( auto const& theImageEltWrap : *itListRange )
        {
            //auto const& theImageElt = boost::unwrap_ref(*it);
            auto const& theImageElt = boost::unwrap_ref(theImageEltWrap);
            for ( uint16_type iloc = 0; iloc < nLocalDofInDualImageElt; ++iloc )
                {
                    for ( uint16_type comp = 0; comp < image_basis_type::nComponents; ++comp )
                        {
                            const size_type gdof = imagedof->localToGlobal( theImageElt, iloc, comp ).index();
                            if (!dof_done[gdof])
                                {
                                    //------------------------
                                    // get the graph row
                                    const auto ig1 = imagedof->mapGlobalProcessToGlobalCluster()[gdof];
                                    const auto theproc = imagedof->procOnGlobalCluster( ig1 );
                                    auto& row = sparsity_graph->row(ig1);
                                    row.template get<0>() = theproc;
                                    row.template get<1>() = gdof;
                                    //------------------------
                                    // the dof point
                                    ublas::column(ptsReal,0 ) = boost::get<0>(imagedof->dofPoint(gdof));
                                    //------------------------
                                    // localisation process
                                    //if (notUseOptLocTest) eltIdLocalised=invalid_v<size_type>;
                                    auto resLocalisation = locTool->run_analysis(ptsReal,eltIdLocalised,theImageElt.vertices()/*theImageElt.G()*/,mpl::int_<interpolation_type::isConforming()>());
                                    for ( bool hasFindPtLocalised : resLocalisation.template get<0>()  )
                                         LOG_IF(ERROR, !hasFindPtLocalised ) << "OperatorInterpolation::updateNoRelationMesh : point localisation fail!\n";
                                    if ( !notUseOptLocTest )
                                        eltIdLocalised = resLocalisation.template get<1>();
                                    //------------------------
                                    // for each localised points
                                    itanal = locTool->result_analysis_begin();
                                    itanal_end = locTool->result_analysis_end();
                                    for ( ;itanal!=itanal_end;++itanal)
                                        {
                                            itL=itanal->second.begin();
                                            ublas::column( ptsRef, 0 ) = boost::get<1>( *itL );

                                            MlocEval = domainbasis->evaluate( ptsRef );

                                            for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                {
                                                    //get global dof
                                                    size_type j =  domaindof->localToGlobal( itanal->first,jloc,comp ).index();
                                                    value_type v = MlocEval( domain_basis_type::nComponents1*jloc
                                                                             + comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof
                                                                             + comp,
                                                                             0 );
                                                    row.template get<2>().insert( domaindof->mapGlobalProcessToGlobalCluster()[j] );
                                                    memory_valueInMatrix[gdof].push_back( std::make_pair( j,v ) );
                                                }
                                        }
                                    dof_done[gdof]=true;
                                } // if (!dof_done[gdof])
                        } //  for ( uint16_type comp = 0; comp < image_basis_type::nComponents; ++comp )
                } // for ( uint16_type iloc = 0; iloc < nLocalDofInDualImageElt; ++iloc )
        } // for( ; it != en; ++ it )
    } // for ( ; itListRange!=enListRange ; ++itListRange)


    //-----------------------------------------
    // compute graph
    sparsity_graph->close(); //sparsity_graph->printPython("mygraphpython.py");
    //-----------------------------------------
    // create matrix
    this->matPtr() = this->backend()->newMatrix( this->domainSpace()->dofOnOff(),
                                                 this->dualImageSpace()->dofOn(),
                                                 sparsity_graph );
    //-----------------------------------------
    // assemble matrix
    for (size_type idx_i=0 ; idx_i<this->dualImageSpace()->nLocalDof() ;++idx_i)
        {
            for (auto it_j=memory_valueInMatrix[idx_i].begin(),en_j=memory_valueInMatrix[idx_i].end() ; it_j!=en_j ; ++it_j)
                {
                    this->matPtr()->set(idx_i,it_j->first,it_j->second);
                }
        }

}

//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::updateNoRelationMeshMPI()
{
    //std::cout << "OperatorInterpolation::updateNoRelationMeshMPI start " << std::endl;

    auto testCommActivities_image=this->dualImageSpace()->worldCommPtr()->hasMultiLocalActivity();

    if (testCommActivities_image.template get<0>())
        {
            //std::cout << "OperatorInterpolation::updateNoRelationMeshMPI hasMultiLocalActivity " << std::endl;
            // save initial activities
            std::vector<int> saveActivities_image = this->dualImageSpace()->worldCommPtr()->activityOnWorld();
            // iterate on each local activity
            const auto colorWhichIsActive = testCommActivities_image.template get<1>();
            auto it_color=colorWhichIsActive.begin();
            auto const en_color=colorWhichIsActive.end();
            for ( ;it_color!=en_color;++it_color )
                {
                    this->dualImageSpace()->worldCommPtr()->applyActivityOnlyOn( *it_color );
                    this->dualImageSpace()->mapOn().worldCommPtr()->applyActivityOnlyOn( *it_color );
                    this->dualImageSpace()->mapOnOff().worldCommPtr()->applyActivityOnlyOn( *it_color );
                    this->updateNoRelationMeshMPI_run(false);
                }
            // revert initial activities
            this->dualImageSpace()->worldCommPtr()->setIsActive(saveActivities_image);
            this->dualImageSpace()->mapOn().worldCommPtr()->setIsActive(saveActivities_image);
            this->dualImageSpace()->mapOnOff().worldCommPtr()->setIsActive(saveActivities_image);
        }
    else
        {
            //std::cout << "OperatorInterpolation::updateNoRelationMeshMPI has One LocalActivity " << std::endl;
            if ( !this->dualImageSpace()->worldCommPtr()->isActive() && !this->domainSpace()->worldCommPtr()->isActive() )
                {
                    this->matPtr() = this->backend()->newZeroMatrix( this->domainSpace()->dofOnOff(),
                                                                     this->dualImageSpace()->dofOn() );
                }
            else
                {
                    this->updateNoRelationMeshMPI_run(true);
                }
        }

}

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::updateNoRelationMeshMPI_run(bool buildNonZeroMatrix)
{

    //std::cout << "OperatorInterpolation::updateNoRelationMeshMPI_run start " << std::endl;

    //-----------------------------------------------------------------------------------------
    // PreProcess : datamap properties and graph
    //-----------------------------------------------------------------------------------------

    const int proc_id = this->dualImageSpace()->worldCommPtr()->localRank();
    const int proc_id_row = this->dualImageSpace()->worldCommPtr()->localRank();
    const int proc_id_col = this->domainSpace()->worldCommPtr()->localRank();
    const int nProc = this->dualImageSpace()->mesh()->worldCommPtr()->size();
    const int nProc_row = this->dualImageSpace()->mesh()->worldCommPtr()->localSize();
    const int nProc_col = this->domainSpace()->mesh()->worldCommPtr()->localSize();
    const int nProc_image = this->dualImageSpace()->mesh()->worldCommPtr()->localSize();
    const int nProc_domain = this->domainSpace()->mesh()->worldCommPtr()->localSize();
    const size_type nrow_dof_on_proc = this->dualImageSpace()->nLocalDof();
    const size_type firstrow_dof_on_proc = this->dualImageSpace()->dof()->firstDofGlobalCluster( proc_id_row );
    const size_type lastrow_dof_on_proc = this->dualImageSpace()->dof()->lastDofGlobalCluster( proc_id_row );
    const size_type firstcol_dof_on_proc = this->domainSpace()->dof()->firstDofGlobalCluster( proc_id_col );
    const size_type lastcol_dof_on_proc = this->domainSpace()->dof()->lastDofGlobalCluster( proc_id_col );


    graph_ptrtype sparsity_graph( new graph_type(this->dualImageSpace()->dof(), this->domainSpace()->dof() ) );


    size_type new_nLocalDofWithoutGhost=this->domainSpace()->nDof()/nProc_row;
    size_type new_nLocalDofWithoutGhost_tempp=this->domainSpace()->nDof()/nProc_row;
    size_type new_nLocalDofWithoutGhostMiss=this->domainSpace()->nDof()%nProc_row;
    if (this->dualImageSpace()->worldCommPtr()->globalSize()==this->domainSpace()->worldCommPtr()->globalSize() )
    {
        new_nLocalDofWithoutGhost = this->domainSpace()->mapOnOff().nLocalDofWithoutGhost();
    }
    else
    {
        if (this->dualImageSpace()->worldCommPtr()->globalRank()==this->dualImageSpace()->worldCommPtr()->masterRank())  new_nLocalDofWithoutGhost+=new_nLocalDofWithoutGhostMiss;
    }

    size_type new_firstdofcol=0,new_lastdofcol=new_nLocalDofWithoutGhost-1;
    bool findMyProc=false;
    int currentProc=0;
    while(!findMyProc)
        {
            if (currentProc==this->dualImageSpace()->worldCommPtr()->globalRank())
                {
                    findMyProc=true;
                }
            else if (currentProc==this->dualImageSpace()->worldCommPtr()->masterRank())
                {
                    new_firstdofcol+=new_nLocalDofWithoutGhost_tempp+new_nLocalDofWithoutGhostMiss;
                    new_lastdofcol+=new_nLocalDofWithoutGhost_tempp+new_nLocalDofWithoutGhostMiss;
                }
            else
                {
                    new_firstdofcol+=new_nLocalDofWithoutGhost_tempp;
                    new_lastdofcol+=new_nLocalDofWithoutGhost_tempp;
                }
            ++currentProc;
        }


    //-----------------------------------------------------------------------------------------
    // relation between domain/image worldcomm and the fusion worldcomm
    //-----------------------------------------------------------------------------------------

    //std::vector<boost::tuple<int,int,int,int> > worldcommFusionProperties;
    boost::tuple<std::vector<rank_type>,std::vector<rank_type>,std::vector<boost::tuple<int,int> > > worldcommFusionProperties;
    if ( this->domainSpace()->worldCommPtr() == this->dualImageSpace()->worldCommPtr() )
    {
        std::vector<rank_type> localMeshRankToWorldCommFusion_domain(nProc_col);
        for( rank_type p = 0 ; p<nProc_col ; ++p )
            localMeshRankToWorldCommFusion_domain[p]=p;
        std::vector<rank_type> localMeshRankToWorldCommFusion_image(nProc_row);
        for( rank_type p = 0 ; p<nProc_row ; ++p )
            localMeshRankToWorldCommFusion_image[p]=p;
        std::vector<boost::tuple<int,int> > procActivitiesOnWorldCommFusion(this->worldCommFusion().globalSize(),boost::make_tuple((int)true,(int)true));
        worldcommFusionProperties.template get<0>() = localMeshRankToWorldCommFusion_domain;
        worldcommFusionProperties.template get<1>() = localMeshRankToWorldCommFusion_image;
        worldcommFusionProperties.template get<2>() = procActivitiesOnWorldCommFusion;
    }
    else if ( this->interpolationType().searchWithCommunication())
    {
        // Attention : marche que si les 2 worldcomms qui s'emboite (mon cas)
        std::vector<rank_type> localMeshRankToWorldCommFusion_domain(nProc_col);
        mpi::all_gather( this->domainSpace()->mesh()->worldCommPtr()->localComm(),
                         this->worldCommFusion().globalRank(),
                         localMeshRankToWorldCommFusion_domain );
        std::vector<rank_type> localMeshRankToWorldCommFusion_image(nProc_row);
        mpi::all_gather( this->dualImageSpace()->mesh()->worldCommPtr()->localComm(),
                         this->worldCommFusion().globalRank(),
                         localMeshRankToWorldCommFusion_image );

        std::vector<boost::tuple<int,int> > procActivitiesOnWorldCommFusion(this->worldCommFusion().globalSize());
        auto dataSendToAllGather = boost::make_tuple( (int)this->domainSpace()->worldCommPtr()->isActive(),(int)this->dualImageSpace()->worldCommPtr()->isActive() );
        mpi::all_gather( this->worldCommFusion().globalComm(),
                         dataSendToAllGather,
                         procActivitiesOnWorldCommFusion );


        //----------------------------------------------//
        // correction to apply if ....
        int firstActiveProc_image=0;
        bool findFirstActive_image=false;
        while (!findFirstActive_image)
        {
            if (procActivitiesOnWorldCommFusion[firstActiveProc_image].template get<1>() )  // if (imageProcIsActive_fusion[firstActiveProc_image])
            {
                findFirstActive_image=true;
            }
            else ++firstActiveProc_image;
        }
        int firstActiveProc_domain=0;
        bool findFirstActive_domain=false;
        while (!findFirstActive_domain)
        {
            if (procActivitiesOnWorldCommFusion[firstActiveProc_domain].template get<0>() ) //if (domainProcIsActive_fusion[firstActiveProc_domain])
            {
                findFirstActive_domain=true;
            }
            else ++firstActiveProc_domain;
        }

        for (int p=0;p<localMeshRankToWorldCommFusion_image.size(); ++p)
        {
            if (!this->dualImageSpace()->worldCommPtr()->isActive())
                localMeshRankToWorldCommFusion_image[p]=p%nProc_image+firstActiveProc_image; // FAIRE COMMMUNICATION!!!!!
        }
        for (int p=0;p<localMeshRankToWorldCommFusion_domain.size(); ++p)
        {
            if (!this->domainSpace()->worldCommPtr()->isActive())
                localMeshRankToWorldCommFusion_domain[p]=p%nProc_domain+firstActiveProc_domain; // FAIRE COMMMUNICATION!!!!!
        }
        //----------------------------------------------//
        // init worldcommFusionProperties
        worldcommFusionProperties.template get<0>() = localMeshRankToWorldCommFusion_domain;
        worldcommFusionProperties.template get<1>() = localMeshRankToWorldCommFusion_image;
        worldcommFusionProperties.template get<2>() = procActivitiesOnWorldCommFusion;

    } // if (this->interpolationType().searchWithCommunication())

    //-----------------------------------------------------------------------------------------
    // determine if ghost dofs must used if the active dof in not present from the range
    //-----------------------------------------------------------------------------------------
    Feel::detail::DomainEltIdFromImageEltId<DomainSpaceType,ImageSpaceType> deifiei( this->domainSpace(), this->dualImageSpace() );
    std::set<size_type> ghostDofUsedToInterpolate = this->defineGhostDofUsedToInterpolate( deifiei );

    //-----------------------------------------------------------------------------------------
    // Start localization process
    //-----------------------------------------------------------------------------------------

    //memory containers
    std::vector< std::list<boost::tuple<int,size_type,double> > > memory_valueInMatrix( this->dualImageSpace()->nLocalDof() );
    for (size_type k=0;k<memory_valueInMatrix.size();++k) memory_valueInMatrix[k].clear();

    std::vector<std::map<size_type,size_type> > memory_col_globalProcessToGlobalCluster(nProc_col);
    std::vector<std::set<size_type> > dof_searchWithProc(this->dualImageSpace()->nLocalDof());

    extrapolation_memory_type dof_extrapolationData(this->dualImageSpace()->nLocalDof());

    //init the localization tool
    //auto locTool = this->domainSpace()->mesh()->tool_localization();
    auto locTool = support(this->domainSpace())->tool_localization();
    bool doExtrapolationAtStart = locTool->doExtrapolation();
    // kdtree parameter
    locTool->kdtree()->nbNearNeighbor(this->interpolationType().nbNearNeighborInKdTree());
    // points in kdtree
    if ( this->interpolationType().onlyLocalizeOnBoundary() ) locTool->updateForUseBoundaryFaces();
    else locTool->updateForUse();
    // no extrapolation in first
    if ( doExtrapolationAtStart && this->interpolationType().searchWithCommunication() ) locTool->setExtrapolation(false);


#if 0
    uint16_type nMPIsearch=15;//5;
#else
    // TODO : put 15 not work when a partition is not gather in one group.
    // we should compute these group and compute barycenter in each of them
    uint16_type nMPIsearch = this->domainSpace()->mesh()->worldCommPtr()->localSize();
#endif
    if( InterpType::isConforming()) nMPIsearch=this->domainSpace()->mesh()->worldCommPtr()->localSize();
    else if (this->domainSpace()->mesh()->worldCommPtr()->localSize()<nMPIsearch) nMPIsearch=this->domainSpace()->mesh()->worldCommPtr()->localSize();
   // only one int this case
   if (!this->interpolationType().searchWithCommunication()) nMPIsearch=1;
   uint16_type counterMPIsearch=1;
   bool FinishMPIsearch=false;
   //if (this->interpolationType().searchWithCommunication()) FinishMPIsearch=true;// not run algo1 !!!!

   boost::mpi::timer mytimer;
   this->worldCommFusion().globalComm().barrier();
   if ( this->worldCommFusion().globalRank() == this->dualImageSpace()->worldCommPtr()->masterRank() )
       LOG(INFO) << " start while " << std::endl;

   size_type nbLocalisationFail=1;
   while(!FinishMPIsearch)
       {
           mytimer.restart();
           auto pointDistribution = this->updateNoRelationMeshMPI_pointDistribution(memory_valueInMatrix,dof_searchWithProc,ghostDofUsedToInterpolate);
           auto memmapGdof = pointDistribution.template get<0>();
           auto memmapComp = pointDistribution.template get<1>();
           auto pointsSearched = pointDistribution.template get<2>();
           auto memmapVertices = pointDistribution.template get<3>();
           //std::cout <<  "proc " << this->worldCommFusion().globalRank() <<  " pointsSearched.size() " << pointsSearched.size() << std::endl;
           double t1 = mytimer.elapsed();
           if ( this->worldCommFusion().globalRank() == this->dualImageSpace()->worldCommPtr()->masterRank() )
               LOG(INFO) << "finish-step1 in " << (boost::format("%1%") % t1).str() << std::endl;
           //this->worldCommFusion().globalComm().barrier();
           mytimer.restart();

           auto memory_localisationFail = this->updateNoRelationMeshMPI_upWithMyWorld( memmapGdof, // input
                                                                                       memmapComp, // input
                                                                                       pointsSearched, // input
                                                                                       memmapVertices, // input
                                                                                       sparsity_graph, // output
                                                                                       memory_valueInMatrix, // output
                                                                                       memory_col_globalProcessToGlobalCluster, // output
                                                                                       dof_searchWithProc, // output,
                                                                                       false,// extrapolation_mode
                                                                                       dof_extrapolationData // empty
                                                                                       );
           //std::cout <<  "proc " << this->worldCommFusion().globalRank() <<  " memory_localisationFail.size() " << memory_localisationFail.size() << std::endl;
           double t2 = mytimer.elapsed();
           if ( this->worldCommFusion().globalRank() == this->dualImageSpace()->worldCommPtr()->masterRank() )
               LOG(INFO) << "finish-step2 in " << (boost::format("%1%") % t2).str() << std::endl;
           //this->worldCommFusion().globalComm().barrier();
           mytimer.restart();


           if (this->interpolationType().searchWithCommunication())
               {
                   auto memory_localisationFail2 = this->updateNoRelationMeshMPI_upWithOtherWorld2( worldcommFusionProperties, //input
                                                                                                    memmapGdof, // input
                                                                                                    memmapComp, // input
                                                                                                    pointsSearched, // input
                                                                                                    memmapVertices, // input
                                                                                                    sparsity_graph, // output
                                                                                                    memory_valueInMatrix, // output
                                                                                                    memory_col_globalProcessToGlobalCluster, // output
                                                                                                    dof_searchWithProc, // output
                                                                                                    false,// extrapolation_mode
                                                                                                    dof_extrapolationData // empty
                                                                                                    );

                   const size_type nbLocalisationFail_loc = memory_localisationFail.size()+memory_localisationFail2.size();
                   mpi::all_reduce( this->worldCommFusion().globalComm(),
                                    nbLocalisationFail_loc,
                                    nbLocalisationFail,
                                    std::plus<double>() );
               }
           else nbLocalisationFail = memory_localisationFail.size();
           //nbLocalisationFail = memory_localisationFail.size()+memory_localisationFail2.size();

           double t3 = mytimer.elapsed();
           if ( this->worldCommFusion().globalRank() == this->dualImageSpace()->worldCommPtr()->masterRank() )
               LOG(INFO) << "finish-step3 in " << (boost::format("%1%") % t3).str() << std::endl;
           mytimer.restart();

           if ( this->worldCommFusion().globalRank() == this->dualImageSpace()->worldCommPtr()->masterRank() )
               LOG(INFO) << " it " << counterMPIsearch << "  nbLocalisationFail " << nbLocalisationFail << std::endl;

           //std::cout <<  "proc " << this->worldCommFusion().globalRank()
           //          << " et " <<nbLocalisationFail << std::endl;
           if (counterMPIsearch<nMPIsearch && nbLocalisationFail>0) ++counterMPIsearch;
           else FinishMPIsearch=true;
       }

   if ( doExtrapolationAtStart && this->interpolationType().searchWithCommunication() ) locTool->setExtrapolation(true);

   if ( doExtrapolationAtStart && nbLocalisationFail>0 )
       {
           LOG(INFO) << " Start Extrapolation" << std::endl;
           std::vector<std::set<size_type> > dof_searchWithProcExtrap(this->dualImageSpace()->nLocalDof());
           //locTool->setExtrapolation(true);
           uint16_type nMPIsearchExtrap=5;
           if (this->domainSpace()->mesh()->worldCommPtr()->localSize()<5) nMPIsearchExtrap=this->domainSpace()->mesh()->worldCommPtr()->localSize();
           // only one int this case
           if (!this->interpolationType().searchWithCommunication()) nMPIsearchExtrap=1;
           uint16_type counterMPIsearchExtrap=1;
           bool FinishMPIsearchExtrap=false;
           // localisation process
           while(!FinishMPIsearchExtrap)
               {
                   auto pointDistribution = this->updateNoRelationMeshMPI_pointDistribution(memory_valueInMatrix,dof_searchWithProcExtrap,ghostDofUsedToInterpolate);
                   auto memmapGdof = pointDistribution.template get<0>();
                   auto memmapComp = pointDistribution.template get<1>();
                   auto pointsSearched = pointDistribution.template get<2>();
                   auto memmapVertices = pointDistribution.template get<3>();
                   //std::cout <<  "proc " << this->worldCommFusion().globalRank() <<  " pointsSearched.size() " << pointsSearched.size() << std::endl;
                   auto memory_localisationFail = this->updateNoRelationMeshMPI_upWithMyWorld( memmapGdof, // input
                                                                                               memmapComp, // input
                                                                                               pointsSearched, // input
                                                                                               memmapVertices, // input
                                                                                               sparsity_graph, // output
                                                                                               memory_valueInMatrix, // output
                                                                                               memory_col_globalProcessToGlobalCluster, // output
                                                                                               dof_searchWithProcExtrap, // output,
                                                                                               true,// extrapolation_mode
                                                                                               dof_extrapolationData // output
                                                                                               );
                   //std::cout <<  "proc " << this->worldCommFusion().globalRank() <<  " memory_localisationFail.size() " << memory_localisationFail.size() << std::endl;
                   if (this->interpolationType().searchWithCommunication())
                       {
                           auto memory_localisationFail2 = this->updateNoRelationMeshMPI_upWithOtherWorld2( worldcommFusionProperties, //input
                                                                                                           memmapGdof, // input
                                                                                                           memmapComp, // input
                                                                                                           pointsSearched, // input
                                                                                                           memmapVertices, // input
                                                                                                           sparsity_graph, // output
                                                                                                           memory_valueInMatrix, // output
                                                                                                           memory_col_globalProcessToGlobalCluster, // output
                                                                                                           dof_searchWithProcExtrap, // output
                                                                                                           true, // extrapolation_mode
                                                                                                           dof_extrapolationData // output
                                                                                                           );
                       }
                   //std::cout <<  "proc " << this->worldCommFusion().globalRank() <<  " memory_localisationFail2.size() " << memory_localisationFail.size() << std::endl;
                   if (counterMPIsearchExtrap<nMPIsearchExtrap) ++counterMPIsearchExtrap;
                   else FinishMPIsearchExtrap=true;
               } // while(!FinishMPIsearchExtrap)

           auto const* imagedof = this->dualImageSpace()->dof().get();
           auto const* domaindof = this->domainSpace()->dof().get();
           auto const* domainbasis = this->domainSpace()->basis().get();
           matrix_node_type ptsRef( domain_mesh_type::nRealDim , 1 );
           matrix_node_type MlocEval(domain_basis_type::nLocalDof*domain_basis_type::nComponents1,1);

           // analysis result
           auto it_extrap = dof_extrapolationData.begin();
           auto const en_extrap = dof_extrapolationData.end();
           for (size_type cpt_gdof=0 ; it_extrap!=en_extrap ; ++it_extrap,++cpt_gdof )
               {
                   auto const gdof = cpt_gdof;
                   auto const pointExtrap = imagedof->dofPoint(gdof).template get<0>();

                   // search nearer
                   int procExtrapoled = 0;
                   double distMin= INT_MAX;
                   auto it_proc = it_extrap->begin();
                   auto const en_proc = it_extrap->end();
                   for ( ; it_proc!=en_proc; ++it_proc)
                       {
                           auto const procCurrent = it_proc->template get<0>();
                           auto const bary = it_proc->template get<1>();
                           // COMPUTE DISTANCE
                           double normDist = 0;
                           for (int q=0;q<image_mesh_type::nRealDim;++q) normDist+=std::pow(bary(q)-pointExtrap(q),2);
                           //std::cout << std::sqrt(normDist) << " "<< bary.size() << " " << pointExtrap.size() << std::endl;
                           if (std::sqrt(normDist)<distMin) { distMin=std::sqrt(normDist);procExtrapoled=procCurrent;}
                       }
                   it_proc = it_extrap->begin();
                   for ( ; it_proc!=en_proc; ++it_proc)
                       {
                           if ( it_proc->template get<0>()==procExtrapoled)
                               {
                                   // get the graph row
                                   auto const ig1 = imagedof->mapGlobalProcessToGlobalCluster()[gdof];
                                   auto const theproc = imagedof->procOnGlobalCluster(ig1);
                                   auto& row = sparsity_graph->row(ig1);
                                   row.template get<0>() = theproc;
                                   row.template get<1>() = ig1 - imagedof->firstDofGlobalCluster(theproc);
                                   // get ref point
                                   ublas::column( ptsRef, 0 ) = it_proc->template get<2>();
                                   // evaluate basis functions for this point
                                   MlocEval = domainbasis->evaluate( ptsRef );
                                   //auto it_jdof = it_proc->template get<3>().begin();
                                   //auto const en_jdof = it_proc->template get<3>().end();
                                   auto const comp = it_proc->template get<4>();
                                   auto const vec_jloc = it_proc->template get<3>();
                                   for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                       //for (uint16_type jloc=0 ; it_jdof != en_jdof ; ++it_jdof,++jloc)
                                       {
                                           //const size_type j_gdof=it_jdof->first;
                                           //const size_type j_gdofGlobalCluster=it_jdof->second;
                                           const size_type j_gdof=vec_jloc[jloc].first;
                                           const size_type j_gdofGlobalCluster=vec_jloc[jloc].second;
                                           row.template get<2>().insert(j_gdofGlobalCluster);//domaindof->mapGlobalProcessToGlobalCluster()[j_gdof]);
                                           // get value
                                           value_type v = MlocEval( domain_basis_type::nComponents1*jloc +
                                                                    comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof +
                                                                    comp, 0 );
                                           // save value
                                           memory_valueInMatrix[gdof].push_back(boost::make_tuple(procExtrapoled,j_gdof,v));
                                           memory_col_globalProcessToGlobalCluster[procExtrapoled][j_gdof]=j_gdofGlobalCluster;//domaindof->mapGlobalProcessToGlobalCluster()[j_gdof];
                                       }
                               } // if ( it_proc->template get<0>()==procExtrapoled)
                       } // for ( ; it_proc!=en_proc; ++it_proc)

               } // for (size_type cpt_gdof=0 ; it_extrap!=en_extrap ; ++it_extrap,++cpt_gdof )
       } // if ( doExtrapolationAtStart )


    //-----------------------------------------------------------------------------------------
    //this->worldCommFusion().barrier();
    //std::cout << "Op---1----- " << std::endl;
    //-----------------------------------------------------------------------------------------
    // compute graph
    sparsity_graph->close();//sparsity_graph->printPython("mygraphpythonMPI.py");
    //-----------------------------------------------------------------------------------------
    //std::cout << "Op---2----- " << std::endl;
    //this->worldCommFusion().barrier();
    //-----------------------------------------------------------------------------------------
    size_type mapCol_nLocalDof = 0;
    for (int p=0;p<nProc_col;++p)
        {
            mapCol_nLocalDof += memory_col_globalProcessToGlobalCluster[p].size();
        }
    //std::cout << "mapCol_nLocalDof " << mapCol_nLocalDof << std::endl;
    std::vector<size_type> mapCol_globalProcessToGlobalCluster(mapCol_nLocalDof);
    //std::vector<size_type> new_mapGlobalClusterToGlobalProcess(new_nLocalDofWithoutGhost);
    std::vector<std::map<size_type,size_type> > mapCol_LocalSpaceDofToLocalInterpDof(nProc_col);
    size_type currentLocalDof=0;
    for (int p = 0 ; p<nProc_col;++p)
        {
            auto it_map = memory_col_globalProcessToGlobalCluster[p].begin();
            auto en_map = memory_col_globalProcessToGlobalCluster[p].end();
            for ( ; it_map!=en_map ; ++it_map)
                {
                    mapCol_globalProcessToGlobalCluster[currentLocalDof]=it_map->second;//memory_col_globalProcessToGlobalCluster[proc_id][i];
                    mapCol_LocalSpaceDofToLocalInterpDof[p][it_map->first]=currentLocalDof;

                    //if ( new_firstdofcol<=it_map->second && new_lastdofcol>=it_map->second)
                    //    new_mapGlobalClusterToGlobalProcess[it_map->second-new_firstdofcol]=currentLocalDof;

                    ++currentLocalDof;
                }
        }

    //-----------------------------------------
    //std::cout << "Op---3----- " << this->worldCommFusion().godRank() << std::endl;
    //this->worldCommFusion().barrier();
    //-----------------------------------------
    // build data map for the columns
    //this->domainSpace()->mapOnOff().showMeMapGlobalProcessToGlobalCluster();
    //this->dualImageSpace()->worldCommPtr()->showMe();
    datamap_ptrtype mapColInterp( new datamap_type(this->dualImageSpace()->worldCommPtr()));// this->domainSpace()->mapOnOff().worldCommPtr());
    mapColInterp->setNDof(this->domainSpace()->mapOnOff().nDof());

    mapColInterp->setNLocalDofWithoutGhost( proc_id, new_nLocalDofWithoutGhost );//  this->domainSpace()->mapOnOff().nLocalDofWithoutGhost() );
    mapColInterp->setNLocalDofWithGhost( proc_id, mapCol_nLocalDof/*this->domainSpace()->mapOnOff().nLocalDofWithGhost()*/ );
    mapColInterp->setFirstDof( proc_id, this->domainSpace()->mapOnOff().firstDof() );
    mapColInterp->setLastDof( proc_id,  this->domainSpace()->mapOnOff().lastDof() );
    mapColInterp->setFirstDofGlobalCluster( proc_id, new_firstdofcol );
    mapColInterp->setLastDofGlobalCluster( proc_id, new_lastdofcol );
    mapColInterp->setMapGlobalProcessToGlobalCluster(mapCol_globalProcessToGlobalCluster);
    //mapColInterp->setMapGlobalClusterToGlobalProcess(new_mapGlobalClusterToGlobalProcess);
    //if ( this->dualImageSpace()->worldCommPtr()->isActive() ) mapColInterp.showMeMapGlobalProcessToGlobalCluster();

    //-----------------------------------------
    //std::cout << "Op---4----- " << this->worldCommFusion().godRank() << " isA " << this->dualImageSpace()->worldCommPtr()->isActive() << std::endl;
    //this->worldCommFusion().barrier();
    //-----------------------------------------
    // create matrix for active process
    if ( this->dualImageSpace()->worldCommPtr()->isActive() )
        {
            this->matPtr() = this->backend()->newMatrix( mapColInterp,//this->domainSpace()->mapOnOff(),
                                                         this->dualImageSpace()->dofOn(),
                                                         sparsity_graph  );
        }
    //-----------------------------------------
    //std::cout << "Op---5----- " << this->worldCommFusion().godRank() << std::endl;
    //this->worldCommFusion().barrier();
    //-----------------------------------------
    // create null matrix for inactive process
    if ( !this->dualImageSpace()->worldCommPtr()->isActive() && buildNonZeroMatrix )
        {
            this->matPtr() = this->backend()->newZeroMatrix( mapColInterp,//this->domainSpace()->mapOnOff(),
                                                             this->dualImageSpace()->dofOn() );
        }
    //-----------------------------------------
    //std::cout << "Op---6----- "  << this->worldCommFusion().godRank() << std::endl;
    //this->worldCommFusion().barrier();
    //-----------------------------------------
    // assemble matrix
    if ( this->dualImageSpace()->worldCommPtr()->isActive() )
        {
            for (size_type idx_i=0 ; idx_i<nrow_dof_on_proc;++idx_i)
                {
                    for (auto it_j=memory_valueInMatrix[idx_i].begin(),en_j=memory_valueInMatrix[idx_i].end() ; it_j!=en_j ; ++it_j)
                        {
                            if (memory_valueInMatrix[idx_i].size()>0)
                            {
                                this->matPtr()->set(idx_i,mapCol_LocalSpaceDofToLocalInterpDof[it_j->template get<0>()][it_j->template get<1>()],it_j->template get<2>());
                            }
                        }
                }
        }
    //-----------------------------------------
    //std::cout << "Op---7----- " << std::endl;
    //this->worldCommFusion().barrier();
    //-----------------------------------------


}

//-----------------------------------------------------------------------------------------------------------------//

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
std::list<boost::tuple<typename OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::size_type,uint16_type> >
OperatorInterpolation<DomainSpaceType,
                      ImageSpaceType,
                      IteratorRange,
                      InterpType>::updateNoRelationMeshMPI_upWithMyWorld(const std::vector< std::vector<size_type> > & memmapGdof,
                                                                         const std::vector< std::vector<uint16_type> > & memmapComp,
                                                                         const std::vector<std::vector<typename image_mesh_type::node_type> > & pointsSearched,
                                                                         const std::vector<std::vector< std::vector<typename image_mesh_type::node_type > > > & memmap_vertices,
                                                                         graph_ptrtype & sparsity_graph,
                                                                         std::vector< std::list<boost::tuple<int,size_type,double> > > & memory_valueInMatrix,
                                                                         std::vector<std::map<size_type,size_type> > & memory_col_globalProcessToGlobalCluster,
                                                                         std::vector<std::set<size_type> > & dof_searchWithProc,
                                                                         bool extrapolation_mode,
                                                                         extrapolation_memory_type & dof_extrapolationData )
{
    const rank_type proc_id = this->domainSpace()->mesh()->worldCommPtr()->localRank();

    std::list<boost::tuple<size_type,uint16_type> > memory_localisationFail;// gdof,comp
    if ( memmapGdof[proc_id].empty() )
        return memory_localisationFail;

    auto const* imagedof = this->dualImageSpace()->dof().get();
    auto const* domaindof = this->domainSpace()->dof().get();
    auto const* domainbasis = this->domainSpace()->basis().get();

    //auto locTool = this->domainSpace()->mesh()->tool_localization();
    auto locTool = support(this->domainSpace())->tool_localization();
    bool notUseOptLocTest = domain_mesh_type::nDim!=domain_mesh_type::nRealDim;
    //if (notUseOptLocTest) locTool->kdtree()->nbNearNeighbor(domain_mesh_type::element_type::numPoints);

    matrix_node_type ptsReal( image_mesh_type::nRealDim, 1 );
    matrix_node_type ptsRef( domain_mesh_type::nDim , 1 );
    matrix_node_type MlocEval(domain_basis_type::nLocalDof*domain_basis_type::nComponents1,1);
    matrix_node_type verticesOfEltSearched;

    size_type eltIdLocalised = invalid_v<size_type>;

    for ( size_type k=0 ; k<memmapGdof[proc_id].size() ; ++k)
        {
            //----------------------------------------------------------------
            if (this->interpolationType().componentsAreSamePoint())
                //------------------------------------------------------------
                {
                    // the searched point
                    ublas::column(ptsReal,0 ) = pointsSearched[proc_id][k];
                    // vertice with conforme case
                    if (InterpType::isConforming())
                        {
                            const uint16_type sizeVertices = memmap_vertices[proc_id][k].size();
                            verticesOfEltSearched.resize( image_mesh_type::nRealDim,sizeVertices);
                            for ( uint16_type v=0;v<sizeVertices;++v)
                                ublas::column(verticesOfEltSearched,v)=memmap_vertices[proc_id][k][v];
                        }

                    // localisation process
                    if (notUseOptLocTest) eltIdLocalised=invalid_v<size_type>;
                    auto resLocalisation = locTool->run_analysis(ptsReal,eltIdLocalised,verticesOfEltSearched,
                                                                 mpl::int_<interpolation_type::isConforming()>());
                    if (!resLocalisation.template get<0>()[0]) // not find
                        {
                            for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                                {
                                    const auto gdof = memmapGdof[proc_id][k+comp];
                                    memory_localisationFail.push_back(boost::make_tuple(gdof,comp) );
                                    dof_searchWithProc[gdof].insert(proc_id);
                                }
                        }
                    else // point found
                        {
                            eltIdLocalised = resLocalisation.template get<1>();

                            if (extrapolation_mode)
                                {
                                    auto const& eltExtrapoled = this->domainSpace()->mesh()->element(eltIdLocalised);
                                    auto const& verticesExtrapoled = eltExtrapoled.vertices();
                                    typename image_mesh_type::node_type bary(verticesExtrapoled.size1());
                                    for (int qi=0;qi<verticesExtrapoled.size1();++qi) bary(qi)=0;//important
                                    for (int qj=0;qj<verticesExtrapoled.size2();++qj)
                                        {
                                            /**/                              bary(0) += ublas::column(verticesExtrapoled,qj)(0);
                                            if (verticesExtrapoled.size1()>1) bary(1) += ublas::column(verticesExtrapoled,qj)(1);
                                            if (verticesExtrapoled.size1()>2) bary(2) += ublas::column(verticesExtrapoled,qj)(2);
                                        }
                                    /**/                              bary(0) /= verticesExtrapoled.size2();
                                    if (verticesExtrapoled.size1()>1) bary(1) /= verticesExtrapoled.size2();
                                    if (verticesExtrapoled.size1()>2) bary(2) /= verticesExtrapoled.size2();
                                    typename image_mesh_type::node_type theRefPtExtrap = locTool->result_analysis().begin()->second.begin()->template get<1>();
                                    std::vector<std::pair<size_type,size_type> > j_gdofs(domain_basis_type::nLocalDof);
                                    for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                                        {
                                            const auto gdof = memmapGdof[proc_id][k+comp];
                                            for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                {
                                                    const size_type j_gdof =  domaindof->localToGlobal( eltIdLocalised,jloc,comp ).index();
                                                    j_gdofs[jloc]=std::make_pair(j_gdof,domaindof->mapGlobalProcessToGlobalCluster()[j_gdof]);
                                                }
                                            dof_extrapolationData[gdof].push_back(boost::make_tuple(proc_id,bary,theRefPtExtrap,j_gdofs,comp));
                                        }
                                }
                            else // no extrapolation
                                {
                                    for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                                        {
                                            const auto gdof = memmapGdof[proc_id][k+comp];
                                            // get the graph row
                                            auto const ig1 = imagedof->mapGlobalProcessToGlobalCluster()[gdof];
                                            auto const theproc = imagedof->procOnGlobalCluster(ig1);
                                            auto& row = sparsity_graph->row(ig1);
                                            row.template get<0>() = theproc;
                                            row.template get<1>() = ig1 - imagedof->firstDofGlobalCluster(theproc);
                                            // for each localised points
                                            auto itanal = locTool->result_analysis().begin();//  result_analysis_begin();
                                            auto const itanal_end = locTool->result_analysis().end();//result_analysis_end();
                                            for ( ;itanal!=itanal_end;++itanal)
                                                {
                                                    const auto itL=itanal->second.begin();
                                                    ublas::column( ptsRef, 0 ) = boost::get<1>(*itL);
                                                    // evaluate basis functions for this point
                                                    MlocEval = domainbasis->evaluate( ptsRef );
                                                    for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                        {
                                                            //get global dof
                                                            const size_type j_gdof =  domaindof->localToGlobal( itanal->first,jloc,comp ).index();
                                                            // up graph
                                                            row.template get<2>().insert(domaindof->mapGlobalProcessToGlobalCluster()[j_gdof]);
                                                            // get value
                                                            const value_type v = MlocEval( domain_basis_type::nComponents1*jloc +
                                                                                           comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof +
                                                                                           comp, 0 );
                                                            // save value
                                                            memory_valueInMatrix[gdof].push_back(boost::make_tuple(proc_id,j_gdof,v));
                                                            memory_col_globalProcessToGlobalCluster[proc_id][j_gdof]=domaindof->mapGlobalProcessToGlobalCluster()[j_gdof];
                                                        }
                                                }
                                            // dof ok : not anymore localise
                                            dof_searchWithProc[gdof].insert(proc_id);
                                        } // comp
                                } // no extrapolation
                        } // else // point found
                    // change k in for
                    k=k+image_basis_type::nComponents-1;
                } // optimization : this->interpolationType().componentsAreSamePoint()
            //----------------------------------------------------
            else // component optimization
                //------------------------------------------------
                {
                    const auto gdof = memmapGdof[proc_id][k];
                    const auto comp = memmapComp[proc_id][k];
                    ublas::column(ptsReal,0 ) = imagedof->dofPoint(gdof).template get<0>();
                    // vertice with conforme case
                    if (InterpType::isConforming())
                        {
                            const uint16_type sizeVertices = memmap_vertices[proc_id][k].size();
                            verticesOfEltSearched.resize( image_mesh_type::nRealDim,sizeVertices);
                            for ( uint16_type v=0;v<sizeVertices;++v)
                                ublas::column(verticesOfEltSearched,v)=memmap_vertices[proc_id][k][v];
                        }

                    // localisation process
                    if (notUseOptLocTest) eltIdLocalised=invalid_v<size_type>;
                    auto resLocalisation = locTool->run_analysis(ptsReal, eltIdLocalised, verticesOfEltSearched, mpl::int_<interpolation_type::isConforming()>());
                    if (!resLocalisation.template get<0>()[0]) // not find
                        {
                            memory_localisationFail.push_back(boost::make_tuple(gdof,comp) );
                            dof_searchWithProc[gdof].insert(proc_id);
                        }
                    else // point found
                        {
                            eltIdLocalised = resLocalisation.template get<1>();
                            if (extrapolation_mode)
                                {
                                    auto const& eltExtrapoled = this->domainSpace()->mesh()->element(eltIdLocalised);
                                    auto const& verticesExtrapoled = eltExtrapoled.vertices();
                                    typename image_mesh_type::node_type bary(verticesExtrapoled.size1());
                                    for (int qi=0;qi<verticesExtrapoled.size1();++qi) bary(qi)=0;//important
                                    for (int qj=0;qj<verticesExtrapoled.size2();++qj)
                                        {
                                            /**/                              bary(0) += ublas::column(verticesExtrapoled,qj)(0);
                                            if (verticesExtrapoled.size1()>1) bary(1) += ublas::column(verticesExtrapoled,qj)(1);
                                            if (verticesExtrapoled.size1()>2) bary(2) += ublas::column(verticesExtrapoled,qj)(2);
                                        }
                                    /**/                              bary(0) /= verticesExtrapoled.size2();
                                    if (verticesExtrapoled.size1()>1) bary(1) /= verticesExtrapoled.size2();
                                    if (verticesExtrapoled.size1()>2) bary(2) /= verticesExtrapoled.size2();
                                    typename image_mesh_type::node_type theRefPtExtrap = locTool->result_analysis().begin()->second.begin()->template get<1>();
                                    std::vector<std::pair<size_type,size_type> > j_gdofs(domain_basis_type::nLocalDof);
                                    for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                        {
                                            const size_type j_gdof =  domaindof->localToGlobal( eltIdLocalised,jloc,comp ).index();
                                            j_gdofs[jloc]=std::make_pair(j_gdof,domaindof->mapGlobalProcessToGlobalCluster()[j_gdof]);
                                        }
                                    dof_extrapolationData[gdof].push_back(boost::make_tuple(proc_id,bary,theRefPtExtrap,j_gdofs,comp));
                                }
                            else
                                {
                                    // get the graph row
                                    auto const ig1 = imagedof->mapGlobalProcessToGlobalCluster()[gdof];
                                    auto const theproc = imagedof->procOnGlobalCluster(ig1);
                                    auto& row = sparsity_graph->row(ig1);
                                    row.template get<0>() = theproc;
                                    row.template get<1>() = ig1 - imagedof->firstDofGlobalCluster(theproc);
                                    // for each localised points
                                    auto itanal = locTool->result_analysis().begin();//  result_analysis_begin();
                                    auto const itanal_end = locTool->result_analysis().end();//result_analysis_end();
                                    for ( ;itanal!=itanal_end;++itanal)
                                        {
                                            const auto itL=itanal->second.begin();
                                            ublas::column( ptsRef, 0 ) = boost::get<1>(*itL);
                                            // evaluate basis functions for this point
                                            MlocEval = domainbasis->evaluate( ptsRef );
                                            for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                {
                                                    //get global dof
                                                    const size_type j_gdof =  domaindof->localToGlobal( itanal->first,jloc,comp ).index();
                                                    // up graph
                                                    row.template get<2>().insert(domaindof->mapGlobalProcessToGlobalCluster()[j_gdof]);
                                                    // get value
                                                    const value_type v = MlocEval( domain_basis_type::nComponents1*jloc +
                                                                                   comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof +
                                                                                   comp, 0 );
                                                    // save value
                                                    memory_valueInMatrix[gdof].push_back(boost::make_tuple(proc_id,j_gdof,v));
                                                    memory_col_globalProcessToGlobalCluster[proc_id][j_gdof]=domaindof->mapGlobalProcessToGlobalCluster()[j_gdof];
                                                }
                                        }
                                    // dof ok : not anymore localise
                                    dof_searchWithProc[gdof].insert(proc_id);
                                }
                        } // else point found
                }// no optimization
        } // for ( size_type k=0 ; k<memmapGdof[proc_id].size() ; ++k)

    return memory_localisationFail;
}

//-----------------------------------------------------------------------------------------------------------------//

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
std::list<boost::tuple<typename OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::size_type,uint16_type> >
OperatorInterpolation<DomainSpaceType, ImageSpaceType,
                      IteratorRange,InterpType>::updateNoRelationMeshMPI_upWithOtherWorld( boost::tuple<std::vector<rank_type>,std::vector<rank_type>,std::vector<boost::tuple<int,int> > > const& worldcommFusionProperties,
                                                                                           std::vector< std::vector<size_type> > const& memmapGdof,
                                                                                           std::vector< std::vector<uint16_type> > const& memmapComp,
                                                                                           std::vector<std::vector<typename image_mesh_type::node_type> > const& pointsSearched,
                                                                                           std::vector<std::vector< std::vector<typename image_mesh_type::node_type > > > const & memmap_vertices,
                                                                                           graph_ptrtype & sparsity_graph,
                                                                                           std::vector< std::list<boost::tuple<int,size_type,double> > > & memory_valueInMatrix,
                                                                                           std::vector<std::map<size_type,size_type> > & memory_col_globalProcessToGlobalCluster,
                                                                                           std::vector<std::set<size_type> > & dof_searchWithProc,
                                                                                           bool extrapolation_mode,
                                                                                           extrapolation_memory_type & dof_extrapolationData )
{
    std::list<boost::tuple<size_type,uint16_type> > memory_localisationFail;// gdof,comp

    auto const* imagedof = this->dualImageSpace()->dof().get();
    auto const* domaindof = this->domainSpace()->dof().get();
    auto const* domainbasis = this->domainSpace()->basis().get();

    const size_type proc_id = this->dualImageSpace()->worldCommPtr()/*worldsComm()[0]*/.localRank();
    const size_type proc_id_image = this->dualImageSpace()->mesh()->worldCommPtr()->localRank();
    const size_type proc_id_domain = this->domainSpace()->mesh()->worldCommPtr()->localRank();
    const size_type nProc = this->dualImageSpace()->mesh()->worldCommPtr()->size();
    const size_type nProc_row = this->dualImageSpace()->mesh()->worldCommPtr()->localSize();
    const size_type nProc_col = this->domainSpace()->mesh()->worldCommPtr()->localSize();
    const size_type nProc_image = this->dualImageSpace()->mesh()->worldCommPtr()->localSize();
    const size_type nProc_domain = this->domainSpace()->mesh()->worldCommPtr()->localSize();

    // localisation tool with matrix node
    //auto locTool = this->domainSpace()->mesh()->tool_localization();
    auto locTool = support(this->domainSpace())->tool_localization();
    bool notUseOptLocTest = domain_mesh_type::nDim!=domain_mesh_type::nRealDim;
    //if (notUseOptLocTest) locTool->kdtree()->nbNearNeighbor(domain_mesh_type::element_type::numPoints);

    matrix_node_type ptsReal( image_mesh_type::nRealDim, 1 );
    matrix_node_type ptsRef( domain_mesh_type::nDim , 1 );
    matrix_node_type MlocEval(domain_basis_type::nLocalDof*domain_basis_type::nComponents1,1);
    matrix_node_type verticesOfEltSearched;

    // random (just to start)
    //size_type eltIdLocalised = this->domainSpace()->mesh()->beginElementWithId(this->domainSpace()->mesh()->worldCommPtr()->localRank())->id();
    // size_type eltIdLocalised = this->domainSpace()->mesh()->beginElementWithProcessId(this->domainSpace()->mesh()->worldCommPtr()->localRank())->id();
    auto const& eltRandom = this->domainSpace()->mesh()->firstElementIteratorWithProcessId()->second;
    size_type eltIdLocalised = eltRandom.id();

    std::vector<bool> dof_done( this->dualImageSpace()->nLocalDof(), false);

    // usefull container
    std::vector<size_type> pointsSearchedSizeWorld(this->dualImageSpace()->mesh()->worldCommPtr()->localComm().size());
    std::vector<typename image_mesh_type::node_type> dataToRecv(1);
    std::vector<uint16_type> dataToRecv_Comp(1,0);
    std::vector< std::vector< typename image_mesh_type::node_type > > dataToRecv_Vertices(1);
    std::vector<typename image_mesh_type::node_type> pointsRefFinded(1);
    std::vector<typename image_mesh_type::node_type> pointsBaryFinded(1);// extrapolation only
    std::vector<bool> pointsRefIsFinded(1,false);
    std::vector<int> pointsIdEltFinded(1,0);
    std::vector<std::vector<int> > pointsDofsColFinded(1,std::vector<int>(1,0));
    std::vector<std::vector<int> > pointsDofsGlobalClusterColFinded(1,std::vector<int>(1,0));
    std::vector<uint16_type> pointsComp(1,0);

#if 0
    // Attention : marche que si les 2 worldcomms qui s'emboite (mon cas)
    std::vector<int> localMeshRankToWorldCommFusion_domain(nProc_col);
    mpi::all_gather( this->domainSpace()->mesh()->worldCommPtr()->localComm(),
                     this->worldCommFusion().globalRank(),
                     localMeshRankToWorldCommFusion_domain );
    std::vector<int> localMeshRankToWorldCommFusion_image(nProc_row);
    mpi::all_gather( this->dualImageSpace()->mesh()->worldCommPtr()->localComm(),
                     this->worldCommFusion().globalRank(),
                     localMeshRankToWorldCommFusion_image );

    std::vector<int> domainProcIsActive_fusion(this->worldCommFusion().globalSize());
    mpi::all_gather( this->worldCommFusion().globalComm(),
                     (int)this->domainSpace()->worldCommPtr()->isActive(),
                     domainProcIsActive_fusion );
    std::vector<int> imageProcIsActive_fusion(this->worldCommFusion().globalSize());
    mpi::all_gather( this->worldCommFusion().globalComm(),
                     (int)this->dualImageSpace()->worldCommPtr()->isActive(),
                     imageProcIsActive_fusion );

    int firstActiveProc_image=0;
    bool findFirstActive_image=false;
    while (!findFirstActive_image)
        {
            if (imageProcIsActive_fusion[firstActiveProc_image])
                {
                    findFirstActive_image=true;
                }
            else ++firstActiveProc_image;
        }
    int firstActiveProc_domain=0;
    bool findFirstActive_domain=false;
    while (!findFirstActive_domain)
        {
            if (domainProcIsActive_fusion[firstActiveProc_domain])
                {
                    findFirstActive_domain=true;
                }
            else ++firstActiveProc_domain;
        }

    for (int p=0;p<localMeshRankToWorldCommFusion_image.size(); ++p)
        {
            if (!this->dualImageSpace()->worldCommPtr()->isActive()) localMeshRankToWorldCommFusion_image[p]=p%nProc_image+firstActiveProc_image; // FAIRE COMMMUNICATION!!!!!
        }
    for (int p=0;p<localMeshRankToWorldCommFusion_domain.size(); ++p)
        {
            if (!this->domainSpace()->worldCommPtr()->isActive()) localMeshRankToWorldCommFusion_domain[p]=p%nProc_domain+firstActiveProc_domain; // FAIRE COMMMUNICATION!!!!!
        }
#else

    auto const& localMeshRankToWorldCommFusion_domain = worldcommFusionProperties.template get<0>();
    auto const& localMeshRankToWorldCommFusion_image = worldcommFusionProperties.template get<1>();
    auto const& activitiesOnWorldCommFusion = worldcommFusionProperties.template get<2>();

    std::vector<int> domainProcIsActive_fusion(this->worldCommFusion().globalSize());
    std::vector<int> imageProcIsActive_fusion(this->worldCommFusion().globalSize());
    for (int p=0 ; p<this->worldCommFusion().globalSize() ; ++p)
    {
        domainProcIsActive_fusion[p] = activitiesOnWorldCommFusion[p].template get<0>();
        imageProcIsActive_fusion[p] = activitiesOnWorldCommFusion[p].template get<1>();
    }

#endif



#if 1
    std::vector<std::vector<int> > searchDistribution(nProc);
    for (int p=0;p<nProc_image;++p)
        {
            searchDistribution[p].clear();
            if ( proc_id_image == p && imageProcIsActive_fusion[this->worldCommFusion().globalRank()] )
                {
                    for (int q=0;q<nProc_domain;++q)
                        {
                            if( pointsSearched[q].size()>0 && localMeshRankToWorldCommFusion_image[p]!=localMeshRankToWorldCommFusion_domain[q] )
                                {
                                    searchDistribution[p].push_back(q);
                                }
                        }
                }

            mpi::broadcast( this->worldCommFusion().globalComm(), searchDistribution[p], localMeshRankToWorldCommFusion_image[p] );
            //mpi::broadcast( this->worldCommFusion().globalComm(), searchDistribution2, localMeshRankToWorldCommFusion_image[proc_id_image] );
        }

#else //OLD
    // searchDistribution (no comm with ourself)
    std::vector<std::list<int> > searchDistribution(nProc);
    for (int p=0;p<nProc_image;++p)
        {
            searchDistribution[p].clear();
            for (int q=0;q<nProc_domain;++q)
                {
                    //if (q!=p)
                    if( (localMeshRankToWorldCommFusion_image[p])!=localMeshRankToWorldCommFusion_domain[q] )
                        {
                            searchDistribution[p].push_back(q);
                        }
                }
        }
#endif


#if 0
    this->worldCommFusion().barrier();
    for (int p=0;p<this->worldCommFusion().globalSize();++p)
        {
            if (p==this->worldCommFusion().globalRank())
                {
                    std::cout << "I am proc " << p << "   ";// << "\n";
                    for (int q=0;q<nProc_image;++q)
                        {
                            std::cout << "["<< q << "] : ";
                            auto it_list = searchDistribution[q].begin();
                            auto en_list = searchDistribution[q].end();
                            for ( ; it_list!=en_list;++it_list)
                                {
                                    std::cout << *it_list <<" ";
                                }
                            std::cout << std::endl;
                        }
                }
            this->worldCommFusion().barrier();
        }
#endif





    // Tag for mpi msg
    const int tag_X = 0, tag_Y = 1, tag_Z = 2, tag_IsFind = 3, tag_IdElt = 4, tag_DofsCol = 5, tag_DofsColGC = 6, tag_Comp = 7, tag_Points=8, tag_Bary=9, tag_Vertices=10;
    //------------------------------
    // proc after proc
    for (int proc=0;proc<nProc_row;++proc)
        {
            if ( proc_id_image == proc && imageProcIsActive_fusion[this->worldCommFusion().globalRank()] )  // send info to rankLocalization
                {
                    for (auto it_rankLocalization=searchDistribution[proc].begin(),en_rankLocalization=searchDistribution[proc].end();
                         it_rankLocalization!=en_rankLocalization;++it_rankLocalization)
                        {
                            const int rankLocalization = *it_rankLocalization;
                            const int rankToSend = localMeshRankToWorldCommFusion_domain[rankLocalization];
                            this->worldCommFusion().globalComm().send(rankToSend,tag_Points,pointsSearched[rankLocalization]);
                            this->worldCommFusion().globalComm().send(rankToSend,tag_Comp,memmapComp[rankLocalization]);
                            if (InterpType::isConforming())
                                this->worldCommFusion().globalComm().send(rankToSend,tag_Vertices,memmap_vertices[rankLocalization]);
                        }
                }
            else
                {
                    for (auto it_rankLocalization=searchDistribution[proc].begin(),en_rankLocalization=searchDistribution[proc].end();
                         it_rankLocalization!=en_rankLocalization;++it_rankLocalization)
                        {
                            const int rankLocalization = *it_rankLocalization;
                            if ( proc_id_domain == rankLocalization && domainProcIsActive_fusion[this->worldCommFusion().globalRank()] ) // get request of proc
                                {
                                    const int rankToRecv = localMeshRankToWorldCommFusion_image[proc];
                                    this->worldCommFusion().globalComm().recv(rankToRecv,tag_Points,dataToRecv);
                                    this->worldCommFusion().globalComm().recv(rankToRecv,tag_Comp,dataToRecv_Comp);
                                    if (InterpType::isConforming())
                                        this->worldCommFusion().globalComm().recv(rankToRecv,tag_Vertices,dataToRecv_Vertices);

                                    const size_type nDataRecv = dataToRecv.size();
                                    // init container
                                    pointsRefFinded.resize(nDataRecv);
                                    pointsRefIsFinded.resize(nDataRecv);std::fill(pointsRefIsFinded.begin(),pointsRefIsFinded.end(),false);
                                    pointsIdEltFinded.resize(nDataRecv);
                                    if (extrapolation_mode) pointsBaryFinded.resize(nDataRecv);
                                    pointsDofsColFinded.resize(nDataRecv,std::vector<int>(1,0));
                                    pointsDofsGlobalClusterColFinded.resize(nDataRecv,std::vector<int>(1,0));
                                    // iterate on points
                                    for (size_type k=0;k<nDataRecv;++k)
                                        {
                                            // get real point to search
                                            ublas::column(ptsReal,0) = dataToRecv[k];
                                            // vertice with conforme case
                                            if (InterpType::isConforming())
                                                {
                                                    const uint16_type sizeVertices = dataToRecv_Vertices[k].size();
                                                    verticesOfEltSearched.resize( image_mesh_type::nRealDim,sizeVertices);
                                                    for ( uint16_type v=0;v<sizeVertices;++v)
                                                        ublas::column(verticesOfEltSearched,v)=dataToRecv_Vertices[k][v];
                                                }
                                            else // random
                                                verticesOfEltSearched = eltRandom.vertices();
                                            // search process
                                            if (notUseOptLocTest) eltIdLocalised=invalid_v<size_type>;
                                            auto resLocalisation = locTool->run_analysis(ptsReal,eltIdLocalised,verticesOfEltSearched,mpl::int_<interpolation_type::value>());
                                            if (resLocalisation.template get<0>()[0]) // is find
                                                {
                                                    eltIdLocalised = resLocalisation.template get<1>();
                                                    //-------------------------------------------------------------------------
                                                    if (!this->interpolationType().componentsAreSamePoint() )// all component
                                                    //-------------------------------------------------------------------------
                                                        {
                                                            const uint16_type comp=dataToRecv_Comp[k];
                                                            pointsRefIsFinded[k]=true;
                                                            pointsIdEltFinded[k]=eltIdLocalised;
                                                            // get point in reference element
                                                            auto itanal = locTool->result_analysis_begin();
                                                            auto const itanal_end = locTool->result_analysis_end();
                                                            for ( ;itanal!=itanal_end;++itanal)
                                                                {
                                                                    auto const itL=itanal->second.begin();
                                                                    //ublas::column( ptsRef, 0 ) = boost::get<1>(*itL);
                                                                    pointsRefFinded[k] = boost::get<1>(*itL);
                                                                }
                                                            // get global process dof and global cluster dof : column data map
                                                            pointsDofsColFinded[k].resize(domain_basis_type::nLocalDof);
                                                            pointsDofsGlobalClusterColFinded[k].resize(domain_basis_type::nLocalDof);
                                                            for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                                {
                                                                    const auto j_gdof = domaindof->localToGlobal( eltIdLocalised,jloc,comp ).index();
                                                                    pointsDofsColFinded[k][jloc] = j_gdof;
                                                                    pointsDofsGlobalClusterColFinded[k][jloc] = domaindof->mapGlobalProcessToGlobalCluster()[j_gdof];
                                                                }
                                                            if (extrapolation_mode)
                                                                {
                                                                    auto const& eltExtrapoled = this->domainSpace()->mesh()->element(eltIdLocalised);
                                                                    auto const verticesExtrapoled = eltExtrapoled.vertices();
                                                                    typename image_mesh_type::node_type bary(verticesExtrapoled.size1());
                                                                    for (int qi=0;qi<verticesExtrapoled.size1();++qi) bary(qi)=0;//important
                                                                    for (int qj=0;qj<verticesExtrapoled.size2();++qj)
                                                                        {
                                                                            /**/                              bary(0) += ublas::column(verticesExtrapoled,qj)(0);
                                                                            if (verticesExtrapoled.size1()>1) bary(1) += ublas::column(verticesExtrapoled,qj)(1);
                                                                            if (verticesExtrapoled.size1()>2) bary(2) += ublas::column(verticesExtrapoled,qj)(2);
                                                                        }
                                                                    /**/                              bary(0) /= verticesExtrapoled.size2();
                                                                    if (verticesExtrapoled.size1()>1) bary(1) /= verticesExtrapoled.size2();
                                                                    if (verticesExtrapoled.size1()>2) bary(2) /= verticesExtrapoled.size2();
                                                                    pointsBaryFinded[k]=bary;
                                                                }
                                                        }
                                                    //-------------------------------------------------------------------------
                                                    else // components optimization
                                                    //-------------------------------------------------------------------------
                                                        {
                                                            auto const ptRefFinded = locTool->result_analysis_begin()->second.begin()->template get<1>();
                                                            for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                                                                {
                                                                    pointsRefIsFinded[k+comp]=true;
                                                                    pointsIdEltFinded[k+comp]=eltIdLocalised;
                                                                    // get point in reference element
                                                                    pointsRefFinded[k+comp] = ptRefFinded;

                                                                    pointsDofsColFinded[k+comp].resize(domain_basis_type::nLocalDof);
                                                                    pointsDofsGlobalClusterColFinded[k+comp].resize(domain_basis_type::nLocalDof);
                                                                    for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                                        {
                                                                            const auto j_gdof = domaindof->localToGlobal( eltIdLocalised,jloc,comp ).index();
                                                                            pointsDofsColFinded[k+comp][jloc] = j_gdof;
                                                                            pointsDofsGlobalClusterColFinded[k+comp][jloc] = domaindof->mapGlobalProcessToGlobalCluster()[j_gdof];
                                                                        }
                                                                }

                                                                    if (extrapolation_mode)
                                                                        {
                                                                            auto const& eltExtrapoled = this->domainSpace()->mesh()->element(eltIdLocalised);
                                                                            auto const verticesExtrapoled = eltExtrapoled.vertices();
                                                                            typename image_mesh_type::node_type bary(verticesExtrapoled.size1());
                                                                            for (int qi=0;qi<verticesExtrapoled.size1();++qi) bary(qi)=0;//important
                                                                            for (int qj=0;qj<verticesExtrapoled.size2();++qj)
                                                                                {
                                                                                    /**/                              bary(0) += ublas::column(verticesExtrapoled,qj)(0);
                                                                                    if (verticesExtrapoled.size1()>1) bary(1) += ublas::column(verticesExtrapoled,qj)(1);
                                                                                    if (verticesExtrapoled.size1()>2) bary(2) += ublas::column(verticesExtrapoled,qj)(2);
                                                                                }
                                                                            /**/                              bary(0) /= verticesExtrapoled.size2();
                                                                            if (verticesExtrapoled.size1()>1) bary(1) /= verticesExtrapoled.size2();
                                                                            if (verticesExtrapoled.size1()>2) bary(2) /= verticesExtrapoled.size2();
                                                                            for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                                                                                pointsBaryFinded[k+comp]=bary;
                                                                        }
                                                                    // change increment in for
                                                                    k=k+image_basis_type::nComponents-1;
                                                        }

                                                    //std::cout << "F";
                                                }
                                            else // Not Find!
                                                {
                                                    //memory_localisationFail
                                                    //std::cout << "NOT FIND"<<std::endl;
                                                }

                                        } // for (size_type k=0;k<nDataRecv;++k)

                                    const int rankToSend = localMeshRankToWorldCommFusion_image[proc];
                                    this->worldCommFusion().globalComm().send(rankToSend,tag_Points,pointsRefFinded);
                                    this->worldCommFusion().globalComm().send(rankToSend,tag_IsFind,pointsRefIsFinded);
                                    this->worldCommFusion().globalComm().send(rankToSend,tag_IdElt,pointsIdEltFinded);
                                    this->worldCommFusion().globalComm().send(rankToSend,tag_DofsCol,pointsDofsColFinded);
                                    this->worldCommFusion().globalComm().send(rankToSend,tag_DofsColGC,pointsDofsGlobalClusterColFinded);
                                    if (extrapolation_mode)
                                        this->worldCommFusion().globalComm().send(rankToSend,tag_Bary,pointsBaryFinded);

                                }
                        } // for (auto it_rankLocalization=...
                }




            if ( proc_id_image == proc && imageProcIsActive_fusion[this->worldCommFusion().globalRank()] )
                {
                    for (auto it_rankLocalization=searchDistribution[proc].begin(),en_rankLocalization=searchDistribution[proc].end();
                         it_rankLocalization!=en_rankLocalization;++it_rankLocalization)
                        {
                            const int rankLocalization = *it_rankLocalization;
                            const int rankToRecv = localMeshRankToWorldCommFusion_domain[rankLocalization];
                            this->worldCommFusion().globalComm().recv(rankToRecv,tag_Points,pointsRefFinded);
                            this->worldCommFusion().globalComm().recv(rankToRecv,tag_IsFind,pointsRefIsFinded);
                            this->worldCommFusion().globalComm().recv(rankToRecv,tag_IdElt,pointsIdEltFinded);
                            this->worldCommFusion().globalComm().recv(rankToRecv,tag_DofsCol,pointsDofsColFinded);
                            this->worldCommFusion().globalComm().recv(rankToRecv,tag_DofsColGC,pointsDofsGlobalClusterColFinded);
                            if (extrapolation_mode)
                                this->worldCommFusion().globalComm().recv(rankToRecv,tag_Bary,pointsBaryFinded);

                            const int rankRecv=rankLocalization;
                            for ( int k=0;k<pointsRefFinded.size();++k)
                                {
                                    if (!pointsRefIsFinded[k])
                                        {
                                            memory_localisationFail.push_back(boost::make_tuple(memmapGdof[rankLocalization][k],memmapComp[rankLocalization][k]));
                                            const auto i_gdof = memmapGdof[rankLocalization][k];
                                            dof_searchWithProc[i_gdof].insert(rankRecv);
                                        }
                                    else
                                        {
                                            const auto i_gdof = memmapGdof[rankLocalization][k];
                                            if (!dof_done[i_gdof])
                                                {

                                                    if (extrapolation_mode)
                                                        {
                                                            auto const bary = pointsBaryFinded[k];
                                                            auto const theRefPtExtrap = pointsRefFinded[k];
                                                            const auto comp = memmapComp[rankLocalization][k];
                                                            std::vector<std::pair<size_type,size_type> > j_gdofs(domain_basis_type::nLocalDof);
                                                            for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                                {
                                                                    const size_type j_gdof =  pointsDofsColFinded[k][jloc];
                                                                    const size_type j_gdof_gc = pointsDofsGlobalClusterColFinded[k][jloc];
                                                                    j_gdofs[jloc]=std::make_pair(j_gdof,j_gdof_gc);
                                                                }
                                                            dof_extrapolationData[i_gdof].push_back(boost::make_tuple(rankLocalization,bary,theRefPtExtrap,j_gdofs,comp));
                                                        }
                                                    else
                                                        {
                                                            //std::cout << "T";
                                                            ublas::column( ptsRef, 0 ) = pointsRefFinded[k];
                                                            //evalute point on the reference element
                                                            MlocEval = domainbasis->evaluate( ptsRef );

                                                            const auto comp = memmapComp[rankLocalization][k];
                                                            const size_type myidElt = pointsIdEltFinded[k];
                                                            const auto ig1 = imagedof->mapGlobalProcessToGlobalCluster()[i_gdof];
                                                            const auto theproc = imagedof->procOnGlobalCluster(ig1);
                                                            auto& row = sparsity_graph->row(ig1);
                                                            row.template get<0>() = theproc;
                                                            row.template get<1>() = ig1 - imagedof->firstDofGlobalCluster(theproc);

                                                            for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                                {
                                                                    //get global process dof
                                                                    const size_type j_gdof =  pointsDofsColFinded[k][jloc];
                                                                    //get global cluster dof
                                                                    const size_type j_gdof_gc = pointsDofsGlobalClusterColFinded[k][jloc];
                                                                    // up graph
                                                                    row.template get<2>().insert(j_gdof_gc);
                                                                    // get value
                                                                    const auto v = MlocEval( domain_basis_type::nComponents1*jloc +
                                                                                             comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof +
                                                                                             comp, 0 );
#if 1
                                                                    // save value
                                                                    memory_valueInMatrix[i_gdof].push_back(boost::make_tuple(rankRecv,j_gdof,v));
                                                                    // usefull to build datamap
                                                                    memory_col_globalProcessToGlobalCluster[rankRecv][j_gdof]=j_gdof_gc;
#endif
                                                                }
                                                            // dof ok : not anymore localise
                                                            dof_searchWithProc[i_gdof].insert(rankRecv);
                                                        }
                                                    dof_done[i_gdof]=true;
                                                } // if (!dof_done[i_gdof])
                                        }
                                }
                        }
                } // if ( proc_id_image == proc && imageProcIsActive_fusion[this->worldCommFusion().globalRank()] )
        }




    return memory_localisationFail;
} // version1

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
std::list<boost::tuple<typename OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::size_type,uint16_type> >
OperatorInterpolation<DomainSpaceType, ImageSpaceType,
                      IteratorRange,InterpType>::updateNoRelationMeshMPI_upWithOtherWorld2( boost::tuple<std::vector<rank_type>,std::vector<rank_type>,std::vector<boost::tuple<int,int> > > const& worldcommFusionProperties,
                                                                                           std::vector< std::vector<size_type> > const& memmapGdof,
                                                                                           std::vector< std::vector<uint16_type> > const& memmapComp,
                                                                                           std::vector<std::vector<typename image_mesh_type::node_type> > const& pointsSearched,
                                                                                           std::vector<std::vector< std::vector<typename image_mesh_type::node_type > > > const & memmap_vertices,
                                                                                           graph_ptrtype & sparsity_graph,
                                                                                           std::vector< std::list<boost::tuple<int,size_type,double> > > & memory_valueInMatrix,
                                                                                           std::vector<std::map<size_type,size_type> > & memory_col_globalProcessToGlobalCluster,
                                                                                           std::vector<std::set<size_type> > & dof_searchWithProc,
                                                                                           bool extrapolation_mode,
                                                                                           extrapolation_memory_type & dof_extrapolationData )
{
    std::list<boost::tuple<size_type,uint16_type> > memory_localisationFail;// gdof,comp

    auto const* imagedof = this->dualImageSpace()->dof().get();
    auto const* domaindof = this->domainSpace()->dof().get();
    auto const* domainbasis = this->domainSpace()->basis().get();

    const size_type proc_id = this->dualImageSpace()->worldCommPtr()/*worldsComm()[0]*/->localRank();
    const size_type proc_id_image = this->dualImageSpace()->mesh()->worldCommPtr()->localRank();
    const size_type proc_id_domain = this->domainSpace()->mesh()->worldCommPtr()->localRank();
    const size_type nProc = this->dualImageSpace()->mesh()->worldCommPtr()->size();
    const size_type nProc_row = this->dualImageSpace()->mesh()->worldCommPtr()->localSize();
    const size_type nProc_col = this->domainSpace()->mesh()->worldCommPtr()->localSize();
    const size_type nProc_image = this->dualImageSpace()->mesh()->worldCommPtr()->localSize();
    const size_type nProc_domain = this->domainSpace()->mesh()->worldCommPtr()->localSize();

    // localisation tool with matrix node
    //auto locTool = this->domainSpace()->mesh()->tool_localization();
    auto locTool = support(this->domainSpace())->tool_localization();
    bool notUseOptLocTest = domain_mesh_type::nDim!=domain_mesh_type::nRealDim;
    //if (notUseOptLocTest) locTool->kdtree()->nbNearNeighbor(domain_mesh_type::element_type::numPoints);

    matrix_node_type ptsReal( image_mesh_type::nRealDim, 1 );
    matrix_node_type ptsRef( domain_mesh_type::nDim , 1 );
    matrix_node_type MlocEval(domain_basis_type::nLocalDof*domain_basis_type::nComponents1,1);
    matrix_node_type verticesOfEltSearched;

    size_type eltIdLocalised = invalid_v<size_type>;

    std::vector<bool> dof_done( this->dualImageSpace()->nLocalDof(), false);




    auto const& localMeshRankToWorldCommFusion_domain = worldcommFusionProperties.template get<0>();
    auto const& localMeshRankToWorldCommFusion_image = worldcommFusionProperties.template get<1>();
    auto const& activitiesOnWorldCommFusion = worldcommFusionProperties.template get<2>();

    std::vector<int> domainProcIsActive_fusion(this->worldCommFusion().globalSize());
    std::vector<int> imageProcIsActive_fusion(this->worldCommFusion().globalSize());
    for (int p=0 ; p<this->worldCommFusion().globalSize() ; ++p)
    {
        domainProcIsActive_fusion[p] = activitiesOnWorldCommFusion[p].template get<0>();
        imageProcIsActive_fusion[p] = activitiesOnWorldCommFusion[p].template get<1>();
    }



    // collective mpi comm to known who comunincates whom
    std::vector<std::vector<rank_type> > searchDistribution(nProc);
    std::vector<rank_type> distributionToSend;
    for (rank_type q=0;q<nProc_domain;++q)
    {
        if ( pointsSearched[q].size()>0 && localMeshRankToWorldCommFusion_image[proc_id_image]!=localMeshRankToWorldCommFusion_domain[q] )
            distributionToSend.push_back(q);
    }
    mpi::all_gather( this->worldCommFusion().globalComm(),
                     distributionToSend,
                     searchDistribution );


#if 0
    this->worldCommFusion().barrier();
    for (int p=0;p<this->worldCommFusion().globalSize();++p)
        {
            if (p==this->worldCommFusion().globalRank())
                {
                    std::cout << "I am proc " << p << "\n";
                    for (int q=0;q<nProc_image;++q)
                        {
                            std::cout << "["<< q << "] : ";
                            auto it_list = searchDistribution[q].begin();
                            auto en_list = searchDistribution[q].end();
                            for ( ; it_list!=en_list;++it_list)
                                {
                                    std::cout << *it_list <<" ";
                                }
                            std::cout << std::endl;
                        }
                }
            this->worldCommFusion().barrier();
        }
#endif


    // containers used in send/recv
    std::map<rank_type, boost::tuple< std::vector<typename image_mesh_type::node_type>,
                                      std::vector<uint16_type>,
                                      std::vector< std::vector<typename image_mesh_type::node_type > >
                                      > > dataToRecv, dataToSend;

    std::map<rank_type, boost::tuple< std::vector<typename image_mesh_type::node_type>,
                                      std::vector<bool>,
                                      std::vector<int>,
                                      std::vector<std::vector<int> >,
                                      std::vector<std::vector<int> >,
                                      std::vector<typename image_mesh_type::node_type> > > dataReRecv, dataReSend;


    // compute number of request isend/irecv
    int nbRequest=searchDistribution[proc_id_image].size();
    for (int proc=0;proc<nProc_row;++proc)
    {
        if ( proc_id_image != proc &&
             std::find(searchDistribution[proc].begin(),searchDistribution[proc].end(),proc_id_image) != searchDistribution[proc].end() )
            ++ nbRequest;
    }

    // build container to send
    for ( rank_type rankLocalization : searchDistribution[proc_id_image] )
    {
        const rank_type rankToSend = localMeshRankToWorldCommFusion_domain[rankLocalization];
        //CHECK( rankToSend == rankLocalization ) << "strange " << rankToSend << " " << rankLocalization << "\n";
        if (InterpType::isConforming()==0)
            CHECK( memmap_vertices[rankLocalization].size() == 0 ) << " memmap_vertices must be empty in conformal case\n";

        dataToSend[rankLocalization] = boost::make_tuple( pointsSearched[rankLocalization], memmapComp[rankLocalization], memmap_vertices[rankLocalization] );
    }

    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;
    // send info
    for ( rank_type rankLocalization : searchDistribution[proc_id_image] )
    {
        const rank_type rankToSend = localMeshRankToWorldCommFusion_domain[rankLocalization];
        reqs[cptRequest] = this->worldCommFusion().globalComm().isend( rankToSend , 0, dataToSend[rankLocalization] );
        ++cptRequest;
    }
    // recv info
    for (rank_type proc=0;proc<nProc_row;++proc)
    {
        if ( proc_id_image != proc &&
             std::find(searchDistribution[proc].begin(),searchDistribution[proc].end(),proc_id_image) != searchDistribution[proc].end() )
        {
            const rank_type rankToRecv = localMeshRankToWorldCommFusion_image[proc];
            reqs[cptRequest] = this->worldCommFusion().globalComm().irecv( rankToRecv , 0, dataToRecv[proc] );
            ++cptRequest;
        }
    }

    // wait all requests
    CHECK( cptRequest == nbRequest ) << cptRequest << " vs " <<  nbRequest << "\n" ;
    mpi::wait_all(reqs, reqs + nbRequest);


    cptRequest=0;
    for ( auto const& dataToTreatBase : dataToRecv )
    {
        //const rank_type rankToRecv = localMeshRankToWorldCommFusion_image[dataToTreatNew->first];
        rank_type theproc = dataToTreatBase.first;
        auto const& dataToTreat = dataToTreatBase.second;
        const size_type nDataRecv = dataToTreat.template get<0>().size();
        auto const& dataToRecv_ptsSearch = dataToTreat.template get<0>();
        auto const& dataToRecv_Comp = dataToTreat.template get<1>();
        auto const& dataToRecv_Vertices = dataToTreat.template get<2>();

        // init container
        dataReSend[theproc];
        auto & pointsRefFinded = dataReSend[theproc].template get<0>();
        auto & pointsRefIsFinded = dataReSend[theproc].template get<1>();
        auto & pointsIdEltFinded = dataReSend[theproc].template get<2>();
        auto & pointsDofsColFinded = dataReSend[theproc].template get<3>();
        auto & pointsDofsGlobalClusterColFinded = dataReSend[theproc].template get<4>();
        auto & pointsBaryFinded = dataReSend[theproc].template get<5>();
        pointsRefFinded.resize(nDataRecv);
        pointsRefIsFinded.resize(nDataRecv);std::fill(pointsRefIsFinded.begin(),pointsRefIsFinded.end(),false);
        pointsIdEltFinded.resize(nDataRecv);
        pointsDofsColFinded.resize(nDataRecv,std::vector<int>(1,0));
        pointsDofsGlobalClusterColFinded.resize(nDataRecv,std::vector<int>(1,0));
        if (extrapolation_mode) pointsBaryFinded.resize(nDataRecv);
        else pointsBaryFinded.clear();

        // iterate on points
        for (size_type k=0;k<nDataRecv;++k)
        {
            // get real point to search
            ublas::column(ptsReal,0) = dataToRecv_ptsSearch[k];
            // vertice with conforme case
            if (InterpType::isConforming())
            {
                const uint16_type sizeVertices = dataToRecv_Vertices[k].size();
                verticesOfEltSearched.resize( image_mesh_type::nRealDim,sizeVertices);
                for ( uint16_type v=0;v<sizeVertices;++v)
                    ublas::column(verticesOfEltSearched,v)=dataToRecv_Vertices[k][v];
            }

            // search process
            if (notUseOptLocTest) eltIdLocalised=invalid_v<size_type>;
            auto resLocalisation = locTool->run_analysis(ptsReal,eltIdLocalised,verticesOfEltSearched,mpl::int_<interpolation_type::isConforming()>());
            if (resLocalisation.template get<0>()[0]) // is find
            {
                eltIdLocalised = resLocalisation.template get<1>();
                //-------------------------------------------------------------------------
                if (!this->interpolationType().componentsAreSamePoint() )// all component
                    //-------------------------------------------------------------------------
                {
                    const uint16_type comp=dataToRecv_Comp[k];
                    pointsRefIsFinded[k]=true;
                    pointsIdEltFinded[k]=eltIdLocalised;
                    // get point in reference element
                    auto itanal = locTool->result_analysis_begin();
                    auto const itanal_end = locTool->result_analysis_end();
                    for ( ;itanal!=itanal_end;++itanal)
                    {
                        auto const itL=itanal->second.begin();
                        pointsRefFinded[k] = boost::get<1>(*itL);
                    }
                    // get global process dof and global cluster dof : column data map
                    pointsDofsColFinded[k].resize(domain_basis_type::nLocalDof);
                    pointsDofsGlobalClusterColFinded[k].resize(domain_basis_type::nLocalDof);
                    for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                    {
                        const auto j_gdof = domaindof->localToGlobal( eltIdLocalised,jloc,comp ).index();
                        pointsDofsColFinded[k][jloc] = j_gdof;
                        pointsDofsGlobalClusterColFinded[k][jloc] = domaindof->mapGlobalProcessToGlobalCluster()[j_gdof];
                    }
                    if (extrapolation_mode)
                    {
                        auto const& eltExtrapoled = this->domainSpace()->mesh()->element(eltIdLocalised);
                        auto const verticesExtrapoled = eltExtrapoled.vertices();
                        typename image_mesh_type::node_type bary(verticesExtrapoled.size1());
                        for (int qi=0;qi<verticesExtrapoled.size1();++qi) bary(qi)=0;//important
                        for (int qj=0;qj<verticesExtrapoled.size2();++qj)
                        {
                            /**/                              bary(0) += ublas::column(verticesExtrapoled,qj)(0);
                            if (verticesExtrapoled.size1()>1) bary(1) += ublas::column(verticesExtrapoled,qj)(1);
                            if (verticesExtrapoled.size1()>2) bary(2) += ublas::column(verticesExtrapoled,qj)(2);
                        }
                        /**/                              bary(0) /= verticesExtrapoled.size2();
                        if (verticesExtrapoled.size1()>1) bary(1) /= verticesExtrapoled.size2();
                        if (verticesExtrapoled.size1()>2) bary(2) /= verticesExtrapoled.size2();
                        pointsBaryFinded[k]=bary;
                    }
                }
                //-------------------------------------------------------------------------
                else // components optimization
                    //-------------------------------------------------------------------------
                {
                    auto const ptRefFinded = locTool->result_analysis_begin()->second.begin()->template get<1>();
                    for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                    {
                        pointsRefIsFinded[k+comp]=true;
                        pointsIdEltFinded[k+comp]=eltIdLocalised;
                        // get point in reference element
                        pointsRefFinded[k+comp] = ptRefFinded;

                        pointsDofsColFinded[k+comp].resize(domain_basis_type::nLocalDof);
                        pointsDofsGlobalClusterColFinded[k+comp].resize(domain_basis_type::nLocalDof);
                        for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                        {
                            const auto j_gdof = domaindof->localToGlobal( eltIdLocalised,jloc,comp ).index();
                            pointsDofsColFinded[k+comp][jloc] = j_gdof;
                            pointsDofsGlobalClusterColFinded[k+comp][jloc] = domaindof->mapGlobalProcessToGlobalCluster()[j_gdof];
                        }
                    }

                    if (extrapolation_mode)
                    {
                        auto const& eltExtrapoled = this->domainSpace()->mesh()->element(eltIdLocalised);
                        auto const verticesExtrapoled = eltExtrapoled.vertices();
                        typename image_mesh_type::node_type bary(verticesExtrapoled.size1());
                        for (int qi=0;qi<verticesExtrapoled.size1();++qi) bary(qi)=0;//important
                        for (int qj=0;qj<verticesExtrapoled.size2();++qj)
                        {
                            /**/                              bary(0) += ublas::column(verticesExtrapoled,qj)(0);
                            if (verticesExtrapoled.size1()>1) bary(1) += ublas::column(verticesExtrapoled,qj)(1);
                            if (verticesExtrapoled.size1()>2) bary(2) += ublas::column(verticesExtrapoled,qj)(2);
                        }
                        /**/                              bary(0) /= verticesExtrapoled.size2();
                        if (verticesExtrapoled.size1()>1) bary(1) /= verticesExtrapoled.size2();
                        if (verticesExtrapoled.size1()>2) bary(2) /= verticesExtrapoled.size2();
                        for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                            pointsBaryFinded[k+comp]=bary;
                    }
                    // change increment in for
                    k=k+image_basis_type::nComponents-1;
                }
                //std::cout << "F";
            }
            else // Not Find!
            {
                //memory_localisationFail
                //std::cout << "NOT FIND"<<std::endl;
            }

        } // for (size_type k=0;k<nDataRecv;++k)

        // send results
        const rank_type rankToReSend = localMeshRankToWorldCommFusion_image[theproc];
        reqs[cptRequest] = this->worldCommFusion().globalComm().isend( rankToReSend , 0, dataReSend[theproc] );
        ++cptRequest;
    }

    // recv results
    for ( rank_type rankLocalization : searchDistribution[proc_id_image] )
    {
        const rank_type rankToReRecv = localMeshRankToWorldCommFusion_domain[rankLocalization];
        reqs[cptRequest] = this->worldCommFusion().globalComm().irecv( rankToReRecv , 0, dataReRecv[rankLocalization] );
        ++cptRequest;
    }

    // wait all requests
    CHECK( cptRequest == nbRequest ) << cptRequest << " vs " <<  nbRequest << "\n" ;
    mpi::wait_all(reqs, reqs + nbRequest);
    delete [] reqs;

    for ( auto const& dataToTreatFinalBase : dataReRecv )
    {
        rank_type theproc = dataToTreatFinalBase.first;
        auto const& dataToTreat = dataToTreatFinalBase.second;
        auto const& dataPointsRefFinded = dataToTreat.template get<0>();
        auto const& dataPointsRefIsFinded = dataToTreat.template get<1>();
        auto const& dataPointsIdEltFinded = dataToTreat.template get<2>();
        auto const& dataPointsDofsColFinded = dataToTreat.template get<3>();
        auto const& dataPointsDofsGlobalClusterColFinded  = dataToTreat.template get<4>();
        auto const& dataPointsBaryFinded = dataToTreat.template get<5>();

        //const int rankRecv=rankLocalization;
        for ( int k=0;k<dataPointsRefFinded.size();++k)
        {
            if (!dataPointsRefIsFinded[k])
            {
                memory_localisationFail.push_back(boost::make_tuple(memmapGdof[theproc][k],memmapComp[theproc][k]));
                const auto i_gdof = memmapGdof[theproc][k];
                dof_searchWithProc[i_gdof].insert(theproc);
            }
            else
            {
                const auto i_gdof = memmapGdof[theproc][k];
                if (!dof_done[i_gdof])
                {

                    if (extrapolation_mode)
                    {
                        auto const bary = dataPointsBaryFinded[k];
                        auto const theRefPtExtrap = dataPointsRefFinded[k];
                        const auto comp = memmapComp[theproc][k];
                        std::vector<std::pair<size_type,size_type> > j_gdofs(domain_basis_type::nLocalDof);
                        for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                        {
                            const size_type j_gdof =  dataPointsDofsColFinded[k][jloc];
                            const size_type j_gdof_gc = dataPointsDofsGlobalClusterColFinded[k][jloc];
                            j_gdofs[jloc]=std::make_pair(j_gdof,j_gdof_gc);
                        }
                        dof_extrapolationData[i_gdof].push_back(boost::make_tuple(theproc,bary,theRefPtExtrap,j_gdofs,comp));
                    }
                    else
                    {
                        ublas::column( ptsRef, 0 ) = dataPointsRefFinded[k];
                        //evalute point on the reference element
                        MlocEval = domainbasis->evaluate( ptsRef );

                        const auto comp = memmapComp[theproc][k];
                        const size_type myidElt = dataPointsIdEltFinded[k];
                        const auto ig1 = imagedof->mapGlobalProcessToGlobalCluster()[i_gdof];
                        const auto theRealProc = imagedof->procOnGlobalCluster(ig1);

                        auto& row = sparsity_graph->row(ig1);
                        row.template get<0>() = theRealProc;
                        row.template get<1>() = ig1 - imagedof->firstDofGlobalCluster(theRealProc);

                        for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                        {
                            //get global process dof
                            const size_type j_gdof =  dataPointsDofsColFinded[k][jloc];
                            //get global cluster dof
                            const size_type j_gdof_gc = dataPointsDofsGlobalClusterColFinded[k][jloc];
                            // up graph
                            row.template get<2>().insert(j_gdof_gc);
                            // get value
                            const auto v = MlocEval( domain_basis_type::nComponents1*jloc +
                                                     comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof +
                                                     comp, 0 );

                            // save value
                            memory_valueInMatrix[i_gdof].push_back(boost::make_tuple(theproc,j_gdof,v));
                            // usefull to build datamap
                            memory_col_globalProcessToGlobalCluster[theproc][j_gdof]=j_gdof_gc;
                        }
                        // dof ok : not anymore localise
                        dof_searchWithProc[i_gdof].insert(theproc);
                    }
                    dof_done[i_gdof]=true;
                } // if (!dof_done[i_gdof])
            }
        }


    } // for ( auto const& dataToTreatFinalBase : dataReRecv )

    return memory_localisationFail;
} // version2


template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
boost::tuple<std::vector< std::vector<typename OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::size_type> >, std::vector< std::vector<uint16_type> >,
             std::vector<std::vector<typename OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::image_mesh_type::node_type> >,
             std::vector<std::vector< std::vector<typename OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::image_mesh_type::node_type > > > >
OperatorInterpolation<DomainSpaceType, ImageSpaceType,
                      IteratorRange,InterpType>::updateNoRelationMeshMPI_pointDistribution(const std::vector< std::list<boost::tuple<int,size_type,double> > > & memory_valueInMatrix,
                                                                                           std::vector<std::set<size_type> > const& dof_searchWithProc,
                                                                                           std::set<size_type> const& ghostDofUsedToInterpolate )
{
    //std::cout << " pointDistribution--1--- " << this->domainSpace()->mesh()->worldCommPtr()->godRank() << std::endl;
    //const size_type proc_id = this->dualImageSpace()->worldsComm()[0].localRank();
    //const size_type nProc = this->dualImageSpace()->mesh()->worldCommPtr()->size();
    //const size_type nProc_image = this->dualImageSpace()->mesh()->worldCommPtr()->localSize();
    const int nProc_domain = this->domainSpace()->mesh()->worldCommPtr()->localSize();

    auto const* imagedof = this->dualImageSpace()->dof().get();
    //iterator_type it, en;

    //auto locTool = this->domainSpace()->mesh()->tool_localization();
    auto locTool = support(this->domainSpace())->tool_localization();

    std::vector<bool> dof_done( this->dualImageSpace()->nLocalDof(), false);
    std::vector< std::list<boost::tuple<size_type,uint16_type> > > memSetGdofAndComp( nProc_domain );
    std::vector< std::list<matrix_node_type> > memSetVertices_conformeInterp( nProc_domain );

#if 0
    // Warning communication!!
    std::vector<typename image_mesh_type::node_type> vecBarycenter(nProc_domain);
    mpi::all_gather( this->domainSpace()->mesh()->worldCommPtr()->localComm(),
                     locTool->barycenter(),
                     vecBarycenter );
#else
    // compute vector of barycenter if necessary and not done : warning communication (mpi::all_gather)!
    if ( this->interpolationType().searchWithCommunication() && !locTool->hasComputedBarycentersWorld() )
        locTool->computeBarycentersWorld();

    // build the vector of barycenter (computed from kdtree point) only if search with comm
    std::vector<boost::tuple<bool,typename image_mesh_type::node_type> > vecBarycenter(nProc_domain);
    if (this->interpolationType().searchWithCommunication())
        vecBarycenter = locTool->barycentersWorld();

    //auto const& vecBarycenter = locTool->barycentersWorld();
#endif
    /*std::cout << " proc " << this->domainSpace()->mesh()->worldCommPtr()->localRank()
              << "  procFuion " << this->worldCommFusion().globalRank()
              << " bary " << locTool->barycenter()
              << std::endl;*/


    double distanceMin=0,distance=0,distanceSquare=0;
    int procForPt=0;

    if ( this->dualImageSpace()->worldCommPtr()->isActive() )
    {
        auto itListRange = M_listRange.begin();
        auto const enListRange = M_listRange.end();
        for ( ; itListRange!=enListRange ; ++itListRange)
        {
            //boost::tie( boost::tuples::ignore, it, en ) = *itListRange;
            //for ( ; it!=en;++it )
            for( auto const& theImageEltWrap : *itListRange )
            {
                //auto const& theImageElt = boost::unwrap_ref(*it);
                auto const& theImageElt = boost::unwrap_ref(theImageEltWrap);
                for ( uint16_type iloc = 0; iloc < nLocalDofInDualImageElt; ++iloc )
                {
                    for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                    {

                        const auto gdof =  imagedof->localToGlobal( theImageElt, iloc, comp ).index();

                        const size_type gdofCluster = imagedof->mapGlobalProcessToGlobalCluster()[gdof];

                        if ( imagedof->dofGlobalProcessIsGhost( gdof ) &&
                             ( ghostDofUsedToInterpolate.find(gdofCluster) == ghostDofUsedToInterpolate.end() ) ) continue;


                        if (!dof_done[gdof] && memory_valueInMatrix[gdof].size()==0)
                        {
                            // the dof point
                            const auto imagePoint = imagedof->dofPoint(gdof).template get<0>();

                            if (this->interpolationType().searchWithCommunication()) // mpi communication
                            {
                                distanceMin=INT_MAX;
                                for ( int proc=0 ; proc<nProc_domain; ++proc)
                                {
                                    // if no point in kdtree, ignore this process
                                    if ( !vecBarycenter[proc].template get<0>() ) continue;

                                    auto const& bary = vecBarycenter[proc].template get<1>();
                                    /**/               distanceSquare  = std::pow(imagePoint(0)-bary(0),2);
                                    if (bary.size()>1) distanceSquare += std::pow(imagePoint(1)-bary(1),2);
                                    if (bary.size()>2) distanceSquare += std::pow(imagePoint(2)-bary(2),2);
                                    distance = std::sqrt( distanceSquare );
                                    if (distance<distanceMin && dof_searchWithProc[gdof].find(proc)==dof_searchWithProc[gdof].end() )
                                    {
                                        procForPt = proc;
                                        distanceMin=distance;
                                    }
                                }
                                memSetGdofAndComp[procForPt].push_back(boost::make_tuple(gdof,comp));
                                if (InterpType::isConforming())  // conformal case
                                    memSetVertices_conformeInterp[procForPt].push_back(theImageElt.vertices());
                            }
                            else // only with myself
                            {
                                memSetGdofAndComp[this->domainSpace()->worldCommPtr()->globalRank()].push_back(boost::make_tuple(gdof,comp));
                                if (InterpType::isConforming()) // conformal case
                                    memSetVertices_conformeInterp[this->domainSpace()->worldCommPtr()->globalRank()].push_back(theImageElt.vertices());
                            }

                            dof_done[gdof]=true;
                        }
                    } // for ( uint16_type comp ... )
                } // for ( uint16_type iloc ... )
            } // for ( ; it!=en;++it )
        } //for ( ; itListRange!=enListRange ; ++itListRange)
    } // isActive

    // memory map (loc index pt) -> global dofs
    std::vector< std::vector<size_type> > memmapGdof( nProc_domain );
    // memory map (loc index pt) -> comp
    std::vector< std::vector<uint16_type> > memmapComp( nProc_domain );
    // points to lacalize
    std::vector<std::vector<typename image_mesh_type::node_type> > pointsSearched( nProc_domain );

    std::vector<std::vector< std::vector<typename image_mesh_type::node_type > > > memmap_vertices( nProc_domain );

    for (int proc=0; proc<nProc_domain;++proc)
        {
            const size_type nData = memSetGdofAndComp[proc].size();
            memmapGdof[proc].resize(nData);
            memmapComp[proc].resize(nData);
            pointsSearched[proc].resize(nData);
            // conforme case
            if(InterpType::isConforming()) memmap_vertices[proc].resize(nData);

            auto it_GdofAndComp = memSetGdofAndComp[proc].begin();
            auto it_vertices = memSetVertices_conformeInterp[proc].begin();
            for (int k=0 ; k<nData ; ++k, ++it_GdofAndComp)
                {
                    memmapGdof[proc][k]=it_GdofAndComp->template get<0>();//gdof;
                    memmapComp[proc][k]=it_GdofAndComp->template get<1>();//comp
                    pointsSearched[proc][k]=imagedof->dofPoint(it_GdofAndComp->template get<0>()).template get<0>();//node
                    if(InterpType::isConforming()) // conforme case
                        {
                            memmap_vertices[proc][k].resize(it_vertices->size2());
                            for (uint16_type v=0;v<it_vertices->size2();++v)
                                memmap_vertices[proc][k][v]=ublas::column(*it_vertices,v);
                            ++it_vertices;
                        }
                }
        }

    return boost::make_tuple(memmapGdof,memmapComp,pointsSearched,memmap_vertices);
}



//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//


namespace detail
{
template <typename RangeType>
struct opinterprangetype
{
    typedef range_t<RangeType> type;
};
}


/**
 * get the type of an OperatorInterpolation
 * \code
 * operator_interpolation_t<space_1_type,space_2_type,elements_pid_t<typename space_2_type::mesh_type>, >
 * \endcode
 */
template<typename DomainSpaceType, typename ImageSpaceType, typename IteratorRange= elements_pid_t<typename ImageSpaceType::mesh_type>, typename InterpType = InterpolationNonConforming >
using operator_interpolation_t =
OperatorInterpolation<DomainSpaceType,
                      ImageSpaceType,
                      range_t<IteratorRange>,
                      InterpType>;

template<typename DomainSpaceType,
         typename ImageSpaceType,
         typename IteratorRange= elements_pid_t<typename ImageSpaceType::mesh_type>,
         typename InterpType = InterpolationNonConforming >
using I_t = operator_interpolation_t<DomainSpaceType,ImageSpaceType,IteratorRange,InterpType>;

template<typename DomainSpaceType,
         typename ImageSpaceType,
         typename IteratorRange= elements_pid_t<typename ImageSpaceType::mesh_type>,
         typename InterpType =InterpolationNonConforming>
using I_ptr_t = std::shared_ptr<I_t<DomainSpaceType,ImageSpaceType,IteratorRange,InterpType>>;


template<typename DomainSpaceType,
         typename ImageSpaceType,
         typename IteratorRange = elements_pid_t<typename ImageSpaceType::mesh_type>,
         typename InterpType = InterpolationGradient<nonconforming_t>>
using Grad_t =
OperatorInterpolation<DomainSpaceType,
                      ImageSpaceType,
                      range_t<IteratorRange>,
                      InterpType>;

template<typename DomainSpaceType,
         typename ImageSpaceType,
         typename IteratorRange = elements_pid_t<typename ImageSpaceType::mesh_type>,
         typename InterpType = InterpolationGradient<nonconforming_t>>
using Grad_ptr_t = std::shared_ptr<Grad_t<DomainSpaceType,
                                            ImageSpaceType,
                                            IteratorRange,
                                            InterpType>>;

template<typename DomainSpaceType,
         typename ImageSpaceType,
         typename IteratorRange = elements_pid_t<typename ImageSpaceType::mesh_type>,
         typename InterpType = InterpolationCurl<nonconforming_t>>
using Curl_t =
OperatorInterpolation<DomainSpaceType,
                      ImageSpaceType,
                      range_t<IteratorRange>,
                      InterpType>;

template<typename DomainSpaceType,
         typename ImageSpaceType,
         typename IteratorRange = elements_pid_t<typename ImageSpaceType::mesh_type>,
         typename InterpType = InterpolationCurl<nonconforming_t>>
using Curl_ptr_t = std::shared_ptr<Curl_t<DomainSpaceType,
                                            ImageSpaceType,
                                            IteratorRange,
                                            InterpType>>;

template<typename DomainSpaceType,
         typename ImageSpaceType,
         typename IteratorRange = elements_pid_t<typename ImageSpaceType::mesh_type>,
         typename InterpType = InterpolationDiv<nonconforming_t>>
using Div_t =
OperatorInterpolation<DomainSpaceType,
                      ImageSpaceType,
                      range_t<IteratorRange>,
                      InterpType>;

template<typename DomainSpaceType,
         typename ImageSpaceType,
         typename IteratorRange = elements_pid_t<typename ImageSpaceType::mesh_type>,
         typename InterpType = InterpolationDiv<nonconforming_t>>
using Div_ptr_t = std::shared_ptr<Div_t<DomainSpaceType,
                                          ImageSpaceType,
                                          IteratorRange,
                                          InterpType>>;
#if 0
template<typename DomainSpaceType, typename ImageSpaceType, typename IteratorRange, typename InterpType >
decltype(auto)
opInterp( std::shared_ptr<DomainSpaceType> const& domainspace,
          std::shared_ptr<ImageSpaceType> const& imagespace,
          IteratorRange const& r,
          typename OperatorInterpolation<DomainSpaceType, ImageSpaceType,typename Feel::detail::opinterprangetype<IteratorRange>::type,InterpType>::backend_ptrtype const& backend,
          InterpType const& interptype,
          bool ddmethod )
{
    typedef OperatorInterpolation<DomainSpaceType,
                                  ImageSpaceType,
                                  typename Feel::detail::opinterprangetype<IteratorRange>::type,
                                  InterpType> operatorinterpolation_type;

    operatorinterpolation_type o ( domainspace,
                                   imagespace,
                                   r,
                                   backend,
                                   interptype,
                                   ddmethod );
    return o;
}
#endif
template<typename DomainSpaceType, typename ImageSpaceType, typename IteratorRange, typename InterpType >
decltype(auto)
opInterpPtr( std::shared_ptr<DomainSpaceType> const& domainspace,
             std::shared_ptr<ImageSpaceType> const& imagespace,
             IteratorRange const& r,
             typename OperatorInterpolation<DomainSpaceType, ImageSpaceType,typename Feel::detail::opinterprangetype<IteratorRange>::type,InterpType>::backend_ptrtype const& backend,
             InterpType const& interptype,
             bool ddmethod,
             OperatorInterpolationMatrixSetup<typename DomainSpaceType::value_type> const& matSetup = OperatorInterpolationMatrixSetup<typename DomainSpaceType::value_type>{} )
{
    typedef OperatorInterpolation<DomainSpaceType,
                                  ImageSpaceType,
                                  typename Feel::detail::opinterprangetype<IteratorRange>::type,
                                  InterpType> operatorinterpolation_type;

    auto opI = std::make_shared<operatorinterpolation_type>( domainspace,
                                                               imagespace,
                                                               r,
                                                               backend,
                                                               interptype,
                                                             ddmethod,
                                                             matSetup );

    return opI;
}


template <typename ... Ts>
auto opInterpolation( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && domainSpace = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_domainSpace);
    auto && imageSpace = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_imageSpace);
    auto && range = args.get_else_invocable(_range, [&imageSpace]() { return elements(support(imageSpace)); } );
    using domain_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(domainSpace)>>>;
    using image_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(imageSpace)>>>;
    auto && backend = args.get_else_invocable(_backend, [](){ return Feel::backend(); } );//Backend<typename domain_space_type::value_type>::build( soption( _name="backend" ) ) );
    auto && type = args.get_else(_type,InterpolationNonConforming() );
    bool ddmethod = args.get_else(_ddmethod,false);
    auto && matrix = args.get_else(_matrix, OperatorInterpolationMatrixSetup<typename image_space_type::value_type>{});

    std::decay_t<decltype(type)> t( type );
    return opInterpPtr( domainSpace,imageSpace,range,backend,std::move(t),ddmethod,matrix );
}

template <typename ... Ts>
auto IPtr( Ts && ... v )
{
    return opInterpolation( std::forward<Ts>(v)... );
}
template <typename ... Ts>
auto I( Ts && ... v )
{
    return *opInterpolation( std::forward<Ts>(v)... );
}

template <typename ... Ts>
auto Grad( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && domainSpace = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_domainSpace);
    auto && imageSpace = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_imageSpace);
    auto && range = args.get_else_invocable(_range, [&imageSpace]() { return elements(support(imageSpace)); } );
    using domain_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(domainSpace)>>>;
    using image_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(imageSpace)>>>;
    auto && backend = args.get_else_invocable(_backend, [](){ return Feel::backend(); } );
    auto && type = args.get_else(_type, makeGradientInterpolation<nonconforming_t>( nonconforming_t() ) );
    bool ddmethod = args.get_else(_ddmethod,false);
    auto && matrix = args.get_else(_matrix, OperatorInterpolationMatrixSetup<typename image_space_type::value_type>{});

    std::decay_t<decltype(type)> t( type );
    return *opInterpPtr( domainSpace,imageSpace,range,backend,std::move(t),ddmethod,matrix );
}

template <typename ... Ts>
auto Curl( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && domainSpace = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_domainSpace);
    auto && imageSpace = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_imageSpace);
    auto && range = args.get_else_invocable(_range, [&imageSpace]() { return elements(support(imageSpace)); } );
    using domain_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(domainSpace)>>>;
    using image_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(imageSpace)>>>;
    auto && backend = args.get_else_invocable(_backend, [](){ return Feel::backend(); } );
    auto && type = args.get_else(_type, makeCurlInterpolation<nonconforming_t>( nonconforming_t() ) );
    bool ddmethod = args.get_else(_ddmethod,false);
    auto && matrix = args.get_else(_matrix, OperatorInterpolationMatrixSetup<typename image_space_type::value_type>{});

    std::decay_t<decltype(type)> t( type );
    return *opInterpPtr( domainSpace,imageSpace,range,backend,std::move(t),ddmethod,matrix );
}

template <typename ... Ts>
auto Div( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && domainSpace = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_domainSpace);
    auto && imageSpace = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_imageSpace);
    auto && range = args.get_else_invocable(_range, [&imageSpace]() { return elements(support(imageSpace)); } );
    using domain_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(domainSpace)>>>;
    using image_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(imageSpace)>>>;
    auto && backend = args.get_else_invocable(_backend, [](){ return Feel::backend(); } );
    auto && type = args.get_else(_type, makeDivInterpolation<nonconforming_t>( nonconforming_t() ) );
    bool ddmethod = args.get_else(_ddmethod,false);
    auto && matrix = args.get_else(_matrix, OperatorInterpolationMatrixSetup<typename image_space_type::value_type>{});

    std::decay_t<decltype(type)> t( type );
    return *opInterpPtr( domainSpace,imageSpace,range,backend,std::move(t),ddmethod,matrix );
}


} // Feel
#endif /* FEELPP_OPERATORINTERPOLATION_HPP */
