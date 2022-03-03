/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 This file is part of the Feel library

 Copyright (C) 2010 University of Coimbra

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
 \file ale.hpp
 \author Goncalo Pena <gpena@mat.uc.pt>
 \date 2010-10-12
 */

#ifndef FEELPP_MODELS_ALE_H
#define FEELPP_MODELS_ALE_H 1

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelmodels/modelcore/modelbase.hpp>
#include <feel/feelmodels/modelmesh/metricmeshadaptation.hpp>

namespace Feel
{
namespace FeelModels
{

template< class Convex, int Order = 1 >
class ALE : public ModelBase
{
public :
    typedef ModelBase super_type;

    typedef ALE< Convex,Order> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef Convex convex_type;
    static const uint16_type Dim = convex_type::nDim;
    static const uint16_type Order_low = convex_type::nOrder;
    typedef Mesh< convex_type > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    using range_elements_type = elements_reference_wrapper_t<mesh_type>;
    using size_type = typename mesh_type::size_type;
    using bc_to_markers_type = std::map< std::string, std::set<std::string> >;

protected :

    template< int i >
    struct MyReferenceFunctionSpace
    {
        typedef bases<Lagrange<i, Vectorial> > basis_type;
        typedef FunctionSpace< mesh_type, basis_type > type;
        typedef std::shared_ptr<type> ptrtype;
        typedef typename type::element_type elt_type;
        typedef std::shared_ptr<elt_type> elt_ptrtype;

    };

    static const int orderUse = mpl::if_< mpl::greater<mpl::int_<Order>,mpl::int_<Order_low> >,
                                          mpl::int_<Order>,
                                          mpl::int_<Order_low> >::type::value ;

public:

    typedef typename MyReferenceFunctionSpace<orderUse>::basis_type ale_map_basis_type;
    typedef typename MyReferenceFunctionSpace<orderUse>::type ale_map_functionspace_type;
    typedef typename MyReferenceFunctionSpace<orderUse>::ptrtype ale_map_functionspace_ptrtype;
    typedef typename MyReferenceFunctionSpace<orderUse>::elt_type ale_map_element_type;


    typedef MetricMeshAdaptation<convex_type> metricmeshadaptation_type;
    typedef std::shared_ptr<metricmeshadaptation_type> metricmeshadaptation_ptrtype;

    /**
     * constructor,copy,desctructor
     */
    explicit ALE( std::string const& prefix="",
                  worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(),
                  ModelBaseRepository const& modelRep = ModelBaseRepository() );
    ALE( ALE const& tc ) = default;
    ALE( ALE && tc ) = default;
    //~ALE();

    /**
     * static builder
     */
    static self_ptrtype build(mesh_ptrtype mesh, std::string const& prefix="",
                              worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(),
                              ModelBaseRepository const& modelRep = ModelBaseRepository() );
    static self_ptrtype build(mesh_ptrtype mesh, range_elements_type const& rangeElt, std::string const& prefix="",
                              worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(),
                              ModelBaseRepository const& modelRep = ModelBaseRepository() );

    /**
     * Add the set of flags that mark the boundary
     */
    void addMarkerInBoundaryCondition( std::string const& bc, std::string const& marker );

    template <typename MarkersType,std::enable_if_t< is_iterable_v<MarkersType>, bool> = true>
    void addMarkersInBoundaryCondition( std::string const& bc, MarkersType const& markers )
        {
            for ( std::string const& marker : markers )
                this->addMarkerInBoundaryCondition( bc, marker );
        }

    bc_to_markers_type const& bcToMarkers() const { return M_bcToMarkers; }
    std::set<std::string> const& markers( std::string const& bc ) const
        {
            auto itFind = M_bcToMarkers.find( bc );
            CHECK( itFind != M_bcToMarkers.end() ) << "no valid bc type : " << bc;
            return itFind->second;
        }

    bool hasBoundaryCondition( std::string const& bc ) const { return M_bcToMarkers.find( bc ) != M_bcToMarkers.end(); }

    virtual void init() = 0;

    virtual ale_map_functionspace_ptrtype const& functionSpace() const = 0;
    virtual ale_map_element_type const& displacement() const = 0;
    virtual ale_map_element_type const& map() const = 0;

    virtual void generateMap( ale_map_element_type const & dispOnBoundary,
                              ale_map_element_type const & oldDisp ) = 0;

    template < typename ExprT >
    void updateMetricMeshAdaptation( Expr<ExprT> const& e )
        {
            if ( !M_metricMeshAdaptation )
                this->initMetricMeshAdaptation();
            this->updateMetricMeshAdaptation( *M_metricMeshAdaptation->feProjection( e ) );
        }
    virtual void initMetricMeshAdaptation() = 0;
    virtual void updateMetricMeshAdaptation( Expr<GinacExVF<2>> const& e ) = 0;
    virtual void updateMetricMeshAdaptation( typename metricmeshadaptation_type::element_scalar_type const& u ) = 0;

protected :
    metricmeshadaptation_ptrtype M_metricMeshAdaptation;
private :
    bc_to_markers_type M_bcToMarkers;
};
} // namespace FeelModels
} // namespace Feel
#endif // FEELPP_MODELS_ALE_H
