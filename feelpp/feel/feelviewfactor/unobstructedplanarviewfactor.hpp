/**
 * @file unobstructedplanarviewfactor.hpp
 * @author Christophe Prud'homme
 * @brief
 * @version 0.1
 * @date 2022-07-22
 *
 * @copyright Copyright (c) 2022 Feel++ Consortium
 * @copyright Copyright (c) 2022 Université de Strasbourg
 *
 */

#pragma once

#include <feel/feelcore/enumerate.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feeldiscr/context.hpp>
#include <feel/feelviewfactor/viewfactorbase.hpp>

namespace Feel {

/**
 * @brief  Computes the view factors for planar faces in unobstructed radiative heat transfer
 * @ingroup ViewFactor
 * @tparam MeshType the mesh type
 * 
 * @code
 * nl::json j{"viewfactor":{}};
 * UnobstructedPlanarViewFactor<Mesh<Simplex<2>>> vf(mesh, specs);
 * vf.compute();
 * @endcode
 */
template<typename MeshType>
class UnobstructedPlanarViewFactor : public ViewFactorBase<MeshType>
{
public:
    using super = ViewFactorBase<MeshType>;
    using value_type = double;
    using mesh_t = typename super::mesh_t;
    using mesh_ptrtype = typename super::mesh_ptrtype;
    using mesh_ptr_t = typename super::mesh_ptr_t;
    UnobstructedPlanarViewFactor() = default;
    UnobstructedPlanarViewFactor( mesh_ptr_t mesh, nl::json const& specs )
        : 
        super( mesh, specs )
        {}
    UnobstructedPlanarViewFactor( const UnobstructedPlanarViewFactor& ) = default;
    UnobstructedPlanarViewFactor( UnobstructedPlanarViewFactor&& ) = default;
    UnobstructedPlanarViewFactor& operator=( const UnobstructedPlanarViewFactor& ) = default;
    UnobstructedPlanarViewFactor& operator=( UnobstructedPlanarViewFactor&& ) = default;
    ~UnobstructedPlanarViewFactor() = default;
//    void init( std::vector<std::string> const& list_of_bdys ) { ViewFactorBase::init( list_of_bdys ); }

    void compute();
};

template<typename MeshType>
void 
UnobstructedPlanarViewFactor<MeshType>::compute()
{
    auto the_im = im( this->mesh_, 2 );
    auto current_pts = [&the_im]( auto f, auto p )
    {
        return the_im.fpoints( f, p.value() );
    };
    for(auto const& [current_index,current_side] : enumerate(this->list_of_bdys_) )
    {
        if ( !this->mesh_->hasMarker( current_side ) )
        {
            throw std::logic_error( "boundary marker " + current_side + " does not exist in mesh" );
        }
        auto current_range = markedfaces( this->mesh_, current_side );
        if ( begin(current_range) != end(current_range) )
        {
            auto current_ctx = context( _element = boost::unwrap_ref( *begin( current_range ) ), _type = on_facets_t(), _geomap = this->mesh_->gm(), _pointset = current_pts );

            for ( auto const& [remote_index, remote_side] : enumerate( this->list_of_bdys_ ) )
            {
                if ( !this->mesh_->hasMarker( current_side ) )
                {
                    throw std::logic_error( "remote boundary marker " + current_side + " does not exist in mesh" );
                }

                // planar surface don't see themselves, hence the view factor is 0
                if ( current_index == remote_index )
                    this->vf_( current_index, remote_index ) = 0.0;
                else
                {
                    auto remote_range = markedfaces( this->mesh_, remote_side );
                    if ( begin(remote_range) != end(remote_range ) )
                    {
                        auto remote_ctx = context( _element = boost::unwrap_ref( *begin( remote_range ) ), _type = on_facets_t(), _geomap = this->mesh_->gm(), _pointset = current_pts );
                        for( auto const& current_wface : markedfaces( this->mesh_, current_side ) )
                        {
                            auto const& current_face = boost::unwrap_ref( current_wface );
                            // get the area of the face
                            this->areas_[current_index] += current_face.measure();
                            auto const& current_pts = emap<value_type>( current_ctx->xReal() );
                            for( auto const& remote_wface : remote_range )
                            {
                                auto const& remote_face = boost::unwrap_ref( remote_wface );
                                // get the area of the face
                                //areas_[remote_index] += remote_face.measure();
                                // compute the view factor
                                
                                auto const& remote_pts = emap<value_type>( remote_ctx->xReal() );

                                for( auto const& [current_index,current_pt]: enumerate(current_pts.colwise()))
                                {
                                    for ( auto const& [remote_index, remote_pt] : enumerate( remote_pts.colwise() ) )
                                    {
                                        auto p2p = (current_pt - remote_pt);
                                        value_type dist = p2p.norm();
                                        // get angles between the two faces
                                        value_type cos1 = p2p.dot( current_ctx->normal(current_index) ) / dist; 
                                        value_type cos2 = p2p.dot( remote_ctx->normal(remote_index) ) / dist;
                                        // add contribution to integrals
                                        this->vf_( current_index, remote_index ) +=
                                                    current_ctx->J( current_index )  * remote_ctx->J( remote_index ) * 
                                                    std::abs( cos1 ) * std::abs( cos2 ) / std::pow( dist, this->exponent_ );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
} // namespace Feel