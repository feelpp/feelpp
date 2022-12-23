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
#include <feel/feelvf/vf.hpp>

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
    auto the_im = im( this->mesh_, this->j_["viewfactor"]["quadrature_order"] );
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
            // constexpr size_type context_v = vm::POINT|vm::NORMAL;
            auto current_ctx = context( _element = boost::unwrap_ref( *begin( current_range ) ), _type = on_facets_t(), _geomap = this->mesh_->gm(),/* _context_at_compile_time = context_v, */_pointset = current_pts );
            for ( auto const& [remote_index, remote_side] : enumerate( this->list_of_bdys_ ) )
            {
                if ( !this->mesh_->hasMarker( remote_side ) )
                {
                    throw std::logic_error( "remote boundary marker " + remote_side + " does not exist in mesh" );
                }

                // planar surface don't see themselves, hence the view factor is 0
                if ( current_index == remote_index && this->j_["viewfactor"]["algorithm"] == "DoubleAreaIntegration")
                    this->vf_( current_index, remote_index ) = 0.0;                                        
                else
                {
                    this->vf_( current_index, remote_index ) = 0.0;
                    this->areas_[current_index] = 0.0;
                    auto remote_range = markedfaces( this->mesh_, remote_side );
                    if ( begin(remote_range) != end(remote_range ) )
                    {                        
                        auto remote_ctx = context( _element = boost::unwrap_ref( *begin( remote_range ) ), _type = on_facets_t(), _geomap = this->mesh_->gm(), /* _context_at_compile_time = context_v,*/ _pointset = current_pts );
                        for( auto const& current_wface : markedfaces( this->mesh_, current_side ) )
                        {
                            auto const& current_face = boost::unwrap_ref( current_wface );
                            // get the area of the face
                            this->areas_[current_index] += current_face.measure();
                            double face_area = current_face.measure();
                            
                            current_ctx->template update<vm::POINT|vm::NORMAL|vm::JACOBIAN>(current_face.element0(),current_face.idInElement0());
                            auto const& current_pts = emap<value_type>( current_ctx->xReal() );
                            for( auto const& remote_wface : remote_range )
                            {
                                auto const& remote_face = boost::unwrap_ref( remote_wface );
                                // get the area of the face
                                //areas_[remote_index] += remote_face.measure();
                                // compute the view factor
                                remote_ctx->template update<vm::POINT|vm::NORMAL|vm::JACOBIAN>(remote_face.element0(),remote_face.idInElement0());
                                auto const& remote_pts = emap<value_type>( remote_ctx->xReal() );

                                for( auto const& [current_index_1,current_pt]: enumerate(current_pts.colwise()))
                                {
                                    if(this->j_["viewfactor"]["algorithm"] == "DoubleAreaIntegration")
                                    {
                                        // Method 2AI - numerical integration of the double integral
                                        for ( auto const& [remote_index_1, remote_pt] : enumerate( remote_pts.colwise() ) )
                                        {                                        
                                            auto p2p = (-current_pts.col(current_index_1) + remote_pts.col(remote_index_1));                                                                             
                                            value_type dist = p2p.norm();
                                            // get angles between the two faces                                                                                
                                            value_type cos1 = p2p.dot( current_ctx->normal(current_index_1) ) / dist;                                         
                                            value_type cos2 = -p2p.dot( remote_ctx->normal(remote_index_1) ) / dist;
                                            // add contribution to integrals 
                                            // w_near * w_far * Jac_near * Jac_far * cos1 * cos2 / (r^exponent * divisor)
                                            this->vf_( current_index, remote_index ) +=
                                                        the_im.weight(current_face.idInElement0(),current_index_1) * the_im.weight(remote_face.idInElement0(),remote_index_1) * 
                                                        current_ctx->J( current_index_1 )  * remote_ctx->J( remote_index_1 ) * 
                                                        cos1 *  cos2  / (std::pow( dist, this->exponent_ )*(this->divisor_));                                                                              
                                        }   
                                    }    
                                    else if (this->j_["viewfactor"]["algorithm"] == "SingleAreaIntegration" && this->mesh_->dimension()==3)
                                    {
                                        // Method 1AI - Hottel and Sarofim, 1967. Radiative Transfer p.48
                                        // F_12 = 1/(2pi A_1) \int_{A_2} \sum_i^{edges A_2} (g_i \cdot n_1)                                                                        
                                        // Only for 3D
                                        for (int i = 0; i< remote_face.nEdges() ; i++)
                                        {
                                            Eigen::VectorXd p0(3);
                                            p0 << remote_face.point(i%3).node()[0],remote_face.point(i%3).node()[1],remote_face.point(i%3).node()[2];
                                            Eigen::VectorXd p1(3);
                                            p1 << remote_face.point((i+1)%3).node()[0],remote_face.point((i+1)%3).node()[1],remote_face.point((i+1)%3).node()[2];
                                            auto p = current_pts.col(current_index_1);                                                                         

                                            auto dot_p = (p0-p).dot(p1-p);
                                            auto cross_p = ((p0-p).template head<3>()).cross((p1-p).template head<3>());
                                            //auto cross_p = cross_product((p0-p),(p1-p));
                                            double norm_cross_p = cross_p.norm();
                                            Eigen::VectorXd normal_vec(3);
                                            if(this->mesh_->dimension()==3)
                                                normal_vec << current_ctx->normal(current_index_1)[0], current_ctx->normal(current_index_1)[1], current_ctx->normal(current_index_1)[2];
                                            auto scalar_prod = cross_p.dot(normal_vec)/norm_cross_p;

                                            this->vf_( current_index, remote_index ) += the_im.weight(current_face.idInElement0(),current_index_1) * 
                                            current_ctx->J( current_index_1 )  * scalar_prod * (pi/2 - math::atan(dot_p/norm_cross_p))/ (2*pi);
                                        }
                                    } 
                                    else 
                                    {
                                        throw std::logic_error( "Integration method not specified. Choose one between the available ones." );
                                    }                                  
                                }
                            }
                        }
                    }
                }                  
                this->vf_( current_index, remote_index ) /= integrate(_range = current_range, _expr = cst(1.0) ).evaluate()(0,0);              
                if(this->j_["viewfactor"]["algorithm"] == "SingleAreaIntegration" && current_index==remote_index)
                {
                    this->vf_( current_index, remote_index ) = math::abs(1 + this->vf_( current_index, remote_index )) ;
                }
            }            
        }        
    }
}
} // namespace Feel