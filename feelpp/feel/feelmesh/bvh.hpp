// Bounding volume hierarchy

#ifndef FEELPP_MESH_BVH_HPP
#define FEELPP_MESH_BVH_HPP

#include <vector>

#include <bvh/v2/bvh.h>
#include <bvh/v2/default_builder.h>
#include <bvh/v2/stack.h>
#include <bvh/v2/tri.h>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feeldiscr/mesh.hpp>

namespace Feel
{

#if 0
// https://en.wikipedia.org/wiki/Orthant
// https://github.com/madmann91/bvh/blob/master/src/bvh/v2/ray.h
struct Orthant
{
    std::uint32_t value = 0;
    static constexpr std::size_t max_dim = sizeof(value) * CHAR_BIT;
    std::uint32_t operator [] (std::size_t i) const { return (value >> i) & std::uint32_t{1}; }
};
#endif

enum class BvhIntersectContext{ anyHint=0, closest, all };

struct BVHEnum
{
    enum class Quality { Low, Medium, High };
};

template <int RealDim>
class BVHRay
{
public:
    using vec_t = eigen_vector_type<RealDim>;
    BVHRay(vec_t const& orig, vec_t const& dir,
           double dmin = 0, double dmax = std::numeric_limits<double>::max() )
        :
        M_origin( orig ),
        M_dir( dir ),
        M_distanceMin( dmin ),
        M_distanceMax( dmax )
        {}
    BVHRay() : BVHRay(vec_t::Zero(),vec_t::Zero()) {}

    BVHRay( BVHRay const& ) = default;
    BVHRay( BVHRay &&) = default;
    BVHRay& operator=( BVHRay && ) = default;
    BVHRay& operator=( BVHRay const& ) = default;

    vec_t const& origin() const noexcept { return M_origin; }
    vec_t const& dir() const noexcept { return M_dir; }
    double distanceMin() const { return M_distanceMin; }
    double distanceMax() const { return M_distanceMax; }

#if 0
    Orthant orthant() const {
        static_assert(RealDim <= Orthant::max_dim);
        Orthant orthant;
        for (int i=0;i<RealDim;++i)
            orthant.value |= std::signbit(M_dir[i]) * (std::uint32_t{1} << i);
        return orthant;
    }
#endif
private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize( Archive& ar, const unsigned int version )
        {
            ar & M_origin;
            ar & M_dir;
            ar & M_distanceMin;
            ar & M_distanceMax;
        }
private:
    vec_t M_origin, M_dir;  // ray origin and dir
    double M_distanceMin, M_distanceMax;
};


template <int RealDim>
class BVHRaysDistributed
{
public:
    using ray_type = BVHRay<RealDim>;
    BVHRaysDistributed() = default;
    BVHRaysDistributed( BVHRaysDistributed && ) = default;
    std::vector<ray_type> const& rays() const { return M_rays; }
    //! return number of local ray
    std::size_t numberOfLocalRay() const { return M_rays.size(); }

    template <typename T>
    void push_back( T && ray ) { M_rays.push_back( std::forward<T>( ray ) ); }
private:
    std::vector<ray_type> M_rays;
};

//! @brief BVH base class
template <typename MeshEntityType>
class BVH : public CommObject
{
public:
    using mesh_entity_type = std::decay_t<MeshEntityType>;
    using index_type = typename mesh_entity_type::index_type;
    static constexpr uint16_type nDim = mesh_entity_type::nDim;
    static constexpr uint16_type nRealDim = mesh_entity_type::nRealDim;
    using value_type = double;
    using vector_realdim_type = Eigen::Matrix<value_type,nRealDim,1>;
    using ray_type = BVHRay<nRealDim>;


    //! @brief Information on the primitive (mesh entity, bounding box, centroid)
    struct BVHPrimitiveInfo
    {
        BVHPrimitiveInfo( mesh_entity_type const& meshEntity )
            :
            M_meshEntity( meshEntity )
            {
                auto verticesUblas = meshEntity.vertices();
                auto G = em_cmatrix_col_type<double>( verticesUblas.data().begin(), nRealDim, mesh_entity_type::numVertices );
                M_bound_min = G.rowwise().minCoeff();
                M_bound_max = G.rowwise().maxCoeff();
                M_bound_min.array() -= 2*FLT_MIN;
                M_bound_max.array() += 2*FLT_MIN;
                //M_centroid = ( M_bound_min + M_bound_max ) * 0.5;
                auto bary = meshEntity.barycenter();
                M_centroid = Eigen::Map<Eigen::Matrix<double,nRealDim,1>>( bary.data().begin() );
            }
        BVHPrimitiveInfo( BVHPrimitiveInfo && ) = default;
        BVHPrimitiveInfo( BVHPrimitiveInfo const& ) = default;
        BVHPrimitiveInfo& operator=( BVHPrimitiveInfo && ) = default;
        BVHPrimitiveInfo& operator=( BVHPrimitiveInfo const& ) = default;

        mesh_entity_type const& meshEntity() const { return M_meshEntity.get(); }
        vector_realdim_type const& boundMin() const noexcept { return M_bound_min; }
        vector_realdim_type const& boundMax() const noexcept { return M_bound_max; }
        vector_realdim_type const& centroid() const noexcept { return M_centroid; }

    private:
        vector_realdim_type M_bound_min;
        vector_realdim_type M_bound_max;
        vector_realdim_type M_centroid;
        std::reference_wrapper<mesh_entity_type const> M_meshEntity;
    };
    using primitiveinfo_type = BVHPrimitiveInfo;

    //! @brief Data returned after apply an intersection with a ray
    struct BVHRayIntersectionResult
    {
        BVHRayIntersectionResult() = default;
        BVHRayIntersectionResult( rank_type processId, index_type primitiveId, double dist )
            :
            M_processId( processId ),
            M_primitiveId( primitiveId ),
            M_distance( dist )
            {}
        BVHRayIntersectionResult( BVHRayIntersectionResult && ) = default;
        BVHRayIntersectionResult( BVHRayIntersectionResult const& ) = default;
        BVHRayIntersectionResult& operator=( BVHRayIntersectionResult && ) = default;
        BVHRayIntersectionResult& operator=( BVHRayIntersectionResult const& ) = default;

        rank_type processId() const noexcept { return M_processId; }
        index_type primitiveId() const noexcept { return M_primitiveId; }
        double distance() const noexcept { return M_distance; }

        template <typename T>
        void setCoordinates( T && coord ) { M_coordinates = std::forward<T>( coord ); M_hasCoordinates = true; }
        // bool hasCoordinates() const noexcept { return M_coordinates.has_value(); }
        //   vector_realdim_type const& coordinates() const noexcept { return *M_coordinates; }
        bool hasCoordinates() const noexcept { return M_hasCoordinates; }
        vector_realdim_type const& coordinates() const noexcept { return M_coordinates; }

    private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize( Archive& ar, const unsigned int version )
            {
                ar & M_processId;
                ar & M_primitiveId;
                ar & M_distance;
                ar & M_hasCoordinates;
                ar & M_coordinates;
#if 0
                if constexpr ( Archive::is_saving::value )
                {
                    bool hasCoordinates = this->hasCoordinates();
                    ar & boost::serialization::make_nvp( "hasCoordinates", hasCoordinates );
                    if ( hasCoordinates )
                        ar & boost::serialization::make_nvp( "coordinates", this->coordinates() );
                }
                else if constexpr ( Archive::is_loading::value )
                {
                    bool hasCoordinates = false;
                    ar & boost::serialization::make_nvp( "hasCoordinates", hasCoordinates );
                    if ( hasCoordinates )
                    {
                        vector_realdim_type coord;
                        ar & boost::serialization::make_nvp( "coordinates", coord );
                        this->setCoordinates( std::move( coord ) );
                    }
                }
#endif
            }
    private:
        rank_type M_processId = invalid_v<rank_type>;
        index_type M_primitiveId = invalid_v<index_type>;
        double M_distance = std::numeric_limits<double>::max();
        //std::optional<vector_realdim_type> M_coordinates;
        bool M_hasCoordinates = false;
        vector_realdim_type M_coordinates = vector_realdim_type::Zero();
    };
    using rayintersection_result_type = BVHRayIntersectionResult;

    //enum class IntersectContext{ anyHint=0, closest, all };

    using IntersectContext = BvhIntersectContext;

    BVH( BVHEnum::Quality quality = BVHEnum::Quality::High, worldcomm_ptr_t worldComm = Environment::worldCommPtr() )
        :
        CommObject( worldComm ),
        M_quality( quality )
        {}
    BVH( BVH && ) = default;
    BVH( BVH const& ) = default;
    virtual ~BVH() {}

    //! return all primitive info
    std::vector<BVHPrimitiveInfo> const& primitiveInfo() const noexcept { return M_primitiveInfo; }

    //! return primitive info at index i
    BVHPrimitiveInfo const& primitiveInfo( index_type i ) const { return M_primitiveInfo.at( i ); }

    //! compute intersection(s) with a ray from the BVH built and return a vector of intersection result
    template<typename... Ts>
    auto intersect( Ts && ... v )
        {
            auto args = NA::make_arguments( std::forward<Ts>(v)... );
            auto && ray = args.get(_ray);
            bool useRobustTraversal = args.get_else(_robust,true);
            IntersectContext ctx = args.get_else(_context,IntersectContext::closest);
            bool parallel = args.get_else(_parallel,this->worldComm().size() > 1);
            value_type tolerance = args.get_else(_tolerance,std::numeric_limits<value_type>::epsilon());

            //bool closestOnly = ctx == IntersectContext::closest;
            using napp_ray_type = std::decay_t<decltype(ray)>;


            std::vector<std::vector<rayintersection_result_type>> ret;
            switch ( ctx )
            {
            case IntersectContext::closest:
            {
                if ( useRobustTraversal )
                    ret = intersectGenericImpl<IntersectContext::closest,true>( ray, tolerance, parallel );
                else
                    ret = intersectGenericImpl<IntersectContext::closest,false>( ray, tolerance, parallel);
                break;
            }
            case IntersectContext::anyHint:
            {
                if ( useRobustTraversal )
                    ret = intersectGenericImpl<IntersectContext::anyHint,true>( ray, tolerance, parallel );
                else
                    ret = intersectGenericImpl<IntersectContext::anyHint,false>( ray, tolerance, parallel );
                break;
            }
            case IntersectContext::all:
            {
                if ( useRobustTraversal )
                    ret = intersectGenericImpl<IntersectContext::all,true>( ray, tolerance, parallel );
                else
                    ret = intersectGenericImpl<IntersectContext::all,false>( ray, tolerance, parallel );
                break;
            }
            }

            // If the input ray is a single ray, return a single vector of intersection results
            if constexpr ( std::is_same_v<napp_ray_type,ray_type> )
                return ret.front();
            else
                return ret;
        }


    template <IntersectContext Ctx,bool useRobustTraversal,typename RayType>
    std::vector<std::vector<rayintersection_result_type>>
    intersectGenericImpl( RayType const& ray, value_type tolerance, bool parallel );

protected:

    struct TlasHit {
        rank_type rank;
        value_type distance;
    };

    //! Returns a list of MPI ranks whose bounding boxes intersect the ray (using TLAS)
    virtual std::vector<TlasHit> getHitPartitions(ray_type const& ray) const
        {
            // Sequential/default fallback: everyone is queried
            std::vector<TlasHit> hitPartitions(this->worldComm().size());
            std::for_each( hitPartitions.begin(), hitPartitions.end(), [this](TlasHit& hit) { hit.rank = this->worldComm().rank(); hit.distance = 0; } );
            return hitPartitions;
        }

    virtual std::vector<rayintersection_result_type> intersectSequential( ray_type const& rayon, value_type tolerance, bool useRobustTraversal = true ) = 0;

    template <typename RangeType>
    void
    updateForUse( RangeType const& range )
        {
            // From the mesh, build the bounding box info for each element and store it in
            // the structure BVHPrimitiveInfo
            M_primitiveInfo.clear();
            M_primitiveInfo.reserve( nelements(range) );
            for ( auto const& eltWrap : range )
            {
                auto const& e = unwrap_ref( eltWrap );
                M_primitiveInfo.push_back( BVHPrimitiveInfo{e} );
            }
        }

protected:
    std::vector<BVHPrimitiveInfo> M_primitiveInfo;
    BVHEnum::Quality M_quality = BVHEnum::Quality::High;
};



//! @brief implementation of BVH tool with an external third party
template <typename MeshEntityType>
class BVH_ThirdParty : public BVH<MeshEntityType>
{
    using super_type = BVH<MeshEntityType>;
    using mesh_entity_type = typename super_type::mesh_entity_type;
    using vector_realdim_type = typename super_type::vector_realdim_type;
    static constexpr uint16_type nRealDim = super_type::nRealDim;

    using value_type = typename super_type::value_type;
    using node_type = bvh::v2::Node<value_type, nRealDim>;
    using backend_bvh_type = bvh::v2::Bvh<node_type>;
    using backend_vector_realdim_type = bvh::v2::Vec<value_type, nRealDim>;
    using backend_precompute_triangle_type = bvh::v2::PrecomputedTri<value_type>;
    using backend_bbox_type = bvh::v2::BBox<value_type, nRealDim>;

    using TlasHit = typename super_type::TlasHit;
public:
    using ray_type = typename super_type::ray_type;
    using rayintersection_result_type = typename super_type::rayintersection_result_type;

    BVH_ThirdParty( BVHEnum::Quality quality, worldcomm_ptr_t worldComm ) : super_type( quality,worldComm ) {}
    BVH_ThirdParty( BVH_ThirdParty && ) = default;

    template <typename RangeType>
    void
    updateForUse( RangeType const& range )
        {
            // up primitiveinfos
            super_type::updateForUse( range );

            // init bvh backend
            using BBox = backend_bbox_type;//bvh::v2::BBox<value_type, nRealDim>;
            std::vector<BBox> bboxes;
            std::vector<backend_vector_realdim_type> centers;
            bboxes.reserve( this->M_primitiveInfo.size() );
            centers.reserve( this->M_primitiveInfo.size() );
            for ( auto const& primInfo : this->M_primitiveInfo )
            {
                bboxes.push_back( BBox{
                        backend_vector_realdim_type::generate([&primInfo] (std::size_t i) { return primInfo.boundMin()[i]; }),
                            backend_vector_realdim_type::generate([&primInfo] (std::size_t i) { return primInfo.boundMax()[i]; })
                            });

                auto const& centroid = primInfo.centroid();
                centers.push_back( backend_vector_realdim_type::generate([&centroid] (std::size_t i) { return centroid[i]; }) );
            }

            typename bvh::v2::DefaultBuilder<node_type>::Config config;
            switch ( this->M_quality )
            {
            default:
            case BVHEnum::Quality::High: config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::High; break;
            case BVHEnum::Quality::Medium: config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::Medium; break;
            case BVHEnum::Quality::Low: config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::Low; break;
            }
            if ( !bboxes.empty() )
                M_bvh = std::make_unique<backend_bvh_type>( bvh::v2::DefaultBuilder<node_type>::build(/*thread_pool,*/ bboxes, centers, config) );
            else
                M_bvh.reset();

            // Permuting the primitive data allows to remove indirections during traversal, which makes it faster.
            static constexpr bool should_permute = true;

            if constexpr ( nRealDim == 3 )
                M_precomputeTriangle.resize( this->M_primitiveInfo.size() );

            for ( std::size_t i = 0; i < this->M_primitiveInfo.size(); ++i )
            {
                auto j = should_permute ? M_bvh->prim_ids[i] : i;
                auto const& primInfo = this->M_primitiveInfo[j];
                auto const& meshEntity = primInfo.meshEntity();
                if constexpr ( nRealDim == 3 )
                {
                    auto const& pt0 = meshEntity.point(0);
                    auto const& pt1 = meshEntity.point(1);
                    auto const& pt2 = meshEntity.point(2);
                    M_precomputeTriangle[i] = backend_precompute_triangle_type{
                        backend_vector_realdim_type::generate([&pt0] (std::size_t i) { return pt0[i]; }),
                        backend_vector_realdim_type::generate([&pt1] (std::size_t i) { return pt1[i]; }),
                        backend_vector_realdim_type::generate([&pt2] (std::size_t i) { return pt2[i]; })
                    };
                }
            }

            // in //, we use TLAS/blas strategy to cull the processes that cannot be hit by the ray,
            // so we need to init the TLAS object that contains the bounding boxes of each process.
            if ( this->worldComm().size() > 1 )
            {
                vector_realdim_type localBoundMin = vector_realdim_type::Constant( std::numeric_limits<value_type>::max() );
                vector_realdim_type localBoundMax = vector_realdim_type::Constant( std::numeric_limits<value_type>::lowest() );

                std::vector<value_type> localBBoxValues( 2 * nRealDim );
                if ( this->M_primitiveInfo.empty() )
                {
                    // fill as an invalid bbox to avoid degenerate TLAS node (can happen if some process have no element)
                    for ( int k = 0; k < nRealDim; ++k )
                    {
                        localBBoxValues[k] = 1;
                        localBBoxValues[nRealDim + k] = -1;
                    }
                }
                else
                {
                    // init TLAS by using the BBoxes of the primitives that are already secured
                    for ( auto const& primInfo : this->M_primitiveInfo )
                    {
                        for ( int k = 0; k < nRealDim; ++k )
                        {
                            localBoundMin[k] = std::min( localBoundMin[k], primInfo.boundMin()[k] );
                            localBoundMax[k] = std::max( localBoundMax[k], primInfo.boundMax()[k] );
                        }
                    }
                    double diag = (localBoundMax - localBoundMin).norm();
                    double epsilon = std::max(1e-5, diag * 1e-4); // 0.01% de la taille locale
                    //double epsilon = 2*FLT_MIN;
                    for ( int k = 0; k < nRealDim; ++k )
                    {
                        localBBoxValues[k] = localBoundMin[k] - epsilon;
                        localBBoxValues[nRealDim + k] = localBoundMax[k] + epsilon;
                    }
                }
                std::vector<value_type> globalBBoxValues( 2 * nRealDim * this->worldComm().size() );
                mpi::all_gather( this->worldComm(), localBBoxValues.data(), localBBoxValues.size(), globalBBoxValues );
                std::vector<BBox> tlasBBoxes;
                tlasBBoxes.reserve( this->worldComm().size() );
                std::vector<backend_vector_realdim_type> tlasCenters;
                tlasCenters.reserve( this->worldComm().size() );
                //std::vector<rank_type> primitiveInfoTlas;
                std::vector<std::tuple<rank_type,backend_bbox_type>> primitiveInfoTlas;
                primitiveInfoTlas.reserve( this->worldComm().size() );
                for ( rank_type p = 0; p < this->worldComm().size(); ++p )
                {
                    // ignore empty bbox (can happen if some process have no element) to avoid degenerate TLAS node
                    auto bboxValues = std::span<value_type>( globalBBoxValues.data() + p * 2 * nRealDim, 2 * nRealDim );
                    if ( bboxValues[0] > bboxValues[nRealDim] ) // invalid bbox, fill as an invalid bbox to avoid degenerate TLAS node
                    {
                        // std::cout << "Warning: process " << p << " has an empty bounding box " << bboxValues[0] << " vs " << bboxValues[nRealDim]
                        //           << " for dim 0, " << bboxValues[1] << " vs " << bboxValues[nRealDim + 1] << " for dim 1, " << bboxValues[2] << " vs " << bboxValues[nRealDim + 2] << " for dim 2" << std::endl;
                        continue;
                    }
                    tlasBBoxes.push_back( BBox{
                            backend_vector_realdim_type::generate([&bboxValues] (std::size_t i) { return bboxValues[i]; }),
                                backend_vector_realdim_type::generate([&bboxValues] (std::size_t i) { return bboxValues[nRealDim + i]; })
                                });
                    auto & bbox = tlasBBoxes.back();
                    tlasCenters.push_back( backend_vector_realdim_type::generate([&bbox] (std::size_t i) { return 0.5 * ( bbox.min[i] + bbox.max[i] ); }) );
                    //auto Eigen::Matrix<value_type,nRealDim,1>;
                    primitiveInfoTlas.push_back( std::make_tuple( p, bbox ) );
                }

                if ( this->worldComm().isMasterRank() && Environment::logVerbosityLevel() > 1 )
                {
                  for ( int i = 0; i < tlasBBoxes.size(); ++i )
                    {
                        auto const& bbox = tlasBBoxes[i];
                        //LOG(INFO) << "global bbox: min: " << bbox.min << " max: " << bbox.max;
                        std::cout << "global bbox: min: " << bbox.min[0] << " " << bbox.min[1] << " " << bbox.min[2]
                                  << " max: " << bbox.max[0] << " " << bbox.max[1] << " " << bbox.max[2]
                                  << " center: " << tlasCenters[i][0] << " " << tlasCenters[i][1] << " " << tlasCenters[i][2]
                                  << " volume: " << (bbox.max[0]-bbox.min[0])*(bbox.max[1]-bbox.min[1])*(bbox.max[2]-bbox.min[2])
                                  << std::endl;
                    }
                }

                typename bvh::v2::DefaultBuilder<node_type>::Config tlasConfig;
                tlasConfig.quality = bvh::v2::DefaultBuilder<node_type>::Quality::High;
                tlasConfig.max_leaf_size = 1; // <-- AJOUT CRUCIAL : force 1 processus par feuille
                if ( tlasBBoxes.size() > 0 )
                    M_bvhTlas = std::make_unique<backend_bvh_type>( bvh::v2::DefaultBuilder<node_type>::build(/*thread_pool,*/ tlasBBoxes, tlasCenters, tlasConfig) );
                else
                    M_bvhTlas.reset();

                // apply permutation to primitiveInfoTlas to be able to directly access it during TLAS traversal (without indirection)
                M_primitiveInfoTlas.resize( primitiveInfoTlas.size() );
                for ( std::size_t k = 0; k < primitiveInfoTlas.size(); ++k )
                {
                    M_primitiveInfoTlas[k] = primitiveInfoTlas[M_bvhTlas->prim_ids[k]];
                }

            } // if world size > 1
        }






    template <typename BVH<MeshEntityType>::IntersectContext Ctx,bool useRobustTraversal,typename RayType>
    std::vector<std::vector<rayintersection_result_type>>
    intersectFullImpl( RayType const& ray, value_type tolerance, bool parallel )
        {
            static constexpr bool closestOnly = Ctx == super_type::IntersectContext::closest;
            static constexpr bool isAnyHint = Ctx == super_type::IntersectContext::anyHint;

#if 1
            if constexpr( std::is_same_v<BVHRaysDistributed<nRealDim>,RayType/*napp_ray_type*/> ) // case rays distributed on process
            {
                auto const& localRays = ray.rays();
                std::vector<std::vector<rayintersection_result_type>> final_results(localRays.size());

                if ( !parallel || this->worldComm().size() <= 1 )
                {
                    for(size_t i = 0; i < localRays.size(); ++i) {
                        final_results[i] = this->intersectSequentialImpl<useRobustTraversal,isAnyHint>( localRays[i], tolerance );
                        if constexpr ( closestOnly )
                            if( final_results[i].size() > 1)
                                final_results[i].resize(1);
                    }
                    return final_results;
                }

                struct Target { int rank; };
                struct RayState {
                    std::vector<Target> targets;
                    int current_idx = 0;
                    bool is_waiting = false;
                    bool active = false;
                };

                std::vector<RayState> states(localRays.size());
                int my_active_rays = 0;
                std::queue<int> ready_rays; // Work Queue
                tic();
                // Initialization: TLAS and queue filling
                for (size_t i = 0; i < localRays.size(); ++i)
                {
                    std::vector<TlasHit> hit_ranks = this->getHitPartitions(localRays[i]);

                    for (auto const& hit : hit_ranks)
                        states[i].targets.push_back({hit.rank});

                    if (!states[i].targets.empty())
                    {
                        states[i].active = true;
                        my_active_rays++;
                        ready_rays.push((int)i); // The ray is ready for processing
                    }
                }
                toc("BVH Intersect - Initialization",Environment::logVerbosityLevel() > 0);
                tic();
                int global_active_rays = 0;
                mpi::all_reduce(this->worldComm(), my_active_rays, global_active_rays, std::plus<int>());
                toc("BVH Intersect - Global Active Rays Count",Environment::logVerbosityLevel() > 0);
                tic();
                // some useful data structures and constants for the main loop
                const int TAG_REQ = 100;
                const int TAG_RES = 101;
                const size_t BATCH_SIZE = 8192;//256;
                const int CHUNK_SIZE = 1024;//8192;//4096;//64;
                const int BATCH_THRESHOLD = 8192;//1024;

                MPI_Comm raw_comm = (MPI_Comm)this->worldComm();
                MPI_Request term_req = MPI_REQUEST_NULL;
                bool checking_term = false;
                int snapshot_my_active_rays = 0;
                int temp_global_active_rays = global_active_rays;
                int loop_counter = 0;

                struct NetworkReq {
                    int ray_id;
                    ray_type ray;
                };
                struct NetworkRes {
                    int ray_id;
                    rayintersection_result_type result;
                };
                struct PendingSend {
                    std::vector<NetworkReq> req_data;
                    std::vector<NetworkRes> res_data;
                    MPI_Request mpi_req;
                };
                std::list<PendingSend> active_sends;

                struct PendingRecv {
                    std::vector<NetworkReq> buffer;
                    MPI_Request req;
                    int source_rank;
                };
                std::list<PendingRecv> active_recvs;


                struct RemoteTask {
                    int source_rank;
                    NetworkReq req;
                };
                std::queue<RemoteTask> pending_remote_compute;
                std::queue<int> pending_local_compute;
                std::unordered_map<int, std::vector<NetworkRes>> pending_outgoing_resps;

                // map rank to requests to send to that rank
                std::map<int, std::vector<NetworkReq>> out_buffers;
                // map rank to requests received from that rank but not yet processed
                std::map<int, std::vector<NetworkReq>> pending_incoming_reqs;
                // Reception buffer for responses
                std::vector<NetworkRes> incoming_resps;

                double total_network_time = 0.0;
                double total_compute_time = 0.0;
                double total_sleep_time = 0.0;

                // Main loop
                while (global_active_rays > 0)
                {
                    // if ( this->worldComm().isMasterRank() )
                    //   std::cout << "start loop iteration, global_active_rays = " << global_active_rays << " counter:"<< loop_counter << std::endl;
                    // tic();
                    // On part du principe qu'on ne va rien faire ce tour-ci
                    bool did_work = false;
                    loop_counter++;
                    bool received_something = false;
                    pending_incoming_reqs.clear(); // The worklist is cleared at the start of the round

                    // ==========================================================
                    // PHASE 1: PURE NETWORK READING & PROCESSING
                    // ==========================================================
                    int flag; MPI_Status status;
                    // Get all incoming requests and put them in pending_incoming_reqs (one batch per sender)

                    double start_network = MPI_Wtime();
                    // 1.A - Get new incoming requests and start Irecv for them
                    MPI_Iprobe(MPI_ANY_SOURCE, TAG_REQ, raw_comm, &flag, &status);
                    while (flag)
                    {
                        did_work = true;
                        received_something = true;
                        int count_bytes;
                        MPI_Get_count(&status, MPI_BYTE, &count_bytes);
                        int num_items = count_bytes / sizeof(NetworkReq);

                        PendingRecv pr;
                        pr.buffer.resize(num_items);
                        pr.source_rank = status.MPI_SOURCE;
                        MPI_Irecv(pr.buffer.data(), count_bytes, MPI_BYTE, pr.source_rank, TAG_REQ, raw_comm, &pr.req);
                        active_recvs.push_back(std::move(pr));

                        MPI_Iprobe(MPI_ANY_SOURCE, TAG_REQ, raw_comm, &flag, &status);
                    }

                    // 1.B - Process completed Irecv requests and put them in pending_remote_compute
                    for (auto it = active_recvs.begin(); it != active_recvs.end(); )
                    {
                        int is_done = 0;
                        MPI_Test(&it->req, &is_done, MPI_STATUS_IGNORE); // Est-ce que le Irecv a fini ?

                        if (is_done)
                        {
                            did_work = true;
                            for (auto const& req : it->buffer)
                            {
                                pending_remote_compute.push({it->source_rank, req});
                            }
                            // The work is done; we're removing the request from the waiting list
                            it = active_recvs.erase(it);
                        }
                        else
                        {
                            // It's not finished downloading yet, so let's move on to the next one
                            ++it;
                        }
                    }

                    double end_network = MPI_Wtime();
                    total_network_time += (end_network - start_network);



                    // get all incoming responses and directly process them (update ray states, final_results, and ready_rays)
                    MPI_Iprobe(MPI_ANY_SOURCE, TAG_RES, raw_comm, &flag, &status);
                    while (flag)
                    {
                        received_something = true;
                        // Get the exact size of the message in bytes
                        int count_bytes = 0;
                        MPI_Get_count(&status, MPI_BYTE, &count_bytes);
                        int num_items = count_bytes / sizeof(NetworkRes);

                        // Prepare the buffer
                        if ( incoming_resps.size() < num_items )
                            incoming_resps.resize(num_items);

                        // Recv message
                        MPI_Recv(incoming_resps.data(), count_bytes, MPI_BYTE, status.MPI_SOURCE, TAG_RES, raw_comm, MPI_STATUS_IGNORE);

                        // treat responses one by one and update states, final_results and ready_rays accordingly
                        for ( std::size_t k = 0; k < num_items; ++k ) {
                            auto const& res = incoming_resps[k];
                            int ray_id = res.ray_id;
                            states[ray_id].is_waiting = false;
                            states[ray_id].current_idx++;

                            bool hasIntersection = res.result.processId() != invalid_v<rank_type>;
                            // store the result if it's an intersection
                            if ( hasIntersection )
                            {
                                if constexpr ( isAnyHint )
                                {
                                    final_results[ray_id] = { res.result };
                                }
                                else if constexpr ( closestOnly )
                                {
                                    if (final_results[ray_id].empty() || res.result.distance() < final_results[ray_id].front().distance())
                                        final_results[ray_id] = { res.result };
                                }
                                else
                                {
                                    final_results[ray_id].push_back( res.result );
                                }
                            }

                            // Put it back in the job queue or mark as completed
                            if constexpr ( isAnyHint )
                            {
                                if ( hasIntersection || states[ray_id].current_idx >= states[ray_id].targets.size() )
                                {
                                    states[ray_id].active = false;
                                    my_active_rays--;
                                }
                                else
                                    ready_rays.push(ray_id);
                            }
                            else
                            {
                                if ( states[ray_id].current_idx < states[ray_id].targets.size() ) {
                                    ready_rays.push(ray_id);
                                }
                                else {
                                    states[ray_id].active = false;
                                    my_active_rays--;
                                }
                            }
                        }
                        MPI_Iprobe(MPI_ANY_SOURCE, TAG_RES, raw_comm, &flag, &status);
                    }

                    // ==========================================================
                    // PHASE 2: LOCAL ROUTING AND BUFFER FILLING
                    // ==========================================================
                    // We unstack ready_rays.
                    // - If it's for another process -> we put it in out_buffers.
                    // - If it's local -> store it in a “local work to do” list
                    //   to compute it WHILE the network is running.
                    while (!ready_rays.empty())
                    {
                        int ray_id = ready_rays.front();
                        ready_rays.pop();
                        int target_rank = states[ray_id].targets[states[ray_id].current_idx].rank;

                        if (target_rank == this->worldComm().rank())
                            pending_local_compute.push(ray_id);
                        else {
                            out_buffers[target_rank].push_back({ ray_id, localRays[ray_id] });
                            states[ray_id].is_waiting = true;
                        }
                    }


                    // ==========================================================
                    // PHASE 3 : SEND INTO MPI NETWORK
                    // ==========================================================
                    for (auto& [target_rank, batch] : out_buffers)
                    {
                        if ( batch.empty() )
                            continue;
                        if ( batch.size() >= BATCH_THRESHOLD/*BATCH_SIZE*/ || !received_something ) {
                            active_sends.push_back({std::move(batch), {}, MPI_REQUEST_NULL});
                            MPI_Isend(active_sends.back().req_data.data(),
                                      active_sends.back().req_data.size() * sizeof(NetworkReq),
                                      MPI_BYTE,
                                      target_rank, TAG_REQ, raw_comm,
                                      &active_sends.back().mpi_req);

                            batch.clear();
                        }
                    }

                    active_sends.remove_if([](PendingSend& ps) {
                                               int is_done = 0;
                                               MPI_Test(&ps.mpi_req, &is_done, MPI_STATUS_IGNORE);
                                               return is_done;
                                           });

                    // ==========================================================
                    // PHASE 4 : HEAVY CPU COMPUTATION (While the network is active)
                    // ==========================================================
                    int chunk_count = 0;
                    double start_compute = MPI_Wtime();
                    if ( !pending_remote_compute.empty() || !pending_local_compute.empty() )
                        did_work = true;

                    // 1. Requests from other processes are processed first
                    while (!pending_remote_compute.empty() && chunk_count < CHUNK_SIZE)
                    {
                        auto task = pending_remote_compute.front();
                        pending_remote_compute.pop();

                        auto local_res = this->intersectSequentialImpl<useRobustTraversal,isAnyHint>(task.req.ray, tolerance);
                        if (!local_res.empty())
                            pending_outgoing_resps[task.source_rank].push_back({task.req.ray_id, local_res.front()});
                        else
                            pending_outgoing_resps[task.source_rank].push_back({task.req.ray_id, rayintersection_result_type(-1, -1, std::numeric_limits<double>::max())});

                        chunk_count++;
                    }

                    // 2. If we haven't reached the CHUNK limit, we calculate our own radii
                    while (!pending_local_compute.empty() && chunk_count < CHUNK_SIZE)
                    {
                        int ray_id = pending_local_compute.front();
                        pending_local_compute.pop();

                        // Heavy computation: intersect the ray with the local BVH
                        auto local_res = this->intersectSequentialImpl<useRobustTraversal,isAnyHint>( localRays[ray_id], tolerance );
                        states[ray_id].current_idx++;

                        bool hasIntersection = !local_res.empty();

                        // store the result if it's an intersection
                        if ( hasIntersection )
                        {
                            if constexpr ( isAnyHint )
                            {
                                final_results[ray_id] = { local_res.front() };
                            }
                            else if constexpr ( closestOnly )
                            {
                                // we guess that local_res is sorted by distance
                                if (final_results[ray_id].empty() || local_res.front().distance() < final_results[ray_id].front().distance())
                                    final_results[ray_id] = { local_res.front() };
                            }
                            else
                            {
                                final_results[ray_id].insert(final_results[ray_id].end(), local_res.begin(), local_res.end());
                            }
                        }

                        // Re-queue or terminate
                        if constexpr ( isAnyHint )
                        {
                            if ( hasIntersection || states[ray_id].current_idx >= states[ray_id].targets.size() )
                            {
                                states[ray_id].active = false;
                                my_active_rays--;
                            }
                            else
                            {
                                ready_rays.push(ray_id);
                            }
                        }
                        else
                        {
                            auto& state = states[ray_id];
                            if (state.current_idx < state.targets.size())
                                ready_rays.push(ray_id);
                            else {
                                state.active = false;
                                my_active_rays--;
                            }
                        }

                        chunk_count++;
                    }

                    double end_compute = MPI_Wtime();
                    total_compute_time += (end_compute - start_compute);

                    if (!did_work)
                    {
                        // If no messages have been received, no buffer has finished downloading,
                        // and there are no ray intersection to calculate... We put the processor to sleep for 100 microseconds.
                        double start_sleep = MPI_Wtime();
                        std::this_thread::sleep_for(std::chrono::microseconds(100));
                        double end_sleep = MPI_Wtime();
                        total_sleep_time += (end_sleep - start_sleep);
                    }

                    // send responses to other processes if we have enough or if we have nothing else to do
                    for (auto& [target_rank, resps] : pending_outgoing_resps)
                    {
                        if (resps.size() >= BATCH_THRESHOLD || (pending_remote_compute.empty() && pending_local_compute.empty()))
                        {
                            if (!resps.empty()) {
                                active_sends.push_back({{}, std::move(resps), MPI_REQUEST_NULL});
                                MPI_Isend(active_sends.back().res_data.data(),
                                          active_sends.back().res_data.size() * sizeof(NetworkRes),
                                          MPI_BYTE, target_rank, TAG_RES, raw_comm, &active_sends.back().mpi_req);
                                // resps has been “moved”; we no longer need to empty it manually
                            }
                        }
                    }

                    // ==========================================================
                    // PHASE 5 : Asynchronous Termination Detection
                    // ==========================================================
                    if (!checking_term)
                    {
                        if (loop_counter % 64 == 0)
                        {
                            snapshot_my_active_rays = my_active_rays;
                            MPI_Iallreduce(&snapshot_my_active_rays, &temp_global_active_rays, 1, MPI_INT, MPI_SUM, raw_comm, &term_req);
                            checking_term = true;
                        }
                    }
                    else
                    {
                        int reduce_done = 0;
                        MPI_Test(&term_req, &reduce_done, MPI_STATUS_IGNORE);
                        if (reduce_done)
                        {
                            // if ( this->worldComm().isMasterRank() )
                            //   std::cout << "["<< this->worldComm().rank() << "] reduction done, global_active_rays = " << temp_global_active_rays << " counter:"<< loop_counter << std::endl;
                            checking_term = false;
                            global_active_rays = temp_global_active_rays;
                        }
                    }




                } // loop while


                if ( Environment::logVerbosityLevel() > 0 )
                    std::cout << "Rang " << this->worldComm().rank() << " - Boucle principale terminée après " << loop_counter << " itérations."
                              << " Temps de calcul pur : " << total_compute_time << " secondes, "
                              << "Temps gestion réseau : " << total_network_time << " secondes, "
                              << "Temps de sommeil : " << total_sleep_time << " secondes."
                              << std::endl;
                // Final cleaning
                for (auto& ps : active_sends)
                    MPI_Wait(&ps.mpi_req, MPI_STATUS_IGNORE);

                toc("BVH Intersect - Main Loop",Environment::logVerbosityLevel() > 0);
                return final_results;
            }
#else
            if constexpr( std::is_same_v<BVHRaysDistributed<nRealDim>,napp_ray_type> ) // case rays distributed on process
            {
                // WARNING: this algo is not good (all_gather of rays then all run bvh), just a quick version for test
                auto const& localRays = ray.rays();
                std::vector<int> resLocalSize( this->worldComm().size() );
                mpi::all_gather( this->worldComm(), (int)localRays.size(), resLocalSize );

                std::vector<ray_type> raysGathered;
                if ( this->worldComm().isMasterRank() )
                {
                    int gatherRaySize = std::accumulate( resLocalSize.begin(), resLocalSize.end(), 0 );
                    raysGathered.resize( gatherRaySize );
                }
                mpi::gatherv( this->worldComm(), localRays, raysGathered.data(), resLocalSize, this->worldComm().masterRank() );
                mpi::broadcast( this->worldComm(), raysGathered, this->worldComm().masterRank() );

                auto intersectGlobal = this->intersect(_ray=raysGathered,_robust=useRobustTraversal,_context=ctx,_parallel=true,_tolerance=tolerance);

                std::vector<std::vector<rayintersection_result_type>> res;
                res.resize( ray.numberOfLocalRay() );
                std::size_t startRayIndexInThisProcess = 0;
                for ( int p=0;p<this->worldComm().rank();++p )
                    startRayIndexInThisProcess += resLocalSize[p];
                std::copy_n(intersectGlobal.cbegin()+startRayIndexInThisProcess, localRays.size(), res.begin());
                return res;
            }
#endif
            else if constexpr ( is_iterable_v<std::decay_t<decltype(ray)>> ) // case rays container are identical all on process (TODO: internal case)
            {
                std::vector<std::vector<rayintersection_result_type>> resSeq;
                resSeq.reserve( ray.size() );
                for ( auto const& currentRay : ray )
                {
                    auto currentResSeq = this->intersectSequential( currentRay, tolerance, useRobustTraversal );
                    if constexpr ( closestOnly )
                        if ( currentResSeq.size() > 1 )
                            currentResSeq.resize(1);
                    resSeq.push_back( std::move( currentResSeq ) );
                }
                if ( !parallel )
                    return resSeq;

                mpi::all_reduce( this->worldComm(), mpi::inplace( resSeq ), [](auto const& x, auto const& y) -> std::vector<std::vector<rayintersection_result_type>> {
                        std::size_t retSize = x.size();
                        std::vector<std::vector<rayintersection_result_type>> ret;
                        DCHECK( x.size() == y.size() ) << "not same size x:" << x.size() << " y:" << y.size();
                        ret.reserve( retSize );
                        for ( int k = 0; k < retSize ; ++k )
                        {
                            auto const& a = x[k];
                            auto const& b = y[k];
                            if ( a.empty() )
                                ret.push_back( b );
                            else if ( b.empty() )
                                ret.push_back( a );
                            else
                            {
                                // WARNING only return one intersection (closest or anyhint)
                                if ( a.front().distance() < b.front().distance() )
                                    ret.push_back( a );
                                else
                                    ret.push_back( b );
                            }
                        }
                        return ret;
                    } );
                return resSeq;
            }
            else // only one ray (all process should have the same ray if parallel=true)
            {
                auto resSeq = this->intersectSequential( ray, tolerance, useRobustTraversal );
                if ( closestOnly && resSeq.size() > 1 )
                    resSeq.resize(1);
                if ( !parallel )
                    return { resSeq };

#if 1
                mpi::all_reduce( this->worldComm(), mpi::inplace( resSeq ), [](auto const& a, auto const& b) -> std::vector<rayintersection_result_type> {
                        if ( a.empty() )
                            return b;
                        else if ( b.empty() )
                            return a;
                        else
                        {
                            // WARNING only return one intersection (closest or anyhint)
                            if ( a.front().distance() < b.front().distance() )
                                return { a.front() };
                            else
                                return { b.front() };
                        }
                    } );
                return { resSeq };
#else
                std::vector<int> resLocalSize( this->worldComm().size() );
                mpi::gather( this->worldComm(), (int)resSeq.size(), resLocalSize, this->worldComm().masterRank() );
                std::vector<rayintersection_result_type> resPar;
                if ( this->worldComm().isMasterRank() )
                {
                    int gatherOutputSize = std::accumulate( resLocalSize.begin(), resLocalSize.end(), 0 );
                    resPar.resize( gatherOutputSize );
                }
                mpi::gatherv( this->worldComm(), resSeq, resPar.data(), resLocalSize, this->worldComm().masterRank() );
                if ( this->worldComm().isMasterRank() )
                {
                    std::sort( resPar.begin(), resPar.end(), [](auto const& res0,auto const& res1){ return res0.distance() < res1.distance(); } );
                    if ( closestOnly && resPar.size() > 1 )
                        resPar.resize(1);
                }
                mpi::broadcast( this->worldComm(), resPar, this->worldComm().masterRank() );
                return resPar;
#endif
            }
        }


private:
    std::vector<TlasHit> getHitPartitions(ray_type const& ray) const override
        {
            auto intersect_bbox = [](const auto& ray, const auto& bbox) {
                                      // Initialize with infinity
                                      value_type t_enter = -std::numeric_limits<value_type>::infinity();
                                      value_type t_exit  =  std::numeric_limits<value_type>::infinity();

                                      for (int i = 0; i < 3; ++i) {
                                          auto invD = 1.0 / ray.dir[i];
                                          auto t0 = (bbox.min[i] - ray.org[i]) * invD;
                                          auto t1 = (bbox.max[i] - ray.org[i]) * invD;

                                          if (invD < 0.0) {
                                              std::swap(t0, t1);
                                          }

                                          t_enter = std::max(t_enter, t0);
                                          t_exit  = std::min(t_exit, t1);

                                          // If the ray completely misses the bounding box
                                          if (t_exit < t_enter) {
                                              return std::make_pair(1.0, -1.0); // Miss
                                          }
                                      }

                                      // Check against the ray's bounds (tmin / tmax)
                                      // If the box is completely behind the ray, or beyond the maximum distance
                                      if (t_exit < ray.tmin || t_enter > ray.tmax)
                                          return std::make_pair(1.0, -1.0); // Miss

                                      // Returns the actual distance of the input (which will indeed be negative if the origin is within it)
                                      return std::make_pair(t_enter, t_exit);
                                  };

            std::vector<TlasHit> hit_ranks;
            if ( !M_bvhTlas )
                return hit_ranks;

            auto rayBackend = bvh::v2::Ray<value_type, nRealDim>{
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.origin()[i]; }),
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.dir()[i]; }),
                ray.distanceMin(), ray.distanceMax()
            };

            bvh::v2::SmallStack<typename backend_bvh_type::Index, 64> stack;
            M_bvhTlas->template intersect<false, /*true*/false>( rayBackend, M_bvhTlas->get_root().index, stack,
                                                                 [this,&hit_ranks,&rayBackend,&intersect_bbox] (std::size_t begin, std::size_t end) {
                                                                     for (std::size_t i = begin; i < end; ++i) {
                                                                         auto const& [rank, bbox] = M_primitiveInfoTlas[i];
                                                                         auto hit = intersect_bbox(rayBackend, bbox);
                                                                         // hit.first = t_enter, hit.second = t_exit
                                                                         // If the test is valid (Hit), add to the result
                                                                         if (hit.first <= hit.second) {
                                                                             hit_ranks.push_back({rank, hit.first});
                                                                         }
                                                                     }
                                                                     return false; // We're not stopping—we want all the intersections
                                                                 }
                                                                 );
            //std::cout << "Rang " << this->worldComm().rank() << " - getHitPartitions: " << hit_ranks.size() << " hits." << std::endl;
            std::sort(hit_ranks.begin(), hit_ranks.end(), [](const TlasHit& a, const TlasHit& b) {
                                                              return a.distance < b.distance; // Sort by entry distance
                                                          });
            return hit_ranks;
        }
    std::vector<rayintersection_result_type> intersectSequential( ray_type const& ray, value_type tolerance, bool useRobustTraversal = true ) override
        {
            auto rayBackend = bvh::v2::Ray<value_type,nRealDim>{
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.origin()[i]; }),
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.dir()[i]; }),
                ray.distanceMin(), ray.distanceMax()
            };
            if ( useRobustTraversal )
                return this->intersectSequentialImpl<true,false>( rayBackend, tolerance );
            else
                return this->intersectSequentialImpl<false,false>( rayBackend, tolerance );
        };

    template <bool UseRobustTraversal,bool IsAnyHit>
    std::vector<rayintersection_result_type> intersectSequentialImpl( ray_type const& ray, value_type tolerance )
        {
            auto rayBackend = bvh::v2::Ray<value_type,nRealDim>{
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.origin()[i]; }),
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.dir()[i]; }),
                ray.distanceMin(), ray.distanceMax()
            };
            return this->intersectSequentialImpl<UseRobustTraversal,IsAnyHit>( rayBackend, tolerance );
        }

    template <bool UseRobustTraversal,bool IsAnyHit>
    std::vector<rayintersection_result_type> intersectSequentialImpl( bvh::v2::Ray<value_type,nRealDim> & rayBackend, value_type tolerance )
        {
            if (  !M_bvh )
                return {};
            static constexpr size_t stack_size = 64;
            static constexpr bool should_permute = true;
            //static constexpr bool isAnyHit = false;
            // Traverse the BVH and get the u, v coordinates of the closest intersection.
            bvh::v2::SmallStack<typename backend_bvh_type::Index, stack_size> stack;
            std::vector<rayintersection_result_type> res;
            M_bvh->template intersect<IsAnyHit, UseRobustTraversal>( rayBackend, M_bvh->get_root().index, stack,
                                                                     [this,&res,&rayBackend,tolerance] (std::size_t begin, std::size_t end) {
                                                                         std::size_t previousResultSize = res.size();
                                                                         for (std::size_t i = begin; i < end; ++i)
                                                                         {
                                                                             std::size_t j = should_permute ? i : M_bvh->prim_ids[i];
                                                                             if constexpr ( nRealDim == 2 )
                                                                             {
                                                                                 CHECK( false ) << "TODO";
                                                                             }
                                                                             else if constexpr ( nRealDim == 3 )
                                                                             {
                                                                                 // NOTE: we apply minus with tolerance because positive tolerance means
                                                                                 // that we can accept intersection outside of triangle at distance given by tolerance
                                                                                 if ( auto hit = M_precomputeTriangle[j].intersect( rayBackend, -tolerance ) )
                                                                                 {
                                                                                     //std::tie(u, v) = *hit;
                                                                                     res.push_back( rayintersection_result_type(this->worldComm().rank(), M_bvh->prim_ids[i], rayBackend.tmax) );
                                                                                     res.back().setCoordinates( this->barycentricToCartesianCoordinates( M_precomputeTriangle[j].convert_to_tri(),
                                                                                                                                                         hit->first, hit->second ) );
                                                                                     if constexpr ( IsAnyHit )
                                                                                         return true;
                                                                                 }
                                                                             }
                                                                         }
                                                                         return res.size() > previousResultSize;
                                                                     });
            //! sort all intersection from the distance (closer to far)
            std::sort( res.begin(), res.end(), [](auto const& res0,auto const& res1){ return res0.distance() < res1.distance(); } );
            return res;
        }

    template <typename TriType>
    vector_realdim_type barycentricToCartesianCoordinates( TriType const& tri, double u, double v ) const {
        auto const& pt0 = tri.p1;
        auto const& pt1 = tri.p2;
        auto const& pt2 = tri.p0;
        return vector_realdim_type{{
                u*pt0[0]+v*pt1[0]+(1-u-v)*pt2[0],
                u*pt0[1]+v*pt1[1]+(1-u-v)*pt2[1],
                u*pt0[2]+v*pt1[2]+(1-u-v)*pt2[2],
            }};
    }

private:
    std::unique_ptr<backend_bvh_type> M_bvh;
    std::vector<backend_precompute_triangle_type> M_precomputeTriangle;
    std::unique_ptr<backend_bvh_type> M_bvhTlas; // top level acceleration structure for parallel case
    std::vector<std::tuple<rank_type,backend_bbox_type> > M_primitiveInfoTlas;
};


//! @brief in house implementation of BVH tool
template <typename MeshEntityType>
class BVH_InHouse : public BVH<MeshEntityType>
{
    using super_type = BVH<MeshEntityType>;
    using self_type = BVH_InHouse<MeshEntityType>;
    using mesh_entity_type = typename super_type::mesh_entity_type;
    using vector_realdim_type = typename super_type::vector_realdim_type;
    static constexpr uint16_type nRealDim = super_type::nRealDim;
    using value_type = typename super_type::value_type;
    using primitiveinfo_type = typename super_type::primitiveinfo_type;
public:
    using ray_type = typename super_type::ray_type;
    using rayintersection_result_type = typename super_type::rayintersection_result_type;

    class BVHNode
    {
        friend class BVH_InHouse<mesh_entity_type>;
    public:
        BVHNode() = default;

        //! return the parent of this node
        BVHNode * parent() const { return M_parent; }

        vector_realdim_type const& boundMin() const noexcept { return M_bounds_min; }
        vector_realdim_type const& boundMax() const noexcept { return M_bounds_max; }
        vector_realdim_type centroid() const { return 0.5*(M_bounds_min + M_bounds_max); }
        int splitAxis() const noexcept { return M_splitAxis; }
        int nPrimitives() const noexcept { return M_nPrimitives; }
        int firstPrimOffset() const noexcept { return M_firstPrimOffset; }

        BVHNode* child( int k ) const { return M_children[k].get(); }

        bool isLeaf() const { return !M_children[0] && !M_children[1]; }


        BVHNode * nearChild( ray_type const& ray ) const
            {
                if( ray.dir()(this->splitAxis()) > 0 )
                    return this->child(0);
                else
                    return this->child(1);
            }

        BVH_InHouse::BVHNode * siblingNode() const
            {
                if ( !M_parent )
                    return nullptr;
                return M_parent->child( this == M_parent->child(0)? 1 : 0 );
            }

        bool checkIntersection(ray_type const& rayon)
            {
                double tmin = 0.0;
                double tmax = FLT_MAX;

                for(int i=0; i<nRealDim; i++)
                {
                    double ratio = 1.0/(rayon.dir()[i]+2*FLT_MIN);
                    double t1 = (M_bounds_min[i]-rayon.origin()[i]) * ratio;
                    double t2 = (M_bounds_max[i]-rayon.origin()[i]) * ratio;
                    if (t1 > t2)
                    {
                        double tTemp = t1;
                        t1 = t2;
                        t2 = tTemp;
                    }
                    if ( t1 > tmin)
                        tmin = t1;
                    if (t2 > tmax)
                        tmax = t2;
                    if (tmin > tmax)
                        return false;
                }

                return true;
            }

        std::pair<bool,double> checkIntersectionWithSegment( ray_type const& ray, std::vector<primitiveinfo_type> const& primitiveInfo ) const
            {
                auto const& meshElt = primitiveInfo[this->firstPrimOffset()].meshEntity();
                auto p1 = Eigen::Map<const Eigen::Matrix<double,nRealDim,1>>( meshElt.point(0).node().data().begin() );
                auto p2 = Eigen::Map<const Eigen::Matrix<double,nRealDim,1>>( meshElt.point(1).node().data().begin() );

                auto const& origin = ray.origin();
                auto const& direction = ray.dir();

                vector_realdim_type v1 = origin - p1;
                vector_realdim_type v2 = p2 - p1;
                vector_realdim_type v3{ -direction[1], direction[0] };

                double dot = v2.dot(v3);
                if (math::abs(dot) < 1e-6)
                    return std::make_pair(false,0);

                double t1 = (v2[0]*v1[1]-v2[1]*v1[0])/ dot;
                double t2 = v1.dot(v3) / dot;

                if (t1 > 2*FLT_MIN && (t2 >= 0.0 && t2 <= 1.0))
                {
#if 0
                    vector_realdim_type w_{
                        origin[0] + direction[0]*t1,
                        origin[1] + direction[1]*t1; };
#endif
                    return std::make_pair(true,t1);
                }
                return std::make_pair(false,t1);
            }

        // Verify if the ray intersects the element
        std::pair<bool,double> checkIntersectionWithTriangle( ray_type const& ray, std::vector<primitiveinfo_type> const& primitiveInfo ) const
            {
                DCHECK( this->isLeaf() ) << "should be a leaf: ";

                auto const& meshElt = primitiveInfo[this->firstPrimOffset()].meshEntity();
                auto p1 = Eigen::Map<const Eigen::Matrix<double,nRealDim,1>>( meshElt.point(0).node().data().begin() );
                auto p2 = Eigen::Map<const Eigen::Matrix<double,nRealDim,1>>( meshElt.point(1).node().data().begin() );
                auto p3 = Eigen::Map<const Eigen::Matrix<double,nRealDim,1>>( meshElt.point(2).node().data().begin() );

                auto const& origin = ray.origin();
                auto const& direction = ray.dir();

                // // normal vector
                auto n1 = (p2-p1).cross(p3-p1);
                n1 = n1/n1.norm();
                double n_dot_dir = direction.dot(n1);
                // Ray is parallel to the triangle's plane
                if (math::abs(n_dot_dir)<1e-6)
                {
                    return std::make_pair(false,0);
                }
                double d = -p1.dot(n1);
                double t_line = -(origin.dot(n1)+d)/n_dot_dir;
                if( t_line <= 1e-10) // intersection not in the same direction as the ray
                    return std::make_pair(false,0);
                // intersection point
                auto w = origin + direction* t_line;

                Eigen::Matrix<double,3,3> m;
                m.col(0) = p2-p1;
                m.col(1) = p3-p1;
                m.col(2) = n1;
                auto w_ = m.inverse()*(w-p1);

                return std::make_pair((w_(0)> 2*FLT_MIN ) && (w_(1)>0) && (w_(0) +  w_(1)<1),t_line);
            }

        std::pair<bool,double> checkLeafIntersection(ray_type const& rayon, std::vector<primitiveinfo_type> const& primitiveInfo)
            {
                if constexpr ( nRealDim == 2 )
                    return checkIntersectionWithSegment( rayon, primitiveInfo );
                else if constexpr ( nRealDim == 3 )
                    return checkIntersectionWithTriangle( rayon, primitiveInfo );
            }

    private:
        BVHNode* setChild( uint16_type k, std::unique_ptr<BVHNode> && childNode )
            {
                if ( childNode->M_parent ) { /*TODO remove child in this parent*/ }

                childNode->M_parent = this;
                M_children[k] = std::move( childNode );
                return M_children[k].get();
            }

        void updateForUse( int firstPrimOffset, int nPrimitives, int splitAxis, Eigen::VectorXd const& bounds_min, Eigen::VectorXd const& bounds_max )
            {
                M_firstPrimOffset = firstPrimOffset;
                M_nPrimitives = nPrimitives;
                M_bounds_min = bounds_min;
                M_bounds_max = bounds_max;
                M_splitAxis = splitAxis;
            }

    private:
        std::array<std::unique_ptr<BVHNode>,2> M_children;
        BVHNode *M_parent = nullptr;
        int M_splitAxis = 0, M_nPrimitives = 0, M_firstPrimOffset = 0;
        vector_realdim_type M_bounds_min, M_bounds_max;
    };

    BVH_InHouse( worldcomm_ptr_t worldComm ) : super_type( BVHEnum::Quality::High, worldComm ) {}

    template <typename RangeType>
    void
    updateForUse( RangeType const& range )
        {
            // up primitiveinfos
            super_type::updateForUse( range );
            // build BVH tree
            this->buildTree();
        }

private:
    // Verify if the ray intersects the whole bounding structure
    // Returns the integer corresponding to the intersected element
    // If no element is intersected, return -1
    std::vector<rayintersection_result_type> intersectSequential( ray_type const& rayon, value_type tolerance, bool useRobustTraversal = true ) override
        {
            M_intersected_leaf = {};
            M_lengths = {};
            if ( !M_rootNode )
                buildTree();
            if ( this->M_primitiveInfo.empty() )
                return {};

            std::vector<rayintersection_result_type> res;
            if ( M_rootNode->checkIntersection(rayon) )
            {
                traverse_stackless( M_rootNode.get(), rayon );
            }
            if ( !M_intersected_leaf.empty() )
            {
                int argmin_lengths = std::distance(M_lengths.begin(), std::min_element(M_lengths.begin(), M_lengths.end()));
                res.push_back( rayintersection_result_type(this->worldComm().rank(), M_intersected_leaf[argmin_lengths], M_lengths[argmin_lengths] ) );
            }
            return res;
        }

    void buildTree()
        {
            if ( M_rootNode )
                return;

            M_rootNode = std::make_unique<BVHNode>();

            std::stack<std::tuple<BVHNode*,int,int,int>> stack;
            stack.push( std::make_tuple(M_rootNode.get(),0,0,this->M_primitiveInfo.size()) );
            // TODO case only one 1 element
            while ( !stack.empty() )
            {
                auto [currentNode,cut_dimension,start_index_primitive,end_index_primitive] = stack.top();
                stack.pop();

                int nPrimitives = end_index_primitive - start_index_primitive;
                auto [bound_min_node,bound_max_node] = nPrimitives > 0 ? this->bounds( start_index_primitive,end_index_primitive ) : std::make_tuple( vector_realdim_type{}, vector_realdim_type{});

                if ( nPrimitives <= 1 )
                {
                    // Create a leaf, since there is only one primitive in the list
                    int firstPrimOffset = M_orderedPrims.size();
                    for (int i = start_index_primitive; i < end_index_primitive; ++i)
                    {
                        int primNum = this->M_primitiveInfo[i].meshEntity().id();
                        M_orderedPrims.push_back(primNum);
                    }
                    currentNode->updateForUse( firstPrimOffset, nPrimitives, -1, bound_min_node, bound_max_node );
                }
                else
                {
                    CHECK( start_index_primitive >=0 && end_index_primitive <= this->M_primitiveInfo.size() ) << start_index_primitive << " " << end_index_primitive;
                    auto mid = (start_index_primitive + end_index_primitive) / 2;
                    std::nth_element(&this->M_primitiveInfo[start_index_primitive], &this->M_primitiveInfo[mid],
                                     &this->M_primitiveInfo[end_index_primitive-1]+1,
                                     [cut_dimension=cut_dimension](primitiveinfo_type const&a, primitiveinfo_type const& b) {
                                         return a.centroid()[cut_dimension] < b.centroid()[cut_dimension];
                                     });

                    int next_cut_dimension=(cut_dimension+1)%nRealDim;
                    auto childNode0 = currentNode->setChild( 0, std::make_unique<BVHNode>() );
                    stack.push( std::make_tuple(childNode0, next_cut_dimension, start_index_primitive, mid) );
                    auto childNode1 = currentNode->setChild( 1, std::make_unique<BVHNode>() );
                    stack.push( std::make_tuple( childNode1, next_cut_dimension, mid, end_index_primitive ) );

                    currentNode->updateForUse( -1, nPrimitives, next_cut_dimension, bound_min_node, bound_max_node );
                }
            }
        }


    std::tuple<vector_realdim_type,vector_realdim_type> bounds( int start_index_primitive, int end_index_primitive ) const
        {
            if ( start_index_primitive >= end_index_primitive )
                throw std::logic_error("Error in BVHNode : compute bounds with no elemnent");

            //vector_realdim_type newBoundsMin, newBoundsMax;
            vector_realdim_type newBoundsMin = this->M_primitiveInfo[start_index_primitive].boundMin();
            vector_realdim_type newBoundsMax = this->M_primitiveInfo[start_index_primitive].boundMax();
            for (int i = start_index_primitive+1; i < end_index_primitive; ++i)
            {
                auto const& primitiveInfo = this->M_primitiveInfo[i];
                for ( uint8_type d=0;d<vector_realdim_type::SizeAtCompileTime;++d )
                {
                    newBoundsMin[d] = std::min( newBoundsMin[d], primitiveInfo.boundMin()[d] );
                    newBoundsMax[d] = std::max( newBoundsMax[d], primitiveInfo.boundMax()[d] );
                }
            }
            return std::make_tuple( std::move(newBoundsMin), std::move(newBoundsMax) );
        }


    void traverse_stackless( BVH_InHouse::BVHNode * tree, ray_type const& rayon )
        {
            auto current_node = M_rootNode->nearChild(rayon);
            if ( !current_node ) // case where root is leaf
            {
                auto [has_intersected_leaf,distance] = M_rootNode->checkLeafIntersection(rayon,this->M_primitiveInfo);
                if ( has_intersected_leaf )
                {
                    M_intersected_leaf.push_back(M_rootNode->firstPrimOffset());
                    M_lengths.push_back(distance);
                }
                return;
            }
            char state = 'P'; // the current node is being traversed from its Parent ('P')

            while(true)
            {
                switch (state)
                {
                case 'C': // the node is being traversed from its child

                    if ( current_node == M_rootNode.get() ) return;

                    if ( current_node == current_node->parent()->nearChild( rayon ) )
                    {
                        current_node = current_node->siblingNode();
                        state = 'S'; // the current node has been accessed from its sibling
                    }
                    else
                    {
                        current_node = current_node->parent();
                        state = 'C'; // the current node has been accessed from its sibling
                    }
                    break;

                case 'S': // the node is being traversed from its Sibling ('S')

                    if ( current_node->checkIntersection( rayon ) == false ) // back to parent
                    {
                        current_node = current_node->parent();
                        state = 'C'; // the current node is being accessed from its child
                    }
                    else if ( current_node->isLeaf() )
                    {
                        auto [has_intersected_leaf,distance] = current_node->checkLeafIntersection(rayon,this->M_primitiveInfo);
                        if ( has_intersected_leaf )
                        {
                            //if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), this->M_primitiveInfo[current_node->firstPrimOffset()].meshEntity().id()) == M_intersected_leaf.end() )
                            if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), current_node->firstPrimOffset() ) == M_intersected_leaf.end() )
                            {
                                //M_intersected_leaf.push_back(this->M_primitiveInfo[current_node->firstPrimOffset()].meshEntity().id());
                                M_intersected_leaf.push_back(current_node->firstPrimOffset());
                                M_lengths.push_back(distance);
                            }
                        }
                        current_node = current_node->parent();
                        state='C'; // the current node is being accessed from its child
                    }
                    else
                    {
                        current_node = current_node->nearChild(rayon);
                        state='P';// the current node has been accessed from its parent
                    }
                    break;

                case 'P':
                    if ( current_node->checkIntersection(rayon) == false )
                    {
                        current_node=current_node->siblingNode();
                        state = 'S'; // the current node has been accessed from its sibling
                    }
                    else if ( current_node->isLeaf() )
                    {
                        auto [has_intersected_leaf,distance] = current_node->checkLeafIntersection(rayon,this->M_primitiveInfo);
                        if ( has_intersected_leaf )
                        {
                            //if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), this->M_primitiveInfo[current_node->firstPrimOffset()].meshEntity().id()) == M_intersected_leaf.end() )
                            if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), current_node->firstPrimOffset()) == M_intersected_leaf.end() )
                            {
                                //M_intersected_leaf.push_back(this->M_primitiveInfo[current_node->firstPrimOffset()].meshEntity().id());
                                M_intersected_leaf.push_back(current_node->firstPrimOffset());
                                M_lengths.push_back(distance);
                            }
                        }
                        current_node = current_node->siblingNode();
                        state = 'S'; // the current node has been accessed from its sibling
                    }
                    else
                    {
                        current_node = current_node->nearChild( rayon );
                        state = 'P'; // the current node has been accessed from its parent
                    }
                    break;

                default:

                    LOG(ERROR) << "ERROR: None of the previous cases has been traversed";

                    throw std::logic_error("Error in BVH traversal: none of the previous cases has been traversed.");

                    break;
                }
            }
        }

private:
    std::unique_ptr<BVHNode> M_rootNode;

    thread_local static inline std::vector<int> M_intersected_leaf;
    thread_local static inline std::vector<double> M_lengths;

    std::vector<int> M_orderedPrims; // order of traversed primitives for depth-first search
};


template <typename MeshEntityType>
template <typename BVH<MeshEntityType>::IntersectContext Ctx,bool useRobustTraversal,typename RayType>
std::vector<std::vector<typename BVH<MeshEntityType>::rayintersection_result_type>>
BVH<MeshEntityType>::intersectGenericImpl( RayType const& ray, value_type tolerance, bool parallel )
{
    auto bvh = dynamic_cast<BVH_ThirdParty<mesh_entity_type>*>( this );
    if (bvh)
        return bvh->template intersectFullImpl<Ctx,useRobustTraversal,RayType>( ray, tolerance, parallel );
    else
        throw std::logic_error("intersectGenericImpl should be called on BVH_ThirdParty");
    return {};
}


template<typename... Ts>
auto boundingVolumeHierarchy( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    using mesh_entity_type = std::remove_const_t<entity_range_t<std::decay_t<decltype(range)>>>;
    std::string const& kind = args.get_else(_kind, mesh_entity_type::nRealDim == 3 ? "third-party" : "in-house");
    BVHEnum::Quality quality = args.get_else(_quality, BVHEnum::Quality::High );
    worldcomm_ptr_t worldcomm = args.get_else(_worldcomm,Environment::worldCommPtr()); // TODO : use default worldcomm from range

    using bvh_type = BVH<mesh_entity_type>;
    std::unique_ptr<bvh_type> bvh;

    if ( kind == "in-house" )
    {
        using bvh_inhouse_type = BVH_InHouse<mesh_entity_type>;
        auto bvhInHouse = std::make_unique<bvh_inhouse_type>(worldcomm);
        bvhInHouse->updateForUse(range);
        bvh = std::move( bvhInHouse );
    }
    else if ( kind == "third-party" )
    {
        if constexpr ( mesh_entity_type::nRealDim != 3 )
            throw std::invalid_argument("third-party only implement with triangle in 3D");
        auto bvhThirdParty = std::make_unique<BVH_ThirdParty<mesh_entity_type>>( quality, worldcomm );
        bvhThirdParty->updateForUse(range);
        bvh = std::move( bvhThirdParty );
    }
    else
        throw std::invalid_argument(fmt::format("invalid bvh arg kind {} (should be third-party or in-house)",kind ));

    return bvh;
}

} // Feel

BOOST_IS_BITWISE_SERIALIZABLE( Feel::BVHRay<2> )
BOOST_IS_BITWISE_SERIALIZABLE( Feel::BVHRay<3> )

#endif /* FEELPP_MESH_BVH_HPP */
