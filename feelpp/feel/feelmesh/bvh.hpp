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


template <int RealDim>
class BVHRay
{
public:
    using vec_t = eigen_vector_type<RealDim>;
    BVHRay(vec_t const& orig, vec_t const& dir)
        :
        M_origin( orig ),
        M_dir( dir )
        {}
    BVHRay( BVHRay const& ) = default;
    BVHRay( BVHRay &&) = default;

    vec_t const& origin() const noexcept { return M_origin; }
    vec_t const& dir() const noexcept { return M_dir; }
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
    vec_t M_origin, M_dir;  // ray origin and dir
};

template <typename MeshEntityType>
class BVH
{
public:
    using mesh_entity_type = std::decay_t<MeshEntityType>;
    static constexpr uint16_type nDim = mesh_entity_type::nDim;
    static constexpr uint16_type nRealDim = mesh_entity_type::nRealDim;
    using vector_realdim_type = Eigen::Matrix<double,nRealDim,1>;
    using ray_type = BVHRay<nRealDim>;

    //! @brief Information on the primitive (bounding box, index, centroid)
    struct BVHPrimitiveInfo
    {
        BVHPrimitiveInfo( mesh_entity_type const& meshElt )
            :
            M_meshElement( meshElt )
            {
                auto verticesUblas = meshElt.vertices();
                auto G = em_cmatrix_col_type<double>( verticesUblas.data().begin(), nRealDim, mesh_entity_type::numVertices );
                M_bound_min = G.rowwise().minCoeff();
                M_bound_max = G.rowwise().maxCoeff();
                M_bound_min.array() -= 2*FLT_MIN;
                M_bound_max.array() += 2*FLT_MIN;
                M_centroid = ( M_bound_min + M_bound_max ) * 0.5;
            }
        BVHPrimitiveInfo( BVHPrimitiveInfo && ) = default;
        BVHPrimitiveInfo( BVHPrimitiveInfo const& ) = default;
        BVHPrimitiveInfo& operator=( BVHPrimitiveInfo && ) = default;
        BVHPrimitiveInfo& operator=( BVHPrimitiveInfo const& ) = default;

        mesh_entity_type const& meshElement() const { return M_meshElement.get(); }
        vector_realdim_type const& boundMin() const noexcept { return M_bound_min; }
        vector_realdim_type const& boundMax() const noexcept { return M_bound_max; }
        vector_realdim_type const& centroid() const noexcept { return M_centroid; }

    private:
        vector_realdim_type M_bound_min;
        vector_realdim_type M_bound_max;
        vector_realdim_type M_centroid;
        std::reference_wrapper<mesh_entity_type const> M_meshElement;
    };
    using primitiveinfo_type = BVHPrimitiveInfo;

    struct BVHRayIntersectionResult
    {
        BVHRayIntersectionResult( BVHPrimitiveInfo const& primitive, double dist )
            :
            M_primitive( primitive ),
            M_distance( dist )
            {}
        BVHRayIntersectionResult( BVHRayIntersectionResult && ) = default;
        BVHRayIntersectionResult( BVHRayIntersectionResult const& ) = default;
        BVHRayIntersectionResult& operator=( BVHRayIntersectionResult && ) = default;
        BVHRayIntersectionResult& operator=( BVHRayIntersectionResult const& ) = default;

        BVHPrimitiveInfo const& primitive() const { return M_primitive.get(); }
    private:
        std::reference_wrapper<BVHPrimitiveInfo const> M_primitive;
        double M_distance;
        //std::array<double,2> M_coordinates;
    };
    using rayintersection_result_type = BVHRayIntersectionResult;

    BVH() = default;
    BVH( BVH && ) = default;
    BVH( BVH const& ) = default;

    virtual std::vector<rayintersection_result_type> intersect( ray_type const& rayon ) = 0;
protected:
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
};

template <typename MeshEntityType>
class BVH_ThirdParty : public BVH<MeshEntityType>
{
    using super_type = BVH<MeshEntityType>;
    using mesh_entity_type = typename super_type::mesh_entity_type;
    using vector_realdim_type = typename super_type::vector_realdim_type;
    static constexpr uint16_type nRealDim = super_type::nRealDim;

    using value_type = double;
    using node_type = bvh::v2::Node<value_type, nRealDim>;
    using backend_bvh_type = bvh::v2::Bvh<node_type>;
    using backend_vector_realdim_type = bvh::v2::Vec<value_type, nRealDim>;
    using backend_precompute_triangle_type = bvh::v2::PrecomputedTri<value_type>;
public:
    using ray_type = typename super_type::ray_type;
    using rayintersection_result_type = typename super_type::rayintersection_result_type;

    BVH_ThirdParty() = default;
    BVH_ThirdParty( BVH_ThirdParty && ) = default;

    template <typename RangeType>
    void
    updateForUse( RangeType const& range )
        {
            // up primitiveinfos
            super_type::updateForUse( range );

            // init bvh backend
            using BBox = bvh::v2::BBox<value_type, nRealDim>;
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

                auto const& meshElt = primInfo.meshElement();
                auto bary = meshElt.barycenter();
                centers.push_back( backend_vector_realdim_type::generate([&bary] (std::size_t i) { return bary[i]; }) );
            }

            typename bvh::v2::DefaultBuilder<node_type>::Config config;
            config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::High;
            M_bvh = std::make_unique<backend_bvh_type>( bvh::v2::DefaultBuilder<node_type>::build(/*thread_pool,*/ bboxes, centers, config) );


            // Permuting the primitive data allows to remove indirections during traversal, which makes it faster.
            static constexpr bool should_permute = true;

            if constexpr ( nRealDim == 3 )
                             M_precomputeTriangle.resize( this->M_primitiveInfo.size() );

            for ( std::size_t i = 0; i < this->M_primitiveInfo.size(); ++i )// auto const& primInfo : this->M_primitiveInfo )
            {
                auto j = should_permute ? M_bvh->prim_ids[i] : i;
                auto const& primInfo = this->M_primitiveInfo[j];
                auto const& meshElt = primInfo.meshElement();
                if constexpr ( nRealDim == 3 )
                             {
                                 auto const& pt0 = meshElt.point(0);
                                 auto const& pt1 = meshElt.point(1);
                                 auto const& pt2 = meshElt.point(2);
                                 M_precomputeTriangle[i] = backend_precompute_triangle_type{
                                         backend_vector_realdim_type::generate([&pt0] (std::size_t i) { return pt0[i]; }),
                                         backend_vector_realdim_type::generate([&pt1] (std::size_t i) { return pt1[i]; }),
                                         backend_vector_realdim_type::generate([&pt2] (std::size_t i) { return pt2[i]; })
                                     };
                             }
            }
        }


    std::vector<rayintersection_result_type> intersect( ray_type const& ray ) override {

        auto rayBackend = bvh::v2::Ray<value_type,nRealDim>{
            backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.origin()[i]; }),
            backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.dir()[i]; }),
            1e-8 // tmin
        };

        return this->intersectImpl<true>( rayBackend );
    };
private:
    template <bool UseRobustTraversal>
    std::vector<rayintersection_result_type> intersectImpl( bvh::v2::Ray<value_type,nRealDim> & rayBackend ) {
        static constexpr size_t stack_size = 64;
        static constexpr bool should_permute = true;
        // Traverse the BVH and get the u, v coordinates of the closest intersection.
        bvh::v2::SmallStack<typename backend_bvh_type::Index, stack_size> stack;
        std::vector<rayintersection_result_type> res;
        M_bvh->template intersect<true/*false*/, UseRobustTraversal>(rayBackend, M_bvh->get_root().index, stack,
                                                     [&] (std::size_t begin, std::size_t end) {
                                                         for (std::size_t i = begin; i < end; ++i) {
                                                             std::size_t j = should_permute ? i : M_bvh->prim_ids[i];
                                                             if (auto hit = M_precomputeTriangle[j].intersect(rayBackend)) {
                                                                 //prim_id = i;
                                                                 res.push_back( rayintersection_result_type{this->M_primitiveInfo[M_bvh->prim_ids[i]], rayBackend.tmax} );
                                                                 return true;
                                                                 //std::tie(u, v) = *hit;
                                                             }
                                                         }
                                                         //return prim_id != invalid_id;
                                                         //return prim_id >= 0;
                                                         return false;
                                                   });

        return res;
    }
private:
    std::unique_ptr<backend_bvh_type> M_bvh;
    std::vector<backend_precompute_triangle_type> M_precomputeTriangle;
    std::vector<int> M_primInv;
};


template <typename MeshEntityType>
class BVHTree : public BVH<MeshEntityType>
{
    using super_type = BVH<MeshEntityType>;
    using self_type = BVHTree<MeshEntityType>;
    using mesh_entity_type = typename super_type::mesh_entity_type;
    using vector_realdim_type = typename super_type::vector_realdim_type;
    static constexpr uint16_type nRealDim = super_type::nRealDim;
    using primitiveinfo_type = typename super_type::primitiveinfo_type;
public:
    using ray_type = typename super_type::ray_type;
    using rayintersection_result_type = typename super_type::rayintersection_result_type;

    class BVHNode
    {
        friend class BVHTree<mesh_entity_type>;
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

        BVHTree::BVHNode * siblingNode() const
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
                auto const& meshElt = primitiveInfo[this->firstPrimOffset()].meshElement();
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

                auto const& meshElt = primitiveInfo[this->firstPrimOffset()].meshElement();
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

    BVHTree() = default;

    template <typename RangeType>
    void
    updateForUse( RangeType const& range )
        {
            // up primitiveinfos
            super_type::updateForUse( range );
            // build BVH tree
            this->buildTree();
        }


    // Verify if the ray intersects the whole bounding structure
    // Returns the integer corresponding to the intersected element
    // If no element is intersected, return -1
    std::vector<rayintersection_result_type> intersect( ray_type const& rayon ) override
        {
            M_intersected_leaf = {};
            M_lengths = {};
            if ( !M_rootNode )
                buildTree();

            std::vector<rayintersection_result_type> res;
            if ( M_rootNode->checkIntersection(rayon) )
            {
                traverse_stackless( M_rootNode.get(), rayon );
            }
            if ( !M_intersected_leaf.empty() )
            {
                int argmin_lengths = std::distance(M_lengths.begin(), std::min_element(M_lengths.begin(), M_lengths.end()));
                //int closer_intersection_element = M_intersected_leaf[argmin_lengths];
                res.push_back( rayintersection_result_type{this->M_primitiveInfo[ M_intersected_leaf[argmin_lengths] ], M_lengths[argmin_lengths] } );
                //return closer_intersection_element;
            }
            return res;
        }
private:

    void buildTree()
        {
            if ( M_rootNode )
                return;

            M_rootNode = std::make_unique<BVHNode>();

            std::stack<std::tuple<BVHNode*,int,int,int>> stack;
            stack.push( std::make_tuple(M_rootNode.get(),0,0,this->M_primitiveInfo.size()) );
            while ( !stack.empty() )
            {
                auto [currentNode,cut_dimension,start_index_primitive,end_index_primitive] = stack.top();
                stack.pop();

                int nPrimitives = end_index_primitive - start_index_primitive;
                DCHECK( nPrimitives>0 ) << nPrimitives;

                auto [bound_min_node,bound_max_node] = this->bounds( start_index_primitive,end_index_primitive );

                if ( nPrimitives == 1 )
                {
                    // Create a leaf, since there is only one primitive in the list
                    int firstPrimOffset = M_orderedPrims.size();
                    for (int i = start_index_primitive; i < end_index_primitive; ++i)
                    {
                        int primNum = this->M_primitiveInfo[i].meshElement().id();
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


    void traverse_stackless( BVHTree::BVHNode * tree, ray_type const& rayon )
        {
            auto current_node = M_rootNode->nearChild(rayon);
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
                            //if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), this->M_primitiveInfo[current_node->firstPrimOffset()].meshElement().id()) == M_intersected_leaf.end() )
                            if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), current_node->firstPrimOffset() ) == M_intersected_leaf.end() )
                            {
                                //M_intersected_leaf.push_back(this->M_primitiveInfo[current_node->firstPrimOffset()].meshElement().id());
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
                            //if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), this->M_primitiveInfo[current_node->firstPrimOffset()].meshElement().id()) == M_intersected_leaf.end() )
                            if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), current_node->firstPrimOffset()) == M_intersected_leaf.end() )
                            {
                                //M_intersected_leaf.push_back(this->M_primitiveInfo[current_node->firstPrimOffset()].meshElement().id());
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


template<typename... Ts>
auto boundingVolumeHierarchy( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    using mesh_entity_type = std::remove_const_t<entity_range_t<std::decay_t<decltype(range)>>>;

    std::string const& kind = args.get_else(_kind, mesh_entity_type::nRealDim == 3 ? "third-party" : "in-house");

    using bvh_type = BVH<mesh_entity_type>;
    std::unique_ptr<bvh_type> bvh;

    if ( kind == "in-house" )
    {
        using bvh_inhouse_type = BVHTree<mesh_entity_type>;
        auto bvhInHouse = std::make_unique<bvh_inhouse_type>();
        bvhInHouse->updateForUse(range);
        bvh = std::move( bvhInHouse );
    }
    else if ( kind == "third-party" )
    {
        if constexpr ( mesh_entity_type::nRealDim != 3 )
            throw std::invalid_argument("third-party only implement with triangle in 3D");
        auto bvhThirdParty = std::make_unique<BVH_ThirdParty<mesh_entity_type>>();
        bvhThirdParty->updateForUse(range);
        bvh = std::move( bvhThirdParty );
    }
    else
        throw std::invalid_argument(fmt::format("invalid bvh arg kind {} (should be third-party or in-house)",kind ));

    return bvh;
}

} // Feel
#endif /* FEELPP_MESH_BVH_HPP */
