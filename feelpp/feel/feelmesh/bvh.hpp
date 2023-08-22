// Bounding volume hierarchy

#ifndef FEELPP_MESH_BVH_HPP
#define FEELPP_MESH_BVH_HPP


#include <vector>
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feeldiscr/mesh.hpp>

namespace Feel
{
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
private:
    vec_t M_origin, M_dir;  // ray origin and dir
};
#if 0
inline constexpr float robust_const(int n)
{
    double epsil = std::numeric_limits<double>::epsilon() * 0.5;
    return (n * epsil) / (1 - n * epsil);
}
#endif
template <typename MeshEntityType>
class BVHTree
{
    using self_type = BVHTree<MeshEntityType>;
    using mesh_entity_type = std::decay_t<MeshEntityType>;
    static constexpr uint16_type nDim = mesh_entity_type::nDim;
    static constexpr uint16_type nRealDim = mesh_entity_type::nRealDim;
    // using VectorRealDim = typename  std::conditional< nRealDim==3,
    //                                                   Eigen::Vector3d,
    //                                                   Eigen::Vector2d>::type;
    using vector_realdim_type = Eigen::Matrix<double,nRealDim,1>;
    typedef typename matrix_node<double>::type matrix_node_type;

    // Information on the primitive (bounding box, index, centroid)
public:
    using ray_type = BVHRay<nRealDim>;

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
#if 0
        bool checkIntersectionWithSegment(BVHRay const& ray)
            {
                Eigen::VectorXd p1(2),p2(2),v1(2),v2(2),v3(2),w_(2);

                auto nodes =  M_mesh->element( M_primitiveInfo[this->firstPrimOffset].M_primitiveNumber).vertices();

                for ( int i = 0; i < nodes.size1(); ++i )
                {
                    p1( i ) = nodes(i,0);
                    p2( i ) = nodes(i,1);
                }
                Eigen::VectorXd origin(2),direction(2);
                direction << ray.dir[0],ray.dir[1];
                origin << ray.origin[0],ray.origin[1];
                v1 = origin - p1;
                v2 = p2 - p1;
                v3 <<-direction[1], direction[0];

                double dot = v2.dot(v3);
                if (math::abs(dot) < 1e-6)
                    return false;

                double t1 = (v2[0]*v1[1]-v2[1]*v1[0])/ dot;
                double t2 = v1.dot(v3) / dot;

                if (t1 > 2*FLT_MIN && (t2 >= 0.0 && t2 <= 1.0))
                {
                    w_[0] =  ray.origin[0] + ray.dir[0]*t1;
                    w_[1] =  ray.origin[1] + ray.dir[1]*t1;
                    return true;
                }

                return false;

            }
#endif

        // TODO VINCENT optimize here
        std::pair<bool,double> checkIntersectionWithSegment(matrix_node_type const& nodes, ray_type const& ray)
            {
                vector_realdim_type p1,p2,v1,v2,v3,w_;

                for ( int i = 0; i < nodes.size1(); ++i )
                {
                    p1( i ) = nodes(i,0);
                    p2( i ) = nodes(i,1);
                }
                vector_realdim_type origin,direction;
                direction << ray.dir()[0],ray.dir()[1];
                origin << ray.origin()[0],ray.origin()[1];
                v1 = origin - p1;
                v2 = p2 - p1;
                v3 <<-direction[1], direction[0];

                double dot = v2.dot(v3);
                if (math::abs(dot) < 1e-6)
                    return std::make_pair(false,0);

                double t1 = (v2[0]*v1[1]-v2[1]*v1[0])/ dot;
                double t2 = v1.dot(v3) / dot;

                if (t1 > 2*FLT_MIN && (t2 >= 0.0 && t2 <= 1.0))
                {
                    w_[0] =  ray.origin()[0] + ray.dir()[0]*t1;
                    w_[1] =  ray.origin()[1] + ray.dir()[1]*t1;
                    return std::make_pair(true,t1);
                }

                return std::make_pair(false,t1);

            }

        // Verify if the ray intersects the element
        std::pair<bool,double> checkIntersectionWithTriangle( ray_type const& ray, std::vector<BVHPrimitiveInfo> const& primitiveInfo ) const
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

        std::pair<bool,double> checkLeafIntersection(ray_type const& rayon, std::vector<BVHPrimitiveInfo> const& primitiveInfo)
            {
                if constexpr ( nRealDim == 2 )
                {
                    auto nodes = primitiveInfo[this->firstPrimOffset()].meshElement().vertices();
                    return checkIntersectionWithSegment( nodes, rayon );
                }
                else if constexpr ( nRealDim == 3 ) //if ( nRealDim==3)
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
            // From the mesh, build the bounding box info for each element and store it in
            // the structure BVHPrimitiveInfo

            M_primitiveInfo.clear();
            M_primitiveInfo.reserve( nelements(range) );
            for ( auto const& eltWrap : range )
            {
                auto const& e = unwrap_ref( eltWrap );
                M_primitiveInfo.push_back( BVHPrimitiveInfo{e} );
            }
            this->buildRootTree();
        }


    // Verify if the ray intersects the whole bounding structure
    // Returns the integer corresponding to the intersected element
    // If no element is intersected, return -1
    int raySearch( ray_type const& rayon )
        {
            M_intersected_leaf = {};
            M_lengths = {};
            if ( !M_rootNode )
                buildRootTree();

            if ( M_rootNode->checkIntersection(rayon) )
            {
                traverse_stackless( M_rootNode.get(), rayon );
            }
            if ( !M_intersected_leaf.empty() )
            {
                int argmin_lengths = std::distance(M_lengths.begin(), std::min_element(M_lengths.begin(), M_lengths.end()));
                int closer_intersection_element = M_intersected_leaf[argmin_lengths];
                return closer_intersection_element;
            }
            else
            {
                return -1;
            }
        }
private:
#if 0
    // From the mesh, build the bounding box info for each element and store it in
    // the structure BVHPrimitiveInfo
    void buildPrimitivesInfo(trace_mesh_ptrtype const& mesh)
        {
            auto rangeElem = elements(mesh);
            M_primitiveInfo.clear();
            M_primitiveInfo.reserve( nelements(rangeElem) );
            for ( auto const& eltWrap : rangeElem )
            {
                auto const& e = unwrap_ref( eltWrap );
                M_primitiveInfo.push_back( BVHPrimitiveInfo{e} );
            }
        }
#endif

    void buildRootTree()
        {
            if ( M_rootNode )
                return;

            M_rootNode = std::make_unique<BVHNode>();

            std::stack<std::tuple<BVHNode*,int,int,int>> stack;
            stack.push( std::make_tuple(M_rootNode.get(),0,0,M_primitiveInfo.size()) );
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
                        int primNum = M_primitiveInfo[i].meshElement().id();
                        M_orderedPrims.push_back(primNum);
                    }
                    currentNode->updateForUse( firstPrimOffset, nPrimitives, -1, bound_min_node, bound_max_node );
                }
                else
                {
                    CHECK( start_index_primitive >=0 && end_index_primitive <= M_primitiveInfo.size() ) << start_index_primitive << " " << end_index_primitive;
                    auto mid = (start_index_primitive + end_index_primitive) / 2;
                    std::nth_element(&M_primitiveInfo[start_index_primitive], &M_primitiveInfo[mid],
                                     &M_primitiveInfo[end_index_primitive-1]+1,
                                     [cut_dimension=cut_dimension](const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b) {
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


    // bool check_intersection_square(BVHRay const& ray)
    // {
    //     Eigen::VectorXd p1(3),p2(3),p3(3),n1(3);

    //     Eigen::VectorXd origin(3),direction(3);
    //     direction << ray.dir[0],ray.dir[1],ray.dir[2];
    //     origin << ray.origin[0],ray.origin[1],ray.origin[2];
    //     n1 << 0.,0.,-1.;
    //     p1 << 0.,0.,1.;
    //     p2 << 1.,0.,1.;
    //     p3 << 0.,1.,1.;

    //     auto t = n1.dot(p1-origin)/n1.dot(direction);
    //     LOG(INFO) << fmt::format("Origin {}, direction {}, normal {}, t {}, p1 {}",origin, direction, n1, t,p1);

    //     if(t<0)
    //         return false;

    //     auto intersection = origin + t*direction;
    //     auto v=intersection-p1;

    //     double width=1;
    //     double height=1;

    //     double proj1=v.dot(p2-p1);
    //     double proj2=v.dot(p3-p1);

    //     if(proj1 < width && proj1>0 && proj2>0 && proj2<height)
    //         return true;

    //     return false;
    // }
    // void loop_over_primitives(BVHTree::BVHNode * tree, BVHRay const& rayon)
    // {
    //     for(auto & primitive : M_primitiveInfo )
    //     {
    //         if (tree->checkIntersectionWithTriangle(rayon,primitive.M_primitiveNumber))
    //             M_intersected_leaf.push_back(primitive.M_primitiveNumber);
    //     }
    // }


    std::tuple<vector_realdim_type,vector_realdim_type> bounds( int start_index_primitive, int end_index_primitive ) const
        {
            if ( start_index_primitive >= end_index_primitive )
                throw std::logic_error("Error in BVHNode : compute bounds with no elemnent");

            //vector_realdim_type newBoundsMin, newBoundsMax;
            vector_realdim_type newBoundsMin = M_primitiveInfo[start_index_primitive].boundMin();
            vector_realdim_type newBoundsMax = M_primitiveInfo[start_index_primitive].boundMax();
            for (int i = start_index_primitive+1; i < end_index_primitive; ++i)
            {
                auto const& primitiveInfo = M_primitiveInfo[i];
                for ( uint8_type d=0;d<vector_realdim_type::SizeAtCompileTime;++d )
                {
                    newBoundsMin[d] = std::min( newBoundsMin[d], primitiveInfo.boundMin()[d] );
                    newBoundsMax[d] = std::max( newBoundsMax[d], primitiveInfo.boundMax()[d] );
                }
            }
            return std::make_tuple( std::move(newBoundsMin), std::move(newBoundsMax) );
        }


    void traverse_stackless(BVHTree::BVHNode * tree, ray_type const& rayon)
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
                        auto [has_intersected_leaf,distance] = current_node->checkLeafIntersection(rayon,M_primitiveInfo);
                        if ( has_intersected_leaf )
                        {
                            if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), M_primitiveInfo[current_node->firstPrimOffset()].meshElement().id()) == M_intersected_leaf.end() )
                            {
                                M_intersected_leaf.push_back(M_primitiveInfo[current_node->firstPrimOffset()].meshElement().id());
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
                        auto [has_intersected_leaf,distance] = current_node->checkLeafIntersection(rayon,M_primitiveInfo);
                        if ( has_intersected_leaf )
                        {
                            if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), M_primitiveInfo[current_node->firstPrimOffset()].meshElement().id()) == M_intersected_leaf.end() )
                            {
                                M_intersected_leaf.push_back(M_primitiveInfo[current_node->firstPrimOffset()].meshElement().id());
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

    std::vector<BVHPrimitiveInfo> M_primitiveInfo;
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
    using bvh_type = BVHTree<mesh_entity_type>;
    auto bvh = std::make_unique<bvh_type>();
    bvh->updateForUse(range);
    return bvh;
}

} // Feel
#endif /* FEELPP_MESH_BVH_HPP */
