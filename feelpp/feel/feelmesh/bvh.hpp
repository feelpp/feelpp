// Bounding volume hierarchy

#ifndef FEELPP_BVH_HPP
#define FEELPP_BVH_HPP 1


#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>

#include <feel/feel.hpp>

namespace Feel
{
    class BVHRay
    { 
    public: 
        using vec_t = Eigen::VectorXd;
        BVHRay(const vec_t &orig, const vec_t &dir) : 
            origin(orig),
            dir(dir) 
        {} 
        BVHRay( const BVHRay &r) :
            origin(r.origin), 
            dir(r.dir)
        {}
        BVHRay() : 
            origin(),
            dir() 
        {} 
        vec_t origin, dir;  // ray origin and dir 
    };

    inline constexpr float robust_const(int n) 
    {
        double epsil = std::numeric_limits<double>::epsilon() * 0.5;
        return (n * epsil) / (1 - n * epsil);
    }

    template<int nDim, int realDim>
    class BVHTree
    {
        using self_type = BVHTree<nDim,realDim>;
        using convex_type = Simplex<nDim,1,realDim>;
        using mesh_type = Mesh<convex_type>;     
        using VectorRealDim = typename  std::conditional< realDim==3,
                                                            Eigen::Vector3d,
                                                            Eigen::Vector2d>::type;        
        using trace_mesh_type = typename  std::conditional< nDim==realDim,
                                                    typename mesh_type::trace_mesh_type,
                                                    mesh_type >::type;
        using trace_mesh_ptrtype =  typename  std::conditional< nDim==realDim,
                                                        typename mesh_type::trace_mesh_ptrtype,
                                                        typename mesh_type::ptrtype >::type ;
        // typedef typename mesh_type::trace_mesh_type trace_mesh_type;
        // typedef typename mesh_type::trace_mesh_ptrtype trace_mesh_ptrtype;
        typedef typename matrix_node<double>::type matrix_node_type;

        // Information on the primitive (bounding box, index, centroid)
        public:
        struct BVHPrimitiveInfo
        {
            BVHPrimitiveInfo(int primitiveNumber, VectorRealDim bounds_min, VectorRealDim bounds_max)
            {
                M_primitiveNumber=primitiveNumber;
                M_bound_min = bounds_min;
                M_bound_max = bounds_max;
                M_centroid = ( M_bound_min + M_bound_max ) * 0.5;            
            }
            int M_primitiveNumber;
            VectorRealDim M_bound_min;
            VectorRealDim M_bound_max;
            VectorRealDim M_centroid;
        };               
        std::vector<BVHPrimitiveInfo> M_primitiveInfo;
        trace_mesh_ptrtype M_mesh;
        thread_local static inline std::vector<int> M_intersected_leaf;
        thread_local static inline std::vector<double> M_lengths;

        std::vector<int> orderedPrims; // order of traversed primitives for depth-first search

        // From the mesh, build the bounding box info for each element and store it in 
        // the structure BVHPrimitiveInfo
        void buildPrimitivesInfo(trace_mesh_ptrtype const& mesh)
        {
            M_mesh = mesh;
            auto elem = elements(mesh);
            int index_Primitive=0;
            for(auto &elt : elem)
            {
                auto const& e = unwrap_ref( elt );                
                VectorRealDim M_bound_min,M_bound_max;
                auto v1 = e.point(0).node();
                double* ptr_data = &v1[0];
                M_bound_min = Eigen::Map<VectorRealDim, Eigen::Unaligned>(ptr_data, v1.size());
                M_bound_max = Eigen::Map<VectorRealDim, Eigen::Unaligned>(ptr_data, v1.size());
            
                for(int j=1;j<e.nPoints();j++)
                {
                    for(int k=0;k<mesh->realDimension();k++)
                    {
                        M_bound_min[k] = std::min(M_bound_min[k],e.point(j).node()[k]);
                        M_bound_max[k] = std::max(M_bound_max[k],e.point(j).node()[k]);
                    }
                }
                for(int k=0;k<mesh->realDimension();k++)
                {
                    M_bound_min[k] = M_bound_min[k] - 2*FLT_MIN;
                    M_bound_max[k] = M_bound_max[k] + 2*FLT_MIN;
                }                
                BVHPrimitiveInfo primitive_i(e.id(),M_bound_min,M_bound_max);
                M_primitiveInfo.push_back(primitive_i);
                index_Primitive++;
            }
        }    

        class BVHNode
        {
            public: 
            void buildLeaf(BVHNode * current_parent, int first, int n, Eigen::VectorXd bounds_min, Eigen::VectorXd bounds_max)
            {
                firstPrimOffset = first;
                nPrimitives = n;
                M_bounds_min = bounds_min;
                M_bounds_max = bounds_max;
                children[0]=nullptr;
                children[1]=nullptr;
                parent = current_parent;
                LOG(INFO) <<fmt::format("leaf built: firstPrimOffset {}, nPrimitives, {}",firstPrimOffset,nPrimitives);
            }

            void buildInternalNode(BVHNode * current_parent, int splitaxisIn, BVHNode *child0, BVHNode *child1)
            {
                children[0]= child0;
                children[1]= child1;
                parent = current_parent;
                M_bounds_min = newBoundsMin(child0->M_bounds_min,child1->M_bounds_min);
                M_bounds_max = newBoundsMax(child0->M_bounds_max,child1->M_bounds_max);
                M_centroid = (M_bounds_min + M_bounds_max) *0.5;
                splitaxis = splitaxisIn;
                nPrimitives = 0;
                LOG(INFO) <<fmt::format("Internal node built: M_bounds_min {}, M_bounds_max {}, splitaxis {}",M_bounds_min,M_bounds_max,splitaxis);
            }

            VectorRealDim newBoundsMin(VectorRealDim bounds1,VectorRealDim bounds2)
            {
                VectorRealDim newBoundsMin(bounds1.size());
                for(int i=0;i<bounds1.size();i++)
                    newBoundsMin[i] = std::min(bounds1[i],bounds2[i]);
                return newBoundsMin;
            }
            VectorRealDim newBoundsMax(VectorRealDim bounds1,VectorRealDim bounds2)
            {
                VectorRealDim newBoundsMax(bounds1.size());
                for(int i=0;i<bounds1.size();i++)
                    newBoundsMax[i] = std::max(bounds1[i],bounds2[i]);
                return newBoundsMax;
            }

            bool isLeaf(){return (nPrimitives !=0) ;}
            

            BVHNode * nearChild(BVHRay const& ray)
            {
                
                if(ray.dir(this->splitaxis)>0)
                    return this->children[0];
                else
                    return this->children[1];

            }
            
            BVHTree::BVHNode * otherChild(BVHTree::BVHNode * parent)
            {
                if (this==parent->children[0])
                    return parent->children[1];
                else
                    return parent->children[0];
            }

            bool checkIntersection(BVHRay const& rayon)
            {
                double tmin = 0.0;
                double tmax = FLT_MAX;

                for(int i=0; i<realDim; i++)
                {                    
                    double ratio = 1.0/(rayon.dir[i]+2*FLT_MIN);
                    double t1 = (M_bounds_min[i]-rayon.origin[i]) * ratio;
                    double t2 = (M_bounds_max[i]-rayon.origin[i]) * ratio;
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
            std::pair<bool,double> checkIntersectionWithSegment(matrix_node_type const& nodes, BVHRay const& ray)
            {
                VectorRealDim p1,p2,v1,v2,v3,w_;     

                for ( int i = 0; i < nodes.size1(); ++i )
                {
                    p1( i ) = nodes(i,0);
                    p2( i ) = nodes(i,1);
                }        
                VectorRealDim origin,direction;
                direction << ray.dir[0],ray.dir[1];
                origin << ray.origin[0],ray.origin[1];        
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
                    w_[0] =  ray.origin[0] + ray.dir[0]*t1;
                    w_[1] =  ray.origin[1] + ray.dir[1]*t1;
                    return std::make_pair(true,t1);
                }

                return std::make_pair(false,t1);            
                
            }
            
            // Verify if the ray intersects the element
            std::pair<bool,double> checkIntersectionWithTriangle( BVHRay const& ray, trace_mesh_ptrtype mesh, std::vector<BVHPrimitiveInfo>& primitiveInfo)
            {
                VectorRealDim p1,p2,p3,n1,w,w_;
                Eigen::Matrix3d m(3,3);            
                    
                auto nodes =  mesh->element(primitiveInfo[this->firstPrimOffset].M_primitiveNumber).vertices();

                for ( int i = 0; i < nodes.size1(); ++i )
                {p1( i ) = nodes(i,0);
                    p2( i ) = nodes(i,1);
                    p3( i ) = nodes(i,2);}
                VectorRealDim origin,direction;
                direction << ray.dir[0],ray.dir[1],ray.dir[2];
                origin << ray.origin[0],ray.origin[1],ray.origin[2];

                
                // // normal vector
                n1 = (p2-p1).cross(p3-p1);
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
                w = origin + direction* t_line;             
               
                m.col(0) = p2-p1;
                m.col(1) = p3-p1;
                m.col(2) = n1;    
                w_ = m.inverse()*(w-p1);

                return std::make_pair((w_(0)> 2*FLT_MIN ) && (w_(1)>0) && (w_(0) +  w_(1)<1),t_line);
                
            }

            std::pair<bool,double> checkLeafIntersection(BVHRay const& rayon, trace_mesh_ptrtype mesh, std::vector<BVHPrimitiveInfo>& primitiveInfo)
            {             
                if constexpr (realDim ==2)
                {
                    auto nodes = mesh->element( primitiveInfo[this->firstPrimOffset].M_primitiveNumber ).vertices();                                        
                    return checkIntersectionWithSegment(nodes,rayon);
                }                
                if constexpr (realDim==3) //if ( realDim==3)
                    return checkIntersectionWithTriangle(rayon,mesh,primitiveInfo);

            }

            BVHNode *children[2];
            BVHNode *parent;
            int splitaxis, nPrimitives,firstPrimOffset;
            VectorRealDim M_bounds_min,M_bounds_max,M_centroid;
            
        };    

        BVHNode *  M_root_tree;

        BVHNode * buildRootTree()
        {
            VectorRealDim bound_min_node,bound_max_node;
            M_root_tree = recursiveBuild(M_root_tree,0,0,M_primitiveInfo.size(),orderedPrims,bound_min_node,bound_max_node);
            return M_root_tree;
        }                        

        BVHNode * recursiveBuild(BVHNode * current_parent, int cut_dimension, int start_index_primitive, int end_index_primitive,
                                std::vector<int> &orderedPrims, VectorRealDim &bound_min_node, VectorRealDim &bound_max_node)
        {
            LOG(INFO) <<fmt::format("cut dimension {}, start index primitive {}, end index primitive {}",cut_dimension,start_index_primitive,end_index_primitive);
            
            BVHNode * node = new BVHTree::BVHNode();
            bound_min_node = M_primitiveInfo[start_index_primitive].M_bound_min;
            bound_max_node = M_primitiveInfo[start_index_primitive].M_bound_max;
            for (int i = start_index_primitive+1; i < end_index_primitive; ++i)
            {
                bound_min_node = node->newBoundsMin(bound_min_node,M_primitiveInfo[i].M_bound_min);
                bound_max_node = node->newBoundsMax(bound_max_node,M_primitiveInfo[i].M_bound_max);
            }
            auto mid = (start_index_primitive + end_index_primitive) / 2;
            std::nth_element(&M_primitiveInfo[start_index_primitive], &M_primitiveInfo[mid], 
                    &M_primitiveInfo[end_index_primitive-1]+1,
                        [cut_dimension](const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b) { 
                            return a.M_centroid[cut_dimension] < b.M_centroid[cut_dimension];
                        });
            int nPrimitives = end_index_primitive - start_index_primitive;
            if (nPrimitives == 1) 
            {
                // Create a leaf, since there is only one primitive in the list
                int firstPrimOffset = orderedPrims.size();
                for (int i = start_index_primitive; i < end_index_primitive; ++i) 
                {
                int primNum = M_primitiveInfo[i].M_primitiveNumber;
                orderedPrims.push_back(primNum);
                }
                node->buildLeaf(current_parent,firstPrimOffset, nPrimitives, bound_min_node,bound_max_node);
                return node;
            }
            else{
                int next_cut_dimension=(cut_dimension+1)%realDim;
                // Create a node, since there are at least two primitives in the list
                node->buildInternalNode(current_parent, next_cut_dimension,
                                    recursiveBuild( node, next_cut_dimension, start_index_primitive, mid, orderedPrims, bound_min_node, bound_max_node),
                                    recursiveBuild( node, next_cut_dimension, mid, end_index_primitive, orderedPrims, bound_min_node, bound_max_node));
            }

            return node;
        }

        // Verify if the ray intersects the whole bounding structure
        // Returns the integer corresponding to the intersected element
        // If no element is intersected, return -1
        int raySearch( BVHRay const& rayon,std::string s)
        {   
            M_intersected_leaf = {};
            M_lengths = {};
            if(!M_root_tree)   
                M_root_tree = buildRootTree();
            // compute tmin and tmax
            double tmin, tmax;
           
            const BVHTree::BVHNode *tn = static_cast<const BVHTree::BVHNode*>( M_root_tree );
            
            VectorRealDim mini = tn->M_bounds_min;
            VectorRealDim maxi = tn->M_bounds_max;
                    
            double div1 = 1./rayon.dir[0];
            double div2 = 1./rayon.dir[1];
            double div3 = 1./rayon.dir[2];


            double t1 = (mini(0) - rayon.origin[0]) * div1;
            double t2 = (maxi(0) - rayon.origin[0]) * div1 * (1 + 2 * robust_const(3));
            double t3 = (mini(1) - rayon.origin[1]) * div2;
            double t4 = (maxi(1) - rayon.origin[1]) * div2 * (1 + 2 * robust_const(3));
            double t5 = (mini(2) - rayon.origin[2]) * div3;
            double t6 = (maxi(2) - rayon.origin[2]) * div3 * (1 + 2 * robust_const(3));

            tmin = std::max(std::max(std::min(t1, t2), std::min(t3, t4)), std::min(t5, t6));
            tmax = std::min(std::min(std::max(t1, t2), std::max(t3, t4)), std::max(t5, t6));            
#if 1
            if(M_root_tree->checkIntersection(rayon))
            {
                traverse_stackless(M_root_tree, rayon);
            }
#endif
            if (!M_intersected_leaf.empty())
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

        void traverse_stackless(BVHTree::BVHNode * tree, BVHRay const& rayon)
        {
            auto current_node = M_root_tree->nearChild(rayon);
            char state = 'P'; // the current node is being traversed from its Parent ('P')      

            while(true)
            {
                switch (state)
                {
                case 'C': // the node is being traversed from its child
                    
                    if(current_node==M_root_tree) return;

                    if(current_node==current_node->parent->nearChild(rayon))
                    {
                        current_node=current_node->otherChild(current_node->parent);
                        state='S'; // the current node has been accessed from its sibling                        
                    }
                    else
                    {
                        current_node=current_node->parent;
                        state='C'; // the current node has been accessed from its sibling
                    }
                    break;

                case 'S': // the node is being traversed from its Sibling ('S')
                    
                    if( current_node->checkIntersection(rayon)==false) // back to parent
                    {                        
                        current_node=current_node->parent;
                        state='C'; // the current node is being accessed from its child                     
                    }
                    else if(current_node->isLeaf())
                    {                        
                        auto [has_intersected_leaf,distance] = current_node->checkLeafIntersection(rayon,M_mesh,M_primitiveInfo);
                        if(has_intersected_leaf)
                        {
                            if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), M_primitiveInfo[current_node->firstPrimOffset].M_primitiveNumber) == M_intersected_leaf.end() )
                            {
                                M_intersected_leaf.push_back(M_primitiveInfo[current_node->firstPrimOffset].M_primitiveNumber);
                                M_lengths.push_back(distance);
                            }                            
                        }
                        current_node=current_node->parent;
                        state='C'; // the current node is being accessed from its child
                    }
                    else
                    {                        
                        current_node = current_node->nearChild(rayon);
                        state='P';// the current node has been accessed from its parent                     
                    }
                    break; 

                case 'P':
                    if( current_node->checkIntersection(rayon)==false )
                    {
                        current_node=current_node->otherChild(current_node->parent);
                        state='S'; // the current node has been accessed from its sibling
                    }
                    else if(current_node->isLeaf())
                    {
                        auto [has_intersected_leaf,distance] = current_node->checkLeafIntersection(rayon,M_mesh,M_primitiveInfo);
                        if(has_intersected_leaf)
                        {
                            if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), M_primitiveInfo[current_node->firstPrimOffset].M_primitiveNumber) == M_intersected_leaf.end() )
                            {
                                M_intersected_leaf.push_back(M_primitiveInfo[current_node->firstPrimOffset].M_primitiveNumber);
                                M_lengths.push_back(distance);
                            }                               
                        }
                        current_node=current_node->otherChild(current_node->parent);
                        state='S'; // the current node has been accessed from its sibling                        
                    }
                    else
                    {                        
                        current_node = current_node->nearChild(rayon);
                        state='P'; // the current node has been accessed from its parent                     
                    }
                    break;

                default:

                    LOG(ERROR) << "ERROR: None of the previous cases has been traversed";
                    
                    throw std::logic_error("Error in BVH traversal: none of the previous cases has been traversed.");

                    break;
                }
            }
        }
    };

} // Feel
#endif /* __BVHTree_H */