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
    class Ray_bvh
    { 
    public: 
        using vec_t = Eigen::VectorXd;
        Ray_bvh(const vec_t &orig, const vec_t &dir) : 
            origin(orig),
            dir(dir) 
        {} 
        Ray_bvh( const Ray_bvh &r) :
            origin(r.origin), 
            dir(r.dir)
        {}
        Ray_bvh() : 
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

    template<int nDim>
    class BVH_tree
    {
        typedef Simplex<nDim,1> convex_type;
        typedef Mesh<convex_type> mesh_type;        
        typedef typename mesh_type::trace_mesh_type trace_mesh_type;
        typedef typename mesh_type::trace_mesh_ptrtype trace_mesh_ptrtype;
        typedef typename matrix_node<double>::type matrix_node_type;

        // Information on the primitive (bounding box, index, centroid)
        public:
        struct BVHPrimitiveInfo
        {
            BVHPrimitiveInfo(int primitiveNumber, Eigen::VectorXd bounds_min, Eigen::VectorXd bounds_max)
            {
                M_primitiveNumber=primitiveNumber;
                M_bound_min = bounds_min;
                M_bound_max = bounds_max;
                M_centroid = ( M_bound_min + M_bound_max ) * 0.5;
            }
            int M_primitiveNumber;
            Eigen::VectorXd M_bound_min;
            Eigen::VectorXd M_bound_max;
            Eigen::VectorXd M_centroid;
        };

        static inline std::vector<BVHPrimitiveInfo> M_primitiveInfo;
        static inline trace_mesh_ptrtype M_mesh;   
        static std::vector<std::pair<int,double>> M_distance_points;
        static inline std::map<int,int> M_primitiveN_to_elId;
        static inline std::vector<int> M_intersected_leaf;
        static inline std::vector<double> M_lengths;

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
                Eigen::VectorXd M_bound_min(mesh->realDimension()),M_bound_max(mesh->realDimension());
                auto v1 = e.point(0).node();
                double* ptr_data = &v1[0];
                M_bound_min = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptr_data, v1.size());
                M_bound_max = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptr_data, v1.size());
            
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
                M_primitiveN_to_elId.insert(std::make_pair(index_Primitive,e.id()));
                M_primitiveInfo.push_back(primitive_i);
                index_Primitive++;
            }
        }    

        struct BVH_node
        {
            void buildLeaf(BVH_node * current_parent, int first, int n, Eigen::VectorXd bounds_min, Eigen::VectorXd bounds_max)
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

            void buildInternalNode(BVH_node * current_parent, int splitaxisIn, BVH_node *child0, BVH_node *child1)
            {
                children[0]= child0;
                children[1]= child1;
                parent = current_parent;
                M_bounds_min = new_Bounds_min(child0->M_bounds_min,child1->M_bounds_min);
                M_bounds_max = new_Bounds_max(child0->M_bounds_max,child1->M_bounds_max);
                M_centroid = (M_bounds_min + M_bounds_max) *0.5;
                splitaxis = splitaxisIn;
                nPrimitives = 0;
                LOG(INFO) <<fmt::format("Internal node built: M_bounds_min {}, M_bounds_max {}, splitaxis {}",M_bounds_min,M_bounds_max,splitaxis);
            }

            Eigen::VectorXd new_Bounds_min(Eigen::VectorXd bounds1,Eigen::VectorXd bounds2)
            {
                Eigen::VectorXd newBoundsMin(bounds1.size());
                for(int i=0;i<bounds1.size();i++)
                    newBoundsMin[i] = std::min(bounds1[i],bounds2[i]);
                return newBoundsMin;
            }
            Eigen::VectorXd new_Bounds_max(Eigen::VectorXd bounds1,Eigen::VectorXd bounds2)
            {
                Eigen::VectorXd newBoundsMax(bounds1.size());
                for(int i=0;i<bounds1.size();i++)
                    newBoundsMax[i] = std::max(bounds1[i],bounds2[i]);
                return newBoundsMax;
            }

            bool isLeaf(){return (nPrimitives !=0) ;}
            

            BVH_node * nearChild(Ray_bvh const& ray)
            {
                
                if(ray.dir(this->splitaxis)>0)
                    return this->children[0];
                else
                    return this->children[1];

            }
            
            BVH_tree::BVH_node * otherChild(BVH_tree::BVH_node * parent)
            {
                if (this==parent->children[0])
                    return parent->children[1];
                else
                    return parent->children[0];
            }

            bool checkIntersection(Ray_bvh const& rayon)
            {
                double tmin = 0.0;
                double tmax = FLT_MAX;

                for(int i=0; i<nDim; i++)
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
            
            bool check_intersection_with_segment(Ray_bvh const& ray)
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

            std::pair<bool,double> check_intersection_with_segment(matrix_node_type const& nodes, Ray_bvh const& ray)
            {
                Eigen::VectorXd p1(2),p2(2),v1(2),v2(2),v3(2),w_(2);     

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
            std::pair<bool,double> check_intersection_with_triangle( Ray_bvh const& ray)
            {
                Eigen::Vector3d p1(3),p2(3),p3(3),n1(3),w(3),w_(3);
                Eigen::Matrix3d m(3,3);
                auto nodes =  M_mesh->element(M_primitiveInfo[this->firstPrimOffset].M_primitiveNumber).vertices();

                for ( int i = 0; i < nodes.size1(); ++i )
                {p1( i ) = nodes(i,0);
                    p2( i ) = nodes(i,1);
                    p3( i ) = nodes(i,2);}
                Eigen::Vector3d origin(3),direction(3);
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

            std::pair<bool,double> checkLeafIntersection(Ray_bvh const& rayon,bool intersected=false)
            {             
                if( nDim ==2)
                {
                    auto nodes = M_mesh->element( M_primitiveInfo[this->firstPrimOffset].M_primitiveNumber ).vertices();                                        
                    return check_intersection_with_segment(nodes,rayon);
                }                
                else //if ( nDim==3)
                    return check_intersection_with_triangle(rayon);

            }

            BVH_node *children[2];
            BVH_node *parent;
            int splitaxis, nPrimitives,firstPrimOffset;
            Eigen::VectorXd M_bounds_min,M_bounds_max,M_centroid;  
        };    

        BVH_node *  M_root_tree;

        BVH_node * buildRootTree()
        {
            M_root_tree = recursiveBuild(M_root_tree,0,0,M_primitiveInfo.size(),orderedPrims);
            return M_root_tree;
        }                        

        BVH_node * recursiveBuild(BVH_node * current_parent, int cut_dimension, int start_index_primitive, int end_index_primitive,
                                std::vector<int> &orderedPrims)
        {
            LOG(INFO) <<fmt::format("cut dimension {}, start index primitive {}, end index primitive {}",cut_dimension,start_index_primitive,end_index_primitive);
            Eigen::VectorXd M_bound_min_node(nDim),M_bound_max_node(nDim);
            BVH_node * node = new BVH_tree::BVH_node();
            M_bound_min_node = M_primitiveInfo[start_index_primitive].M_bound_min;
            M_bound_max_node = M_primitiveInfo[start_index_primitive].M_bound_max;
            for (int i = start_index_primitive+1; i < end_index_primitive; ++i)
            {
                M_bound_min_node = node->new_Bounds_min(M_bound_min_node,M_primitiveInfo[i].M_bound_min);
                M_bound_max_node = node->new_Bounds_max(M_bound_max_node,M_primitiveInfo[i].M_bound_max);
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
                node->buildLeaf(current_parent,firstPrimOffset, nPrimitives, M_bound_min_node,M_bound_max_node);
                return node;
            }
            else{
                // Create a node, since there are at least two primitives in the list
                node->buildInternalNode(current_parent,(cut_dimension+1)%nDim,
                                    recursiveBuild( node, (cut_dimension+1)%nDim, start_index_primitive, mid, orderedPrims),
                                    recursiveBuild( node, (cut_dimension+1)%nDim, mid, end_index_primitive, orderedPrims));
            }

            return node;
        }

        // Verify if the ray intersects the whole bounding structure

        void ray_search( Ray_bvh const& rayon,std::string s)
        {   
            M_intersected_leaf.clear();
            M_lengths.clear();
            if(!M_root_tree)   
                M_root_tree = buildRootTree();
            // compute tmin and tmax
            double tmin, tmax;
           
            const BVH_tree::BVH_node *tn = static_cast<const BVH_tree::BVH_node*>( M_root_tree );
            
            Eigen::VectorXd mini = tn->M_bounds_min;
            Eigen::VectorXd maxi = tn->M_bounds_max;
                    
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
        }

        // bool check_intersection_square(Ray_bvh const& ray)        
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
        // void loop_over_primitives(BVH_tree::BVH_node * tree, Ray_bvh const& rayon)
        // {
        //     for(auto & primitive : M_primitiveInfo )
        //     {
        //         if (tree->check_intersection_with_triangle(rayon,primitive.M_primitiveNumber))
        //             M_intersected_leaf.push_back(primitive.M_primitiveNumber);
        //     }
        // }

        void traverse_stackless(BVH_tree::BVH_node * tree, Ray_bvh const& rayon)
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
                        auto [has_intersected_leaf,distance] = current_node->checkLeafIntersection(rayon);
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
                        auto [has_intersected_leaf,distance] = current_node->checkLeafIntersection(rayon);
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
                    std::cout << "ERROR: None of the previous cases has been traversed" <<std::endl;
                    break;
                }
            }
        }
    };

} // Feel
#endif /* __BVH_tree_H */