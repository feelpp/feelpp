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

        std::vector<int> orderedPrims; // order of traversed primitives for depth-first search

        // From the mesh, build the bounding box info for each element and store it in 
        // the structure BVHPrimitiveInfo
        void buildPrimitivesInfo(trace_mesh_ptrtype mesh)
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
                    // std::cout << e.point(j)<<std::endl;
                    for(int k=0;k<mesh->realDimension();k++)
                    {
                        M_bound_min[k] = std::min(M_bound_min[k],e.point(j).node()[k])-2*FLT_MIN;
                        M_bound_max[k] = std::max(M_bound_max[k],e.point(j).node()[k])+2*FLT_MIN;
                    }
                }
                // std::cout << "ID=" << e.id() << "VERTICES=" << e.vertices() <<std::endl;
                // std::cout << "Bound min " << M_bound_min << "Bound max" << M_bound_max <<std::endl;
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

            void buildInternalNode(BVH_node * current_parent, int splitaxis, BVH_node *child0, BVH_node *child1)
            {
                children[0]= child0;
                children[1]= child1;
                parent = current_parent;
                M_bounds_min = new_Bounds_min(child0->M_bounds_min,child1->M_bounds_min);
                M_bounds_max = new_Bounds_max(child0->M_bounds_max,child1->M_bounds_max);
                M_centroid = (M_bounds_min + M_bounds_max) *0.5;
                splitaxis = splitaxis;
                nPrimitives = 0;
                LOG(INFO) <<fmt::format("Internal node built: M_bounds_min {}, M_bounds_max {}, splitaxis {}",M_bounds_min,M_bounds_max,splitaxis);
            }

            Eigen::VectorXd new_Bounds_min(Eigen::VectorXd bounds1,Eigen::VectorXd bounds2)
            {
                Eigen::VectorXd newBoundsMin(bounds1.size());
                for(auto i=0;i<bounds1.size();i++)
                    newBoundsMin[i] = std::min(bounds1[i],bounds2[i]);
                return newBoundsMin;
            }
            Eigen::VectorXd new_Bounds_max(Eigen::VectorXd bounds1,Eigen::VectorXd bounds2)
            {
                Eigen::VectorXd newBoundsMax(bounds1.size());
                for(auto i=0;i<bounds1.size();i++)
                    newBoundsMax[i] = std::max(bounds1[i],bounds2[i]);
                return newBoundsMax;
            }

            bool isLeaf(){return (nPrimitives !=0) ;}

            double split_v(int axis)
            {
                return M_centroid[axis];
            }

            BVH_node * nearChild(Ray_bvh ray)
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

            bool checkIntersection(Ray_bvh rayon)
            {
                double div1 = 1./(rayon.dir[0]+2*FLT_MIN);
                double div2 = 1./(rayon.dir[1]+2*FLT_MIN);
                double div3 = 1./(rayon.dir[2]+2*FLT_MIN);


                double t1 = (M_bounds_min(0) - rayon.origin[0]) * div1;
                double t2 = (M_bounds_max(0) - rayon.origin[0]) * div1 * (1 + 2 * robust_const(3));
                double t3 = (M_bounds_min(1) - rayon.origin[1]) * div2;
                double t4 = (M_bounds_max(1) - rayon.origin[1]) * div2 * (1 + 2 * robust_const(3));
                double t5 = (M_bounds_min(2) - rayon.origin[2]) * div3;
                double t6 = (M_bounds_max(2) - rayon.origin[2]) * div3 * (1 + 2 * robust_const(3));

                double tmin = std::max(std::max(std::min(t1, t2), std::min(t3, t4)), std::min(t5, t6));
                double tmax = std::min(std::min(std::max(t1, t2), std::max(t3, t4)), std::max(t5, t6));

                if (tmax < 0)
                {
                    return false;
                }
                if (tmin > tmax)
                {                    
                    return false;
                }

                // double tmin = 0.0;
                // double tmax = FLT_MAX;

                // for(int i=0; i<nDim; i++)
                // {
                //     if(math::abs(rayon.dir[i]) < 1e-4)
                //     {
                //         if(rayon.origin[i] < M_bounds_min[i] || rayon.origin[i] > M_bounds_max[i]) 
                //             return false;                        
                //     }
                //     else
                //     {
                //         double ratio = 1.0/rayon.dir[i];
                //         double t1 = (M_bounds_min[i]-rayon.origin[i]) * ratio;
                //         double t2 = (M_bounds_max[i]-rayon.origin[i]) * ratio;
                        // if (t1 > t2) 
                        // {
                        //     double tTemp = t1;
                        //     t1 = t2;
                        //     t2 = tTemp;
                        // }
                        // if ( t1 > tmin) 
                        //     tmin = t1;
                        // if (t2 > tmax)
                        //     tmax = t2;
                        // if (tmin > tmax)
                        //     return false;
                //     }
                // }
                
                return true;         
            }
            
            // Verify if the ray intersects the element
            bool check_intersection_with_triangle( Ray_bvh &ray)
            {
                Eigen::Vector3d p1(3),p2(3),p3(3),n1(3),w(3),w_(3);
                Eigen::Matrix3d m(3,3);

                auto nodes =  M_mesh->element( M_primitiveInfo[this->firstPrimOffset].M_primitiveNumber).vertices();

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
                
                // intersection point
                w = origin + direction* ( (p1-origin).dot(n1))/(direction.dot(n1));            
                if(w==origin)
                    return false;
                m.col(0) = p2-p1;
                m.col(1) = p3-p1;
                m.col(2) = n1;    
                w_ = m.inverse()*(w-p1);    

                // std::cout << "Normal n1 " <<n1 << " distance " << ((p1-origin).dot(n1))/(direction.dot(n1))<< std::endl;
                // std::cout << "Intersection pt " <<w << " w_ " << w_<< std::endl;
      

                // M_distance_points.push_back(std::make_pair(this->M_primitiveNumber,w));

                return (w_(0)>=0) && (w_(1)>=0) && (w_(0) +  w_(1)<=1);
                
            }

            bool checkLeafIntersection(Ray_bvh rayon)
            {             
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

        void ray_search( Ray_bvh &rayon,std::string s)
        {   
      
            if(!M_root_tree)   
                M_root_tree = buildRootTree();
            // compute tmin and tmax
            double tmin, tmax;
           
            const BVH_tree::BVH_node *tn = static_cast<const BVH_tree::BVH_node*>( M_root_tree );
            // ////LOG(INFO) << fmt::format("initial node: {}",*tn);
            Eigen::VectorXd mini = tn->M_bounds_min;
            Eigen::VectorXd maxi = tn->M_bounds_max;
        
            // std::cout << "mini " << mini << " maxi " << maxi << std::endl;
            // spdlog::debug("Check bounding box coordinates : min = ({},{},{}) max = ({},{},{})", mini(0),mini(1),mini(2), maxi(0),maxi(1),maxi(2));
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
                std::cout << "Intersection with BVH is found" <<std::endl;
                traverse_stackless(M_root_tree, rayon);
            }            
            //run_ray_search_Book( M_root_tree, rayon, tmin, tmax, 0 ,s);       
#endif
        }


        void traverse_stackless(BVH_tree::BVH_node * tree, Ray_bvh rayon)
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
                        // std::cout << state << std::endl;
                    }
                    else
                    {
                        current_node=current_node->parent;
                        state='C'; // the current node has been accessed from its sibling
                        // std::cout << state << std::endl;
                    }
                    break;

                case 'S': // the node is being traversed from its Sibling ('S')
                    
                    if( current_node->checkIntersection(rayon)==false) // back to parent
                    {
                        // std::cout << fmt::format("the ray with point {} and dir {} does not intersect the node with min {} max {}", 
                        //         rayon.origin,rayon.dir,current_node->M_bounds_min, current_node->M_bounds_max) <<std::endl;
                        current_node=current_node->parent;
                        state='C'; // the current node is being accessed from its child
                        // std::cout << state << std::endl;
                    }
                    else if(current_node->isLeaf())
                    {
                        std::cout << fmt::format("we are on a leaf and the ray with point {} and dir {} intersects the node with min {} max {}", 
                                rayon.origin,rayon.dir,current_node->M_bounds_min, current_node->M_bounds_max) <<std::endl;
                        bool has_intersected_leaf = current_node->checkLeafIntersection(rayon);
                        if(has_intersected_leaf)
                        {
                            std::cout << "Leaf intersected" <<std::endl;
                        }
                        current_node=current_node->parent;
                        state='C'; // the current node is being accessed from its child
                        // std::cout << state << std::endl;
                    }
                    else
                    {
                        // std::cout << fmt::format("the ray with point {} and dir {} intersects the node with min {} max {}", 
                        //         rayon.origin,rayon.dir,current_node->M_bounds_min, current_node->M_bounds_max) <<std::endl;
                        current_node = current_node->nearChild(rayon);
                        state='P';// the current node has been accessed from its parent
                        // std::cout << state << std::endl;
                    }
                    break; 

                case 'P':
                    if( current_node->checkIntersection(rayon)==false )
                    {
                        // std::cout << fmt::format("the ray with point {} and dir {} does not intersect the node with min {} max {}", 
                        //         rayon.origin,rayon.dir,current_node->M_bounds_min, current_node->M_bounds_max) <<std::endl;
                        current_node=current_node->otherChild(current_node->parent);
                        state='S'; // the current node has been accessed from its sibling
                        // std::cout << state << std::endl;
                    }
                    else if(current_node->isLeaf())
                    {
                        std::cout << fmt::format("we are on a leaf and the ray with point {} and dir {} intersects the node with min {} max {}", 
                                rayon.origin,rayon.dir,current_node->M_bounds_min, current_node->M_bounds_max) <<std::endl;
                        bool has_intersected_leaf = current_node->checkLeafIntersection(rayon);
                        if(has_intersected_leaf)
                        {
                            std::cout << "Leaf intersected" <<std::endl;
                        }
                        current_node=current_node->otherChild(current_node->parent);
                        state='S'; // the current node has been accessed from its sibling
                        // std::cout << state << std::endl;
                    }
                    else
                    {
                        // std::cout << fmt::format("the ray with point {} and dir {} intersects the node with min {} max {}", 
                        //         rayon.origin,rayon.dir,current_node->M_bounds_min, current_node->M_bounds_max) <<std::endl;
                        current_node = current_node->nearChild(rayon);
                        state='P'; // the current node has been accessed from its parent
                        // std::cout << state << std::endl;
                    }
                    break;

                default:
                    std::cout << "ERROR: None of the previous cases has been traversed" <<std::endl;
                    break;
                }
            }
        }
        void
        run_ray_search_Book( BVH_tree::BVH_node * tree, Ray_bvh rayon, double tmin, double tmax, uint16_type iter,std::string s )
        {
            if(tree==NULL)
            {
                LOG(INFO) << fmt::format("node is null");
                return;
            }
            if ( ! tree->isLeaf() )
            {
                LOG(INFO)<< fmt::format("Origin: {}; direction: {} ",rayon.origin,rayon.dir);
                
                const BVH_tree::BVH_node *tn = static_cast<const BVH_tree::BVH_node*>( tree );
                BVH_tree::BVH_node *nearChild, *farChild;
                
                int splitaxis=0;
                size_type N = tn->M_bounds_min.size();
                if ( ( iter%N )==0 )
                    splitaxis = 0;
                else if ( ( iter%N )==1 )
                    splitaxis = 1;
                else if ( ( iter%N )==2 )
                    splitaxis = 2;

                bool first = rayon.origin[splitaxis] > tn->split_v;                
                LOG(INFO) << fmt::format("tn->split v :{}, split dir {}",tn->split_v,splitaxis);        
                if(first)
                {
                    nearChild = tn->right;
                    farChild = tn->left;            
                    LOG(INFO) << fmt::format("right first");
                }
                else
                {
                    nearChild = tn->left;
                    farChild = tn->right;
                    LOG(INFO) << fmt::format("left first");
                }        

                if (std::abs(rayon.dir[splitaxis])<=1e-8)
                {
                    run_ray_search_Book(nearChild,rayon,tmin,tmax,iter+1,s);
                }
                else
                {
                    double div = 1./ (rayon.dir[splitaxis]);
                    double tsplit = (tn->split_v - rayon.origin[splitaxis]) * div;
                    LOG(INFO) << fmt::format("tmin {} tsplit {} tmax {} ",tmin,tsplit,tmax);
                    if (0. <= tsplit && tsplit < tmax + 1e-8)
                    {
                        run_ray_search_Book(nearChild,rayon,tmin,tsplit,iter+1,s);
                        rayon.origin += tsplit * rayon.dir;

                        run_ray_search_Book(farChild,rayon,0.,tmax-tsplit,iter+1,s);
                    }
                    else
                    {
                        run_ray_search_Book(nearChild,rayon,tmin,tmax,iter+1,s);
                    }

                    
                }
            }
            else       
            {
                // LOG(INFO) << fmt::format("problem here 5");
                const BVH_tree::BVH_node *tl = static_cast<const BVH_tree::BVH_node*>( tree );
                if(s=="print_info")        
                {
                    LOG(INFO) << fmt::format("leaf: {}",*tl);
                    LOG(INFO) << fmt::format("===========");
                }
                for ( size_type i=tl->firstPrimOffset; i<tl->firstPrimOffset+tl->nPrimitives; ++i )
                    ;//triangle_intersection(rayon, i);                
                
                return;
            }



        }

    };

} // Feel
#endif /* __BVH_tree_H */