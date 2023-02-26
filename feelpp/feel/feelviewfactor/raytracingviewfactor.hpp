/**
 * @file raytracingviewfactor.hpp
 * @author Christophe Prud'homme (christophe.prudhomme@cemosis.fr)
 * @brief
 * @version 0.1
 * @date 2022-07-22
 *
 * @copyright Copyright (c) 2022 Feel++ Consortium
 * @copyright Copyright (c) 2022 Université de Strasbourg
 *
 */
#pragma once

#include <feel/feel.hpp>
#include <feel/feelviewfactor/viewfactorbase.hpp>
#include <nanoflann.hpp>
#include <feel/feelviewfactor/kdtreevectorofvectorsadaptor.hpp>
// #include <feel/feeldiscr/createsubmesh.hpp>
using namespace nanoflann;
#include <random>
#include <cmath>
#include <feel/feelmesh/kdtree.hpp>

namespace Feel {
// Compute random direction, uniformly distributed on the sphere or on a circle
std::vector<double> get_random_direction(std::vector<double> &random_direction)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();              

    std::default_random_engine generator1(seed);
    std::default_random_engine generator2(seed1);
    std::uniform_real_distribution<double> xi2(0,1);
    std::uniform_real_distribution<double> xi1(0,1);
    if(random_direction.size()==3)
    {
        double phi = 2.*M_PI*xi1(generator1);
        double theta = math::asin(math::sqrt(xi2(generator2)));
        random_direction[0]=math::sin(theta)*math::cos(phi);
        random_direction[1]=math::sin(theta)*math::sin(phi);
        random_direction[2]=math::cos(theta);
    }
    else if (random_direction.size()==2)
    {
        double phi = 2.*M_PI*xi1(generator1);
        random_direction[0]=math::cos(phi);
        random_direction[1]=math::sin(phi);
    }
    else
    {
        throw std::logic_error( "Wrong dimension " + std::to_string(random_direction.size()) + " for the random direction" );
    }    

    return random_direction;
}   

// 3D case
// Compute the sum of the areas of three triangles
double element_area(Eigen::VectorXd const& point,Eigen::VectorXd const& el_p1,Eigen::VectorXd const& el_p2,Eigen::VectorXd const& el_p3)
{
    double area = 0.;            
    if(point.size()==3)
    {
        Eigen::Vector3d point_3d,el_p1_3d,el_p2_3d,el_p3_3d;
        point_3d << point[0], point[1],point[2];
        el_p1_3d << el_p1[0], el_p1[1],el_p1[2];
        el_p2_3d << el_p2[0], el_p2[1],el_p2[2];
        el_p3_3d << el_p3[0], el_p3[1],el_p3[2];

        auto v1 = (point_3d-el_p1_3d).cross(point_3d-el_p2_3d);
        auto v2 = (point_3d-el_p2_3d).cross(point_3d-el_p3_3d);
        auto v3 = (point_3d-el_p3_3d).cross(point_3d-el_p1_3d);
        area = math::sqrt(v1.dot(v1))/2. + math::sqrt(v2.dot(v2))/2. + math::sqrt(v3.dot(v3))/2. ;            
        return area;
    }
    else
    {
        throw std::logic_error( "Wrong area calculation for the random direction" );
        return -1.;
    }
}

// 2D case
// Compute the sum of the lengths of two segments
double element_area(Eigen::VectorXd const& point,Eigen::VectorXd const& el_p1,Eigen::VectorXd const& el_p2)
{     
    auto vector1 = point-el_p1;
    auto vector2 = point-el_p2;     

    return math::sqrt(vector1.dot(vector1)) + math::sqrt(vector2.dot(vector2));
}
// 3D case
// Compare the area of the 2d simplex V1V2V3 (as sum of 3 subtriangles V_iV_jB) and the 2d triangle
// created by the intersection P of the ray with the plane of  V1V2V3 (as sum of V_iV_jP)
bool isOnSurface(Eigen::VectorXd const &point,Eigen::VectorXd const &el_p1,Eigen::VectorXd const &el_p2,Eigen::VectorXd const &el_p3)
{
    auto c = (el_p1+el_p2+el_p3)/3.;
    auto elem_area = element_area(c, el_p1,el_p2,el_p3); 
    auto area = element_area(point, el_p1,el_p2,el_p3);             
    if ((area-elem_area)<1e-5)
        return true;
    else
        return false;
}

// 2D case, test if a point is on a segment
bool isOnSurface(Eigen::VectorXd const &point,Eigen::VectorXd const &el_p1,Eigen::VectorXd const &el_p2)
{
    auto v1 = (el_p2-el_p1);
    auto v2 = (point-el_p1);

    auto cross = v2[1]*v1[0]-v1[1]*v2[0];

    if(math::abs(cross)>1e-6)
        return false;
    else
    {
        if (math::abs(el_p2[0]-el_p1[0]) >= math::abs(el_p2[1]-el_p1[1]))
        {
            return el_p2[0]-el_p1[0] > 0 ? 
            el_p1[0] <= point[0] && point[0] <= el_p2[0] :
            el_p2[0] <= point[0] && point[0] <= el_p1[0];
        }
        else
        {
            return el_p2[1]-el_p1[1] > 0 ? 
            el_p1[1] <= point[1] && point[1] <= el_p2[1] :
                el_p2[1] <= point[1] && point[1] <= el_p1[1];
        }
    }
}
#if 0
template <typename MeshType>
class RayTracingViewFactor : public ViewFactorBase<MeshType>
{
    typedef typename MeshType::ptrtype mesh_ptrtype;
    typedef typename MeshType::trace_mesh_ptrtype tr_mesh_ptrtype;
    typedef typename MeshType::face_type face_type;
    typedef std::vector<std::vector<double> > my_vector_of_vectors_t;   
    typedef typename matrix_node<double>::type matrix_node_type;
    typedef KDTreeVectorOfVectorsAdaptor< my_vector_of_vectors_t, double >  KdTree_type;

public:
    struct Ray 
        { 
            public: 
                //using vec_t = Eigen::Vector3d;
                using vec_t = Eigen::VectorXd;
                Ray(const vec_t &orig, const vec_t &dir) : 
                    origin(orig),
                    dir(dir) 
                {}     
                vec_t origin, dir;  // ray origin and dir 
        };
    using value_type = double;
    RayTracingViewFactor() = default;
    RayTracingViewFactor( std::vector<std::string> const& list_of_bdys ) 
        : 
        list_of_bdys( list_of_bdys )
        {}
    RayTracingViewFactor( const RayTracingViewFactor& ) = default;
    RayTracingViewFactor( RayTracingViewFactor&& ) = default;
    RayTracingViewFactor& operator=( const RayTracingViewFactor& ) = default;
    RayTracingViewFactor& operator=( RayTracingViewFactor&& ) = default;
    ~RayTracingViewFactor() = default;
    void init( std::vector<std::string> const& list_of_bdys ) { ViewFactorBase<MeshType>::init( list_of_bdys ); }    
    tr_mesh_ptrtype M_submesh;
    mesh_ptrtype M_mesh;
    Eigen::MatrixXd M_view_factors_matrix;
    std::vector<std::string> M_markers_string,list_of_bdys;
    std::vector<int> M_markers_int;
    int M_Nrays;
    std::vector<int> M_point_indices;
    KdTree_type* M_kdtree_boundariesptr;
    Eigen::VectorXd M_view_factor_row;
    bool M_inward;
    KDTree M_kd_tree_feelpp;
    

    RayTracingViewFactor(mesh_ptrtype mesh,int Nrays,std::vector<std::string> surface_markers={},bool inward=true)
    {
        M_mesh = mesh;
        M_Nrays = Nrays;
        M_inward= inward;
        if(surface_markers.empty())
            M_submesh = createSubmesh(_mesh=mesh,_range=boundaryfaces(mesh),_update=MESH_ADD_ELEMENTS_INFO);    
        else
            M_submesh = createSubmesh(_mesh=mesh,_range=markedfaces(mesh,surface_markers),_update=MESH_ADD_ELEMENTS_INFO);

        my_vector_of_vectors_t  boundary_points; 
        int boundary_pts_index = 0;

        boundary_points.resize(M_submesh->numPoints());
        M_point_indices.resize(M_submesh->numPoints());

        for(auto point : M_submesh->points())
        {   
            boundary_points[boundary_pts_index].resize(int(M_mesh->dimension()));		
            for(int i=0;i<int(M_mesh->dimension());i++)
                boundary_points[boundary_pts_index][i] = point.second.node()[i];
            
            M_point_indices[boundary_pts_index] = point.first;
            M_kd_tree_feelpp.addPoint(point.second.node(),point.first);

            boundary_pts_index++;
            
        }    
        //M_kd_tree_feelpp=build_tree(M_kd_tree_feelpp.points().begin(),M_kd_tree_feelpp.points().end(),0,-1);
        std::cout <<   "number of points" << M_submesh->numPoints() << std::endl;                                   
        M_kdtree_boundariesptr = new KDTreeVectorOfVectorsAdaptor(M_mesh->dimension() /*dim*/, boundary_points, 10 /* max leaf */ );                              


        M_view_factors_matrix.resize(M_submesh->markerNames().size(),M_submesh->markerNames().size());
        M_view_factor_row.resize(M_submesh->markerNames().size());
        for(auto &m : M_submesh->markerNames())
        {        
            //std::cout << m.first << m.second[0] << m.second[1] << m.second[2] << m.second[3] << std::endl;
            if(m.second[1]==M_submesh->dimension())
            {    
                M_markers_string.push_back(m.first);
                M_markers_int.push_back(m.second[0]);
            }
        }
        this->init(M_markers_string);
    }


    Eigen::VectorXd get_random_point(matrix_node_type const& element_points)
    {            
        int dimension;

        dimension = column(element_points, 0).size();
        if(dimension==3)
        {
            // Choose points in a parallelogram, uniformly
            Eigen::VectorXd p1(dimension),p2(dimension),p3(dimension),v(dimension),u(dimension),p(dimension);
            for(int i=0;i<3;i++)
            {
                p1(i)=column(element_points, 0)[i];
                p2(i)=column(element_points, 1)[i];
                p3(i)=column(element_points, 2)[i];                
            }
            v = p2-p1;
            u = p3-p1;
            while(true)
            {
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();             
                unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();             
                std::default_random_engine generator1(seed),generator2(seed1);
                std::uniform_real_distribution<double> xi1(0,1),xi2(0,1);
                double s = xi1(generator1);
                double t = xi2(generator2);
                // If the point is on the left of the diagonal, keep it, else take the symmetric one
                bool in_triangle = (s + t <= 1);
                if(in_triangle)
                    p = p1 + s * u + t * v;  
                else 
                    p= p1 + (1 - s) * u + (1 - t) * v;

                if (isOnSurface(p,p1,p2,p3))
                    return p;
            }
        }
        else if(dimension==2)
        {
            Eigen::VectorXd p1(dimension),p2(dimension),v(dimension),p(dimension);
            for(int i=0;i<2;i++)
            {
                p1(i)=column(element_points, 0)[i];
                p2(i)=column(element_points, 1)[i];                
            }
            v = p2-p1;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();              
            std::default_random_engine generator1(seed);
            std::uniform_real_distribution<double> xi1(0,1);
            double s = xi1(generator1);
            p = p1 + s * v;                  
            if (isOnSurface(p,p1,p2))
                return p;            
        }
        else
        {
            Eigen::VectorXd p1(dimension);
            for(int i=0;i<dimension;i++)
            {
                p1(i)=column(element_points, 0)[i];
            }
            throw std::logic_error( "Problem for the computation of the random point" );
            return p1;
        }

    }        
    std::pair<bool,Eigen::VectorXd> check_intersection_with_triangle(matrix_node_type const& nodes, Ray &ray)
    {
        Eigen::VectorXd p1(3),p2(3),p3(3),n1(3),w(3),w_(3);
        Eigen::Matrix3d m(3,3);

        for ( int i = 0; i < nodes.size1(); ++i )
        {p1( i ) = nodes(i,0);
            p2( i ) = nodes(i,1);
            p3( i ) = nodes(i,2);}
        Eigen::VectorXd origin(3),direction(3);
        direction << ray.dir[0],ray.dir[1],ray.dir[2];
        origin << ray.origin[0],ray.origin[1],ray.origin[2];


        // tell the compiler that the VectorXd is actually a 3d vector to allow the cross product to be done
        n1 = ((p2-p1).head<3>()).cross((p3-p1).head<3>());
        //n1_normalised = n1/n1.norm();
        double n_dot_dir = direction.dot(n1);
        // Ray is parallel to the triangle's plane
        if (math::abs(n_dot_dir)<1e-6)
        {
            return std::make_pair(false,origin);
        }
        double d = -p1.dot(n1);
        double t = -(origin.dot(n1)+d)/n_dot_dir;
        if(t<=0) // require that the intersection point does not belong to the plane of the triangle
        {
            return std::make_pair(false,origin);
        }
        w = origin + direction * t;  
        if(isOnSurface(w,p1,p2,p3))
        {
            return std::make_pair(true,w);
        }
        else
        {
            return std::make_pair(false,w);
        }               
        
    }

    std::pair<bool,Eigen::VectorXd> check_intersection_with_segment(matrix_node_type const& nodes, Ray &ray)
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
            return std::make_pair(false,w_);

        double t1 = (v2[0]*v1[1]-v2[1]*v1[0])/ dot;
        double t2 = v1.dot(v3) / dot;
        //std::cout << t1 << " " << t2 << std::endl;
        if (t1 >= 0.0 && (t2 >= 0.0 && t2 <= 1.0))
        {
            w_[0] =  ray.origin[0] + ray.dir[0]*t1;
            w_[1] =  ray.origin[1] + ray.dir[1]*t1;
            return std::make_pair(true,w_);
        }

        return std::make_pair(false,w_);            
        
    }

    std::tuple<double,double,std::vector<double>> check_if_intersection_with_tree(const Ray &rayon)
    {
        // FROM "Real time collision detection", Ericson C, p.180
        double tmin = 0.;
        double tmax = INT_MAX*1;
        std::vector<double> first_intersection_point(3); // the intersection with the first slab of the bounding box
        for (int i = 0; i < 3; i++) {       
            if (std::abs(rayon.dir[i]) < 1e-6) 
            {
                // Ray is parallel to slab. No hit if origin not within slab
                if (rayon.origin[i] < M_kdtree_boundariesptr->index->root_bbox[i].low || rayon.origin[i] > M_kdtree_boundariesptr->index->root_bbox[i].high) 
                    return std::make_tuple(-INT_MAX,INT_MAX,first_intersection_point);
            }
            else 
            {        
                // Compute intersection t value of ray with near and far plane of slab                
                double t1 = (M_kdtree_boundariesptr->index->root_bbox[i].low - rayon.origin[i]) * 1. / rayon.dir[i];
                double t2 = (M_kdtree_boundariesptr->index->root_bbox[i].high - rayon.origin[i]) * 1. / rayon.dir[i];
                // Make t1 be the intersection with the near plane, t2 with the far plane
                if (t1 > t2) std::swap(t1, t2);
                // Compute the intersection of slab intersection intervals
                if (t1 > tmin) tmin = t1;
                if (t2 < tmax) tmax = t2; 
                // Exit with no collision as soon as slab intersection becomes empty if (tmin > tmax) return 0;
                if (tmin > tmax) return std::make_tuple(-INT_MAX,INT_MAX,first_intersection_point);
            }
        }             
        for(int i=0;i<3;i++)
        {
            first_intersection_point[i]= rayon.origin[i] + rayon.dir[i] * tmin;
        } 
        return std::make_tuple(tmin,tmax,first_intersection_point);
    }

    std::vector<int> get_intersection_indices(const Ray &rayon, std::vector<int> &leaf_indices,bool test=false)
    {              
        auto [tmin, tmax, first_intersection_point] = check_if_intersection_with_tree(rayon);
        if(tmin==-INT_MAX && tmax==INT_MAX)
        {
            return leaf_indices;
        }
        // Starting from the root_node, find the leaves where an intersection is still possible                  
        if(!test)
          VisitNodes(M_kdtree_boundariesptr->index->root_node, first_intersection_point, rayon, tmax-tmin,leaf_indices);      
        return leaf_indices;
    }

    mesh_ptrtype mesh() {return M_mesh;}
    std::vector<std::string> markerNames(){return M_markers_string;}
    void compute()
    {
        for(int i=0;i<M_markers_string.size();i++)
            {
                // Surface with marker i launches rays to surface of marker j
                // this->vf_.row(i)=computeViewFactor(M_markers_string[i]); 
                std::cout << computeViewFactor_feel(M_markers_string[i]) << std::endl;
            }                        
    };

    Eigen::VectorXd computeViewFactor(std::string marker)
    {
        int dim = M_mesh->dimension();
        M_view_factor_row.setZero();        
        auto ray_submesh = createSubmesh(_mesh=M_submesh,_range=markedelements(M_submesh,marker));               
        std::vector<double> random_direction(dim);            
        for(auto const &el: ray_submesh->elements())
        {
            for(int i=0;i<M_Nrays;i++)
            {              
                // Construct the ray emitting from a random point of the element
                auto random_origin = get_random_point(el.second.vertices());                   

                get_random_direction(random_direction);  

                Eigen::VectorXd rand_dir(dim); 
                Eigen::VectorXd p1(dim),p2(dim),p3(dim);
                if(dim==3)
                {
                    for(int i=0;i<3;i++)
                    {
                        p1(i)=column(el.second.vertices(), 0)[i];
                        p2(i)=column(el.second.vertices(), 1)[i];
                        p3(i)=column(el.second.vertices(), 2)[i];                        
                        rand_dir(i) = random_direction[i];
                    }                               
                    auto element_normal = ((p2-p1).head<3>()).cross((p3-p1).head<3>());
                    element_normal.normalize();
                    if(rand_dir.dot(element_normal)>0.)
                        rand_dir *= -1.;                        
                }
                else if(dim==2)
                {
                    for(int i=0;i<dim;i++)
                    {
                        p1(i)=column(el.second.vertices(), 0)[i];
                        p2(i)=column(el.second.vertices(), 1)[i];                        
                        rand_dir(i) = random_direction[i];
                    }                               
                    Eigen::VectorXd element_normal(dim); 
                    element_normal << -(p2[1]-p1[1]),p2[0]-p1[0];
                    element_normal.normalize();
                    if(rand_dir.dot(element_normal)>0.)
                        rand_dir *= -1.; 
                }

                Ray ray(random_origin,rand_dir);

                // Find the intersection indices (=point labels) between the tree and the ray
                std::vector<int> leaf_indices={}; 
                get_intersection_indices(ray,leaf_indices);                                            
                
                // Check if the ray intersects an element
                auto intersection= check_intersection_with_elements(ray,leaf_indices,M_point_indices);
                // Find the marker's index corresponding to the element being intersected by the ray
                auto index_view_factor = std::find(M_markers_int.begin(), M_markers_int.end(), intersection.second);
                
                if(intersection.first)
                    M_view_factor_row(std::distance(M_markers_int.begin(),index_view_factor))++;
            }  
        }              
        M_view_factor_row /=(1.*M_Nrays*ray_submesh->numElements());

        std::cout << M_Nrays*ray_submesh->numElements() << std::endl;
        
        return M_view_factor_row;
    }

    Eigen::VectorXd computeViewFactor_feel(std::string marker)
    {
        int dim = M_mesh->dimension();
        M_view_factor_row.setZero();        
        auto ray_submesh = createSubmesh(_mesh=M_submesh,_range=markedelements(M_submesh,marker));               
        std::vector<double> random_direction(dim);            
        for(auto const &el: ray_submesh->elements())
        {
            for(int i=0;i<M_Nrays;i++)
            {              
                // Construct the ray emitting from a random point of the element
                auto random_origin = get_random_point(el.second.vertices());                   

                get_random_direction(random_direction);  

                Eigen::VectorXd rand_dir(dim); 
                Eigen::VectorXd p1(dim),p2(dim),p3(dim);
                if(dim==3)
                {
                    for(int i=0;i<dim;i++)
                    {
                        p1(i)=column(el.second.vertices(), 0)[i];
                        p2(i)=column(el.second.vertices(), 1)[i];
                        p3(i)=column(el.second.vertices(), 2)[i];                        
                        rand_dir(i) = random_direction[i];
                    }                               
                    auto element_normal = ((p2-p1).head<3>()).cross((p3-p1).head<3>());
                    element_normal.normalize();
                    if(rand_dir.dot(element_normal)>0.)
                        rand_dir *= -1.;                        
                }
                else if(dim==2)
                {
                    for(int i=0;i<dim;i++)
                    {
                        p1(i)=column(el.second.vertices(), 0)[i];
                        p2(i)=column(el.second.vertices(), 1)[i];                        
                        rand_dir(i) = random_direction[i];
                    }                               
                    Eigen::VectorXd element_normal(dim); 
                    element_normal << -(p2[1]-p1[1]),p2[0]-p1[0];
                    element_normal.normalize();
                    if(rand_dir.dot(element_normal)>0.)
                        rand_dir *= -1.; 
                }

                Feel::Ray ray(random_origin,rand_dir);
                Ray ray2(random_origin,rand_dir);

                // Find the intersection indices (=point labels) between the tree and the ray
                M_kd_tree_feelpp.ray_search(ray);
                std::vector<int> leaf_indices_feelpp={};

                for(auto pt : M_kd_tree_feelpp.pointsNearNeighbor() )                    
                    leaf_indices_feelpp.push_back(boost::get<3>(pt));
                
                std::cout << leaf_indices_feelpp << std::endl;                                            
                
                // Check if the ray intersects an element
                auto [intersection,marker_index]= check_intersection_with_elements_feel(ray2,leaf_indices_feelpp,M_point_indices);
                // Find the marker's index corresponding to the element being intersected by the ray
                auto index_view_factor = std::find(M_markers_int.begin(), M_markers_int.end(), marker_index);
                
                if(intersection)
                    M_view_factor_row(std::distance(M_markers_int.begin(),index_view_factor))++;
            }  
        }              
        M_view_factor_row /=(1.*M_Nrays*ray_submesh->numElements());

        std::cout << M_Nrays*ray_submesh->numElements() << std::endl;
        
        return M_view_factor_row;
    }

    // Check intersection with elements of the mesh, by looking at those connected to each point
    std::pair<bool,int> 
    check_intersection_with_elements(Ray &ray, std::vector<int> &leaf_indices, std::vector<int> &point_indices)
    {
        bool intersection_ = false;

        for(int j=0;j<leaf_indices.size();j++)
        {
            auto current_mesh_node = M_kdtree_boundariesptr->index->vAcc[leaf_indices.data()[j]];        
            // find current_mesh_node in mesh
            auto mesh_node = M_submesh->point(point_indices[current_mesh_node]);         

            for(auto &el : mesh_node.elements())
            {                               
                auto [intersection_happened,intersection_point] =check_intersection_with_triangle(M_submesh->element(el.first).vertices(),ray);
                if (intersection_happened)
                {
                    intersection_ =  true;                        
                    auto marker_id = M_submesh->element(el.first).marker().value();                        
                    std::cout << "intersection_point" << std::endl;
                    std::cout << intersection_point << std::endl;
                    //std::cout << M_markers_string[marker_id] << std::endl;
                    return std::make_pair(intersection_,marker_id);
                }             
            }  
        }
        std::cout << "no intersection for Ray " << ray.origin.data()[0] << " "<< ray.origin.data()[1] << " "<< ray.origin.data()[2] 
        << " " << ray.dir.data()[0] <<" " << ray.dir.data()[1] << " " << ray.dir.data()[2] << std::endl;
        return std::make_pair(intersection_,0);
    }

    std::pair<bool,int> 
    check_intersection_with_elements_feel(Ray &ray, std::vector<int> &leaf_indices, std::vector<int> &point_indices)
    {
        bool intersection_ = false;

        for(int j=0;j<leaf_indices.size();j++)
        {            
            // find current_mesh_node in mesh
            auto mesh_node = M_submesh->point(leaf_indices[j]);   
            // std::cout << "Submesh point "<< M_submesh->point(leaf_indices[j]) << "Submesh point indices" << mesh_node << std::endl;

            for(auto &el : mesh_node.elements())
            {                               
                auto [intersection_happened,intersection_point] =check_intersection_with_triangle(M_submesh->element(el.first).vertices(),ray);
                if (intersection_happened)
                {
                    intersection_ =  true;                        
                    auto marker_id = M_submesh->element(el.first).marker().value();                        
                    std::cout << "intersection_point" << std::endl;
                    std::cout << intersection_point << std::endl;
                    // std::cout << M_markers_string[marker_id] << std::endl;
                    return std::make_pair(intersection_,marker_id);
                }             
            }  
        }
        std::cout << "no intersection for Ray " << ray.origin.data()[0] << " "<< ray.origin.data()[1] << " "<< ray.origin.data()[2] 
        << " " << ray.dir.data()[0] <<" " << ray.dir.data()[1] << " " << ray.dir.data()[2] << std::endl;
        return std::make_pair(intersection_,0);
    }

    // // Visit the kdtree and look for the leaves closer to the intersection point of the ray
    // // FROM "Real time collision detection", Ericson C, p.323
    void VisitNodes2(KdTree_type::index_t::NodePtr &node, std::vector<double> &point,Ray const& rayon, double tmax,std::vector<int> &leaf_indices)
    {               
        if(node == nullptr)
            return;
        /* If "node" is a leaf node, collect the indices of the associated points. */
        if ((node->child1 == nullptr) && (node->child2 == nullptr))
        {               
            std::cout << "Enter in return" << std::endl;                 
            for(int j=node->node_type.lr.left; j<node->node_type.lr.right;j++)        
                leaf_indices.push_back(j);
            // std::cout << leaf_indices << std::endl;
            return;
        }
            
        /* Choose the child branch to be taken first */
        int idx = node->node_type.sub.divfeat; // The splitting dimension x=0, y=1, z=2
        
        double val  = point[idx]; // The [idx] coordinate of the intersection point "point"                

        std::vector<KdTree_type::index_t::NodePtr>    bestChild;
        KdTree_type::index_t::NodePtr    otherChild;

        int first = (val > node->node_type.sub.divhigh);
        if(first) // if this is true, choose the right child node
        {
            std::cout << "Enter in first" << std::endl;
            std::cout << "idx " << idx << std::endl;
            bestChild  = {node->child2};
            otherChild = {node->child1};  
            
        }
        else if((val < node->node_type.sub.divlow))
        {
            std::cout << "Enter in else first" << std::endl;
            std::cout << "idx " << idx << std::endl;
            bestChild  = {node->child1};
            otherChild = {node->child2};  
        }
        else
        {
            std::cout << "Enter in else else first" << std::endl;
            std::cout << "idx " << idx << std::endl;
            bestChild  = {node->child1,node->child2};
            otherChild = {};
        }

        if (std::abs(rayon.dir[idx]) <= 1e-6) {

            for(int j=0;j<bestChild.size();j++)
            {
                std::cout << "Enter in rayon dir first for" << std::endl;
                std::cout << "idx " << idx << std::endl;
                VisitNodes(bestChild[j], point, rayon, tmax,leaf_indices);  
            }        
                        
        }
            else {
            // Find t value for intersection between segment and rightmost split plane
            // (nanoflann uses two split planes)
            double t = (node->node_type.sub.divhigh - val) / rayon.dir[idx];                
            if ((t>=tmax))//  && rayon.dir[idx]>0) || (rayon.dir[idx]<0 && t<tmax)) // if the segment cuts the right split plane outside the bounding box, try with the left plane
            {
                t = (node->node_type.sub.divlow - val) / rayon.dir[idx];
                bestChild  = {node->child1,node->child2};
                otherChild = {};
            }
            // Test if line segment intersects splitting plane
            if (0.0 <= t && t < tmax) {
                // If it does, traverse the near child first, then the far one
                for(int j=0;j<bestChild.size();j++)
                {
                    std::cout << "Enter in after" << std::endl;
                    std::cout << "idx " << idx << std::endl;
                    VisitNodes(bestChild[j], point, rayon, t,leaf_indices);
                }

            std::vector<double> first_intersection_point(3);
            for(int i=0;i<3;i++)
                    first_intersection_point[i]= point[i] + rayon.dir[i] * t;
            std::cout << "Enter in other child" << std::endl;
            std::cout << "idx " << idx << std::endl;
            VisitNodes(otherChild, first_intersection_point, rayon, tmax - t,leaf_indices); } 
            else {                    
                // If not, just traverse the near child
            std::cout << "end with this" << std::endl;
            std::cout << "idx " << idx << std::endl;
            for(int j=0;j<bestChild.size();j++)
                VisitNodes(bestChild[j], point, rayon, tmax,leaf_indices); 
            }
                    
        }

    }
    void VisitNodes(KdTree_type::index_t::NodePtr &node, std::vector<double> &point,Ray const& rayon, double tmax,std::vector<int> &leaf_indices)
    {               
        if(node == nullptr)
            return;
        /* If "node" is a leaf node, collect the indices of the associated points. */
        if ((node->child1 == nullptr) && (node->child2 == nullptr))
        {               
            // std::cout << "Enter in return" << std::endl;                 
            for(int j=node->node_type.lr.left; j<node->node_type.lr.right;j++)        
                leaf_indices.push_back(j);
            // std::cout << leaf_indices << std::endl;
            return;
        }
            
        /* Choose the child branch to be taken first */
        int idx = node->node_type.sub.divfeat; // The splitting dimension x=0, y=1, z=2
        
        double val  = point[idx]; // The [idx] coordinate of the intersection point "point"                

        std::vector<KdTree_type::index_t::NodePtr>    bestChild;
        KdTree_type::index_t::NodePtr    otherChild;

        int first = (val > (node->node_type.sub.divhigh+node->node_type.sub.divlow)*0.5);
        if(first) // if this is true, choose the right child node
        {
            // std::cout << "Enter in first" << std::endl;
            // std::cout << "idx " << idx << std::endl;
            bestChild  = {node->child2};
            otherChild = {node->child1};  
            
        }
        else
        {
            // std::cout << "Enter in else first" << std::endl;
            // std::cout << "idx " << idx << std::endl;
            bestChild  = {node->child1};
            otherChild = {node->child2};  
        }
        
        if (std::abs(rayon.dir[idx]) <= 1e-6) {

            for(int j=0;j<bestChild.size();j++)
            {
                // std::cout << "Enter in rayon dir first for" << std::endl;
                // std::cout << "idx " << idx << std::endl;
                VisitNodes(bestChild[j], point, rayon, tmax,leaf_indices);  
            }        
                        
        }
            else {
            // Find t value for intersection between segment and rightmost split plane
            // (nanoflann uses two split planes)
            double t = (0.5*(node->node_type.sub.divhigh+node->node_type.sub.divlow) - val) / rayon.dir[idx];                  
            // Test if line segment intersects splitting plane
            if (0.0 <= t && t < tmax) {
                // If it does, traverse the near child first, then the far one
                for(int j=0;j<bestChild.size();j++)
                {
                    // std::cout << "Enter in after" << std::endl;
                    // std::cout << "idx " << idx << std::endl;
                    VisitNodes(bestChild[j], point, rayon, t,leaf_indices);
                }
            std::vector<double> first_intersection_point(3);
            for(int i=0;i<3;i++)
                    first_intersection_point[i]= point[i] + rayon.dir[i] * t;
            // std::cout << "Enter in other child" << std::endl;
            // std::cout << "idx " << idx << std::endl;
            VisitNodes(otherChild, first_intersection_point, rayon, tmax - t,leaf_indices); } 
            else {                    
                // If not, just traverse the near child
            // std::cout << "end with this" << std::endl;
            // std::cout << "idx " << idx << std::endl;
            for(int j=0;j<bestChild.size();j++)
                VisitNodes(bestChild[j], point, rayon, tmax,leaf_indices); 
            }
                    
        }

    }

};  
#endif
} // namespace Feel
