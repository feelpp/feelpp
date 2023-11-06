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

#include <random>
#include <cmath>
#include <iostream>
#include <future>

#include <feel/feelviewfactor/viewfactorbase.hpp>
#include <feel/feelmesh/bvh.hpp>
#include <feel/feelvf/vf.hpp>
#include <nanoflann.hpp>
#include <feel/feelviewfactor/kdtreevectorofvectorsadaptor.hpp>

using namespace nanoflann;


namespace Feel {

// Compute random direction, uniformly distributed on the sphere or on a circle
void getRandomDirection(std::vector<double> &random_direction, std::mt19937 & M_gen, std::mt19937 & M_gen2,Eigen::VectorXd normal)
{
    std::uniform_real_distribution<double> xi2(0,1);
    std::uniform_real_distribution<double> xi1(0,1);

    int size = random_direction.size();
    Eigen::VectorXd z_axis(size), crossProd1(size), crossProd2(size),crossProd3(size), direction(size);
    Eigen::MatrixXd matrix1(size,size), matrix2(size,size);
    z_axis << 0,0,1;

    if(random_direction.size()==3)
    {
        double phi = 2.*M_PI*xi1(M_gen);
        double theta = math::asin(math::sqrt( xi2(M_gen2)));
        random_direction[0]=math::sin(theta)*math::cos(phi);
        random_direction[1]=math::sin(theta)*math::sin(phi);
        random_direction[2]=math::cos(theta);

        if( math::abs(normal.dot(z_axis)) < 1)  // if the normal is not aligned with z-axis,
                                                // rotate random_direction according to the rotation bringing normal to z-axis
        {
            direction << random_direction[0],random_direction[1],random_direction[2];
            crossProd1 = ((z_axis).head<3>()).cross((normal).head<3>());
            crossProd2 = ((crossProd1).head<3>()).cross((z_axis).head<3>());
            crossProd3 = ((crossProd1).head<3>()).cross((normal).head<3>());

            matrix1.col(0) = z_axis;
            matrix1.col(1) = crossProd2;
            matrix1.col(2) = crossProd1;

            matrix2.row(0) = normal;
            matrix2.row(1) = crossProd3;
            matrix2.row(2) = crossProd1;


            direction = matrix1 * (matrix2 * direction);

            random_direction.resize(direction.size());
            Eigen::VectorXd::Map(&random_direction[0], direction.size()) = direction;
        }
    }
    else if (random_direction.size()==2)
    {
        double phi = M_PI*xi1(M_gen);
        random_direction[0]=math::cos(phi);
        random_direction[1]=math::sin(phi);
    }
    else
    {
        throw std::logic_error( "Wrong dimension " + std::to_string(random_direction.size()) + " for the random direction" );
    }

    //return random_direction;
}

// 3D case
// Compute the sum of the areas of three triangles
double elementArea(Eigen::VectorXd const& point,Eigen::VectorXd const& el_p1,Eigen::VectorXd const& el_p2,Eigen::VectorXd const& el_p3)
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
double elementArea(Eigen::VectorXd const& point,Eigen::VectorXd const& el_p1,Eigen::VectorXd const& el_p2)
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
    auto elem_area = elementArea(c, el_p1,el_p2,el_p3);
    auto area = elementArea(point, el_p1,el_p2,el_p3);
    if (math::abs(area-elem_area)<1e-5)
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
template <typename MeshType>
class RayTracingViewFactor : public ViewFactorBase<MeshType>
{
    using super = ViewFactorBase<MeshType>;
    typedef typename MeshType::ptrtype mesh_ptrtype;

    using tr_mesh_ptrtype = trace_mesh_ptr_t<MeshType>;
    using tr_mesh_t = trace_mesh_t<MeshType>;

    typedef typename MeshType::face_type face_type;
    typedef typename matrix_node<double>::type matrix_node_type;
    using bvh_type = BVH<element_t<tr_mesh_t>>;
    using bvh_ray_type = typename bvh_type::ray_type;
public:
    using value_type = double;
    RayTracingViewFactor() = default;
    RayTracingViewFactor( std::vector<std::string> const& list_of_bdys )
        :
        list_of_bdys( list_of_bdys )
        {}
    RayTracingViewFactor(mesh_ptrtype mesh, nl::json const& specs )
        :
        super(mesh,specs)
        {
            M_Nrays = specs["viewfactor"]["Nrays"].get<double>() ;
            M_Nthreads = specs["viewfactor"]["Nthreads"].get<double>() ;
            if(this->list_of_bdys_.empty())
                M_submesh = createSubmesh(_mesh=mesh,_range=boundaryfaces(mesh),_update=MESH_ADD_ELEMENTS_INFO);
            else
                M_submesh = createSubmesh(_mesh=mesh,_range=markedfaces(mesh,this->list_of_bdys_),_update=MESH_ADD_ELEMENTS_INFO);

            M_bvh_tree = boundingVolumeHierarchy(_range=elements(M_submesh));

            M_view_factor_row.resize(this->list_of_bdys_.size());

            for(auto &m : this->list_of_bdys_)
            {
                M_markers_string.push_back(m);
                for(auto &m_mesh : M_submesh->markerNames())
                {
                    // std::cout << m_mesh.first << " "<< m_mesh.second[0] << " "<< m_mesh.second[1] << " "<< m_mesh.second[2] << " "<< m_mesh.second[3] << std::endl;
                    if(m_mesh.second[1]==M_submesh->dimension() && m_mesh.first==m)
                    {
                        M_markers_int.push_back(m_mesh.second[0]);
                        break;
                    }
                }
            }
            // std::cout <<      M_view_factor_row << std::endl;
            std::cout <<      M_markers_string << std::endl;
            std::cout <<      M_markers_int << std::endl;

            std::random_device rd;  // Will be used to obtain a seed for the random number engine
            std::random_device rd2;  // Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
            std::mt19937 gen2(rd2()); // Standard mersenne_twister_engine seeded with rd()
            gen.seed(std::chrono::high_resolution_clock::now()
                                   .time_since_epoch()
                                   .count());
            gen2.seed(std::chrono::high_resolution_clock::now()
                                   .time_since_epoch()
                                   .count());

            M_gen=gen;
            M_gen2=gen2;
        }
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
    int M_Nthreads;
    std::vector<int> M_point_indices;
    Eigen::VectorXd M_view_factor_row;
    std::shared_ptr<bvh_type> M_bvh_tree;
    std::random_device M_rd;  // Will be used to obtain a seed for the random number engine
    std::random_device M_rd2;  // Will be used to obtain a seed for the random number engine
    std::mt19937 M_gen; // Standard mersenne_twister_engine seeded with rd()
    std::mt19937 M_gen2; // Standard mersenne_twister_engine seeded with rd2()

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
                unsigned seed2 = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                unsigned seed3 = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                std::default_random_engine generator3(seed2),generator4(seed3);
                std::uniform_real_distribution<double> xi1(0,1),xi2(0,1);
                double s = xi1(generator3);
                double t = xi2(generator4);
                // If the point is on the left of the diagonal, keep it, else take the symmetric one
                bool in_triangle = (s + t <= 1);
                if(in_triangle)
                    p = p1 + s * u + t * v;
                else
                    p= p1 + (1 - s) * u + (1 - t) * v;

                if (isOnSurface(p,p1,p2,p3))
                    return p;
                else
                {
                    throw std::logic_error("Point not on triangle, but it must be");
                    return p1;
                }
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
            std::default_random_engine generator3(seed);
            std::uniform_real_distribution<double> xi1(0,1);
            double s = xi1(generator3);
            p = p1 + s * v;
            if (isOnSurface(p,p1,p2))
                return p;
            else
            {
                throw std::logic_error("Point not on segment, but it must be");
                return p1;
            }
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

    mesh_ptrtype mesh() {return M_mesh;}
    std::vector<std::string> markerNames(){return M_markers_string;}
    void compute(bool elementwise=false)
    {
        if(this->j_["viewfactor"]["type"]=="Raytracing")
        {
            for(int i=0;i<this->list_of_bdys_.size();i++)
            {
                // Surface with marker i launches rays to surface of marker j
                this->vf_.row(i)=computeViewFactor_bvh(this->list_of_bdys_[i]);
                std::cout << this->vf_.row(i) << std::endl;
            }
        }
        else
        {
            std::cout<< "Viewfactor not computed, error in the algorithm choice" <<std::endl;
        }
    };

    Eigen::VectorXd computeViewFactor_bvh(std::string const& marker)
    {
        int dim = M_mesh->dimension();
        M_view_factor_row.setZero();

        auto ray_submesh = createSubmesh(_mesh=M_submesh,_range=markedelements(M_submesh,marker));
        std::vector<double> random_direction(dim);


        for(auto const &el : ray_submesh->elements()) // from each element of the submesh, launch M_Nrays randomly oriented
        {
            auto rays_from_element = [&](int n_rays_thread){
                Eigen::VectorXd row_vf_matrix(M_view_factor_row.size());
                row_vf_matrix.setZero();
                // for(int i=0;i<M_Nrays;i++)
                for(int i=0;i<n_rays_thread;i++)
                {

                    // Construct the ray emitting from a random point of the element
                    auto random_origin = get_random_point(el.second.vertices());

                    Eigen::VectorXd rand_dir(dim);
                    Eigen::VectorXd p1(dim),p2(dim),p3(dim),origin(3);

                    if(dim==3)
                    {
                        for(int i=0;i<dim;i++)
                        {
                            p1(i)=column(el.second.vertices(), 0)[i];
                            p2(i)=column(el.second.vertices(), 1)[i];
                            p3(i)=column(el.second.vertices(), 2)[i];
                            origin(i) = random_origin[i];
                            rand_dir(i) = random_direction[i];
                        }
                        auto element_normal = ((p3-p1).head<3>()).cross((p2-p1).head<3>());
                        element_normal.normalize();
                        getRandomDirection(random_direction,M_gen,M_gen2,element_normal);
                        for(int i=0;i<dim;i++)
                        {
                            rand_dir(i) = random_direction[i];
                        }
                        if(rand_dir.dot(element_normal)<0.) // if the ray direction is inward to the emitting surface, invert its direction
                        {
                            rand_dir(0) = rand_dir(0) * (-1.);
                            rand_dir(1) = rand_dir(1) * (-1.);
                            rand_dir(2) = rand_dir(2) * (-1.);
                        }

                    }
                    Eigen::VectorXd pQ1(3),pQ2(3),pQ3(3),pQ4(3),normal(3);

                    bvh_ray_type ray(origin,rand_dir,1e-8); // warning, put a minimal distance > 0

                    auto rayIntersectionResult = M_bvh_tree->intersect(ray) ;
                    if ( !rayIntersectionResult.empty() )
                    {
                        auto marker_index = rayIntersectionResult.front().primitive().meshEntity().marker().value();

                        // Find the marker's index corresponding to the element being intersected by the ray
                        auto index_view_factor = std::find(M_markers_int.begin(), M_markers_int.end(), marker_index);

                        // M_view_factor_row(std::distance(M_markers_int.begin(),index_view_factor))++;
                        row_vf_matrix(std::distance(M_markers_int.begin(),index_view_factor))++;
                    }
                    else{
                        std::cout << "RAY NOT INTERSECTING: Rand origin" << ray.origin() << "rand_dir" <<  ray.dir() <<std::endl;
                    }
                }
                return row_vf_matrix;
            };
                std::vector<int> n_rays_thread;
                n_rays_thread.push_back(M_Nrays - (M_Nthreads-1) * (int)(M_Nrays / M_Nthreads));
                for(int t= 1; t < M_Nthreads; ++t){
                   n_rays_thread.push_back( M_Nrays / M_Nthreads);
                }
                // Used to store the future results
                std::vector<std::future<Eigen::VectorXd>> futures;
                for(int t = 0; t < M_Nthreads; ++t){
                    // Start a new asynchronous task
                    // std::cout << n_rays_thread[t] << std::endl;
                    futures.emplace_back(std::async(std::launch::async, rays_from_element, n_rays_thread[t]));
                }


                for(std::future<Eigen::VectorXd>& f : futures){
                    // Wait for the result to be ready
                    M_view_factor_row += f.get();
                }

        }
        M_view_factor_row /=(1.*M_Nrays*ray_submesh->numElements());
        std::cout << M_Nrays*ray_submesh->numElements() << std::endl;

        auto index_marker =std::find(M_markers_string.begin(),M_markers_string.end(),marker);
        this->areas_[std::distance(M_markers_string.begin(),index_marker)] = integrate(_range=elements(ray_submesh),_expr=cst(1.)).evaluate()(0,0);

        std::cout << "Areas_" << marker << " " << this->areas_[std::distance(M_markers_string.begin(),index_marker)] << std::endl;

        return M_view_factor_row;
    }

};
}
