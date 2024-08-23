
namespace Feel
{
        // Choose a random pair of indices in the discretized azimuth and altitude vectors
    template <typename MeshType>
    void
    ShadingMask<MeshType>::getRandomDirectionSM(std::vector<double> &random_direction, std::mt19937 & M_gen, std::mt19937 & M_gen2, int& index_azimuth, int& index_altitude)
    {
        std::uniform_int_distribution<int> dist_azimuth(0,M_azimuthSize-1);
        std::uniform_int_distribution<int> dist_altitude(0,M_altitudeSize-1);

        int size = random_direction.size();

        if(random_direction.size()==3)
        {
            index_azimuth = dist_azimuth(M_gen);
            index_altitude = dist_altitude(M_gen2);
            double phi = -( M_azimuthAngles[index_azimuth] ) + M_PI*0.5 ; // recover spherical coordinate from azimuth angle
            double theta = M_PI*0.5 - M_altitudeAngles[index_altitude]; // recover spherical coordinate from altitude

            random_direction[0]=math::sin(theta)*math::cos(phi);
            random_direction[1]=math::sin(theta)*math::sin(phi);
            random_direction[2]=math::cos(theta);

        }
        else
        {
            throw std::logic_error( "Wrong dimension " + std::to_string(random_direction.size()) + " for the random direction" );
        }

    }

    // Get a random point from the surface of the triangle
    template <typename MeshType>
    Eigen::VectorXd
    ShadingMask<MeshType>::get_random_point(matrix_node_type const& element_points)
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
        else
        {
            Eigen::VectorXd p1(dimension);
            for(int i=0;i<dimension;i++)
            {
                p1(i)=column(element_points, 0)[i];
            }
            throw std::logic_error( "Problem in the computation of the random point" );
            return p1;
        }

    }

    // 3D case
    // Compute the sum of the areas of three subtriangles
    template <typename MeshType>
    double 
    ShadingMask<MeshType>::elementArea(Eigen::VectorXd const& point,Eigen::VectorXd const& el_p1,Eigen::VectorXd const& el_p2,Eigen::VectorXd const& el_p3)
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

    // 3D case
    // Compare the area of the 2d simplex V1V2V3 (as sum of 3 subtriangles V_iV_jB) and the 2d triangle
    // created by the intersection P of the ray with the plane of  V1V2V3 (as sum of V_iV_jP)
    template <typename MeshType>
    bool 
    ShadingMask<MeshType>::isOnSurface(Eigen::VectorXd const &point,Eigen::VectorXd const &el_p1,Eigen::VectorXd const &el_p2,Eigen::VectorXd const &el_p3)
    {
        auto c = (el_p1+el_p2+el_p3)/3.;

        auto elem_area = elementArea(c, el_p1,el_p2,el_p3);
        auto area = elementArea(point, el_p1,el_p2,el_p3);
        if (math::abs(area-elem_area)/area<1e-5)
            return true;
        else
            return false;
    }

}