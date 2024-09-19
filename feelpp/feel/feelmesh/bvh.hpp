// Bounding volume hierarchy

#ifndef FEELPP_MESH_BVH_HPP
#define FEELPP_MESH_BVH_HPP


#include <vector>

#include <bvh/v2/bvh.h>
#include <bvh/v2/default_builder.h>
#include <bvh/v2/stack.h>
#include <bvh/v2/tri.h>

//Addendum for GPU part
#include <bvh/v2/stream.h>
#include <bvh/v2/node.h>


//#include <bvh/v2/vec.h>
//#include <bvh/v2/ray.h>
//#include <bvh/v2/thread_pool.h>
//#include <bvh/v2/executor.h>


#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feeldiscr/mesh.hpp>

//#include <feel/feeltask/taskpu.hpp>
#include "hip/hip_runtime.h"
#include "hip/hip_runtime_api.h"


#include <assert.h>
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <optional>
#include <random>
#include <cfloat>


/*
#include <thrust/device_vector.h> 
#include <thrust/transform.h> 
#include <thrust/functional.h> 
#include <thrust/execution_policy.h>
#include <thrust/random.h>
*/


#include "thrust/device_vector.h"
#include "thrust/transform.h"
#include "thrust/functional.h"
#include "thrust/execution_policy.h"
#include "thrust/random.h"
#include "thrust/host_vector.h"
#include "thrust/device_vector.h"
#include "thrust/sort.h"

#include "thrust/generate.h"
#include "thrust/sort.h"
#include "thrust/copy.h"
#include "thrust/count.h"



#include <limits>
#include <climits>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <atomic>




/*
namespace bvhhip
{

#define SWAP(T, a, b) do { T tmp = a; a = b; b = tmp; } while (0)
	

struct Vec3 {
    float x, y, z;

    __host__ __device__
    Vec3 operator-(const Vec3& v) const {
        return {x - v.x, y - v.y, z - v.z};
    }

    __host__ __device__
    Vec3 operator+(const Vec3& v) const {
        return {x + v.x, y + v.y, z + v.z};
    }

    __host__ __device__
    Vec3 operator*(float scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

     __host__ __device__
        Vec3 operator*(const Vec3& other) const {
        return Vec3(x * other.x, y * other.y, z * other.z);
    }

    __host__ __device__
    float dot(const Vec3& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    __host__ __device__
    Vec3 cross(const Vec3& v) const {
        return {
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        };
    }

    __host__ __device__
    Vec3 () : x(0), y(0), z(0) {} 

    __host__ __device__
    Vec3 (float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}  

    __host__ __device__
    Vec3 init(float x, float y, float z) {
        Vec3 v;
        v.x = x;
        v.y = y;
        v.z = z;
        return v;
    }

    __host__ __device__
    Vec3 float3ToVec3(const float3& f) {
        return Vec3(f.x, f.y, f.z);
    }
};


struct F3Triangle {
    float3 v0, v1, v2;
};


struct Triangle {
    Vec3 v0, v1, v2;

    __host__ __device__
    bool intersect(const Vec3& rayOrigin, const Vec3& rayDir, float& t) const {
        Vec3 edge1 = v1 - v0;
        Vec3 edge2 = v2 - v0;
        Vec3 h = rayDir.cross(edge2);
        float a = edge1.dot(h);

        if (a > -1e-8 && a < 1e-8) return false; // Ray is parallel to triangle

        float f = 1.0f / a;
        Vec3 s = rayOrigin - v0;
        float u = f * s.dot(h);
        if (u < 0.0f || u > 1.0f) return false;

        Vec3 q = s.cross(edge1);
        float v = f * rayDir.dot(q);
        if (v < 0.0f || u + v > 1.0f) return false;

        t = f * edge2.dot(q);
        return t > 1e-8; // Intersection occurs
    }
};


struct Box {
    Vec3 min;
    Vec3 max;

    __host__ __device__
    bool intersect(const Vec3& rayOrigin, const Vec3& rayDir, float& t) const {
        float tMin = (min.x - rayOrigin.x) / rayDir.x;
        float tMax = (max.x - rayOrigin.x) / rayDir.x;

        if (tMin > tMax) SWAP(float,tMin, tMax);

        float tyMin = (min.y - rayOrigin.y) / rayDir.y;
        float tyMax = (max.y - rayOrigin.y) / rayDir.y;

        if (tyMin > tyMax) SWAP(float,tyMin, tyMax);

        if ((tMin > tyMax) || (tyMin > tMax))
            return false;

        if (tyMin > tMin)
            tMin = tyMin;

        if (tyMax < tMax)
            tMax = tyMax;

        float tzMin = (min.z - rayOrigin.z) / rayDir.z;
        float tzMax = (max.z - rayOrigin.z) / rayDir.z;

        if (tzMin > tzMax) SWAP(float,tzMin, tzMax);

        if ((tMin > tzMax) || (tzMin > tMax))
            return false;

        t = tMin;
        return true;
    }
};


struct BoxCentroid {
    float3 centroid;
    int index;
};

struct TriangleCentroid {
    float3 centroid;
    int index;
};

struct Rectangle {
    Vec3 v0, v1, v2, v3;
};

struct F3Ray {
    float3 origin;
    float3 direction;
};

struct Ray {
    Vec3 origin, direction;
};


struct BVHNode {
    float3 min, max; // Min and max coordinates of the bounding box
    int leftChild, rightChild; // Indices of child nodes (-1 for leaves)
    int triangleIndex;
    int splitAxis; //add
    int boxIndex; //add
   
};

float3 toFloat3(const Vec3& v) { return {v.x, v.y, v.z}; }

std::vector<F3Triangle> boxToTriangles(const Box& box) {
    std::vector<F3Triangle> triangles;
    triangles.reserve(12);

    Vec3 vertices[8] = {
        {box.min.x, box.min.y, box.min.z},
        {box.max.x, box.min.y, box.min.z},
        {box.max.x, box.max.y, box.min.z},
        {box.min.x, box.max.y, box.min.z},
        {box.min.x, box.min.y, box.max.z},
        {box.max.x, box.min.y, box.max.z},
        {box.max.x, box.max.y, box.max.z},
        {box.min.x, box.max.y, box.max.z}
    };

    // Definition of the 12 triangles (2 per face)
    int indices[12][3] = {
        {0, 1, 2}, {0, 2, 3}, // Front face
        {1, 5, 6}, {1, 6, 2}, // Right face
        {5, 4, 7}, {5, 7, 6}, // Back face
        {4, 0, 3}, {4, 3, 7}, // Left face
        {3, 2, 6}, {3, 6, 7}, // Top face
        {4, 5, 1}, {4, 1, 0}  // Bottom face
    };

    for (int i = 0; i < 12; ++i) {
        F3Triangle tri;
        tri.v0 = toFloat3(vertices[indices[i][0]]);
        tri.v1 = toFloat3(vertices[indices[i][1]]);
        tri.v2 = toFloat3(vertices[indices[i][2]]);
        triangles.push_back(tri);
    }
    return triangles;
}


// Intersection function between a ray and a plan
__host__ __device__ std::optional<Vec3> rayPlaneIntersect(const Ray& ray, const Vec3& planePoint, const Vec3& planeNormal) {
    constexpr float epsilon = 1e-6f;

    float denom = planeNormal.dot(ray.direction);

    // Vérification si le rayon est parallèle au plan
    if (fabs(denom) < epsilon) {
        return {}; // Pas d'intersection
    }

    Vec3 p0l0 = planePoint - ray.origin;
    float t = p0l0.dot(planeNormal) / denom;

    // Check if the ray is parallel to the plane
    if (t < 0) {
        return {}; // The plan is behind the ray
    }

    // Calculate the intersection point
    return { ray.origin + ray.direction * t };
}

// Intersection function between a ray and a rectangle
__host__ __device__ std::optional<Vec3> rayRectangleIntersect(const Ray& ray, const Rectangle& rect) {
    // We must define the plane of the rectangle
    Vec3 edge1 = rect.v1 - rect.v0;
    Vec3 edge2 = rect.v3 - rect.v0;
    Vec3 normal = edge1.cross(edge2); // Normal of the rectangle

    constexpr float epsilon = 1e-6f;
    float denom = normal.dot(ray.direction);

    // Check if the ray is parallel to the plane
    if (fabs(denom) < epsilon) {
        return {};// No intersection
    }

    // Calculation of the distance t at which the ray intersects the plane
    Vec3 p0l0 = rect.v0 - ray.origin;
    float t = p0l0.dot(normal) / denom;

    // Check if the intersection is in front of the ray
    if (t < 0) {
        return {};// The rectangle is behind the ray
    }

    // The rectangle is behind the ray
    Vec3 intersectionPoint = ray.origin + ray.direction * t;

    // Check if the intersection point is inside the rectangle
    Vec3 c;

    // Check for the first side
    Vec3 edge00 = rect.v1 - rect.v0;
    Vec3 vp0 = intersectionPoint - rect.v0;
    c = edge00.cross(vp0);
    if (normal.dot(c) < 0) return {}; // The point is outside

    // Checking for the second side
    Vec3 edge01 = rect.v2 - rect.v1;
    Vec3 vp1 = intersectionPoint - rect.v1;
    c = edge01.cross(vp1);
    if (normal.dot(c) < 0) return {}; // The point is outside

    // Checking for the third side
    Vec3 edge02 = rect.v3 - rect.v2;
    Vec3 vp2 = intersectionPoint - rect.v2;
    c = edge02.cross(vp2);
    if (normal.dot(c) < 0) return {}; // The point is outside

    // Checking for the fourth side
    Vec3 edge03 = rect.v0 - rect.v3;
    Vec3 vp3 = intersectionPoint - rect.v3;
    c = edge03.cross(vp3);
    if (normal.dot(c) < 0) return {}; // The point is outside

    return intersectionPoint;
}

// Function to calculate the bounding box of a triangle
__host__ __device__
void calculateBoundingBox(const F3Triangle& triangle, float3& min, float3& max)
{
    min = make_float3(fminf(fminf(triangle.v0.x, triangle.v1.x), triangle.v2.x),
                      fminf(fminf(triangle.v0.y, triangle.v1.y), triangle.v2.y),
                      fminf(fminf(triangle.v0.z, triangle.v1.z), triangle.v2.z));
    max = make_float3(fmaxf(fmaxf(triangle.v0.x, triangle.v1.x), triangle.v2.x),
                      fmaxf(fmaxf(triangle.v0.y, triangle.v1.y), triangle.v2.y),
                      fmaxf(fmaxf(triangle.v0.z, triangle.v1.z), triangle.v2.z));
}


// Function to build a simple BVH (medium construction method)
void buildBVHWithTriangleVersion1(thrust::device_vector<F3Triangle>& triangles, thrust::device_vector<BVHNode>& nodes)
{
    int numTriangles = triangles.size();
    nodes.resize(2 * numTriangles - 1);

    // Initialize the sheets
    for (int i = 0; i < numTriangles; ++i) {
        BVHNode* raw_ptr = thrust::raw_pointer_cast(nodes.data());
        BVHNode& node = raw_ptr[numTriangles - 1 + i];

        calculateBoundingBox(triangles[i], node.min, node.max);
        node.triangleIndex = i;
        node.leftChild = node.rightChild = -1;
    }

    // Build the internal nodes
    for (int i = numTriangles - 2; i >= 0; --i) {
        BVHNode* raw_ptr = thrust::raw_pointer_cast(nodes.data());
        BVHNode& node = raw_ptr[i];
        int leftChild = 2 * i + 1;
        int rightChild = 2 * i + 2;

        node.leftChild = leftChild;
        node.rightChild = rightChild;
        node.triangleIndex = -1;

        BVHNode leftNode = nodes[leftChild];
        BVHNode rightNode = nodes[rightChild];
        node.min = make_float3(fminf(leftNode.min.x, rightNode.min.x),
                               fminf(leftNode.min.y, rightNode.min.y),
                               fminf(leftNode.min.z, rightNode.min.z));


        node.max = make_float3(fmaxf(leftNode.max.x, rightNode.max.x),
                               fmaxf(leftNode.max.y, rightNode.max.y),
                               fmaxf(leftNode.max.z, rightNode.max.z));
    }
}

// Function to make a dot
__host__ __device__
float dot(const float3& a, const float3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Function to make a cross
__host__ __device__
float3 cross(const float3& a, const float3& b) {
    return float3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

// Function to return a length
__host__ __device__
float length(const float3& v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Function to normalize
__host__ __device__
float3 normalize(const float3& v) {
    float len = length(v);
    if (len > 0) {
        return float3(v.x / len, v.y / len, v.z / len);
    }
    return v;
}

// Function to write a float3
__host__ __device__
void print_float3(const float3& v) {
    printf("%f %f %f\n",v.x,v.y,v.z);
}

__device__ bool rayTriangleIntersect(const F3Ray& ray, const F3Triangle& triangle, float& t, float3& intersectionPoint) {
    float3 edge1 = triangle.v1 - triangle.v0;
    float3 edge2 = triangle.v2 - triangle.v0;
    float3 h = cross(ray.direction, edge2);
    float a = dot(edge1, h);

    if (a > -1e-6 && a < 1e-6) return false;

    float f = 1.0f / a;
    float3 s = ray.origin - triangle.v0;
    float u = f * dot(s, h);

    if (u < 0.0f || u > 1.0f) return false;

    float3 q = cross(s, edge1);
    float v = f * dot(ray.direction, q);

    if (v < 0.0f || u + v > 1.0f) return false;

    t = f * dot(edge2, q);

     // Calculate the intersection point
    if (t > 1e-6) {
        intersectionPoint = ray.origin + t * ray.direction;
        //printf("%f %f %f\n",intersectionPoint.x,intersectionPoint.y,intersectionPoint.z); OK
    }
    else
    {
        intersectionPoint =make_float3(INFINITY, INFINITY, INFINITY);
    }

    return (t > 1e-6);
}

__global__ void rayTracingKernel(BVHNode* nodes, F3Triangle* triangles, F3Ray* rays, int* hitResults, float* distance,float3* intersectionPoint, int numRays) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numRays) return;

    F3Ray ray = rays[idx];
    int stack[64];
    int stackPtr = 0;
    stack[stackPtr++] = 0;

    float closestT = INFINITY;
    int closestTriangle = -1;
    float3 closestIntersectionPoint=make_float3(INFINITY, INFINITY, INFINITY);
    bool isView=false;

    while (stackPtr > 0) {
        int nodeIdx = stack[--stackPtr];
        BVHNode& node = nodes[nodeIdx];

        // Ray-box intersection test
        float tmin = (node.min.x - ray.origin.x) / ray.direction.x;
        float tmax = (node.max.x - ray.origin.x) / ray.direction.x;
        if (tmin > tmax) SWAP(float,tmin, tmax);

        float tymin = (node.min.y - ray.origin.y) / ray.direction.y;
        float tymax = (node.max.y - ray.origin.y) / ray.direction.y;
        if (tymin > tymax) SWAP(float,tymin, tymax);

        if ((tmin > tymax) || (tymin > tmax)) continue;

        if (tymin > tmin) tmin = tymin;
        if (tymax < tmax) tmax = tymax;

        float tzmin = (node.min.z - ray.origin.z) / ray.direction.z;
        float tzmax = (node.max.z - ray.origin.z) / ray.direction.z;
        if (tzmin > tzmax) SWAP(float,tzmin, tzmax);

        if ((tmin > tzmax) || (tzmin > tmax)) continue;

        if (tzmin > tmin) tmin = tzmin;
        if (tzmax < tmax) tmax = tzmax;

        if (tmax < 0) continue;

        if (node.triangleIndex != -1) {
            // Sheet: test the intersection with the triangle
            float t;
            float3 intersectionPointT;
            if (rayTriangleIntersect(ray, triangles[node.triangleIndex], t,intersectionPointT)) {

                if (isView) printf("[%i] <%f %f %f>\n",idx,intersectionPointT.x,intersectionPointT.y,intersectionPointT.z);

                if (t < closestT) {
                    closestT = t;
                    closestTriangle = node.triangleIndex;
                    closestIntersectionPoint=intersectionPointT;
                }
            }
        } else {
            // Internal node: add children to the stack
            stack[stackPtr++] = node.leftChild;
            stack[stackPtr++] = node.rightChild;
        }
    }

    hitResults[idx]        = closestTriangle;
    distance[idx]          = closestT;
    intersectionPoint[idx] = closestIntersectionPoint;
    //if (closestTriangle!=-1) { printf("t=%f\n",closestT); } OK
}
}
*/



__global__ void kernel(float* x,float* y,int n){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    //int idx = hipblockIdx.x * hipblockDim.x + hipthreadIdx.x;
    if (idx < n) y[idx] += 1;
}




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

struct BVHEnum
{
    enum class Quality { Low, Medium, High };
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
    using vector_realdim_type = Eigen::Matrix<double,nRealDim,1>;
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
        void setCoordinates( T && coord ) { M_coordinates = std::forward<T>( coord ); }
        bool hasCoordinates() const noexcept { return M_coordinates.has_value(); }
        vector_realdim_type const& coordinates() const noexcept { return *M_coordinates; }

    private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize( Archive& ar, const unsigned int version )
            {
                ar & M_processId;
                ar & M_primitiveId;
                ar & M_distance;

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
            }
    private:
        rank_type M_processId = invalid_v<rank_type>;
        index_type M_primitiveId = invalid_v<index_type>;
        double M_distance = std::numeric_limits<double>::max();
        std::optional<vector_realdim_type> M_coordinates;
    };
    using rayintersection_result_type = BVHRayIntersectionResult;

    enum class IntersectContext{ anyHint=0, closest, all };

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

            bool closestOnly = ctx == IntersectContext::closest;
            using napp_ray_type = std::decay_t<decltype(ray)>;
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

                auto intersectGlobal = this->intersect(_ray=raysGathered,_robust=useRobustTraversal,_context=ctx,_parallel=true);

                std::vector<std::vector<rayintersection_result_type>> res;
                res.resize( ray.numberOfLocalRay() );
                std::size_t startRayIndexInThisProcess = 0;
                for ( int p=0;p<this->worldComm().rank();++p )
                    startRayIndexInThisProcess += resLocalSize[p];
                std::copy_n(intersectGlobal.cbegin()+startRayIndexInThisProcess, localRays.size(), res.begin());
                return res;
            }
            else if constexpr ( is_iterable_v<std::decay_t<decltype(ray)>> ) // case rays container are identical all on process (TODO: internal case)
            {
                std::vector<std::vector<rayintersection_result_type>> resSeq;
                resSeq.reserve( ray.size() );
				
				//A modifier ici pour le ray tracing HIP AMD GPU THRUST
                for ( auto const& currentRay : ray )
                {
                    auto currentResSeq = this->intersectSequential( currentRay,useRobustTraversal );
                    if ( closestOnly && currentResSeq.size() > 1 )
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
                auto resSeq = this->intersectSequential( ray,useRobustTraversal );
                if ( closestOnly && resSeq.size() > 1 )
                    resSeq.resize(1);
                if ( !parallel )
                    return resSeq;

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
                return resSeq;
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
protected:

    virtual std::vector<rayintersection_result_type> intersectSequential( ray_type const& rayon, bool useRobustTraversal = true ) = 0;

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

    using value_type = double;
    using node_type = bvh::v2::Node<value_type, nRealDim>;
    using backend_bvh_type = bvh::v2::Bvh<node_type>;
    using backend_vector_realdim_type = bvh::v2::Vec<value_type, nRealDim>;
    using backend_precompute_triangle_type = bvh::v2::PrecomputedTri<value_type>;
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
            M_bvh = std::make_unique<backend_bvh_type>( bvh::v2::DefaultBuilder<node_type>::build(/*thread_pool,*/ bboxes, centers, config) );


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
        }

private:

    std::vector<rayintersection_result_type> intersectSequential( ray_type const& ray, bool useRobustTraversal = true ) override
        {
            auto rayBackend = bvh::v2::Ray<value_type,nRealDim>{
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.origin()[i]; }),
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.dir()[i]; }),
                ray.distanceMin(), ray.distanceMax()
            };
            if ( useRobustTraversal )
                return this->intersectImpl<true>( rayBackend );
            else
                return this->intersectImpl<false>( rayBackend );
        };

    template <bool UseRobustTraversal>
    std::vector<rayintersection_result_type> intersectImpl( bvh::v2::Ray<value_type,nRealDim> & rayBackend )
        {
            static constexpr size_t stack_size = 64;
            static constexpr bool should_permute = true;
            static constexpr bool isAnyHit = false;
            // Traverse the BVH and get the u, v coordinates of the closest intersection.
            bvh::v2::SmallStack<typename backend_bvh_type::Index, stack_size> stack;
            std::vector<rayintersection_result_type> res;
            M_bvh->template intersect<isAnyHit, UseRobustTraversal>( rayBackend, M_bvh->get_root().index, stack,
                                                                     [this,&res,&rayBackend] (std::size_t begin, std::size_t end) {
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
                                                                                 if (auto hit = M_precomputeTriangle[j].intersect(rayBackend))
                                                                                 {
                                                                                     //std::tie(u, v) = *hit;
                                                                                     res.push_back( rayintersection_result_type(this->worldComm().rank(), M_bvh->prim_ids[i], rayBackend.tmax) );
                                                                                     res.back().setCoordinates( this->barycentricToCartesianCoordinates( M_precomputeTriangle[j].convert_to_tri(), hit->first, hit->second ) );
                                                                                     if constexpr ( isAnyHit )
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
};



//! @brief implementation of BVH tool with an external third party in GPU
template <typename MeshEntityType>
class BVH_GPUParty : public BVH<MeshEntityType>
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

    bool isSaveBVH;
    bool isLoadBVH;
    std::string M_FileName;


public:
    using ray_type = typename super_type::ray_type;
    using rayintersection_result_type = typename super_type::rayintersection_result_type;

    BVH_GPUParty( BVHEnum::Quality quality, worldcomm_ptr_t worldComm ) : super_type( quality,worldComm ) 
    { 
        M_FileName="bvh";
    }

    BVH_GPUParty( BVH_GPUParty && ) = default;

    static std::optional<backend_bvh_type> load_bvh(const std::string& file_name) {
            std::ifstream in(file_name, std::ofstream::binary);
            if (!in) return std::nullopt;
            bvh::v2::StdInputStream stream(in);
        return std::make_optional(backend_bvh_type::deserialize(stream));
    }   

    static bool save_bvh(const backend_bvh_type& bvh, const std::string& file_name) {
            std::ofstream out(file_name, std::ofstream::binary);
            if (!out) return false;
            bvh::v2::StdOutputStream stream(out);
            bvh.serialize(stream);
        return true;
    }

    void setFileName(std::string s) { M_FileName=s; } 




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


            std::cout << "Size primitiveInfo="<<this->M_primitiveInfo.size()<< std::endl;
            std::cout << "Size bboxes="<<bboxes.size()<< std::endl;
            std::cout << "bboxes sizeof="<<sizeof(bboxes)<< std::endl;

            for ( auto const& primInfo : this->M_primitiveInfo )
            {
                bboxes.push_back( BBox{
                        backend_vector_realdim_type::generate([&primInfo] (std::size_t i) { return primInfo.boundMin()[i]; }),
                        backend_vector_realdim_type::generate([&primInfo] (std::size_t i) { return primInfo.boundMax()[i]; })
                    });

                auto const& centroid = primInfo.centroid();
                centers.push_back( backend_vector_realdim_type::generate([&centroid] (std::size_t i) { return centroid[i]; }) );

            }

            std::cout << "PRIMITIVES\n";
            std::cout << "Size primitiveInfo="<<this->M_primitiveInfo.size()<< std::endl;
            for (int k=0; k<centers.size();++k)
            {
                std::cout<<"["<<k<<"]=<"<<centers[k].values[0]<<","<<centers[k].values[1]<<","<<centers[k].values[2]<<">\n";
            };

            std::cout << "BBOXES\n";
            std::cout << "Size of bboxes="<<bboxes.size()<< std::endl;
            for (int k=0; k<bboxes.size();++k)
            {
                std::cout<<"["<<k<<"] min=<"<<bboxes[k].min.values[0]<<","<<bboxes[k].min.values[1]<<","<<bboxes[k].min.values[2]<<"> ";
                std::cout<<"max=<"<<bboxes[k].max.values[0]<<","<<bboxes[k].max.values[1]<<","<<bboxes[k].max.values[2]<<"> ";
                std::cout<<"Centroid=<"<<centers[k].values[0]<<","<<centers[k].values[1]<<","<<centers[k].values[2]<<">\n";
            };

            for (int k=0; k<bboxes.size();++k)
            {
                std::cout<<"["<<k<<"] primInfo="<<this->M_primitiveInfo[k].boundMin()[0]<<","
                <<this->M_primitiveInfo[k].boundMin()[1]<<","
                <<this->M_primitiveInfo[k].boundMin()[2]<<"\n";
            };


        

            typename bvh::v2::DefaultBuilder<node_type>::Config config;
            switch ( this->M_quality )
            {
            default:
            case BVHEnum::Quality::High: config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::High; break;
            case BVHEnum::Quality::Medium: config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::Medium; break;
            case BVHEnum::Quality::Low: config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::Low; break;
            }

            //M_bvh = std::make_unique<backend_bvh_type>( bvh::v2::DefaultBuilder<node_type>::build(/*thread_pool,*/ bboxes, centers, config) );
            auto bvh_v2 = bvh::v2::DefaultBuilder<node_type>::build(/*thread_pool,*/ bboxes, centers, config);

            
            isSaveBVH=true;
            if (isSaveBVH)
            {
                std::cout << "Save BVH" << std::endl;
                save_bvh(bvh_v2, M_FileName+".bin");
            }

            isLoadBVH=false;
            if (isLoadBVH)
            {
                std::cout << "Load BVH" << std::endl;
                auto other_bvh = load_bvh(M_FileName+".bin");
                if (bvh_v2 == other_bvh) { std::cout << "The deserialized BVH is the same as the original one" << std::endl; }
                else { std::cerr << "The deserialized BVH does not match the original one" << std::endl; }
            }


            M_bvh = std::make_unique<backend_bvh_type>(std::move(bvh_v2));

            std::cout << "Size of M_bvh="<<sizeof(M_bvh)<< std::endl;

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


        }

private:

    std::vector<rayintersection_result_type> intersectSequential( ray_type const& ray, bool useRobustTraversal = true ) override
        {
            auto rayBackend = bvh::v2::Ray<value_type,nRealDim>{
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.origin()[i]; }),
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.dir()[i]; }),
                ray.distanceMin(), ray.distanceMax()
            };
            if ( useRobustTraversal )
                return this->intersectImpl<true>( rayBackend );
            else
                return this->intersectImpl<false>( rayBackend );
        };

    template <bool UseRobustTraversal>
    std::vector<rayintersection_result_type> intersectImpl( bvh::v2::Ray<value_type,nRealDim> & rayBackend )
        {
            static constexpr size_t stack_size = 64;
            static constexpr bool should_permute = true;
            static constexpr bool isAnyHit = false;
            // Traverse the BVH and get the u, v coordinates of the closest intersection.
            bvh::v2::SmallStack<typename backend_bvh_type::Index, stack_size> stack;
            std::vector<rayintersection_result_type> res;
            M_bvh->template intersect<isAnyHit, UseRobustTraversal>( rayBackend, M_bvh->get_root().index, stack,
                                                                     [this,&res,&rayBackend] (std::size_t begin, std::size_t end) {
                                                                         std::size_t previousResultSize = res.size();

                                                                         std::cout << "end="<<end<< std::endl;

                                                                         for (std::size_t i = begin; i < end; ++i)
                                                                         {
                                                                             std::size_t j = should_permute ? i : M_bvh->prim_ids[i];
                                                                             if constexpr ( nRealDim == 2 )
                                                                             {
                                                                                 CHECK( false ) << "TODO";
                                                                             }
                                                                             else if constexpr ( nRealDim == 3 )
                                                                             {
                                                                                 if (auto hit = M_precomputeTriangle[j].intersect(rayBackend))
                                                                                 {
                                                                                     //std::tie(u, v) = *hit;
                                                                                     res.push_back( rayintersection_result_type(this->worldComm().rank(), M_bvh->prim_ids[i], rayBackend.tmax) );
                                                                                     res.back().setCoordinates( this->barycentricToCartesianCoordinates( M_precomputeTriangle[j].convert_to_tri(), hit->first, hit->second ) );
                                                                                     if constexpr ( isAnyHit )
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
};



/***************************************************************************************************************************************************/






//! @brief implementation of BVH tool with an external third party in GPU
template <typename MeshEntityType>
class BVH_HIP_Party : public BVH<MeshEntityType>
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

    BVH_HIP_Party( BVHEnum::Quality quality, worldcomm_ptr_t worldComm ) : super_type( quality,worldComm ) 
    { }

    BVH_HIP_Party( BVH_HIP_Party && ) = default;

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


            std::cout << "Size primitiveInfo="<<this->M_primitiveInfo.size()<< std::endl;
            std::cout << "Size bboxes="<<bboxes.size()<< std::endl;
            std::cout << "bboxes sizeof="<<sizeof(bboxes)<< std::endl;

            for ( auto const& primInfo : this->M_primitiveInfo )
            {
                bboxes.push_back( BBox{
                        backend_vector_realdim_type::generate([&primInfo] (std::size_t i) { return primInfo.boundMin()[i]; }),
                        backend_vector_realdim_type::generate([&primInfo] (std::size_t i) { return primInfo.boundMax()[i]; })
                    });

                auto const& centroid = primInfo.centroid();
                centers.push_back( backend_vector_realdim_type::generate([&centroid] (std::size_t i) { return centroid[i]; }) );

            }
			
			/*
				// Load the mesh
				std::vector<F3Triangle> hostTriangles;	
				...
				// Transfer triangles to GPU
				thrust::device_vector<F3Triangle> deviceTriangles = hostTriangles;	

				// Building the BVH
				thrust::device_vector<BVHNode> deviceNodes;
				buildBVHWithTriangleVersion1(deviceTriangles, deviceNodes);	
                std::cout<<"[INFO]: BVH built with " << deviceNodes.size() << " nodes" << std::endl;				
			*/

            std::cout << "PRIMITIVES\n";
            std::cout << "Size primitiveInfo="<<this->M_primitiveInfo.size()<< std::endl;
            for (int k=0; k<centers.size();++k)
            {
                std::cout<<"["<<k<<"]=<"<<centers[k].values[0]<<","<<centers[k].values[1]<<","<<centers[k].values[2]<<">\n";
            };

            std::cout << "BBOXES\n";
            std::cout << "Size of bboxes="<<bboxes.size()<< std::endl;
            for (int k=0; k<bboxes.size();++k)
            {
                std::cout<<"["<<k<<"] min=<"<<bboxes[k].min.values[0]<<","<<bboxes[k].min.values[1]<<","<<bboxes[k].min.values[2]<<"> ";
                std::cout<<"max=<"<<bboxes[k].max.values[0]<<","<<bboxes[k].max.values[1]<<","<<bboxes[k].max.values[2]<<"> ";
                std::cout<<"Centroid=<"<<centers[k].values[0]<<","<<centers[k].values[1]<<","<<centers[k].values[2]<<">\n";
            };

            for (int k=0; k<bboxes.size();++k)
            {
                std::cout<<"["<<k<<"] primInfo="<<this->M_primitiveInfo[k].boundMin()[0]<<","
                <<this->M_primitiveInfo[k].boundMin()[1]<<","
                <<this->M_primitiveInfo[k].boundMin()[2]<<"\n";
            };


        
			/*
            typename bvh::v2::DefaultBuilder<node_type>::Config config;
            switch ( this->M_quality )
            {
            default:
				case BVHEnum::Quality::High: config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::High; break;
				case BVHEnum::Quality::Medium: config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::Medium; break;
				case BVHEnum::Quality::Low: config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::Low; break;
            }

            //M_bvh = std::make_unique<backend_bvh_type>( bvh::v2::DefaultBuilder<node_type>::build(bboxes, centers, config) );
            auto bvh_v2 = bvh::v2::DefaultBuilder<node_type>::build(bboxes, centers, config);

            M_bvh = std::make_unique<backend_bvh_type>(std::move(bvh_v2));


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
			*/


        }

private:

    std::vector<rayintersection_result_type> intersectSequential( ray_type const& ray, bool useRobustTraversal = true ) override
        {
            auto rayBackend = bvh::v2::Ray<value_type,nRealDim>{
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.origin()[i]; }),
                backend_vector_realdim_type::generate([&ray] (std::size_t i) { return ray.dir()[i]; }),
                ray.distanceMin(), ray.distanceMax()
            };
            if ( useRobustTraversal )
                return this->intersectImpl<true>( rayBackend );
            else
                return this->intersectImpl<false>( rayBackend );
        };

    template <bool UseRobustTraversal>
    std::vector<rayintersection_result_type> intersectImpl( bvh::v2::Ray<value_type,nRealDim> & rayBackend )
        {
            static constexpr size_t stack_size = 64;
            static constexpr bool should_permute = true;
            static constexpr bool isAnyHit = false;
            // Traverse the BVH and get the u, v coordinates of the closest intersection.
            bvh::v2::SmallStack<typename backend_bvh_type::Index, stack_size> stack;
            std::vector<rayintersection_result_type> res;
            M_bvh->template intersect<isAnyHit, UseRobustTraversal>( rayBackend, M_bvh->get_root().index, stack,
                                                                     [this,&res,&rayBackend] (std::size_t begin, std::size_t end) {
                                                                         std::size_t previousResultSize = res.size();

                                                                         std::cout << "end="<<end<< std::endl;

                                                                         for (std::size_t i = begin; i < end; ++i)
                                                                         {
                                                                             std::size_t j = should_permute ? i : M_bvh->prim_ids[i];
                                                                             if constexpr ( nRealDim == 2 )
                                                                             {
                                                                                 CHECK( false ) << "TODO";
                                                                             }
                                                                             else if constexpr ( nRealDim == 3 )
                                                                             {
                                                                                 if (auto hit = M_precomputeTriangle[j].intersect(rayBackend))
                                                                                 {
                                                                                     //std::tie(u, v) = *hit;
                                                                                     res.push_back( rayintersection_result_type(this->worldComm().rank(), M_bvh->prim_ids[i], rayBackend.tmax) );
                                                                                     res.back().setCoordinates( this->barycentricToCartesianCoordinates( M_precomputeTriangle[j].convert_to_tri(), hit->first, hit->second ) );
                                                                                     if constexpr ( isAnyHit )
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
};







/***************************************************************************************************************************************************/

// Remarks: This model does not work properly. Depending on the complexity of the model, there are errors.
//! @brief in house implementation of BVH tool
template <typename MeshEntityType>
class BVH_InHouse : public BVH<MeshEntityType>
{
    using super_type = BVH<MeshEntityType>;
    using self_type = BVH_InHouse<MeshEntityType>;
    using mesh_entity_type = typename super_type::mesh_entity_type;
    using vector_realdim_type = typename super_type::vector_realdim_type;
    static constexpr uint16_type nRealDim = super_type::nRealDim;
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
    std::vector<rayintersection_result_type> intersectSequential( ray_type const& rayon, bool useRobustTraversal = true ) override
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


template<typename... Ts>
auto boundingVolumeHierarchy( Ts && ... v )
{

/*
    const int N = 100;
    thrust::device_vector<float> d_x(N, 1.0f); 
    thrust::device_vector<float> d_y(N, 2.0f);
*/


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
    else if ( kind == "gpu-party" )
    {
        if constexpr ( mesh_entity_type::nRealDim != 3 )
            throw std::invalid_argument("gpu-party only implement with triangle in 3D");
        auto bvhGPUParty = std::make_unique<BVH_GPUParty<mesh_entity_type>>( quality, worldcomm );
        bvhGPUParty->updateForUse(range);
        bvh = std::move( bvhGPUParty );
    }
	
	else if ( kind == "hip-party" )
    {
        if constexpr ( mesh_entity_type::nRealDim != 3 )
            throw std::invalid_argument("gpu-party only implement with triangle in 3D");
        auto bvhHIPParty = std::make_unique<BVH_HIP_Party<mesh_entity_type>>( quality, worldcomm );
        bvhHIPParty->updateForUse(range);
        bvh = std::move( bvhHIPParty );
    }
    else
        throw std::invalid_argument(fmt::format("invalid bvh arg kind {} (should be third-party or in-house)",kind ));

    return bvh;
}

} // Feel
#endif /* FEELPP_MESH_BVH_HPP */