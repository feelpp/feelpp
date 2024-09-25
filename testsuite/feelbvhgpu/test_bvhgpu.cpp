#define BOOST_TEST_MODULE bvhgpu_tests


/*
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/bvh.hpp>

#include <feel/feeldiscr/pdh.hpp>
*/

#include <feel/feelmesh/ranges.hpp>

#include <feel/feelcore/enumerate.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/bvh.hpp>
#include <fmt/chrono.h>

#include <feel/feelcore/json.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>

#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/partitionmesh.hpp>
#include <feel/feelfilters/partitionio.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelvf/vf.hpp>



//#define COMPILE_WITH_HIP
//#include <feel/feeltask/taskpu.hpp>


//#include <hip/hip_runtime.h>
//#include <hip/hip_runtime_api.h>

#include "hip/hip_runtime.h"
#include "hip/hip_runtime_api.h"
//#include "hipblas.h"
//#include "hipsolver.h"
//#include "hipblas-export.h"

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


#include <feel/feeltask/taskpu.hpp>

using namespace Feel;

#define BLOCK_SIZE 128


//#include "Gu.hpp"
//#include "GuHip.hpp"


#include "bvhHybrid.hpp"


namespace bvhhipTestIfCompile
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


// Function to load a triangle mesh from a simple OBJ file
bool loadOBJTriangle(const std::string& filename, std::vector<F3Triangle>& triangles) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::vector<float3> vertices;
    std::string line;
    bool isView=false;
    while (std::getline(file, line)) {
        if (line[0] == 'v') {
            float x, y, z;
            sscanf(line.c_str(), "v %f %f %f", &x, &y, &z);
            vertices.push_back(make_float3(x, y, z));
            if (isView) std::cout <<"v=<" << x <<","<<y<<","<<z<< ">\n";
        } else if (line[0] == 'f') {
            unsigned int i1, i2, i3;
            sscanf(line.c_str(), "f %u %u %u", &i1, &i2, &i3);
            if (isView) std::cout <<"f=<" << i1 <<","<<i2<<","<<i3<< ">\n";
            triangles.push_back({vertices[i1-1], vertices[i2-1], vertices[i3-1]});
        }
    }
    return true;
}


}


//=========================================================================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

__global__ void kernel(float* x,float* y,int n){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    //int idx = hipblockIdx.x * hipblockDim.x + hipthreadIdx.x;
    if (idx < n) y[idx] += 1;
}


void __device__ vector_add_device(const double vecA, const double vecB, double &vecC)
{
    vecC = vecA + vecB;
}

void __global__ vector_add(const double *vecA, const double *vecB, double *vecC, const int nb)
{
    const int i = hipBlockDim_x*hipBlockIdx_x+hipThreadIdx_x;

    if (i < nb)
        vector_add_device(vecA[i], vecB[i], vecC[i]); 
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//=========================================================================================================================


//=========================================================================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

struct saxpy_functor
{
    const float a;
    saxpy_functor(float _a) : a(_a) {}
    __host__ __device__
    float operator()(const float& x, const float& y) const { return a * x + y; }
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//=========================================================================================================================


//=========================================================================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void check_solution(double coeff,double* a_in,double* b_in,double* c_in,int nb)
{
  printf("[INFO]: ");
	int errors = 0;
  	for (int i = 0; i < nb; i++) {
	    if (c_in[i] != coeff*(a_in[i] + b_in[i])) { errors++; }
	}
  	if (errors!=0) { printf("FAILED: %d errors\n",errors); } else { printf ("WELL DONE PASSED! :-)\n"); }
}

void write_vector(std::string ch,double* v,int nb)
{
  std::cout<<"[INFO]: "<<ch<<"> ";
	for (int i = 0; i < nb; i++) { std::cout<<int(v[i]); }
  std::cout<<std::endl;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//=========================================================================================================================




//=========================================================================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Section to use if you want to make meshes, otherwise it doesn't work.
inline
AboutData
makeAbout()
{
    AboutData about( "test_bvhgpu" ,
                     "test_bvhgpu" ,
                     "0.2",
                     "nD(n=2,3) test bvh",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2024 Feel++ Consortium" );

    about.addAuthor( "Noname", "developer", "Noname@cemosis.fr", "" );
    return about;
}

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description opts( "Test Environment options" );
    opts.add_options()
        ( "mesh2D.filename", po::value<std::string>(), "mesh2D.filename" )
        ( "mesh3D.filename", po::value<std::string>(), "mesh3D.filename" )
        ;
    return opts;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//=========================================================================================================================


//=========================================================================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template <typename BvhType,typename RayIntersectionResultType>
void printRayIntersectionResultsBeta( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs )
{
    if ( bvh->worldComm().isMasterRank() )
        BOOST_TEST_MESSAGE( "Number of intersection: " << rirs.size() );
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            BOOST_TEST_MESSAGE( " --  ProcessId: " << rir.processId() );
            BOOST_TEST_MESSAGE( " --  PrimitiveId: " << rir.primitiveId() );
            BOOST_TEST_MESSAGE( " --  Distance: " << rir.distance() );

            auto const& prim = bvh->primitiveInfo( rir.primitiveId() );
            BOOST_TEST_MESSAGE( " --  Mesh entity id: " << prim.meshEntity().id() );
            BOOST_TEST_MESSAGE( " --  Mesh entity barycenter: " << prim.meshEntity().barycenter() );
            BOOST_TEST_MESSAGE( " ---------------------------------" );
        }
        bvh->worldComm().barrier();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}

template <typename BvhType,typename RayIntersectionResultType>
void printRayIntersectionResults( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs )
{
    if ( bvh->worldComm().isMasterRank() )
        std::cout << "Number of intersection: " << rirs.size()<< "\n";
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            std::cout << " --  ProcessId: " << rir.processId()<< "\n";
            std::cout << " --  PrimitiveId: " << rir.primitiveId()<< "\n";
            std::cout << " --  Distance: " << rir.distance()<< "\n";

            auto const& prim = bvh->primitiveInfo( rir.primitiveId() );
            std::cout << " --  Mesh entity id: " << prim.meshEntity().id()<< "\n";
            std::cout << " --  Mesh entity barycenter: " << prim.meshEntity().barycenter()<< "\n";
            //std::cout << " --  Mesh entity barycenter: " << prim.meshEntity()
            std::cout <<" ---------------------------------"<< "\n" ;
        }
        bvh->worldComm().barrier();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}
 

template <typename RangeType>
void test3D( RangeType const& range )
{
    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    Eigen::Vector3d origin1={1000.0,0.0,0.0};
    Eigen::Vector3d direction_perp_1={1.,0.,0.};

    Eigen::Vector3d origin2={-10.0,-0.25,-0.25};
    Eigen::Vector3d direction_perp_2={1.,0.,0.};

    Eigen::Vector3d origin3={-5.0,-0.20,-0.20};
    Eigen::Vector3d direction_perp_3={1.,0.,0.};

    std::vector<bvh_ray_type> rays;

    //rays.push_back( bvh_ray_type(origin1,direction_perp_1) );
    rays.push_back( bvh_ray_type(origin2,direction_perp_2) );
    rays.push_back( bvh_ray_type(origin3,direction_perp_3) );

    BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;
    for (int k=0;k<rays.size();++k)
    {
        raysDistributed.push_back( rays[k] );
    }


    std::cout<<"\n";
    std::cout<<"\n";

    std::cout<<"=====================================================================================+============\n";
    std::cout<<"++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++\n";
    std::cout<<"[INFO]: Bounding Colume Hierarchy\n";

    auto bvhInHouse = boundingVolumeHierarchy(_range=range,_kind="in-house");
    //auto bvhThirdParty = boundingVolumeHierarchy(_range=range,_kind="third-party");
    auto bvhThirdPartyLow = boundingVolumeHierarchy(_range=range,_kind="third-party",_quality=BVHEnum::Quality::High);

    //auto bvhGpuParty = boundingVolumeHierarchy(_range=range,_kind="gpu-party",_quality=BVHEnum::Quality::High);

    std::cout<<"++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++\n";
    std::cout<<"========================================================================================+=========\n";

    std::cout<<"\n";
    std::cout<<"\n";

    /*
    for ( auto bvhCurrent : {bvhInHouse.get(),bvhThirdParty.get(),bvhThirdPartyLow.get()} )
    {
        for ( auto const& ray : rays )
        {
            std::cout << "rays parts"<< "\n";
            auto rayIntersectionResult1 = bvhCurrent->intersect(_ray=ray) ;
            BOOST_CHECK_MESSAGE(!rayIntersectionResult1.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults( bvhCurrent,rayIntersectionResult1 );
        }

    
        auto multiRayIntersectionResult = bvhCurrent->intersect(_ray=rays);
        for ( auto const& rayIntersectionResult : multiRayIntersectionResult)
        {
            std::cout << "multiRayIntersectionResul parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults( bvhCurrent,rayIntersectionResult );
        }

        auto multiRayDistributedIntersectionResult = bvhCurrent->intersect(_ray=raysDistributed);
        for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
        {
            std::cout << "multiRayDistributedIntersectionResult parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults( bvhCurrent,rayIntersectionResult );
        }
    }
    */



    //BVHRay<mesh_entity_type::nRealDim> my_ray(origin2,direction_perp_2, 1e-8 );

    std::cout<<"=================================================================================================\n";
    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";

    for ( auto const& ray : rays )
        {
            std::cout << "rays parts"<< "\n";
            auto rayIntersectionResult1 = bvhInHouse->intersect(_ray=ray) ;
            BOOST_CHECK_MESSAGE(!rayIntersectionResult1.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhInHouse,rayIntersectionResult1 );
        }

    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    std::cout<<"=================================================================================================\n";


    std::cout<<"=================================================================================================\n";
    std::cout<<"*************************************************************************************************\n";

    auto multiRayDistributedIntersectionResult = bvhThirdPartyLow->intersect(_ray=raysDistributed);
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
        {
            std::cout << "multiRayDistributedIntersectionResult parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhThirdPartyLow,rayIntersectionResult);
        }

    std::cout<<"=================================================================================================\n";
    std::cout<<"=================================================================================================\n";

/*
     std::cout<<"=================================================================================================\n";
    std::cout<<"*************************************************************************************************\n";

    auto multiRayDistributedIntersectionGpuResult = bvhGpuParty->intersect(_ray=raysDistributed);
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionGpuResult )
        {
            std::cout << "GPU parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhGpuParty,rayIntersectionResult);
        }

    std::cout<<"=================================================================================================\n";
    std::cout<<"=================================================================================================\n";
*/

    //std::map<std::string,std::unique_ptr<BVH<typename tr_mesh_type::element_type>>> M_bvh_tree_vector;

    /*
    using mesh_type = MeshType;
    typedef typename MeshType::ptrtype mesh_ptrtype;
    using tr_mesh_type = typename std::conditional<MeshType::nDim==MeshType::nRealDim,
                                                typename MeshType::trace_mesh_type,
                                                typename MeshType::type >::type ;


    std::unique_ptr<BVH<typename tr_mesh_type::element_type>> M_bvh;
    */

    //auto rayIntersectionResult =  M_bvh->intersect(ray) ;
    //if ( !rayIntersectionResult.empty() ) {  std::cout<<"[INFO]: Intersect\n"; }    

}




template <typename RangeType>
void test3DWithHybrid( RangeType const& range )
{

    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHHybrid::BVHRay<mesh_entity_type::nRealDim>;

    Eigen::Vector3d origin1={1000.0,0.0,0.0};
    Eigen::Vector3d direction_perp_1={1.,0.,0.};

    Eigen::Vector3d origin2={-10.0,-0.25,-0.25};
    Eigen::Vector3d direction_perp_2={1.,0.,0.};

    Eigen::Vector3d origin3={-5.0,-0.20,-0.20};
    Eigen::Vector3d direction_perp_3={1.,0.,0.};


    std::vector<bvh_ray_type> rays;

    rays.push_back( bvh_ray_type(origin1,direction_perp_1) );
    rays.push_back( bvh_ray_type(origin2,direction_perp_2) );
    rays.push_back( bvh_ray_type(origin3,direction_perp_3) );

    BVHHybrid::BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;

    for (int k=0;k<rays.size();++k)
    {
        raysDistributed.push_back( rays[k] );
    }

    std::cout<<"\n";
    std::cout<<"\n";

    std::cout<<"=====================================================================================+============\n";
    std::cout<<"++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++\n";
    std::cout<<"[INFO]: Bounding Colume Hierarchy\n";

    //auto bvhInHouse        = BVHHybrid::boundingVolumeHierarchy(_range=range,_kind="in-house");
    //auto bvhThirdParty     = BVHHybrid::boundingVolumeHierarchy(_range=range,_kind="third-party");
    //auto bvhThirdPartyLow  = BVHHybrid::boundingVolumeHierarchy(_range=range,_kind="third-party",_quality=BVHEnum::Quality::High);
    //auto bvhGpuParty       = BVHHybrid::boundingVolumeHierarchy(_range=range,_kind="gpu-party",_quality=BVHEnum::Quality::High);
    //auto bvhHIPParty       = BVHHybrid::boundingVolumeHierarchy(_range=range,_kind="hip-party",_quality=BVHEnum::Quality::High);

    auto bvhHIPParty       = BVHHybrid::boundingVolumeHierarchy(_range=range,_kind="hip-party");

    std::cout<<"++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++\n";
    std::cout<<"========================================================================================+=========\n";

    std::cout<<"\n";
    std::cout<<"\n";

/*
    std::cout<<"=================================================================================================\n";
    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";

    for ( auto const& ray : rays )
        {
            std::cout << "rays parts"<< "\n";
            auto rayIntersectionResult1 = bvhInHouse->intersect(_ray=ray) ;
            BOOST_CHECK_MESSAGE(!rayIntersectionResult1.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhInHouse,rayIntersectionResult1 );
        }

    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    std::cout<<"=================================================================================================\n";
*/


/*
    std::cout<<"=================================================================================================\n";
    std::cout<<"*************************************************************************************************\n";

    auto multiRayDistributedIntersectionResult = bvhThirdPartyLow->intersect(_ray=raysDistributed);
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
        {
            std::cout << "multiRayDistributedIntersectionResult parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhThirdPartyLow,rayIntersectionResult);
        }

    std::cout<<"=================================================================================================\n";
    std::cout<<"=================================================================================================\n";

*/

/*
    std::cout<<"=================================================================================================\n";
    std::cout<<"*************************************************************************************************\n";


    auto multiRayDistributedIntersectionGpuResult = bvhGpuParty->intersect(_ray=raysDistributed);
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionGpuResult )
        {
            std::cout << "GPU parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhGpuParty,rayIntersectionResult);
        }

    std::cout<<"=================================================================================================\n";
    std::cout<<"=================================================================================================\n";
*/

    std::cout<<"=================================================================================================\n";
    std::cout<<"* HIP GPU ****************************************************************************************\n";


    auto multiRayDistributedIntersectionHipResult = bvhHIPParty->intersect(_ray=raysDistributed);
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionHipResult )
        {
            std::cout << "Hip parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhHIPParty,rayIntersectionResult);
        }

    std::cout<<"=================================================================================================\n";
    std::cout<<"=================================================================================================\n";
    

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//=========================================================================================================================


BOOST_AUTO_TEST_SUITE( bvh_intersection_gpu_tests )

/*
BOOST_AUTO_TEST_CASE( test_gpu_with_specx )
{
    using namespace Feel;
    std::cout<<std::endl;
    int numDevices=0;
    HIP_CHECK(hipGetDeviceCount(&numDevices));
    std::cout<<"[INFO]: Get numDevice = "<<numDevices<<"\n";
    int deviceID=0;
    HIP_CHECK(hipGetDevice(&deviceID));
    std::cout<<"[INFO]: Get deviceID activated = "<<deviceID<<"\n";

    int nbThreads=3;
    int nbElements=10;
    double c0=2.0; 
    double c1=1.5; 
    int size_array = sizeof(double) * nbElements;

    double *h_vecA = (double *)malloc(size_array);
    double *h_vecB = (double *)malloc(size_array);
    double *h_vecC = (double *)malloc(size_array);

    for (int i = 0; i < nbElements; i++)
    {
        h_vecA[i] = 1;
        h_vecB[i] = 2;
        h_vecC[i] = 0;
    }

    write_vector("Vector A",h_vecA,nbElements);
    write_vector("Vector B",h_vecB,nbElements);


    auto FC1=[size_array,nbElements](const double* vh_vecA,const double* vh_vecB,double* vh_vecC0) {  
          double *d_vecA,*d_vecB,*d_vecC;
          hipMalloc((void **)&d_vecA, size_array);
          hipMalloc((void **)&d_vecB, size_array);
          hipMalloc((void **)&d_vecC, size_array);
          hipMemcpy(d_vecA, vh_vecA, size_array, hipMemcpyHostToDevice);
          hipMemcpy(d_vecB, vh_vecB, size_array, hipMemcpyHostToDevice);
          int grid_size = (nbElements + BLOCK_SIZE - 1) / BLOCK_SIZE;
          hipLaunchKernelGGL(vector_add, grid_size,BLOCK_SIZE, 0, 0, d_vecA, d_vecB, d_vecC, nbElements);
          hipMemcpy(vh_vecC0, d_vecC, size_array, hipMemcpyDeviceToHost);
          hipFree(d_vecA); hipFree(d_vecB); hipFree(d_vecC);
        return true;
    };
    
    int numMode=3;
    //NOTA: normal configuration
        if (numMode==1) {
            FC1(h_vecA,h_vecB,h_vecC);
        }

    //NOTA: we use normal specx without gpu
        if (numMode==2) {
            SpRuntime runtime(nbThreads);
            runtime.task(SpRead(h_vecA),SpRead(h_vecB),SpWrite(h_vecC),FC1);
            runtime.waitAllTasks();
            runtime.stopAllThreads();
            runtime.generateDot("Test_hip_AMD.dot", true);
            runtime.generateTrace("Test_hip_AMD.svg");  
        }


    //NOTA: we use my task class
        if (numMode==3) {
            int numTypeThread=3;  //0:No thread  1:Std::thread 2:std:async 3:Specx  
            Task::Task TsK( nbThreads, numTypeThread );
            TsK.setSave( false );
            TsK.setInfo( false );
            TsK.setFileName( "./Test_hip_AMD" );
            TsK.add( _param(h_vecA,h_vecB,h_vecC), _tasks = FC1 );
            TsK.run();
            TsK.close();
        }



    write_vector("Vector C = (A + B)",h_vecC,nbElements);
    check_solution(1.0,h_vecA,h_vecB,h_vecC,nbElements);
    free(h_vecA); free(h_vecB); free(h_vecC);
    std::cout<<"\n";
}
*/

BOOST_AUTO_TEST_CASE( test_load_mesh3 )
{

    //TEST Thrust
    thrust::device_vector<float> d_x(10, 1.0f);

    using namespace Feel;
     using Feel::cout;   
    //typedef Mesh<Simplex<2> > mesh_type;

    using mesh_type = Mesh<Simplex<3,1,3>>; //<Dim,Order,RDim>


    //auto mesh = loadMesh(_mesh=new  mesh_type,_filename="/nvme0/lemoinep/feelppGPUBeta/feelpp/testsuite/feeltasks/cases/mesh3d/mymesh3d.geo");
    //auto mesh = loadMesh(_mesh=new  mesh_type,_filename="/nvme0/lemoinep/feelppGPUBeta/feelpp/testsuite/feelbvhgpu/cases/mesh3d/mymesh3d.geo");
    //auto mesh = loadMesh(_mesh=new  mesh_type,_filename="/nvme0/lemoinep/feelppGPUBeta/feelpp/testsuite/feelbvhgpu/cubic.geo");
    auto mesh = loadMesh(_mesh=new  mesh_type,_filename="/nvme0/lemoinep/feelppGPUBeta/feelpp/testsuite/feelbvhgpu/cases/mesh/Test.geo");

    std::cout << "[INFO] maxNumElement : "<< mesh->maxNumElements() << std::endl;
    std::cout << "[INFO] maxNumFace    : "<< mesh->maxNumFaces() << std::endl;
    std::cout << "[INFO] maxNumPoints  : "<< mesh->maxNumPoints() << std::endl;
    std::cout << "[INFO] maxNumVerices : "<< mesh->maxNumVertices() << std::endl;

    auto rangeFaces = markedfaces(mesh,{"CavityBottom","CavitySides","CavityTop","Up3","Down3","Back3","Left3","Rigth3"});
    auto submesh = createSubmesh(_mesh=mesh,_range=rangeFaces);

    

    auto Xhd0 = Pdh<0>(mesh);
    auto measures = Xhd0->element();
    measures.on(_range=elements(mesh),_expr=vf::meas());
    double measMin = measures.min();
    double measMax = measures.max();
    size_type nbdyfaces = nelements(boundaryfaces(mesh));

        std::cout << "[INFO]: mesh entities" << std::endl;
        std::cout << "[INFO]:    number of elements : " << mesh->numGlobalElements() << std::endl;
        std::cout << "[INFO]:    number of faces : " << mesh->numGlobalFaces() << std::endl;
        std::cout << "[INFO]:    number of boundary faces : " << nbdyfaces << std::endl;
        std::cout << "[INFO]:    number of points : " << mesh->numGlobalPoints() << std::endl;
        std::cout << "[INFO]:    number of vertices : " << mesh->numGlobalVertices() << std::endl;
        std::cout << "[INFO]: mesh sizes" << std::endl;
        std::cout << "[INFO]:    h max : " << mesh->hMax() << std::endl;
        std::cout << "[INFO]:    h min : " << mesh->hMin() << std::endl;
        std::cout << "[INFO]:    h avg : " << mesh->hAverage() << std::endl;
        std::cout << "[INFO]:    measure : " << mesh->measure() << "\t" << measMin<<" : " << measMax << std::endl;
        std::cout << "[INFO]: Number of Partitions : " << mesh->numberOfPartitions() << std::endl ;
        std::cout << "[INFO]: nbdyfaces : " <<nbdyfaces<<"\n";
        std::cout<<"\n";

    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    test3D( rangeFaces );
    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    //test3D( elements(submesh) );

    std::cout<<"+*-*+*-*+*-*+*+*-*+*-*+*-*+*-*+*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n";
    test3DWithHybrid( rangeFaces );
    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    std::cout<<"\n";
}

/*
BOOST_AUTO_TEST_CASE( test_thrust )
{
    // Taille des vecteurs 
    const int N = 100;
    // Allocation et initialisation des vecteurs sur le device
    thrust::device_vector<float> d_x(N, 1.0f); 
    thrust::device_vector<float> d_y(N, 2.0f);
    // Scalaire a 
    float a = 2.0f;

    for (int i = 0; i < d_x.size(); i++) { std::cout << d_x[i] << " "; }
    std::cout << std::endl;

    for (int i = 0; i < d_y.size(); i++) { std::cout << d_y[i] << " "; }
    std::cout << std::endl;

    // Exécution de l'opération SAXPY
    thrust::transform(thrust::device,
        d_x.begin(), 
        d_x.end(), 
        d_y.begin(), 
        d_y.begin(), 
        saxpy_functor(a));

    // Affichage des résultats
    for (int i = 0; i < d_y.size(); i++) { std::cout << d_y[i] << " "; }
    std::cout << std::endl;

    std::cout << "[INFO]: WELL DONE :-) Test_thrust!"<<"\n";
}
*/

/*
BOOST_AUTO_TEST_CASE( test_thrust_bvh_hip_amd_gpu )
{
    using namespace bvhhipTestIfCompile;
    std::string filename;
    filename = "/nvme0/lemoinep/feelppGPUBeta/feelpp/testsuite/feelbvhgpu/cases/Triangle2Cube.obj";
    std::cout<<"\n";
    std::cout<<"[INFO]: BEGIN BVH WITH TRIANGLE\n";
    std::cout<<"[INFO]: LOAD SCENE >>>"<<filename<<"\n";

    // Load the mesh
    std::vector<F3Triangle> hostTriangles;
    loadOBJTriangle(filename, hostTriangles);
    //getchar();
      

  
    // Transfer triangles to GPU
    thrust::device_vector<F3Triangle> deviceTriangles = hostTriangles;

    // Building the BVH
    thrust::device_vector<BVHNode> deviceNodes;
    buildBVHWithTriangleVersion1(deviceTriangles, deviceNodes);

    std::cout<<"[INFO]: BVH built with " << deviceNodes.size() << " nodes" << std::endl;

    // Generate test several rays
    std::cout<<"[INFO]: GENERATE TEST SEVERAL RAYS\n";
    const int numRays = 1024;

    thrust::host_vector<float3> hostIntersectionPoint(numRays);
    
    thrust::host_vector<F3Ray> hostRays(numRays);
    for (int i = 0; i < numRays; ++i) {
        hostRays[i].origin = make_float3(0, 0, 10);
        hostRays[i].direction = make_float3(
            (float)rand() / RAND_MAX * 2 - 1,
            (float)rand() / RAND_MAX * 2 - 1,
            -1
        );
        hostRays[i].direction = normalize(hostRays[i].direction);
        hostIntersectionPoint[i]=make_float3(INFINITY, INFINITY, INFINITY);
    }
    thrust::device_vector<F3Ray>  deviceRays              = hostRays;
    thrust::device_vector<float3> deviceIntersectionPoint = hostIntersectionPoint;

    // Allocate memory for the results
    thrust::device_vector<int>   deviceHitResults(numRays);
    thrust::device_vector<float> deviceDistanceResults(numRays);

    // Start the ray tracing kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (numRays + threadsPerBlock - 1) / threadsPerBlock;
    rayTracingKernel<<<blocksPerGrid, threadsPerBlock>>>(
        thrust::raw_pointer_cast(deviceNodes.data()),
        thrust::raw_pointer_cast(deviceTriangles.data()),
        thrust::raw_pointer_cast(deviceRays.data()),
        thrust::raw_pointer_cast(deviceHitResults.data()),
        thrust::raw_pointer_cast(deviceDistanceResults.data()),
        thrust::raw_pointer_cast(deviceIntersectionPoint.data()),
        numRays
    );

    // Retrieve the results
    thrust::host_vector<int>   hostHitResults = deviceHitResults;
    thrust::host_vector<float> hostDistanceResults = deviceDistanceResults;
    hostIntersectionPoint=deviceIntersectionPoint;

    // Count intersections
    int hitCount = thrust::count_if(hostHitResults.begin(), hostHitResults.end(),thrust::placeholders::_1 != -1);
    std::cout<<"[INFO]: Number of rays that intersected the mesh : " << hitCount << " / " << numRays << std::endl;
    for (int i=0;i<hostHitResults.size();++i)
    {
        if (hostHitResults[i]!=-1)
        {
            std::cout<<"      Intersection found with Num Ray ["<<i<<"] ori= <"<<hostRays[i].origin.x<<","<<hostRays[i].origin.y<<","<<hostRays[i].origin.z<<"> ";
            std::cout<<" dir= <"<<hostRays[i].direction.x<<","<<hostRays[i].direction.y<<","<<hostRays[i].direction.z<<"> ";
            std::cout<<" Hit= "<<hostHitResults[i]<<" min dist="<<hostDistanceResults[i];
            std::cout<<" IntersectionPoint= <"<<hostIntersectionPoint[i].x<<","<<hostIntersectionPoint[i].y<<","<<hostIntersectionPoint[i].z<<"> "<<"\n";
        }
    }
    std::cout<<"\n";
    std::cout<<"[INFO]: END BVH WITH TRIANGLE\n";
    std::cout<<"\n";


}
*/



BOOST_AUTO_TEST_SUITE_END()




