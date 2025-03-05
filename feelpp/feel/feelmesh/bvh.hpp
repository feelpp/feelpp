
#pragma clang diagnostic ignored "-Wunused-result"

#include <vector>

#include <bvh/v2/bvh.h>
#include <bvh/v2/default_builder.h>
#include <bvh/v2/node.h>
#include <bvh/v2/stack.h>
#include <bvh/v2/stream.h>
#include <bvh/v2/tri.h>

#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>

#include <algorithm>
#include <assert.h>
#include <cfloat>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#if defined( FEELPP_HAS_HIP )
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>

#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/random.h>
#include <thrust/sort.h>
#include <thrust/transform.h>

#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/generate.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/partition.h>
#include <thrust/sort.h>
#include <thrust/system/hip/vector.h>

#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/reduce.h>
#endif // FEELPP_HAS_HIP

#include <atomic>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stack>

// BEGIN:: FOR INFORMATION
// Internet sources of inspiration
//    https://en.wikipedia.org/wiki/Orthant
//    https://github.com/madmann91/bvh/blob/master/src/bvh/v2/ray.h
//    https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html
//    https://github.com/scratchapixel/scratchapixel-code/tree/main/ray-tracing-rendering-a-triangle
//    https://github.com/ToruNiina/lbvh
// END:: FOR INFORMATION

#ifdef FEELPP_HAS_HIP

#define HIP_CHECK( command )                                                                \
    {                                                                                       \
        hipError_t status = command;                                                        \
        if ( status != hipSuccess )                                                         \
        {                                                                                   \
            std::cerr << "Error: HIP reports " << hipGetErrorString( status ) << std::endl; \
            std::abort();                                                                   \
        }                                                                                   \
    }

// Conditional assertion based on debug mode for HIP operations
#ifdef NDEBUG
#define HIP_ASSERT( x ) x // No assertion in release mode
#else
#define HIP_ASSERT( x ) ( assert( ( x ) == hipSuccess ) ) // Assert in debug mode
#endif

#include "lbvh/aabb.cuh"
#include "lbvh/bvh.cuh"
#include "lbvh/morton_code.cuh"
#include "lbvh/predicator.cuh"
#include "lbvh/query.cuh"
#include "lbvh/utility.cuh"

namespace bvhHip
{

struct Vec3
{
    float x, y, z;

    __host__ __device__ __inline__ Vec3()
        : x( 0 ), y( 0 ), z( 0 ) {}
    __host__ __device__ __inline__ Vec3( float a )
        : x( a ), y( a ), z( a ) {}
    __host__ __device__ __inline__ Vec3( float x, float y, float z )
        : x( x ), y( y ), z( z ) {}

    __host__ __device__ __inline__ Vec3 operator+( const Vec3& v ) const { return Vec3( x + v.x, y + v.y, z + v.z ); }
    __host__ __device__ __inline__ Vec3 operator-( const Vec3& v ) const { return Vec3( x - v.x, y - v.y, z - v.z ); }
    __host__ __device__ __inline__ Vec3 operator*( float f ) const { return Vec3( x * f, y * f, z * f ); }
    __host__ __device__ __inline__ Vec3 operator*( const Vec3& v ) const { return Vec3( x * v.x, y * v.y, z * v.z ); }
    __host__ __device__ __inline__ Vec3 operator/( const Vec3& v ) const { return Vec3( x / v.x, y / v.y, z / v.z ); }
    __host__ __device__ __inline__ Vec3 operator/( float f ) const
    {
        float inv = 1.0f / f;
        return Vec3( x * inv, y * inv, z * inv );
    }

    __host__ __device__ __inline__ float& operator[]( size_t i ) { return ( &x )[i]; }
    __host__ __device__ __inline__ const float& operator[]( size_t i ) const { return ( &x )[i]; }
};

__host__ __device__ __inline__ Vec3 min( const Vec3& a, const Vec3& b )
{
    return Vec3( fminf( a.x, b.x ), fminf( a.y, b.y ), fminf( a.z, b.z ) );
}

__host__ __device__ __inline__ Vec3 max( const Vec3& a, const Vec3& b )
{
    return Vec3( fmaxf( a.x, b.x ), fmaxf( a.y, b.y ), fmaxf( a.z, b.z ) );
}

__host__ __device__ __inline__ Vec3 cross( const Vec3& a, const Vec3& b )
{
    return Vec3( a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x );
}

__host__ __device__ __inline__ float dot( const Vec3& a, const Vec3& b )
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

struct Ray
{
    Vec3 origin, direction;
    size_t id;
};

struct Triangle
{
    Vec3 v0, v1, v2;
    size_t id;
};

struct AABB
{
    Vec3 min, max;
};

struct BVHNode
{
    AABB bounds;
    size_t leftChild;
    size_t rightChild;
    size_t firstTriangleIndex;
    size_t triangleCount;

    size_t firstPrimitive;
    size_t primitiveCount;

    size_t triangleIndex;
};

struct Intersection
{
    bool hit;
    float t;
    size_t triangleIndex;
};

struct TriangleInfo
{
    Vec3 centroid;
    size_t index;
};

__host__ __device__ __inline__ float angleScalar( const Vec3 v1, const Vec3 v2 )
{
    float p = ( v1.x ) * ( v2.x ) + ( v1.y ) * ( v2.y ) + ( v1.z ) * ( v2.z );
    float n1 = sqrt( v1.x * v1.x + v1.y * v1.y + v1.z * v1.z );
    float n2 = sqrt( v2.x * v2.x + v2.y * v2.y + v2.z * v2.z );
    float d = n1 * n2;
    float res = 0.0f;
    if ( d > 0.0f )
    {
        float r = p / d;
        if ( r > 1.0f ) r = 1.0f;
        res = acos( r );
    }
    return ( res ); // in radian
}

__host__ __device__ __inline__ float calculateHalfOpeningAngle( const Triangle& triangle, const Vec3& origin )
{
    // This function will be used to speed up the calculations and will adapt the limit angle of the sameDirection function
    Vec3 barycenter = {
        ( triangle.v0.x + triangle.v1.x + triangle.v2.x ) / 3.0f,
        ( triangle.v0.y + triangle.v1.y + triangle.v2.y ) / 3.0f,
        ( triangle.v0.z + triangle.v1.z + triangle.v2.z ) / 3.0f };
    float distance = sqrt( pow( barycenter.x - origin.x, 2 ) +
                           pow( barycenter.y - origin.y, 2 ) +
                           pow( barycenter.z - origin.z, 2 ) );
    Vec3 edge1 = { triangle.v1.x - triangle.v0.x, triangle.v1.y - triangle.v0.y, triangle.v1.z - triangle.v0.z };
    Vec3 edge2 = { triangle.v2.x - triangle.v0.x, triangle.v2.y - triangle.v0.y, triangle.v2.z - triangle.v0.z };
    Vec3 cross = {
        edge1.y * edge2.z - edge1.z * edge2.y,
        edge1.z * edge2.x - edge1.x * edge2.z,
        edge1.x * edge2.y - edge1.y * edge2.x };
    float area = 0.5f * sqrt( cross.x * cross.x + cross.y * cross.y + cross.z * cross.z );
    float solidAngle = area / ( distance * distance );
    float halfOpeningAngle = asin( sqrt( solidAngle / ( 4 * M_PI ) ) );
    return halfOpeningAngle;
}

__host__ __device__ __inline__ bool sameDirection( Triangle& tri, Ray& ray, const float& angleLim )
{ // To be modified soon according to the radius of the triangle object
    Vec3 dT;
    dT.x = ( tri.v0.x + tri.v1.x + tri.v2.x ) / 3.0f - ray.origin.x;
    dT.y = ( tri.v0.y + tri.v1.y + tri.v2.y ) / 3.0f - ray.origin.y;
    dT.z = ( tri.v0.z + tri.v1.z + tri.v2.z ) / 3.0f - ray.origin.z;
    float angle1 = fabs( angleScalar( dT, ray.direction ) );
    float angle2 = calculateHalfOpeningAngle( tri, ray.origin );
    // return ( (angle1 <= angleLim) && (angleLim<=angle2 ) );
    // return ( (angle1 <= angleLim) ); // TODO : to define the conditions well
    return ( angle1 <= angleLim ) && ( angle1 <= angle2 );
}

__host__ __device__ __inline__ bool sameDirectionTest( Triangle& tri, Ray& ray, const float& angleLim )
{
    Vec3 dT0 = tri.v0 - ray.origin;
    Vec3 dT1 = tri.v1 - ray.origin;
    Vec3 dT2 = tri.v2 - ray.origin;
    Vec3 dT3;
    dT3.x = ( tri.v0.x + tri.v1.x + tri.v2.x ) / 3.0f - ray.origin.x;
    dT3.y = ( tri.v0.y + tri.v1.y + tri.v2.y ) / 3.0f - ray.origin.y;
    dT3.z = ( tri.v0.z + tri.v1.z + tri.v2.z ) / 3.0f - ray.origin.z;

    bool b0 = ( fabs( angleScalar( dT0, ray.direction ) ) <= angleLim );
    bool b1 = ( fabs( angleScalar( dT1, ray.direction ) ) <= angleLim );
    bool b2 = ( fabs( angleScalar( dT2, ray.direction ) ) <= angleLim );
    bool b3 = ( fabs( angleScalar( dT3, ray.direction ) ) <= angleLim );

    return ( b0 || b1 || b2 || b3 );
}

__device__ __inline__ void swap( float& a, float& b )
{
    float temp = a;
    a = b;
    b = temp;
}

// BEGIN::RAY TRACING
__device__ __inline__ bool rayTriangleIntersect( const Ray& ray, const Triangle& tri, float& t, Vec3& intersectionPoint )
{
    Vec3 edge1 = tri.v1 - tri.v0;
    Vec3 edge2 = tri.v2 - tri.v0;
    Vec3 h = cross( ray.direction, edge2 );
    float a = dot( edge1, h );
    const float EPSILON = 1e-8f;
    if ( fabs( a ) < EPSILON ) return false;

    float f = 1.0f / a;
    Vec3 s = ray.origin - tri.v0;
    float u = f * dot( s, h );

    if ( u < 0.0f || u > 1.0f ) return false;

    Vec3 q = cross( s, edge1 );
    float v = f * dot( ray.direction, q );

    if ( v < 0.0f || u + v > 1.0f ) return false;

    t = f * dot( edge2, q );

    // If t is negative, the intersection is behind the origin of the ray
    // Which means that the origin is inside the triangle
    if ( t < -EPSILON )
    {
        t = 0.0f;
        intersectionPoint = ray.origin;
        return true;
    }

    // If t is very close to zero, consider that the origin is on the triangle
    if ( fabs( t ) < EPSILON )
    {
        intersectionPoint = ray.origin;
        return true;
    }

    // Normal intersection in front of the ray origin
    if ( t > EPSILON )
    {
        intersectionPoint.x = ray.origin.x + t * ray.direction.x;
        intersectionPoint.y = ray.origin.y + t * ray.direction.y;
        intersectionPoint.z = ray.origin.z + t * ray.direction.z;
        return true;
    }

    intersectionPoint.x = INFINITY;
    intersectionPoint.y = INFINITY;
    intersectionPoint.z = INFINITY;

    return false;
}

__device__ __inline__ bool rayTriangleIntersectSurfaceEdge( const Ray& ray, const Triangle& tri, float& t, Vec3& intersectionPoint, bool* hitEdge, bool* hitVertex )
{
    // This will solve the problem of intersection of radius and vertex or edge of the triangle.
    Vec3 edge1 = tri.v1 - tri.v0;
    Vec3 edge2 = tri.v2 - tri.v0;
    Vec3 h = cross( ray.direction, edge2 );
    float a = dot( edge1, h );

    // Check if the ray is parallel to the triangle
    if ( a > -1e-6f && a < 1e-6f ) return false;

    float f = 1.0f / a;
    Vec3 s = ray.origin - tri.v0;
    float u = f * dot( s, h );

    // Checking barycentric coordinates
    if ( u < -1e-6f || u > 1.0f ) return false;

    Vec3 q = cross( s, edge1 );
    float v = f * dot( ray.direction, q );

    if ( v < -1e-6f || u + v > 1.0f + 1e-6f ) return false;

    // Calculation of t (distance to intersection)
    t = f * dot( edge2, q );

    // Check if the intersection is in front of the ray
    if ( t >= 0 )
    {
        intersectionPoint.x = ray.origin.x + t * ray.direction.x;
        intersectionPoint.y = ray.origin.y + t * ray.direction.y;
        intersectionPoint.z = ray.origin.z + t * ray.direction.z;

        // Check if the intersection occurs on an edge
        *hitEdge = ( u <= 1e-6f || v <= 1e-6f || u + v >= 1.0f - 1e-6f );

        // Check if the intersection occurs on a vertex
        float epsilon = 1e-6f; // Tolerance to determine if we touch a vertex

        auto distanceSquared = []( const Vec3& a, const Vec3& b )
        {
            float dx = a.x - b.x;
            float dy = a.y - b.y;
            float dz = a.z - b.z;
            return dx * dx + dy * dy + dz * dz;
        };

        *hitVertex = ( distanceSquared( intersectionPoint, tri.v0 ) < epsilon * epsilon ||
                       distanceSquared( intersectionPoint, tri.v1 ) < epsilon * epsilon ||
                       distanceSquared( intersectionPoint, tri.v2 ) < epsilon * epsilon );

        return true;
    }
    else
    {
        intersectionPoint.x = INFINITY;
        intersectionPoint.y = INFINITY;
        intersectionPoint.z = INFINITY;
        return false;
    }
}

__device__ bool rayAABBIntersect0( const Ray& ray, const AABB& aabb )
{
    Vec3 invDir = Vec3( 1.0f / ray.direction.x, 1.0f / ray.direction.y, 1.0f / ray.direction.z );
    Vec3 tMin = ( aabb.min - ray.origin ) * invDir;
    Vec3 tMax = ( aabb.max - ray.origin ) * invDir;
    Vec3 t1 = Vec3( fminf( tMin.x, tMax.x ), fminf( tMin.y, tMax.y ), fminf( tMin.z, tMax.z ) );
    Vec3 t2 = Vec3( fmaxf( tMin.x, tMax.x ), fmaxf( tMin.y, tMax.y ), fmaxf( tMin.z, tMax.z ) );
    float tNear = fmaxf( fmaxf( t1.x, t1.y ), t1.z );
    float tFar = fminf( fminf( t2.x, t2.y ), t2.z );
    return tNear <= tFar;
}

__device__ __inline__ bool rayAABBIntersect( const Ray& ray, const AABB& aabb )
{
    Vec3 invDir = Vec3( 1.0f / ray.direction.x, 1.0f / ray.direction.y, 1.0f / ray.direction.z );
    float tx1 = ( aabb.min.x - ray.origin.x ) * invDir.x;
    float tx2 = ( aabb.max.x - ray.origin.x ) * invDir.x;

    float tmin = fminf( tx1, tx2 );
    float tmax = fmaxf( tx1, tx2 );

    float ty1 = ( aabb.min.y - ray.origin.y ) * invDir.y;
    float ty2 = ( aabb.max.y - ray.origin.y ) * invDir.y;

    tmin = fmaxf( tmin, fminf( ty1, ty2 ) );
    tmax = fminf( tmax, fmaxf( ty1, ty2 ) );

    float tz1 = ( aabb.min.z - ray.origin.z ) * invDir.z;
    float tz2 = ( aabb.max.z - ray.origin.z ) * invDir.z;

    tmin = fmaxf( tmin, fminf( tz1, tz2 ) );
    tmax = fminf( tmax, fmaxf( tz1, tz2 ) );

    return tmax >= tmin;
}

__global__ void raytraceKernel(
    Ray* rays,
    size_t numRays,
    BVHNode* bvhNodes,
    Triangle* triangles,
    int* hitTriangles,
    float* distance,
    Vec3* intersectionPoint,
    int* hitId )

{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if ( idx >= numRays ) return;

    Ray ray = rays[idx];
    size_t stack[64];
    size_t stackPtr = 0;
    stack[stackPtr++] = 0;

    float closestT = INFINITY;
    size_t closestTriangle = -1;
    size_t closesIntersectionId = -1;

    Vec3 intersectionPointT;
    intersectionPointT.x = INFINITY;
    intersectionPointT.y = INFINITY;
    intersectionPointT.z = INFINITY;
    Vec3 closestIntersectionPoint;
    closestIntersectionPoint.x = INFINITY;
    closestIntersectionPoint.y = INFINITY;
    closestIntersectionPoint.z = INFINITY;

    const float angleLim = 0.6f;

    bool isView = false; // isView = true;

    while ( stackPtr > 0 )
    {
        size_t nodeIdx = stack[--stackPtr];
        BVHNode& node = bvhNodes[nodeIdx];

        if ( !rayAABBIntersect( ray, node.bounds ) ) continue;

        if ( node.triangleCount > 0 )
        {
            for ( size_t i = 0; i < node.triangleCount; ++i )
            {
                Triangle& tri = triangles[node.firstTriangleIndex + i];
                // if (sameDirection(tri,ray,angleLim))
                if ( sameDirectionTest( tri, ray, angleLim ) )
                {
                    float t;
                    if ( rayTriangleIntersect( ray, tri, t, intersectionPointT ) )
                    {

                        // if ( isView ) printf( "      Num Ray[%i] <%f %f %f>\n", idx, intersectionPointT.x, intersectionPointT.y, intersectionPointT.z );
                        if ( isView ) printf( "      Num Ray[%zi] <%f %f %f>\n", idx, intersectionPointT.x, intersectionPointT.y, intersectionPointT.z );
                        if ( t < closestT )
                        {
                            closestT = t;
                            closestTriangle = node.firstTriangleIndex + i;
                            closestIntersectionPoint = intersectionPointT;
                            // closesIntersectionId = triangles[closestTriangle].id;
                            closesIntersectionId = tri.id;
                        }
                    }
                }
            }
        }
        else
        {
            stack[stackPtr++] = node.leftChild;
            stack[stackPtr++] = node.rightChild;
        }
    }

    hitTriangles[idx] = closestTriangle;
    distance[idx] = closestT;
    intersectionPoint[idx] = closestIntersectionPoint;
    hitId[idx] = closesIntersectionId;
}

__global__ void raytraceKernel2(
    Ray* rays,
    size_t numRays,
    BVHNode* bvhNodes,
    Triangle* triangles,
    int* hitTriangles,
    float* distance,
    Vec3* intersectionPoint,
    int* hitId )

{

    // hipEvent_t start, stop;
    // hipEventCreate(&start);

    // an improved version to avoid errors if the stack size is insufficient and the mesh is large.
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if ( idx >= numRays ) return;

    Ray ray = rays[idx];
    size_t stack[64];
    size_t stackPtr = 0;
    stack[stackPtr++] = 0;

    float closestT = INFINITY;
    size_t closestTriangle = -1;
    size_t closesIntersectionId = -1;

    Vec3 intersectionPointT;
    intersectionPointT.x = INFINITY;
    intersectionPointT.y = INFINITY;
    intersectionPointT.z = INFINITY;
    Vec3 closestIntersectionPoint;
    closestIntersectionPoint.x = INFINITY;
    closestIntersectionPoint.y = INFINITY;
    closestIntersectionPoint.z = INFINITY;

    const float angleLim = 0.6f;

    bool isView = false; // isView = true;

    while ( stackPtr > 0 && stackPtr < 64 )
    {
        size_t nodeIdx = stack[--stackPtr];
        BVHNode& node = bvhNodes[nodeIdx];

        if ( !rayAABBIntersect( ray, node.bounds ) ) continue;

        if ( node.triangleCount > 0 )
        {
            for ( size_t i = 0; i < node.triangleCount; ++i )
            {
                Triangle& tri = triangles[node.firstTriangleIndex + i];

                if ( sameDirectionTest( tri, ray, angleLim ) )
                {
                    float t;
                    if ( rayTriangleIntersect( ray, tri, t, intersectionPointT ) )
                    {

                        // if ( isView ) printf( "      Num Ray[%i] <%f %f %f>\n", idx, intersectionPointT.x, intersectionPointT.y, intersectionPointT.z );
                        if ( isView ) printf( "      Num Ray[%zu] <%f %f %f>\n", idx, intersectionPointT.x, intersectionPointT.y, intersectionPointT.z );
                        if ( t < closestT )
                        {
                            closestT = t;
                            closestTriangle = node.firstTriangleIndex + i;
                            closestIntersectionPoint = intersectionPointT;
                            closesIntersectionId = triangles[closestTriangle].id;
                        }
                    }
                }
            }
        }
        else if ( stackPtr < 62 )
        {
            stack[stackPtr++] = node.leftChild;
            stack[stackPtr++] = node.rightChild;
        }
    }

    hitTriangles[idx] = closestTriangle;
    distance[idx] = closestT;
    intersectionPoint[idx] = closestIntersectionPoint;
    hitId[idx] = closesIntersectionId;
}

__device__ __inline__ bool rayTriangleIntersect4( const Ray& ray, const Triangle& tri,
                                                  float& t, Vec3& intersectionPoint )
{
    Vec3 edge1 = tri.v1 - tri.v0;
    Vec3 edge2 = tri.v2 - tri.v0;
    Vec3 h = cross( ray.direction, edge2 );
    float a = dot( edge1, h );
    const float EPSILON = 1e-8f;
    if ( fabs( a ) < EPSILON )
        return false; // Rayon parallèle au triangle

    float f = 1.0f / a;
    Vec3 s = ray.origin - tri.v0;
    float u = f * dot( s, h );

    if ( u < 0.0f || u > 1.0f )
        return false;

    Vec3 q = cross( s, edge1 );
    float v = f * dot( ray.direction, q );

    if ( v < 0.0f || u + v > 1.0f )
        return false;

    t = f * dot( edge2, q );

    // If t is negative, the intersection is behind the origin of the ray
    // Which means that the origin is inside the triangle
    /*
     if (t < -EPSILON) {
       t = 0.0f;
       intersectionPoint = ray.origin;
       return true;
     }
     */

    // If t is very close to zero, consider that the origin is on the triangle

    /*
    if (fabs(t) < EPSILON) {
      intersectionPoint = ray.origin;
      return true;
    }
    */

    // Normal intersection in front of the ray origin
    // if (t > EPSILON) {
    if ( t >= 0.0f )
    {
        intersectionPoint.x = ray.origin.x + t * ray.direction.x;
        intersectionPoint.y = ray.origin.y + t * ray.direction.y;
        intersectionPoint.z = ray.origin.z + t * ray.direction.z;
        return true;
    }

    intersectionPoint.x = INFINITY;
    intersectionPoint.y = INFINITY;
    intersectionPoint.z = INFINITY;

    return false;
}

__device__ bool rayAABBIntersect4Old( const Ray& ray, const AABB& aabb )
{

    bool isView = false; // isView=true;
    const float EPSILON = 1e-8f;
    Vec3 invDir = Vec3( 1.0f / ray.direction.x, 1.0f / ray.direction.y,
                        1.0f / ray.direction.z );
    Vec3 tMin = ( aabb.min - ray.origin ) * invDir;
    Vec3 tMax = ( aabb.max - ray.origin ) * invDir;
    Vec3 t1 = min( tMin, tMax );
    Vec3 t2 = max( tMin, tMax );
    float tNear = fmax( fmax( t1.x, t1.y ), t1.z );
    float tFar = fmin( fmin( t2.x, t2.y ), t2.z );

    if ( isView )
    {
        printf( "Ray origin: (%f, %f, %f), direction: (%f, %f, %f)\n", ray.origin.x,
                ray.origin.y, ray.origin.z, ray.direction.x, ray.direction.y,
                ray.direction.z );
        printf( "AABB min: (%f, %f, %f), max: (%f, %f, %f)\n", aabb.min.x, aabb.min.y,
                aabb.min.z, aabb.max.x, aabb.max.y, aabb.max.z );
        printf( "tNear: %f, tFar: %f\n", tNear, tFar );
    }

    if ( tFar < 0 )
    {
        return tNear < 0;
    }
    return tNear <= tFar + EPSILON;
}

__device__ __inline__ bool rayAABBIntersect4( const Ray& ray, const AABB& aabb )
{
    const float EPSILON = 1e-8f;
    float3 invDir;
    invDir.x = 1.0f / ray.direction.x;
    invDir.y = 1.0f / ray.direction.y;
    invDir.z = 1.0f / ray.direction.z;

    float3 t1, t2;
    t1.x = ( aabb.min.x - ray.origin.x ) * invDir.x;
    t1.y = ( aabb.min.y - ray.origin.y ) * invDir.y;
    t1.z = ( aabb.min.z - ray.origin.z ) * invDir.z;
    t2.x = ( aabb.max.x - ray.origin.x ) * invDir.x;
    t2.y = ( aabb.max.y - ray.origin.y ) * invDir.y;
    t2.z = ( aabb.max.z - ray.origin.z ) * invDir.z;

    float tNear = fmaxf( fmaxf( fminf( t1.x, t2.x ), fminf( t1.y, t2.y ) ), fminf( t1.z, t2.z ) );
    float tFar = fminf( fminf( fmaxf( t1.x, t2.x ), fmaxf( t1.y, t2.y ) ), fmaxf( t1.z, t2.z ) );

    return tFar >= 0.0f && tNear <= tFar + EPSILON;
}

__global__ void
raytraceKernel_Parallel( Ray* rays, size_t numRays, BVHNode* bvhNodes,
                         Triangle* triangles, int* hitTriangles,
                         float* distance, Vec3* intersectionPoint, int* hitId )

{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if ( idx >= numRays )
        return;

    Ray ray = rays[idx];
    size_t stack[64];
    size_t stackPtr = 0;
    // On commence par la racine du BVH
    stack[stackPtr++] = 0;

    float closestT = INFINITY;
    size_t closestTriangle = -1;
    size_t closesIntersectionId = -1;

    Vec3 intersectionPointT;
    intersectionPointT.x = INFINITY;
    intersectionPointT.y = INFINITY;
    intersectionPointT.z = INFINITY;
    Vec3 closestIntersectionPoint;
    closestIntersectionPoint.x = INFINITY;
    closestIntersectionPoint.y = INFINITY;
    closestIntersectionPoint.z = INFINITY;

    bool isView = false; // isView = true;

    const float angleLim = 0.6f;

    while ( stackPtr > 0 )
    {
        size_t nodeIdx = stack[--stackPtr];
        BVHNode& node = bvhNodes[nodeIdx];

        // if (nodeIdx < 0 || nodeIdx >= numRays) continue;

        if ( nodeIdx < 0 )
            continue;

        // if (!rayAABBIntersect4(ray, node.bounds)) continue;

        if ( node.triangleCount == 1 )
        {
            Triangle& tri = triangles[node.triangleIndex];
            // if (sameDirectionTest(tri,ray,angleLim))
            {
                float t;
                if ( rayTriangleIntersect4( ray, tri, t, intersectionPointT ) )
                {

                    // if ( isView ) printf( "      Num Ray[%i] <%f %f %f> dist=%f\n", idx, intersectionPointT.x,intersectionPointT.y, intersectionPointT.z, t );

                    if ( isView ) printf( "      Num Ray[%zu] <%f %f %f> dist=%f\n", idx, intersectionPointT.x, intersectionPointT.y, intersectionPointT.z, t );

                    if ( t < closestT )
                    {
                        closestT = t;
                        closestTriangle = node.triangleIndex;
                        closestIntersectionPoint = intersectionPointT;
                        closesIntersectionId = triangles[closestTriangle].id;

                        // if ( isView ) printf( "      Close Num Ray[%i] <%f %f %f> dist=%f\n", idx, intersectionPointT.x,intersectionPointT.y, intersectionPointT.z, t );
                        if ( isView ) printf( "      Close Num Ray[%zu] <%f %f %f> dist=%f\n", idx, intersectionPointT.x, intersectionPointT.y, intersectionPointT.z, t );
                    }
                }
            }
        }
        else
        {
            stack[stackPtr++] = node.rightChild;
            stack[stackPtr++] = node.leftChild;
        }
    }

    // if (closesIntersectionId>0) printf("      Num Ray[%i] dist=%f <%f %f
    // %f>\n", idx, closestT,intersectionPointT.x, intersectionPointT.y,
    // intersectionPointT.z);

    // if (closesIntersectionId > 0) printf("      Num Ray[%i] dist=%f\n", idx,
    // closestT);

    hitTriangles[idx] = closestTriangle;
    distance[idx] = closestT;
    intersectionPoint[idx] = closestIntersectionPoint;
    hitId[idx] = closesIntersectionId;
}

__global__ void raytraceKernel_Parallel3( Ray* rays, size_t numRays, BVHNode* bvhNodes,
                                          size_t numNodes,
                                          Triangle* triangles, size_t numTriangles, int* hitTriangles,
                                          float* distance, Vec3* intersectionPoint, int* hitId )
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if ( idx >= numRays || numNodes <= 0 || numTriangles <= 0 )
        return;

    Ray ray = rays[idx];
    constexpr size_t MAX_STACK_SIZE = 64;
    size_t stack[MAX_STACK_SIZE];
    size_t stackPtr = 0;

    stack[stackPtr++] = 0;

    float closestT = INFINITY;
    size_t closestTriangle = -1;
    size_t closestIntersectionId = -1;
    Vec3 intersectionPointT;
    intersectionPointT.x = INFINITY;
    intersectionPointT.y = INFINITY;
    intersectionPointT.z = INFINITY;

    Vec3 closestIntersectionPoint;
    closestIntersectionPoint.x = INFINITY;
    closestIntersectionPoint.y = INFINITY;
    closestIntersectionPoint.z = INFINITY;

    constexpr float angleLim = 0.6f;

    while ( stackPtr > 0 && stackPtr < MAX_STACK_SIZE )
    {
        size_t nodeIdx = stack[--stackPtr];
        if ( nodeIdx < 0 || nodeIdx >= numNodes )
            continue;

        BVHNode& node = bvhNodes[nodeIdx];

        if ( node.triangleCount == 1 && node.triangleIndex >= 0 && node.triangleIndex < numTriangles )
        {

            if ( !rayAABBIntersect4( ray, node.bounds ) ) // add test 1
                continue;

            Triangle& tri = triangles[node.triangleIndex];
            // if (sameDirectionTest(tri, ray, angleLim)) //add test 2
            {
                float t;
                Vec3 intersectionPointT;
                if ( rayTriangleIntersect4( ray, tri, t, intersectionPointT ) )
                {
                    if ( t < closestT )
                    {
                        closestT = t;
                        closestTriangle = node.triangleIndex;
                        closestIntersectionPoint = intersectionPointT;
                        closestIntersectionId = tri.id;
                    }
                }
            }
        }
        else if ( stackPtr < MAX_STACK_SIZE - 1 )
        {
            if ( node.leftChild >= 0 && node.leftChild < numNodes )
                stack[stackPtr++] = node.leftChild;
            if ( node.rightChild >= 0 && node.rightChild < numNodes && stackPtr < MAX_STACK_SIZE )
                stack[stackPtr++] = node.rightChild;
        }
    }
    hitTriangles[idx] = closestTriangle;
    distance[idx] = closestT;
    intersectionPoint[idx] = closestIntersectionPoint;
    hitId[idx] = closestIntersectionId;
}

// END::RAY TRACING

// BEGIN::BVH GPU
__host__ __device__ __inline__ void calculateBoundingBox( const Triangle& triangle, Vec3& min_values, Vec3& max_values )
{
    min_values = min( triangle.v0, min( triangle.v1, triangle.v2 ) );
    max_values = max( triangle.v0, max( triangle.v1, triangle.v2 ) );
}

__global__ void initializeLeaves( Triangle* triangles, BVHNode* nodes, size_t numTriangles )
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if ( idx < numTriangles )
    {
        BVHNode& node = nodes[numTriangles - 1 + idx];
        calculateBoundingBox( triangles[idx], node.bounds.min, node.bounds.max );
        node.triangleIndex = idx;
        node.leftChild = node.rightChild = -1;
        node.firstTriangleIndex = idx;
        node.triangleCount = 1;
    }
}

void buildBVH_GPU_Version2( Triangle* d_triangles, BVHNode* d_nodes,
                            int numTriangles )
{
    int totalNodes = 2 * numTriangles - 1;
    int blockSize = 1024;
    int numBlocks = ( numTriangles + blockSize - 1 ) / blockSize;
    hipLaunchKernelGGL( initializeLeaves, dim3( numBlocks ), dim3( blockSize ), 0, 0,
                        d_triangles, d_nodes, numTriangles );

    BVHNode* h_nodes = new BVHNode[2 * numTriangles - 1];
    HIP_ASSERT( hipMemcpy( h_nodes, d_nodes,
                           ( 2 * numTriangles - 1 ) * sizeof( BVHNode ),
                           hipMemcpyDeviceToHost ) );

    for ( int i = numTriangles - 2; i >= 0; --i )
    {
        BVHNode& node = h_nodes[i];
        int leftChild = 2 * i + 1;
        int rightChild = 2 * i + 2;
        node.leftChild = leftChild;
        node.rightChild = rightChild;
        node.triangleIndex = -1;

        BVHNode& leftNode = h_nodes[leftChild];
        BVHNode& rightNode = h_nodes[rightChild];
        node.bounds.min = min( leftNode.bounds.min, rightNode.bounds.min );
        node.bounds.max = max( leftNode.bounds.max, rightNode.bounds.max );
    }
    HIP_ASSERT( hipMemcpy( d_nodes, h_nodes,
                           ( 2 * numTriangles - 1 ) * sizeof( BVHNode ),
                           hipMemcpyHostToDevice ) );
    delete[] h_nodes;
}

__global__ void buildEvaluationNodes( BVHNode* nodes, int numTriangles )
{
    for ( int i = numTriangles - 2; i >= 0; --i )
    {
        BVHNode& node = nodes[i];
        int leftChild = 2 * i + 1;
        int rightChild = 2 * i + 2;
        node.leftChild = leftChild;
        node.rightChild = rightChild;
        node.triangleIndex = -1;
        BVHNode& leftNode = nodes[leftChild];
        BVHNode& rightNode = nodes[rightChild];
        node.bounds.min = min( leftNode.bounds.min, rightNode.bounds.min );
        node.bounds.max = max( leftNode.bounds.max, rightNode.bounds.max );
    }
    //__syncthreads();
}

void buildBVH_GPU_Version3( Triangle* d_triangles, BVHNode* d_nodes, int numTriangles )
{
    // size_t blockSize = 512;
    // size_t blockSize = 1024;
    // size_t numBlocks = ( numTriangles + blockSize - 1 ) / blockSize;
    int blockSize = 1024;
    int numBlocks = ( numTriangles + blockSize - 1 ) / blockSize;
    hipLaunchKernelGGL( initializeLeaves, dim3( numBlocks ), dim3( blockSize ), 0, 0, d_triangles, d_nodes, numTriangles );
    hipDeviceSynchronize();
    hipLaunchKernelGGL( buildEvaluationNodes, dim3( 1 ), dim3( 1 ), 0, 0, d_nodes, numTriangles );
    hipDeviceSynchronize();
}



__global__ void initializeLeaves4(Triangle* d_triangles, BVHNode* d_nodes, int numTriangles) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < numTriangles) {
        BVHNode& node = d_nodes[idx + numTriangles - 1]; 
        node.bounds.min = d_triangles[idx].v0; 
        node.bounds.max = d_triangles[idx].v0;
        node.bounds.min = min(node.bounds.min, d_triangles[idx].v1);
        node.bounds.max = max(node.bounds.max, d_triangles[idx].v1);
        node.bounds.min = min(node.bounds.min, d_triangles[idx].v2);
        node.bounds.max = max(node.bounds.max, d_triangles[idx].v2);
        node.triangleIndex = idx;
        node.leftChild = -1; 
        node.rightChild = -1;
    }
}


__global__ void buildEvaluationNodes4(BVHNode* nodes, int numTriangles) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < numTriangles - 1) {
        int nodeIndex = numTriangles - 2 - idx;
        BVHNode& node = nodes[nodeIndex];
        int leftChild = 2 * nodeIndex + 1;
        int rightChild = 2 * nodeIndex + 2;
        node.leftChild = leftChild;
        node.rightChild = rightChild;
        node.triangleIndex = -1;
        BVHNode& leftNode = nodes[leftChild];
        BVHNode& rightNode = nodes[rightChild];
        node.bounds.min = min(leftNode.bounds.min, rightNode.bounds.min);
        node.bounds.max = max(leftNode.bounds.max, rightNode.bounds.max);
    }
}


void buildBVH_GPU_Version4(Triangle* d_triangles, BVHNode* d_nodes, int numTriangles) {
    if (numTriangles <= 0) return; 

    int blockSize = 256; 
    int numBlocks = (numTriangles + blockSize - 1) / blockSize;


    hipLaunchKernelGGL(initializeLeaves4, dim3(numBlocks), dim3(blockSize), 0, 0, d_triangles, d_nodes, numTriangles);
    HIP_CHECK(hipGetLastError());

    if (numTriangles > 1) { 
        int evaluationBlockSize = 256;
        int evaluationNumBlocks = (numTriangles - 1 + evaluationBlockSize - 1) / evaluationBlockSize;
        hipLaunchKernelGGL(buildEvaluationNodes4, dim3(evaluationNumBlocks), dim3(evaluationBlockSize), 0, 0, d_nodes, numTriangles);
        HIP_CHECK(hipGetLastError());
    }
    HIP_CHECK(hipDeviceSynchronize()); 
}




// Bellow new versions...

__device__ __inline__ bool compareTriangles( const TriangleInfo& a, const TriangleInfo& b, size_t axis )
{
    return a.centroid[axis] < b.centroid[axis];
}

__global__ void initTriangleInfo( Triangle* triangles, TriangleInfo* triInfo, size_t numTriangles )
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if ( idx < numTriangles )
    {
        Triangle& tri = triangles[idx];
        triInfo[idx].centroid = ( tri.v0 + tri.v1 + tri.v2 ) / 3.0f;
        triInfo[idx].index = idx;
    }
}

__global__ void bitonicSort( TriangleInfo* triInfo, size_t j, size_t k, size_t numTriangles, size_t axis )
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t ixj = i ^ j;

    if ( ( ixj > i ) && ( i < numTriangles ) && ( ixj < numTriangles ) )
    {
        bool ascending = ( ( i & k ) == 0 );
        if ( compareTriangles( triInfo[i], triInfo[ixj], axis ) == ascending )
        {
            TriangleInfo temp = triInfo[i];
            triInfo[i] = triInfo[ixj];
            triInfo[ixj] = temp;
        }
    }
}

__global__ void buildBVHNodes( BVHNode* nodes, TriangleInfo* triInfo, Triangle* triangles, size_t numTriangles )
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t totalNodes = 2 * numTriangles - 1;

    if ( idx < totalNodes )
    {
        BVHNode& node = nodes[idx];

        if ( idx >= numTriangles - 1 )
        { // Leaf node
            size_t triIdx = triInfo[idx - ( numTriangles - 1 )].index;
            node.bounds.min = node.bounds.max = triangles[triIdx].v0;
            node.bounds.min = min( node.bounds.min, min( triangles[triIdx].v1, triangles[triIdx].v2 ) );
            node.bounds.max = max( node.bounds.max, max( triangles[triIdx].v1, triangles[triIdx].v2 ) );
            node.triangleIndex = triIdx;
            node.triangleCount = 1;
            node.leftChild = node.rightChild = -1;
        }
        else
        { // Internal node
            size_t leftChild = 2 * idx + 1;
            size_t rightChild = 2 * idx + 2;
            node.leftChild = leftChild;
            node.rightChild = rightChild;
            node.triangleIndex = -1;

            node.bounds.min = min( nodes[leftChild].bounds.min, nodes[rightChild].bounds.min );
            node.bounds.max = max( nodes[leftChild].bounds.max, nodes[rightChild].bounds.max );
            node.triangleCount = nodes[leftChild].triangleCount + nodes[rightChild].triangleCount;
        }
    }
}

__global__ void buildBVHNodesParallel( BVHNode* nodes, TriangleInfo* triInfo,
                                       Triangle* triangles, size_t numTriangles )
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t totalNodes = 2 * numTriangles - 1;

    if ( idx < totalNodes )
    {
        BVHNode& node = nodes[idx];

        if ( idx >= numTriangles - 1 )
        { // Leaf node
            size_t triIdx = triInfo[idx - ( numTriangles - 1 )].index;
            // Initialization of the node limits with the first top of the triangle
            node.bounds.min = node.bounds.max = triangles[triIdx].v0;

            // Calculation of AABB limits for the triangle
            node.bounds.min =
                min( node.bounds.min, min( triangles[triIdx].v1, triangles[triIdx].v2 ) );
            node.bounds.max =
                max( node.bounds.max, max( triangles[triIdx].v1, triangles[triIdx].v2 ) );
            node.triangleIndex = triIdx;
            node.triangleCount = 1;
            node.leftChild = -1;  // no children for a leaf knot
            node.rightChild = -1; // no children for a leaf knot
        }
        else
        { // Internal node
            size_t leftChild = 2 * idx + 1;
            size_t rightChild = 2 * idx + 2;

            // Verification that children exist before accessing their limits
            if ( leftChild < totalNodes && rightChild < totalNodes )
            {
                node.leftChild = leftChild;
                node.rightChild = rightChild;
                node.triangleIndex = -1; // no associated triangle

                // Calculation of AABB limits encompassing children
                node.bounds.min =
                    min( nodes[leftChild].bounds.min, nodes[rightChild].bounds.min );
                node.bounds.max =
                    max( nodes[leftChild].bounds.max, nodes[rightChild].bounds.max );
                node.triangleCount =
                    nodes[leftChild].triangleCount + nodes[rightChild].triangleCount;
            }
            else
            {
                // If children do not exist (rare case), we can reset the node with
                // default values
                node.bounds.min = Vec3{ 0.0f, 0.0f, 0.0f }; // Default values
                node.bounds.max = Vec3{ 0.0f, 0.0f, 0.0f }; // Default values
                node.triangleIndex = -1;
                node.triangleCount = 0;
                node.leftChild = -1;
                node.rightChild = -1;
            }
        }
    }
}

__global__ void computeExtents( TriangleInfo* triInfo, size_t numTriangles, float* minExtents, float* maxExtents )
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if ( idx < numTriangles )
    {
        for ( size_t axis = 0; axis < 3; ++axis )
        {
            atomicMin( &minExtents[axis], triInfo[idx].centroid[axis] );
            atomicMax( &maxExtents[axis], triInfo[idx].centroid[axis] );
        }
    }
}

void buildBVH_GPU_Parallel( Triangle* d_triangles, BVHNode* d_nodes,
                            size_t numTriangles )
{

    // std::cout << "[INFO]: buildBVH_GPU_Parallel\n";

    size_t totalNodes = 2 * numTriangles - 1;
    // size_t blockSize = 512;
    size_t blockSize = 1024;
    size_t numBlocks = ( numTriangles + blockSize - 1 ) / blockSize;

    TriangleInfo* d_triInfo;
    HIP_ASSERT( hipMalloc( &d_triInfo, numTriangles * sizeof( TriangleInfo ) ) );
    hipLaunchKernelGGL( initTriangleInfo, dim3( numBlocks ), dim3( blockSize ), 0, 0,
                        d_triangles, d_triInfo, numTriangles );

    // Sort by axis X=0
    for ( size_t k = 2; k <= numTriangles; k *= 2 )
    {
        for ( size_t j = k / 2; j > 0; j /= 2 )
        {
            hipLaunchKernelGGL( bitonicSort, dim3( numBlocks ), dim3( blockSize ), 0, 0,
                                d_triInfo, j, k, numTriangles, 0 );
        }
    }

    numBlocks = ( totalNodes + blockSize - 1 ) / blockSize;
    hipLaunchKernelGGL( buildBVHNodes, dim3( numBlocks ), dim3( blockSize ), 0, 0,
                        d_nodes, d_triInfo, d_triangles, numTriangles );

    hipFree( d_triInfo );
}

void buildBVH_GPU_Parallel_Best_Axis( Triangle* d_triangles, BVHNode* d_nodes,
                                      size_t numTriangles )
{
    std::cout << "[INFO]: buildBVH_GPU_Parallel\n";

    size_t totalNodes = 2 * numTriangles - 1;
    // size_t blockSize = 512;
    size_t blockSize = 1024;
    size_t numBlocks = ( numTriangles + blockSize - 1 ) / blockSize;

    TriangleInfo* d_triInfo;
    HIP_ASSERT( hipMalloc( &d_triInfo, numTriangles * sizeof( TriangleInfo ) ) );
    hipLaunchKernelGGL( initTriangleInfo, dim3( numBlocks ), dim3( blockSize ), 0, 0,
                        d_triangles, d_triInfo, numTriangles );

    // Compute extents
    float *d_minExtents, *d_maxExtents;
    HIP_ASSERT( hipMalloc( &d_minExtents, 3 * sizeof( float ) ) );
    HIP_ASSERT( hipMalloc( &d_maxExtents, 3 * sizeof( float ) ) );

    // Initialize extents
    float initMin = std::numeric_limits<float>::max();
    float initMax = std::numeric_limits<float>::lowest();
    HIP_ASSERT( hipMemset( d_minExtents, *reinterpret_cast<int*>( &initMin ),
                           3 * sizeof( float ) ) );
    HIP_ASSERT( hipMemset( d_maxExtents, *reinterpret_cast<int*>( &initMax ),
                           3 * sizeof( float ) ) );

    hipLaunchKernelGGL( computeExtents, dim3( numBlocks ), dim3( blockSize ), 0, 0,
                        d_triInfo, numTriangles, d_minExtents, d_maxExtents );

    // Copy extents back to host
    float h_minExtents[3], h_maxExtents[3];
    HIP_ASSERT( hipMemcpy( h_minExtents, d_minExtents, 3 * sizeof( float ),
                           hipMemcpyDeviceToHost ) );
    HIP_ASSERT( hipMemcpy( h_maxExtents, d_maxExtents, 3 * sizeof( float ),
                           hipMemcpyDeviceToHost ) );

    // Find axis with largest extent
    size_t bestAxis = 0;
    float maxExtent = h_maxExtents[0] - h_minExtents[0];
    for ( size_t axis = 1; axis < 3; ++axis )
    {
        float extent = h_maxExtents[axis] - h_minExtents[axis];
        if ( extent > maxExtent )
        {
            maxExtent = extent;
            bestAxis = axis;
        }
    }

    // Sort using the best axis
    for ( size_t k = 2; k <= numTriangles; k *= 2 )
    {
        for ( size_t j = k / 2; j > 0; j /= 2 )
        {
            hipLaunchKernelGGL( bitonicSort, dim3( numBlocks ), dim3( blockSize ), 0, 0,
                                d_triInfo, j, k, numTriangles, bestAxis );
        }
    }

    numBlocks = ( totalNodes + blockSize - 1 ) / blockSize;
    hipLaunchKernelGGL( buildBVHNodesParallel, dim3( numBlocks ), dim3( blockSize ), 0,
                        0, d_nodes, d_triInfo, d_triangles, numTriangles );

    hipFree( d_triInfo );
    hipFree( d_minExtents );
    hipFree( d_maxExtents );
}

//***********************
//***********************
//***********************

//***********************
//***********************
//***********************
// END::GPU

} // namespace bvhHip

// LBVH Explorer
namespace bvhLinearExplorer
{

struct Ray
{
    float4 origin;
    float4 direction;
    int id; // add
};

struct Triangle
{
    float4 v1, v2, v3;
    int id;
};

struct HitRay
{
    float distanceResults;
    int hitResults;
    int idResults;
    float3 intersectionPoint;
};

struct aabb_getter
{
    __device__ __inline__ lbvh::aabb<float>
    operator()( const Triangle& tri ) const noexcept
    {
        lbvh::aabb<float> retval;
        retval.lower =
            make_float4( fminf( fminf( tri.v1.x, tri.v2.x ), tri.v3.x ),
                         fminf( fminf( tri.v1.y, tri.v2.y ), tri.v3.y ),
                         fminf( fminf( tri.v1.z, tri.v2.z ), tri.v3.z ), 0.0f );
        retval.upper =
            make_float4( fmaxf( fmaxf( tri.v1.x, tri.v2.x ), tri.v3.x ),
                         fmaxf( fmaxf( tri.v1.y, tri.v2.y ), tri.v3.y ),
                         fmaxf( fmaxf( tri.v1.z, tri.v2.z ), tri.v3.z ), 0.0f );
        return retval;
    }
};

struct AABB
{
    float4 min;
    float4 max;
};

struct merge_aabb
{
    __device__ __inline__ lbvh::aabb<float>
    operator()( const lbvh::aabb<float>& a, const lbvh::aabb<float>& b ) const
    {
        return lbvh::merge( a, b );
    }
};

struct CenterGlobalSpaceBox
{
    float4 min;
    float4 max;
    float width;
    float height;
    float depth;
    float radius;
    float volume;
    float4 position;
};

__host__ __device__ __inline__ float length( const float3& v )
{
    return sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
}

__host__ __device__ __inline__ float length( const float4& v )
{
    return sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
}

__host__ __device__ __inline__ float angleScalar( const float4 v1,
                                                  const float4 v2 )
{
    float p = ( v1.x * v2.x ) + ( v1.y * v2.y ) + ( v1.z * v2.z );
    float n1 = sqrt( v1.x * v1.x + v1.y * v1.y + v1.z * v1.z );
    float n2 = sqrt( v2.x * v2.x + v2.y * v2.y + v2.z * v2.z );
    float d = n1 * n2;
    if ( d > 0.0f )
    {
        float r = p / d;
        r = fmaxf( -1.0f, fminf( 1.0f, r ) ); // Clamp r to [-1, 1]
        return acosf( r );
    }
    return 0.0f;
}

__host__ __device__ __inline__ float
calculateHalfOpeningAngle( const Triangle& triangle, const float4& origin )
{
    // This function will be used to speed up the calculations and will adapt the
    // limit angle of the sameDirection function
    float4 barycenter = { ( triangle.v1.x + triangle.v2.x + triangle.v3.x ) / 3.0f,
                          ( triangle.v1.y + triangle.v2.y + triangle.v3.y ) / 3.0f,
                          ( triangle.v1.z + triangle.v2.z + triangle.v3.z ) / 3.0f,
                          0.0f };

    float distance =
        sqrt( pow( barycenter.x - origin.x, 2 ) + pow( barycenter.y - origin.y, 2 ) +
              pow( barycenter.z - origin.z, 2 ) );

    float4 edge1 = { triangle.v2.x - triangle.v1.x, triangle.v2.y - triangle.v1.y,
                     triangle.v2.z - triangle.v1.z, 0.0f };

    float4 edge2 = { triangle.v3.x - triangle.v1.x, triangle.v3.y - triangle.v1.y,
                     triangle.v3.z - triangle.v1.z, 0.0f };

    float4 cross = { edge1.y * edge2.z - edge1.z * edge2.y,
                     edge1.z * edge2.x - edge1.x * edge2.z,
                     edge1.x * edge2.y - edge1.y * edge2.x, 0.0f };

    float area =
        0.5f * sqrt( cross.x * cross.x + cross.y * cross.y + cross.z * cross.z );
    float solidAngle = area / ( distance * distance );
    float halfOpeningAngle = asin( sqrt( solidAngle / ( 4 * M_PI ) ) );

    return halfOpeningAngle;
}

__host__ __device__ __inline__ bool
sameDirection( const Triangle& tri, const Ray& ray,
               const float& angleLim )
{ // To be modified soon according to the
  // radius of the triangle object
    float4 dT;
    dT.x = ( tri.v1.x + tri.v2.x + tri.v3.x ) / 3.0f - ray.origin.x;
    dT.y = ( tri.v1.y + tri.v2.y + tri.v3.y ) / 3.0f - ray.origin.y;
    dT.z = ( tri.v1.z + tri.v2.z + tri.v3.z ) / 3.0f - ray.origin.z;
    float angle1 = angleScalar( dT, ray.direction );
    float angle2 = calculateHalfOpeningAngle( tri, ray.origin );
    // return (angle <= angleLim);
    return ( angle1 <= angleLim ) && ( angle1 <= angle2 );
}

__host__ __device__ __inline__ bool
sameDirectionTest( const Triangle& tri, const Ray& ray, const float& angleLim )
{
    float4 dT1 = tri.v1 - ray.origin;
    float4 dT2 = tri.v2 - ray.origin;
    float4 dT3 = tri.v3 - ray.origin;
    float4 dT4;
    dT3.x = ( tri.v1.x + tri.v2.x + tri.v3.x ) / 3.0f - ray.origin.x;
    dT3.y = ( tri.v1.y + tri.v2.y + tri.v3.y ) / 3.0f - ray.origin.y;
    dT3.z = ( tri.v1.z + tri.v2.z + tri.v3.z ) / 3.0f - ray.origin.z;

    bool b1 = ( fabs( angleScalar( dT1, ray.direction ) ) <= angleLim );
    bool b2 = ( fabs( angleScalar( dT2, ray.direction ) ) <= angleLim );
    bool b3 = ( fabs( angleScalar( dT3, ray.direction ) ) <= angleLim );
    bool b4 = ( fabs( angleScalar( dT4, ray.direction ) ) <= angleLim );
    return ( b1 || b2 || b3 || b4 );
}

struct distance_calculator
{
    __device__ __inline__ float operator()( const float4 point,
                                            const Triangle& tri ) const noexcept
    {
        float4 center = make_float4( ( tri.v1.x + tri.v2.x + tri.v3.x ) / 3.0f,
                                     ( tri.v1.y + tri.v2.y + tri.v3.y ) / 3.0f,
                                     ( tri.v1.z + tri.v2.z + tri.v3.z ) / 3.0f, 0.0f );
        return ( point.x - center.x ) * ( point.x - center.x ) +
               ( point.y - center.y ) * ( point.y - center.y ) +
               ( point.z - center.z ) * ( point.z - center.z );
    }
};

__host__ __device__ __inline__ void normalizeRayDirection( Ray& ray )
{
    float len = sqrtf( ray.direction.x * ray.direction.x +
                       ray.direction.y * ray.direction.y +
                       ray.direction.z * ray.direction.z );
    if ( len > 0 )
    {
        float invLen = 1.0f / len;
        ray.direction.x *= invLen;
        ray.direction.y *= invLen;
        ray.direction.z *= invLen;
    }
}

__device__ __inline__ float3 cross( float3 a, float3 b )
{
    return make_float3( a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x );
}

__device__ __inline__ float dot( float3 a, float3 b )
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ __inline__ float3 float4_to_float3( const float4& v )
{
    return make_float3( v.x, v.y, v.z );
}

__device__ __inline__ bool
raySphereIntersection( const float4& rayOrigin, const float4& rayDirection,
                       const float4& sphereCenter, const float sphereRadius,
                       float4& intersectionPoint, float& distance )
{
    // Calculate oc (rayOrigin - sphereCenter)
    float4 oc =
        make_float4( rayOrigin.x - sphereCenter.x, rayOrigin.y - sphereCenter.y,
                     rayOrigin.z - sphereCenter.z, 0.0f );

    // Calculate a, b, c for the quadratic equation
    float a = oc.x * oc.x + oc.y * oc.y +
              oc.z * oc.z; // Incorrect usage, corrected below
    float a_correct = rayDirection.x * rayDirection.x +
                      rayDirection.y * rayDirection.y +
                      rayDirection.z * rayDirection.z;
    float b = 2.0f * ( oc.x * rayDirection.x + oc.y * rayDirection.y +
                       oc.z * rayDirection.z );
    float c =
        oc.x * oc.x + oc.y * oc.y + oc.z * oc.z - sphereRadius * sphereRadius;

    float discriminant = b * b - 4 * a_correct * c;

    if ( discriminant < 0 )
    {
        return false; // No intersection
    }

    float sqrtDiscriminant = sqrtf( discriminant );
    float t1 = ( -b - sqrtDiscriminant ) / ( 2.0f * a_correct );
    float t2 = ( -b + sqrtDiscriminant ) / ( 2.0f * a_correct );

    // Choose the smallest positive distance
    if ( t1 > 0 && t1 < t2 )
    {
        distance = t1;
    }
    else if ( t2 > 0 )
    {
        distance = t2;
    }
    else
    {
        return false; // Both intersections are behind the ray
    }

    // Calculate the intersection point
    intersectionPoint =
        make_float4( rayOrigin.x + distance * rayDirection.x,
                     rayOrigin.y + distance * rayDirection.y,
                     rayOrigin.z + distance * rayDirection.z, 0.0f );

    return true;
}

__device__ __inline__ bool
rayTriangleIntersect( const Ray& ray, const Triangle& triangle, float& t )
{
    float4 edge1 = triangle.v2 - triangle.v1;
    float4 edge2 = triangle.v3 - triangle.v1;
    float4 h =
        make_float4( ray.direction.y * edge2.z - ray.direction.z * edge2.y,
                     ray.direction.z * edge2.x - ray.direction.x * edge2.z,
                     ray.direction.x * edge2.y - ray.direction.y * edge2.x, 0 );
    float a = edge1.x * h.x + edge1.y * h.y + edge1.z * h.z;
    if ( a > -1e-8 && a < 1e-8 )
        return false;
    float f = 1.0f / a;
    float4 s = ray.origin - triangle.v1;
    float u = f * ( s.x * h.x + s.y * h.y + s.z * h.z );
    if ( u < 0.0 || u > 1.0 )
        return false;
    float4 q =
        make_float4( s.y * edge1.z - s.z * edge1.y, s.z * edge1.x - s.x * edge1.z,
                     s.x * edge1.y - s.y * edge1.x, 0 );
    float v = f * ( ray.direction.x * q.x + ray.direction.y * q.y +
                    ray.direction.z * q.z );
    if ( v < 0.0 || u + v > 1.0 )
        return false;
    t = f * ( edge2.x * q.x + edge2.y * q.y + edge2.z * q.z );
    return ( t > 1e-8 );
}

__device__ __inline__ float dot( const float4& a, const float4& b )
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ __inline__ float4 cross( const float4& a, const float4& b )
{
    return make_float4( a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x, 0.0f );
}

__device__ __inline__ bool pointInTriangle2( const float4& point,
                                             const Triangle& triangle,
                                             float epsilon = 1e-5f )
{
    float4 v0 = triangle.v2 - triangle.v1;
    float4 v1 = triangle.v3 - triangle.v1;
    float4 v2 = point - triangle.v1;
    bool isView = true;
    isView = false;

    // Check if point lies on the plane of the triangle
    float4 normal = cross( v0, v1 );
    float distanceToPlane = dot( normal, v2 ) / length( normal );

    if ( isView )
        printf( "DistanceToPlan: %f\n", distanceToPlane );

    if ( abs( distanceToPlane ) > epsilon )
    {
        return false; // Point does not lie on the plane
    }

    float d00 = dot( v0, v0 );
    float d01 = dot( v0, v1 );
    float d11 = dot( v1, v1 );
    float d20 = dot( v2, v0 );
    float d21 = dot( v2, v1 );

    float denom = d00 * d11 - d01 * d01;
    float v = ( d11 * d20 - d01 * d21 ) / denom;
    float w = ( d00 * d21 - d01 * d20 ) / denom;
    float u = 1.0f - v - w;
    if ( isView )
    {
        printf( "Point: (%f, %f, %f)\n", point.x, point.y, point.z );
        printf( "Triangle: (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n",
                triangle.v1.x, triangle.v1.y, triangle.v1.z, triangle.v2.x,
                triangle.v2.y, triangle.v2.z, triangle.v3.x, triangle.v3.y,
                triangle.v3.z );
        printf( "u: %f, v: %f, w: %f\n", u, v, w );
    }

    bool result = ( u >= -epsilon ) && ( v >= -epsilon ) && ( w >= -epsilon ) &&
                  ( u + v + w <= 1.0f + epsilon );
    if ( isView )
        printf( "Result: %s\n", result ? "true" : "false" );

    return result;
}

__device__ __inline__ float4 normalize( float4 v )
{
    float length = sqrtf( v.x * v.x + v.y * v.y + v.z * v.z );
    if ( length > 0 )
    {
        float invLength = 1.0f / length;
        v.x *= invLength;
        v.y *= invLength;
        v.z *= invLength;
    }
    return v;
}

//**************************************************************************************************************
//--------------------------------------------------------------------------------------------------------------

__host__ __device__ __inline__ void initializeDirections( float4* directions )
{
    const float c1 = 1.0f / sqrt( 3.0f );
    const float4 predefinedDirections[14] = {
        make_float4( 1.0f, 0.0f, 0.0f, 0.0f ), make_float4( -1.0f, 0.0f, 0.0f, 0.0f ),
        make_float4( 0.0f, 1.0f, 0.0f, 0.0f ), make_float4( 0.0f, -1.0f, 0.0f, 0.0f ),
        make_float4( 0.0f, 0.0f, 1.0f, 0.0f ), make_float4( 0.0f, 0.0f, -1.0f, 0.0f ),
        make_float4( c1, c1, c1, 0.0f ), make_float4( c1, c1, -c1, 0.0f ),
        make_float4( -c1, c1, c1, 0.0f ), make_float4( -c1, c1, -c1, 0.0f ),
        make_float4( c1, -c1, c1, 0.0f ), make_float4( c1, -c1, -c1, 0.0f ),
        make_float4( -c1, -c1, c1, 0.0f ), make_float4( -c1, -c1, -c1, 0.0f ) };

    for ( int i = 0; i < 14; ++i )
    {
        directions[i] = predefinedDirections[i];
    }
}

__host__ __device__ __inline__ void
buildInformationGlobalAABB( const lbvh::aabb<float>& aabb,
                            CenterGlobalSpaceBox& gBox )
{
    gBox.min = aabb.lower;
    gBox.max = aabb.upper;
    gBox.width = aabb.upper.x - aabb.lower.x;
    gBox.height = aabb.upper.y - aabb.lower.y;
    gBox.depth = aabb.upper.z - aabb.lower.z;
    gBox.radius = length( aabb.upper - aabb.lower ) * 0.5f;
    gBox.volume = gBox.width * gBox.height * gBox.depth;
    gBox.position = make_float4( ( aabb.upper.x + aabb.lower.x ) * 0.5f,
                                 ( aabb.upper.y + aabb.lower.y ) * 0.5f,
                                 ( aabb.upper.z + aabb.lower.z ) * 0.5f, 0.0f );
}

__global__ void initializeDirectionsKernel( float4* directions )
{
    initializeDirections( directions );
}

template <typename T, typename U>
__global__ void rayTracingKernelExplorationOptimized(
    lbvh::bvh_device<T, U> bvh_dev, Ray* rays, HitRay* d_HitRays, int numRays,
    float4* directions, const CenterGlobalSpaceBox* d_gBox )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if ( idx >= numRays )
        return;

    Ray ray = rays[idx];

    d_HitRays[idx].hitResults = -1;
    d_HitRays[idx].distanceResults = INFINITY; // distance
    d_HitRays[idx].intersectionPoint = make_float3( INFINITY, INFINITY, INFINITY );
    d_HitRays[idx].idResults = -1;

    constexpr float epsilon = 0.001f;
    constexpr float angleLimit = 0.6f;
    constexpr int maxIterations = 50;
    constexpr float angleToTriangleLim = 0.2f;

    float delta = -epsilon * 1.0f; // Initial delta for ray advancement
    bool foundCandidate = false;
    Triangle closestTriangle;
    int closestTriangleId = -1;
    int step = 1;
    float t;
    bool isViewInfo = true;
    isViewInfo = false;

    if ( step == 1 )
    {
        step = 2;
        const float epsilonC = 0.01f;
        for ( int i = 0; i < 14; ++i )
        {
            float4 currentPosition = ray.origin + directions[i] * epsilonC;
            const auto nearestTriangleIndex = lbvh::query_device(
                bvh_dev, lbvh::nearest( currentPosition ), distance_calculator() );

            if ( nearestTriangleIndex.first != 0xFFFFFFFF )
            {
                const Triangle& hitTriangle =
                    bvh_dev.objects[nearestTriangleIndex.first];
                if ( pointInTriangle2( currentPosition, hitTriangle, 0.01f ) )
                {

                    if ( isViewInfo )
                        printf( "in step1-level1\n" );

                    rayTriangleIntersect( rays[idx], hitTriangle, t );
                    {
                        if ( isViewInfo )
                            printf( "in step1-level2 it=%i\n", i );
                        float4 hit_point = ray.origin + ray.direction * t;
                        d_HitRays[idx].hitResults = nearestTriangleIndex.first;
                        d_HitRays[idx].distanceResults = t; // distance
                        d_HitRays[idx].intersectionPoint =
                            make_float3( hit_point.x, hit_point.y, hit_point.z );
                        d_HitRays[idx].idResults = hitTriangle.id;
                        step = -1;
                    }
                    break; // Exit the loop once we find a valid intersection
                }
            }
        }
    }

    //__syncthreads();

    if ( step == 2 )
    {
        if ( isViewInfo )
            printf( "in step2-level1\n" );

        float4 intersectionPointWithSphereScene;
        float distanceSphereScene = INFINITY;
        bool isSceneIntersection = raySphereIntersection(
            ray.origin, ray.direction, d_gBox->position, d_gBox->radius,
            intersectionPointWithSphereScene, distanceSphereScene );

        if ( isViewInfo )
        {
            printf( "isSceneIntersection=%i\n", isSceneIntersection );
            printf( "  dist=%f\n", distanceSphereScene );
            printf( "  position=<%f,%f,%f>\n", intersectionPointWithSphereScene.x,
                    intersectionPointWithSphereScene.y,
                    intersectionPointWithSphereScene.z );
        }

        if ( isSceneIntersection == 0 )
            return; // no objects are in this direction

        float distRayOriTogBoxCenter = length( ray.origin - d_gBox->position );
        float rfar = INFINITY;
        if ( distRayOriTogBoxCenter < 2.0f * d_gBox->radius )
        {
            rfar = 2.0f * d_gBox->radius;
        }
        else
        {
            rfar = 2.0f * distRayOriTogBoxCenter;
        }

        float dp = 0.0f;
        float angleToTriangle = INFINITY;
        float distanceToTriangle = 0.0f;
        float halfOpeningAngle = INFINITY;
        float inActivationAngle = INFINITY;
        float4 currentPosition;
        float4 currentPositionLast = currentPosition;

        for ( int iteration = 0; iteration < maxIterations; ++iteration )
        {
            currentPositionLast = currentPosition;
            currentPosition = ray.origin + ray.direction * delta;

            if ( currentPosition == currentPositionLast )
            {
                delta += distanceToTriangle * 0.15f + epsilon;
                currentPosition = ray.origin + ray.direction * delta;
            }

            auto nearestTriangleIndex = lbvh::query_device(
                bvh_dev, lbvh::nearest( currentPosition ), distance_calculator() );
            if ( nearestTriangleIndex.first != 0xFFFFFFFF )
            {
                const Triangle& hitTriangle =
                    bvh_dev.objects[nearestTriangleIndex.first];
                float4 triangleCenter =
                    ( hitTriangle.v1 + hitTriangle.v2 + hitTriangle.v3 ) / 3.0f;
                float4 directionToTriangle = triangleCenter - ray.origin;

                angleToTriangle = fabs( angleScalar( directionToTriangle, ray.direction ) );
                distanceToTriangle = length( directionToTriangle );
                halfOpeningAngle =
                    calculateHalfOpeningAngle( hitTriangle, currentPosition );
                inActivationAngle =
                    halfOpeningAngle - angleToTriangle - angleToTriangleLim;

                // if (halfOpeningAngle > 0.01f) { // Close object check
                if ( rayTriangleIntersect( ray, hitTriangle, t ) )
                {
                    if ( isViewInfo )
                        printf( "in step2-level2 it=%i\n", iteration );
                    float4 hit_point = ray.origin + ray.direction * t;
                    d_HitRays[idx].hitResults = nearestTriangleIndex.first;
                    d_HitRays[idx].distanceResults = length( ray.origin - hit_point );
                    d_HitRays[idx].intersectionPoint =
                        make_float3( hit_point.x, hit_point.y, hit_point.z );
                    d_HitRays[idx].idResults = hitTriangle.id;
                    return;
                }
                //}
            }

            if ( angleToTriangle > angleToTriangleLim )
            {
                // delta += epsilon * exp(iteration*0.25f);
                delta += distanceToTriangle * 0.75f + epsilon;
            }

            /*
                  if (inActivationAngle > 0.0f) {
                    // delta += epsilon * exp(iteration*0.25f);
                    delta += distanceToTriangle * 0.75f + epsilon;
                  }
            */

        } // END FOR
    }
}

__device__ __inline__ bool
checkOverlap( const float4& observer, const float4& obj1, const float& radius1,
              const float4& obj2, const float& radius2 )
{
    float4 diff1 = make_float4( obj1.x - observer.x, obj1.y - observer.y,
                                obj1.z - observer.z, 0.0f );
    float4 diff2 = make_float4( obj2.x - observer.x, obj2.y - observer.y,
                                obj2.z - observer.z, 0.0f );
    float dist1 =
        sqrtf( diff1.x * diff1.x + diff1.y * diff1.y + diff1.z * diff1.z );
    float dist2 =
        sqrtf( diff2.x * diff2.x + diff2.y * diff2.y + diff2.z * diff2.z );
    float angle1 = atan2f( radius1, dist1 );
    float angle2 = atan2f( radius2, dist2 );
    float dotProduct = diff1.x * diff2.x + diff1.y * diff2.y + diff1.z * diff2.z;
    float angleBetween = acosf( dotProduct / ( dist1 * dist2 ) );
    return angleBetween < ( angle1 + angle2 );
}

template <typename T, typename U>
__global__ void rayTracingKernelExplorationOptimized2(
    lbvh::bvh_device<T, U> bvh_dev, Ray* rays, HitRay* d_HitRays, int numRays,
    float4* directions, const CenterGlobalSpaceBox* d_gBox )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if ( idx >= numRays )
        return;

    Ray ray = rays[idx];
    HitRay& hitRay = d_HitRays[idx];

    // Initialize hit results
    d_HitRays[idx].hitResults = -1;
    d_HitRays[idx].distanceResults = INFINITY; // distance
    d_HitRays[idx].intersectionPoint = make_float3( INFINITY, INFINITY, INFINITY );
    d_HitRays[idx].idResults = -1;

    constexpr float epsilon = 0.001f;
    constexpr int maxIterations = 50;
    constexpr float epsilonC = 0.01f;
    constexpr float angleToTriangleLim = 0.2f;

    bool isViewInfo = true;
    isViewInfo = false;

    float t;

    // Step 1: Quick intersection check
    for ( int i = 0; i < 14; ++i )
    {
        float4 currentPosition = ray.origin + directions[i] * epsilonC;
        auto nearestTriangleIndex = lbvh::query_device(
            bvh_dev, lbvh::nearest( currentPosition ), distance_calculator() );

        if ( nearestTriangleIndex.first != 0xFFFFFFFF )
        {
            const Triangle& hitTriangle = bvh_dev.objects[nearestTriangleIndex.first];
            if ( pointInTriangle2( currentPosition, hitTriangle, 0.01f ) )
            {
                if ( isViewInfo )
                    printf( "in step1-level1\n" );

                rayTriangleIntersect( ray, hitTriangle, t );
                {
                    if ( isViewInfo )
                        printf( "in step1-level2 it=%i\n", i );
                    updateHitResults( hitRay, nearestTriangleIndex.first, t, ray,
                                      hitTriangle );
                    return;
                }
            }
        }
    }

    // Step 2: Iterative ray marching
    if ( isViewInfo )
        printf( "in step2-level1\n" );
    float delta = epsilon;
    float4 currentPosition, lastPosition = ray.origin;

    float4 intersectionPointWithSphereScene;
    float distanceSphereScene = INFINITY;
    bool isSceneIntersection = raySphereIntersection(
        ray.origin, ray.direction, d_gBox->position, d_gBox->radius,
        intersectionPointWithSphereScene, distanceSphereScene );

    if ( isViewInfo )
    {
        printf( "isSceneIntersection=%i\n", isSceneIntersection );
        printf( "  dist=%f\n", distanceSphereScene );
        printf( "  position=<%f,%f,%f>\n", intersectionPointWithSphereScene.x,
                intersectionPointWithSphereScene.y,
                intersectionPointWithSphereScene.z );
    }

    if ( isSceneIntersection == 0 )
        return; // no objects are in this direction

    float distRayOriTogBoxCenter = length( ray.origin - d_gBox->position );
    float rfar = INFINITY;
    if ( distRayOriTogBoxCenter < 2.0f * d_gBox->radius )
    {
        rfar = 2.0f * d_gBox->radius;
    }
    else
    {
        rfar = 2.0f * distRayOriTogBoxCenter;
    }

    float angleToTriangle;
    float distanceToTriangle;
    float halfOpeningAngle;
    float inActivationAngle;

    for ( int iteration = 0; iteration < maxIterations; ++iteration )
    {
        currentPosition = ray.origin + ray.direction * delta;
        if ( currentPosition == lastPosition )
        {
            delta += epsilon;
            continue;
        }
        lastPosition = currentPosition;
        // if (delta>rfar) printf("OUT it=%i\n",iteration);
        if ( delta > rfar )
            return; // we find ourselves outside the encompassing sphere of the scene.
                    // There is nothing more to look for after that.

        auto nearestTriangleIndex = lbvh::query_device(
            bvh_dev, lbvh::nearest( currentPosition ), distance_calculator() );
        if ( nearestTriangleIndex.first != 0xFFFFFFFF )
        {
            const Triangle& hitTriangle = bvh_dev.objects[nearestTriangleIndex.first];

            // if (isViewInfo)
            //     printf("in step2 num triangle=%i\n",nearestTriangleIndex.first);

            float4 positionToTriangle =
                ( hitTriangle.v1 + hitTriangle.v2 + hitTriangle.v3 ) / 3.0f;
            float4 directionToTriangle = positionToTriangle - currentPosition;
            distanceToTriangle = length( directionToTriangle );

            angleToTriangle = fabs( angleScalar( directionToTriangle, ray.direction ) );

            halfOpeningAngle =
                calculateHalfOpeningAngle( hitTriangle, currentPosition );

            inActivationAngle =
                halfOpeningAngle - ( angleToTriangle - angleToTriangleLim );

            // printf("halfOpeningAngle= %f angleToTriangle=%f inActivationAngle=%f
            // num triangle=%i\n", halfOpeningAngle,
            // angleToTriangle,inActivationAngle,nearestTriangleIndex.first);
            float4 positionInDir = ray.direction * ( distanceToTriangle - 0.001f );

            if ( 1 == 0 )
            {
                printf( "\n" );
                printf( "currentPosition     = <%f %f %f>\n", currentPosition.x,
                        currentPosition.y, currentPosition.z );
                printf( "positionToTriangle  = <%f %f %f>\n", positionToTriangle.x,
                        positionToTriangle.y, positionToTriangle.z );
                printf( "positionInDir       = <%f %f %f>\n", positionInDir.x,
                        positionInDir.y, positionInDir.z );
                printf( "directionToTriangle = <%f %f %f>\n", directionToTriangle.x,
                        directionToTriangle.y, directionToTriangle.z );

                if ( checkOverlap( { 0.0f, 0.0f, 0.0f, 0.0f }, directionToTriangle,
                                   halfOpeningAngle, positionInDir, angleToTriangleLim ) )
                {
                    printf( "Overlap TRUE num triangle=%i\n", nearestTriangleIndex.first );
                }
                else
                {
                    printf( "Overlap FALSE num triangle=%i\n", nearestTriangleIndex.first );
                }
            }

            /*
                  if (checkOverlap({0.0f, 0.0f, 0.0f, 0.0f},
                                   {10.0f, 0.0f, 1.0f, 0.0f},
                                   2,
                                   {9.0f, 0.0f, 1.0f, 0.0f},
                                   10))
                  {
                    printf("Overlap TRUE num triangle 2\n");
                  } else {
                    printf("Overlap FALSE num triangle 2\n");
                  }
            */

            if ( rayTriangleIntersect( ray, hitTriangle, t ) )
            {
                if ( isViewInfo )
                    printf( "in step2-level2 it=%i\n", iteration );
                updateHitResults( hitRay, nearestTriangleIndex.first, t, ray,
                                  hitTriangle );
                return;
            }

            delta += ( angleToTriangle > angleToTriangleLim )
                         ? ( distanceToTriangle * 0.75f + epsilon )
                         : epsilon;
        }
        else
        {
            delta += epsilon;
        }
    }
}

} // END namespace bvhLinearExplorer

#endif

using namespace Feel;

namespace Feel
{

template <int RealDim>
class BVHRay
{
  public:
    using vec_t = eigen_vector_type<RealDim>;
    BVHRay( vec_t const& orig, vec_t const& dir,
            double dmin = 0, double dmax = std::numeric_limits<double>::max() )
        : M_origin( orig ),
          M_dir( dir ),
          M_distanceMin( dmin ),
          M_distanceMax( dmax )
    {
    }
    BVHRay()
        : BVHRay( vec_t::Zero(), vec_t::Zero() ) {}

    BVHRay( BVHRay const& ) = default;
    BVHRay( BVHRay&& ) = default;
    BVHRay& operator=( BVHRay&& ) = default;
    BVHRay& operator=( BVHRay const& ) = default;

    vec_t const& origin() const noexcept { return M_origin; }
    vec_t const& dir() const noexcept { return M_dir; }
    double distanceMin() const { return M_distanceMin; }
    double distanceMax() const { return M_distanceMax; }
    size_t id;

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
        ar& M_origin;
        ar& M_dir;
        ar& M_distanceMin;
        ar& M_distanceMax;
    }

  private:
    vec_t M_origin, M_dir; // ray origin and dir
    double M_distanceMin, M_distanceMax;
};

template <int RealDim>
class BVHRaysDistributed
{
  public:
    using ray_type = BVHRay<RealDim>;
    BVHRaysDistributed() = default;
    BVHRaysDistributed( BVHRaysDistributed&& ) = default;
    std::vector<ray_type> const& rays() const { return M_rays; }
    //! return number of local ray
    std::size_t numberOfLocalRay() const { return M_rays.size(); }

    template <typename T>
    void push_back( T&& ray ) { M_rays.push_back( std::forward<T>( ray ) ); }

  private:
    std::vector<ray_type> M_rays;
};

struct BVHEnum
{
    enum class Quality
    {
        Low,
        Medium,
        High
    };
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
    using vector_realdim_type = Eigen::Matrix<double, nRealDim, 1>;
    using ray_type = BVHRay<nRealDim>;

    //! @brief Information on the primitive (mesh entity, bounding box, centroid)
    struct BVHPrimitiveInfo
    {
        BVHPrimitiveInfo( mesh_entity_type const& meshEntity )
            : M_meshEntity( meshEntity )
        {
            auto verticesUblas = meshEntity.vertices();
            auto G = em_cmatrix_col_type<double>( verticesUblas.data().begin(), nRealDim, mesh_entity_type::numVertices );
            M_bound_min = G.rowwise().minCoeff();
            M_bound_max = G.rowwise().maxCoeff();
            M_bound_min.array() -= 2 * FLT_MIN;
            M_bound_max.array() += 2 * FLT_MIN;
            // M_centroid = ( M_bound_min + M_bound_max ) * 0.5;
            auto bary = meshEntity.barycenter();
            M_centroid = Eigen::Map<Eigen::Matrix<double, nRealDim, 1>>( bary.data().begin() );
        }
        BVHPrimitiveInfo( BVHPrimitiveInfo&& ) = default;
        BVHPrimitiveInfo( BVHPrimitiveInfo const& ) = default;
        BVHPrimitiveInfo& operator=( BVHPrimitiveInfo&& ) = default;
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
            : M_processId( processId ),
              M_primitiveId( primitiveId ),
              M_distance( dist )
        {
            M_id = -1;
            // M_id = 0;
        }
        BVHRayIntersectionResult( BVHRayIntersectionResult&& ) = default;
        BVHRayIntersectionResult( BVHRayIntersectionResult const& ) = default;
        BVHRayIntersectionResult& operator=( BVHRayIntersectionResult&& ) = default;
        BVHRayIntersectionResult& operator=( BVHRayIntersectionResult const& ) = default;

        rank_type processId() const noexcept { return M_processId; }
        index_type primitiveId() const noexcept { return M_primitiveId; }
        double distance() const noexcept { return M_distance; }
        size_t get_id() const noexcept { return M_id; }

        template <typename T>
        void setCoordinates( T&& coord ) { M_coordinates = std::forward<T>( coord ); }
        bool hasCoordinates() const noexcept { return M_coordinates.has_value(); }
        vector_realdim_type const& coordinates() const noexcept { return *M_coordinates; }
        size_t M_id;

      private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize( Archive& ar, const unsigned int version )
        {
            ar& M_processId;
            ar& M_primitiveId;
            ar& M_distance;

            if constexpr ( Archive::is_saving::value )
            {
                bool hasCoordinates = this->hasCoordinates();
                ar& boost::serialization::make_nvp( "hasCoordinates", hasCoordinates );
                if ( hasCoordinates )
                    ar& boost::serialization::make_nvp( "coordinates", this->coordinates() );
            }
            else if constexpr ( Archive::is_loading::value )
            {
                bool hasCoordinates = false;
                ar& boost::serialization::make_nvp( "hasCoordinates", hasCoordinates );
                if ( hasCoordinates )
                {
                    vector_realdim_type coord;
                    ar& boost::serialization::make_nvp( "coordinates", coord );
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

    enum class IntersectContext
    {
        anyHint = 0,
        closest,
        all
    };

    BVH( BVHEnum::Quality quality = BVHEnum::Quality::High, worldcomm_ptr_t worldComm = Environment::worldCommPtr() )
        : CommObject( worldComm ),
          M_quality( quality )
    {
    }
    BVH( BVH&& ) = default;
    BVH( BVH const& ) = default;
    virtual ~BVH() {}

    //! return all primitive info
    std::vector<BVHPrimitiveInfo> const& primitiveInfo() const noexcept { return M_primitiveInfo; }

    //! return primitive info at index i
    BVHPrimitiveInfo const& primitiveInfo( index_type i ) const { return M_primitiveInfo.at( i ); }

    //! compute intersection(s) with a ray from the BVH built and return a vector of intersection result
    template <typename... Ts>
    auto intersect( Ts&&... v )
    {

        auto args = NA::make_arguments( std::forward<Ts>( v )... );
        auto&& ray = args.get( _ray );
        bool useRobustTraversal = args.get_else( _robust, true );
        IntersectContext ctx = args.get_else( _context, IntersectContext::closest );
        bool parallel = args.get_else( _parallel, this->worldComm().size() > 1 );

        if ( this->isGPUHip() ) parallel = false;

        bool closestOnly = ctx == IntersectContext::closest;
        using napp_ray_type = std::decay_t<decltype( ray )>;
        if constexpr ( std::is_same_v<BVHRaysDistributed<nRealDim>, napp_ray_type> ) // case rays distributed on process
        {

#if 1
            // Old Version
            int worldSize = this->worldComm().size(); // Nombre de processus
            int worldRank = this->worldComm().rank(); // Rang du processus actuel
            printf( "IIIIIIIIIIIIIIIIIII worldSize=%i worldRank=%i\n", worldSize, worldRank );

            tic();
            // WARNING: this algo is not good (all_gather of rays then all run bvh), just a quick version for test
            tic();
            auto const& localRays = ray.rays();
            std::vector<int> resLocalSize( this->worldComm().size() );
            mpi::all_gather( this->worldComm(), (int)localRays.size(), resLocalSize );
            auto timeDuration_intersect_block1_1 = toc( "timeDuration_intersect_block1_1 : all_gather" );

            tic();
            std::vector<ray_type> raysGathered;
            if ( this->worldComm().isMasterRank() )
            {
                int gatherRaySize = std::accumulate( resLocalSize.begin(), resLocalSize.end(), 0 );
                raysGathered.resize( gatherRaySize );
            }
            auto timeDuration_intersect_block1_2 = toc( "timeDuration_intersect_block1_2 : accumulation" );

            tic();
            mpi::gatherv( this->worldComm(), localRays, raysGathered.data(), resLocalSize, this->worldComm().masterRank() );
            auto timeDuration_intersect_block1_3 = toc( "timeDuration_intersect_block1_3 : gatherv" );

            tic();
            mpi::broadcast( this->worldComm(), raysGathered, this->worldComm().masterRank() );
            auto timeDuration_intersect_block1_4 = toc( "timeDuration_intersect_block1_4 : broacast" );

            tic();
            auto intersectGlobal = this->intersect( _ray = raysGathered, _robust = useRobustTraversal, _context = ctx, _parallel = true );
            // auto intersectGlobal = this->intersect( _ray = raysGathered, _robust = useRobustTraversal, _context = ctx);
            auto timeDuration_intersect_block1_5 = toc( "timeDuration_intersect_block1_5 : intersect" );

            tic();
            std::vector<std::vector<rayintersection_result_type>> res;
            res.resize( ray.numberOfLocalRay() );
            std::size_t startRayIndexInThisProcess = 0;
            for ( int p = 0; p < this->worldComm().rank(); ++p )
                startRayIndexInThisProcess += resLocalSize[p];
            std::copy_n( intersectGlobal.cbegin() + startRayIndexInThisProcess, localRays.size(), res.begin() );
            auto timeDuration_intersect_block1_6 = toc( "timeDuration_intersect_block1_6 : add all intersection" );

            auto timeDuration_intersect_block1 = toc( "timeDuration_intersect_block1 : all block in function" );
            return res;
#endif

#if 0
       
            // Old Version
            int worldSize = this->worldComm().size(); // Nombre de processus
            int worldRank = this->worldComm().rank(); // Rang du processus actuel
            printf( "IIIIIIIIIIIIIIIIIII worldSize=%i worldRank=%i\n", worldSize, worldRank );

            tic();
            // WARNING: this algo is not good (all_gather of rays then all run bvh), just a quick version for test
            tic();
            auto const& localRays = ray.rays();
            std::vector<int> resLocalSize( this->worldComm().size() );
            mpi::all_gather( this->worldComm(), (int)localRays.size(), resLocalSize );
            auto timeDuration_intersect_block1_1 = toc( "timeDuration_intersect_block1_1 : all_gather" );

            tic();
            std::vector<ray_type> raysGathered;
            if ( this->worldComm().isMasterRank() )
            {
                int gatherRaySize = std::accumulate( resLocalSize.begin(), resLocalSize.end(), 0 );
                raysGathered.resize( gatherRaySize );
            }
            auto timeDuration_intersect_block1_2 = toc( "timeDuration_intersect_block1_2 : accumulation" );

            tic();
            mpi::gatherv( this->worldComm(), localRays, raysGathered.data(), resLocalSize, this->worldComm().masterRank() );
            auto timeDuration_intersect_block1_3 = toc( "timeDuration_intersect_block1_3 : gatherv" );

            tic();
            mpi::broadcast( this->worldComm(), raysGathered, this->worldComm().masterRank() );
            auto timeDuration_intersect_block1_4 = toc( "timeDuration_intersect_block1_4 : broacast" );

            tic();
            auto intersectGlobal = this->intersect( _ray = raysGathered, _robust = useRobustTraversal, _context = ctx, _parallel = true );
            // auto intersectGlobal = this->intersect( _ray = raysGathered, _robust = useRobustTraversal, _context = ctx);
            auto timeDuration_intersect_block1_5 = toc( "timeDuration_intersect_block1_5 : intersect" );


            this->worldComm().barrier();

            tic();
            std::vector<std::vector<rayintersection_result_type>> res;
            res.resize( ray.numberOfLocalRay() );
            std::size_t startRayIndexInThisProcess = 0;
            for ( int p = 0; p < this->worldComm().rank(); ++p )
                startRayIndexInThisProcess += resLocalSize[p];
            std::copy_n( intersectGlobal.cbegin() + startRayIndexInThisProcess, localRays.size(), res.begin() );
            auto timeDuration_intersect_block1_6 = toc( "timeDuration_intersect_block1_6 : add all intersection" );

            auto timeDuration_intersect_block1 = toc( "timeDuration_intersect_block1 : all block in function" );
            return res;
#endif
        }
        else if constexpr ( is_iterable_v<std::decay_t<decltype( ray )>> ) // case rays container are identical all on process (TODO: internal case)
        {
            tic();
            std::vector<std::vector<rayintersection_result_type>> resSeq;
            resSeq.reserve( ray.size() );

            auto timeDuration_resSeqreserve = toc( "timeDuration_resSeqreserve" );

            if ( !( this->isGPUHip() ) )
            {
                for ( auto const& currentRay : ray )
                {
                    auto currentResSeq = this->intersectSequential( currentRay, useRobustTraversal );
                    if ( closestOnly && currentResSeq.size() > 1 )
                        currentResSeq.resize( 1 );
                    resSeq.push_back( std::move( currentResSeq ) );
                }
            }
            else
            {
                tic();
                resSeq = this->intersectAllRaysWithGPU( ray );
                auto timeDurationFunctionIntersectAllRaysWithGPU = toc( "timeDurationFunctionIntersectAllRaysWithGPU" );
            }

            if ( !parallel )
                return resSeq;

            tic();
            mpi::all_reduce( this->worldComm(), mpi::inplace( resSeq ), []( auto const& x, auto const& y ) -> std::vector<std::vector<rayintersection_result_type>>
                             {
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
                        return ret; } );
            auto timeDuration_intersect_block2 = toc( "timeDuration_intersect_block2" );

            return resSeq;
        }
        else // only one ray (all process should have the same ray if parallel=true)
        {
            // Before addendum
            /*
            auto resSeq = this->intersectSequential( ray,useRobustTraversal );
            if ( closestOnly && resSeq.size() > 1 )
                resSeq.resize(1);
            */

            // TODO CTRL This part if OK

            tic();
            std::vector<std::vector<rayintersection_result_type>> resSeq;
            if ( !( this->isGPUHip() ) )
            {
                resSeq = this->intersectSequential( ray, useRobustTraversal );
            }
            else
            {
                resSeq = this->intersectAllRaysWithGPU( ray );
            }

            if ( closestOnly && resSeq.size() > 1 )
                resSeq.resize( 1 );

            auto timeDuration_intersect_block3 = toc( "timeDuration_intersect_block3" );

            if ( !parallel )
                return resSeq;

#if 1
            mpi::all_reduce( this->worldComm(), mpi::inplace( resSeq ), []( auto const& a, auto const& b ) -> std::vector<rayintersection_result_type>
                             {
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
                        } } );
            return resSeq;
#else
            tic();
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
                std::sort( resPar.begin(), resPar.end(), []( auto const& res0, auto const& res1 )
                           { return res0.distance() < res1.distance(); } );
                if ( closestOnly && resPar.size() > 1 )
                    resPar.resize( 1 );
            }
            mpi::broadcast( this->worldComm(), resPar, this->worldComm().masterRank() );

            auto timeDuration_intersect_block4 = toc( "timeDuration_intersect_block4" );

            return resPar;
#endif
        }
    }

  protected:
    virtual std::vector<rayintersection_result_type> intersectSequential( ray_type const& rayon, bool useRobustTraversal = true ) = 0;
    virtual std::vector<std::vector<rayintersection_result_type>> intersectAllRaysWithGPU( std::vector<ray_type> const& rayons ) = 0;
    virtual bool isGPUHip() = 0;

    template <typename RangeType>
    void
    updateForUse( RangeType const& range )
    {
        // From the mesh, build the bounding box info for each element and store it in
        // the structure BVHPrimitiveInfo
        M_primitiveInfo.clear();
        M_primitiveInfo.reserve( nelements( range ) );
        for ( auto const& eltWrap : range )
        {
            auto const& e = unwrap_ref( eltWrap );
            M_primitiveInfo.push_back( BVHPrimitiveInfo{ e } );
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

    BVH_ThirdParty( BVHEnum::Quality quality, worldcomm_ptr_t worldComm )
        : super_type( quality, worldComm ) {}
    BVH_ThirdParty( BVH_ThirdParty&& ) = default;

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
                backend_vector_realdim_type::generate( [&primInfo]( std::size_t i )
                                                       { return primInfo.boundMin()[i]; } ),
                backend_vector_realdim_type::generate( [&primInfo]( std::size_t i )
                                                       { return primInfo.boundMax()[i]; } ) } );

            auto const& centroid = primInfo.centroid();
            centers.push_back( backend_vector_realdim_type::generate( [&centroid]( std::size_t i )
                                                                      { return centroid[i]; } ) );
        }

        typename bvh::v2::DefaultBuilder<node_type>::Config config;
        switch ( this->M_quality )
        {
        default:
        case BVHEnum::Quality::High:
            config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::High;
            break;
        case BVHEnum::Quality::Medium:
            config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::Medium;
            break;
        case BVHEnum::Quality::Low:
            config.quality = bvh::v2::DefaultBuilder<node_type>::Quality::Low;
            break;
        }
        M_bvh = std::make_unique<backend_bvh_type>( bvh::v2::DefaultBuilder<node_type>::build( /*thread_pool,*/ bboxes, centers, config ) );

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
                auto const& pt0 = meshEntity.point( 0 );
                auto const& pt1 = meshEntity.point( 1 );
                auto const& pt2 = meshEntity.point( 2 );
                M_precomputeTriangle[i] = backend_precompute_triangle_type{
                    backend_vector_realdim_type::generate( [&pt0]( std::size_t i )
                                                           { return pt0[i]; } ),
                    backend_vector_realdim_type::generate( [&pt1]( std::size_t i )
                                                           { return pt1[i]; } ),
                    backend_vector_realdim_type::generate( [&pt2]( std::size_t i )
                                                           { return pt2[i]; } ) };
            }
        }
    }

  private:
    // Value indicating whether we are in GPU mode
    bool isGPUHip() override
    {
        return ( false );
    }

    std::vector<std::vector<rayintersection_result_type>> intersectAllRaysWithGPU( std::vector<ray_type> const& rayons ) override
    {
        CHECK( false ) << "no implementation";
        std::vector<std::vector<rayintersection_result_type>> res;
        return ( res );
    }

    std::vector<rayintersection_result_type> intersectSequential( ray_type const& ray, bool useRobustTraversal = true ) override
    {
        int idRay = ray.id;
        auto rayBackend = bvh::v2::Ray<value_type, nRealDim>{
            backend_vector_realdim_type::generate( [&ray]( std::size_t i )
                                                   { return ray.origin()[i]; } ),
            backend_vector_realdim_type::generate( [&ray]( std::size_t i )
                                                   { return ray.dir()[i]; } ),
            ray.distanceMin(), ray.distanceMax() };

        if ( useRobustTraversal )
            return this->intersectImpl<true>( rayBackend, idRay );
        else
            return this->intersectImpl<false>( rayBackend, idRay );
    };

    template <bool UseRobustTraversal>
    std::vector<rayintersection_result_type> intersectImpl( bvh::v2::Ray<value_type, nRealDim>& rayBackend, int idRay )
    {
        static constexpr size_t stack_size = 64;
        static constexpr bool should_permute = true;
        static constexpr bool isAnyHit = false;
        // Traverse the BVH and get the u, v coordinates of the closest intersection.
        bvh::v2::SmallStack<typename backend_bvh_type::Index, stack_size> stack;
        std::vector<rayintersection_result_type> res;
        M_bvh->template intersect<isAnyHit, UseRobustTraversal>( rayBackend, M_bvh->get_root().index, stack,
                                                                 [this, &res, &rayBackend, &idRay]( std::size_t begin, std::size_t end )
                                                                 {
                                                                     std::size_t previousResultSize = res.size();
                                                                     for ( std::size_t i = begin; i < end; ++i )
                                                                     {
                                                                         std::size_t j = should_permute ? i : M_bvh->prim_ids[i];
                                                                         if constexpr ( nRealDim == 2 )
                                                                         {
                                                                             CHECK( false ) << "TODO";
                                                                         }
                                                                         else if constexpr ( nRealDim == 3 )
                                                                         {
                                                                             if ( auto hit = M_precomputeTriangle[j].intersect( rayBackend ) )
                                                                             {
                                                                                 // std::tie(u, v) = *hit;
                                                                                 res.push_back( rayintersection_result_type( this->worldComm().rank(), M_bvh->prim_ids[i], rayBackend.tmax ) );
                                                                                 res.back().setCoordinates( this->barycentricToCartesianCoordinates( M_precomputeTriangle[j].convert_to_tri(), hit->first, hit->second ) );
                                                                                 res.back().M_id = idRay; // ray.id()[i];
                                                                                 // std::cout<<"in intersectImpl rank="<<this->worldComm().rank()<<" id="<<idRay<<" "<<M_bvh->prim_ids[i]<<"\n";
                                                                                 if constexpr ( isAnyHit )
                                                                                     return true;
                                                                             }
                                                                         }
                                                                     }
                                                                     return res.size() > previousResultSize;
                                                                 } );
        //! sort all intersection from the distance (closer to far)
        std::sort( res.begin(), res.end(), []( auto const& res0, auto const& res1 )
                   { return res0.distance() < res1.distance(); } );

        return res;
    }

    template <typename TriType>
    vector_realdim_type barycentricToCartesianCoordinates( TriType const& tri, double u, double v ) const
    {
        auto const& pt0 = tri.p1;
        auto const& pt1 = tri.p2;
        auto const& pt2 = tri.p0;
        return vector_realdim_type{ {
            u * pt0[0] + v * pt1[0] + ( 1 - u - v ) * pt2[0],
            u * pt0[1] + v * pt1[1] + ( 1 - u - v ) * pt2[1],
            u * pt0[2] + v * pt1[2] + ( 1 - u - v ) * pt2[2],
        } };
    }

  private:
    std::unique_ptr<backend_bvh_type> M_bvh;
    std::vector<backend_precompute_triangle_type> M_precomputeTriangle;
};

/***************************************************************************************************************************************************/
// Added additional tools to the BVH_ThirdParty class (Saving the BVH, displaying information to control values, etc.)
// (... out )

/***************************************************************************************************************************************************/
// GPU part for BVH and Ray Tracing

//! @brief implementation of BVH tool with an external third party in GPU

#if defined( FEELPP_HAS_HIP )

template <typename MeshEntityType>
class BVH_HIP_Party : public BVH<MeshEntityType>
{
    using super_type = BVH<MeshEntityType>;
    using mesh_entity_type = typename super_type::mesh_entity_type;
    using vector_realdim_type = typename super_type::vector_realdim_type;
    static constexpr uint16_type nRealDim = super_type::nRealDim;

  public:
    using ray_type = typename super_type::ray_type;
    using rayintersection_result_type = typename super_type::rayintersection_result_type;

    lbvh::bvh<float, bvhLinearExplorer::Triangle, bvhLinearExplorer::aabb_getter> bvhl;
    lbvh::bvh_device<float, bvhLinearExplorer::Triangle> bvhl_dev;
    bvhLinearExplorer::Triangle* deviceLBVHTriangles;
    bvhLinearExplorer::CenterGlobalSpaceBox* d_gBox;

    bvhHip::BVHNode* devicebvhHipNodes;
    bvhHip::Triangle* deviceHipTriangles;
    size_t nbTriangles;
    size_t nbNodes;

    int numDevice;
    int numVersion;
    int modeGPU; // 1 - HIP 4 - LBVH Explorer
    bool isUnifiedMemory;
    int rank;
    bool isView;
    bool isViewDataRT;

    BVH_HIP_Party( BVHEnum::Quality quality, worldcomm_ptr_t worldComm )
        : super_type( quality, worldComm )
    {
        int num_gpus = 0;
        hipGetDeviceCount( &num_gpus );
        // numDevice = 0;
        // numDevice = 3;
        rank = worldComm->rank();
        numDevice = rank % num_gpus;
        printf( "[INFO GPU]: constructor rank %i numDevice=%i  num_gpus=%i\n", rank, numDevice, num_gpus );
        numVersion = 2;
        modeGPU = 1;
        isUnifiedMemory = true; // isUnifiedMemory = false;
        isView = false;         // isView = true;
        isViewDataRT = false;   // isViewDataRT = true;

        numVersion = 2;
        // isView = true; isViewDataRT = true; isUnifiedMemory = false; numVersion = 0;
        // isView = true; isViewDataRT = true; isUnifiedMemory = true; numVersion = 1;
        // isView = true; isViewDataRT = true; isUnifiedMemory = true; numVersion = 2;
        // isView = true; isViewDataRT = true;

        //modeGPU = 4;
        //numVersion = 4;
    }

    ~BVH_HIP_Party()
    {
        //  isView=true;
        //  Memory cleaning
        HIP_ASSERT( hipFree( devicebvhHipNodes ) );
        HIP_ASSERT( hipFree( deviceHipTriangles ) );
        if ( isView ) std::cout << "[INFO GPU]: [GPU MEMORY CLEANING DONE]"
                                << "\n";
    }

    BVH_HIP_Party( BVH_HIP_Party&& ) = default;

    template <typename RangeType>
    void
    updateForUse( RangeType const& range )
    {
        // isView = true;
        // up primitiveinfos
        super_type::updateForUse( range );
        // init bvh backend
        if ( isView ) std::cout << "[INFO GPU]: Size primitiveInfo=" << this->M_primitiveInfo.size() << " [" << rank << ":" << numDevice << "]"
                                << "\n";
        //...

        // hip device used
        hipSetDevice( numDevice );
        int numDeviceActivated;
        hipGetDevice( &numDeviceActivated );
        if ( isView ) std::cout << "[INFO GPU]: Num Device Activated=" << numDeviceActivated << " [" << rank << ":" << numDevice << "]"
                                << "\n";

        // Definition of operating modes
        bool isModeBox = true;
        isModeBox = false;
        bool isModeDirectInDevice = false;
        // isModeDirectInDevice = true; // Todo CTRL in infinity and size max of HIP GPU

        std::chrono::steady_clock::time_point t_begin_bvh_gpu, t_end_bvh_gpu;
        long int t_laps_bvh_gpu = 0;

        if ( modeGPU == 1 )
        {
            if ( !isUnifiedMemory )
            {
                tic();
                std::vector<bvhHip::Triangle> HostHipTriangles;

                // Load Mesh  in host
                for ( size_t k = 0; k < this->M_primitiveInfo.size(); ++k )
                {
                    size_t id = this->M_primitiveInfo[k].meshEntity().id();
                    auto const& primInfo = this->M_primitiveInfo[k];
                    auto const& meshEntity = primInfo.meshEntity();
                    auto const& pt0 = meshEntity.point( 0 );
                    auto const& pt1 = meshEntity.point( 1 );
                    auto const& pt2 = meshEntity.point( 2 );
                    bvhHip::Triangle ltri;
                    ltri.v0 = bvhHip::Vec3( pt0[0], pt0[1], pt0[2] );
                    ltri.v1 = bvhHip::Vec3( pt1[0], pt1[1], pt1[2] );
                    ltri.v2 = bvhHip::Vec3( pt2[0], pt2[1], pt2[2] );
                    ltri.id = id;
                    HostHipTriangles.push_back( ltri );
                }

                size_t numTriangles = HostHipTriangles.size();

                HIP_ASSERT( hipMalloc( &deviceHipTriangles, numTriangles * sizeof( bvhHip::Triangle ) ) );
                HIP_ASSERT( hipMemcpy( deviceHipTriangles, HostHipTriangles.data(), numTriangles * sizeof( bvhHip::Triangle ), hipMemcpyHostToDevice ) );

                HIP_ASSERT( hipMalloc( &devicebvhHipNodes, ( 2 * numTriangles - 1 ) * sizeof( bvhHip::BVHNode ) ) );

                auto timeDataTransfertDuration = toc( "timeDataTransfertTriangleDuration" );
                if ( isView ) std::cout << "[INFO GPU]: TransfertDataTriangle=" << numTriangles << " [" << rank << ":" << numDevice << "]"
                                        << "\n";

                t_begin_bvh_gpu = std::chrono::steady_clock::now();
                if ( numVersion == 0 ) bvhHip::buildBVH_GPU_Version2( deviceHipTriangles, devicebvhHipNodes, numTriangles );
                if ( numVersion == 1 ) bvhHip::buildBVH_GPU_Version3( deviceHipTriangles, devicebvhHipNodes, numTriangles );
                if ( numVersion == 4 ) bvhHip::buildBVH_GPU_Version4( deviceHipTriangles, devicebvhHipNodes, numTriangles );
                if ( numVersion == 2 ) bvhHip::buildBVH_GPU_Parallel_Best_Axis( deviceHipTriangles, devicebvhHipNodes, numTriangles );
                t_end_bvh_gpu = std::chrono::steady_clock::now();

                if ( isView ) std::cout << "[INFO GPU]: BVH GPU FINISHED"
                                        << " [" << rank << ":" << numDevice << "]"
                                        << "\n";
                t_laps_bvh_gpu = std::chrono::duration_cast<std::chrono::microseconds>( t_end_bvh_gpu - t_begin_bvh_gpu ).count();
                if ( isView ) std::cout << "[INFO GPU]: Elapsed microseconds inside BVH  : " << t_laps_bvh_gpu << " us"
                                        << " [" << rank << ":" << numDevice << "]"
                                        << "\n";

                nbTriangles = numTriangles;
                nbNodes = 2 * nbTriangles - 1;
                ;
            }
            else
            {
                tic();
                size_t numTriangles = this->M_primitiveInfo.size();

#if 1
                HIP_ASSERT( hipMallocManaged( &deviceHipTriangles, numTriangles * sizeof( bvhHip::Triangle ) ) );
                HIP_ASSERT( hipMallocManaged( &devicebvhHipNodes, ( 2 * numTriangles - 1 ) * sizeof( bvhHip::BVHNode ) ) );
                // Load Mesh  in host-device

                for ( size_t k = 0; k < this->M_primitiveInfo.size(); ++k )
                {
                    size_t id = this->M_primitiveInfo[k].meshEntity().id();
                    auto const& primInfo = this->M_primitiveInfo[k];
                    auto const& meshEntity = primInfo.meshEntity();
                    auto const& pt0 = meshEntity.point( 0 );
                    auto const& pt1 = meshEntity.point( 1 );
                    auto const& pt2 = meshEntity.point( 2 );
                    deviceHipTriangles[k].v0 = bvhHip::Vec3( pt0[0], pt0[1], pt0[2] );
                    deviceHipTriangles[k].v1 = bvhHip::Vec3( pt1[0], pt1[1], pt1[2] );
                    deviceHipTriangles[k].v2 = bvhHip::Vec3( pt2[0], pt2[1], pt2[2] );
                    deviceHipTriangles[k].id = id;
                }
#endif

#if 0
                HIP_ASSERT( hipMallocManaged( &deviceHipTriangles, numTriangles * sizeof( bvhHip::Triangle ) ) );
                HIP_ASSERT( hipMallocManaged( &devicebvhHipNodes, ( 2 * numTriangles - 1 ) * sizeof( bvhHip::BVHNode ) ) );
                auto initializeTrianglesRange = [this](bvhHip::Triangle* deviceHipTriangles,
                                                size_t start, size_t end) {
                for (size_t k = start; k < end; ++k)
                {
                    auto const& primInfo = this->M_primitiveInfo[k];
                    auto const& meshEntity = primInfo.meshEntity();
                    auto const& pt0 = meshEntity.point(0);
                    auto const& pt1 = meshEntity.point(1);
                    auto const& pt2 = meshEntity.point(2);
                    deviceHipTriangles[k].v0 = bvhHip::Vec3(pt0[0], pt0[1], pt0[2]);
                    deviceHipTriangles[k].v1 = bvhHip::Vec3(pt1[0], pt1[1], pt1[2]);
                    deviceHipTriangles[k].v2 = bvhHip::Vec3(pt2[0], pt2[1], pt2[2]);
                    deviceHipTriangles[k].id = meshEntity.id();
                }
                };

                size_t numThreads = std::thread::hardware_concurrency();
                std::vector<std::thread> threads;
                size_t chunkSize = numTriangles / numThreads;
                size_t remainder = numTriangles % numThreads;
                size_t start = 0;
                for (size_t i = 0; i < numThreads; ++i)
                {
                    size_t end = start + chunkSize + (i < remainder ? 1 : 0);
                    threads.emplace_back(initializeTrianglesRange, deviceHipTriangles, start, end);
                    start = end;
                }

                for (auto& thread : threads)
                {
                    thread.join();
                }
#endif

                auto timeDataTransfertDuration = toc( "timeDataTransfertTriangleDuration" );
                if ( isView ) std::cout << "[INFO GPU]: TransfertDataTriangle With Unified Memory=" << numTriangles << " [" << rank << ":" << numDevice << "]"
                                        << "\n";

                t_begin_bvh_gpu = std::chrono::steady_clock::now();
                if ( numVersion == 0 ) bvhHip::buildBVH_GPU_Version2( deviceHipTriangles, devicebvhHipNodes, numTriangles );
                if ( numVersion == 1 ) bvhHip::buildBVH_GPU_Version3( deviceHipTriangles, devicebvhHipNodes, numTriangles );
                if ( numVersion == 4 ) bvhHip::buildBVH_GPU_Version4( deviceHipTriangles, devicebvhHipNodes, numTriangles );
                if ( numVersion == 2 ) bvhHip::buildBVH_GPU_Parallel_Best_Axis( deviceHipTriangles, devicebvhHipNodes, numTriangles );
                t_end_bvh_gpu = std::chrono::steady_clock::now();

                if ( isView ) std::cout << "[INFO GPU]: BVH GPU FINISHED"
                                        << " [" << rank << ":" << numDevice << "]"
                                        << "\n";
                t_laps_bvh_gpu = std::chrono::duration_cast<std::chrono::microseconds>( t_end_bvh_gpu - t_begin_bvh_gpu ).count();
                if ( isView ) std::cout << "[INFO GPU]: Elapsed microseconds inside BVH  : " << t_laps_bvh_gpu << " us"
                                        << " [" << rank << ":" << numDevice << "]"
                                        << "\n";
                nbTriangles = numTriangles;
                nbNodes = 2 * nbTriangles - 1;
            }

        } // END modeGPU==1

        if ( modeGPU == 4 ) // LBVH Explorer
        {
            std::vector<bvhLinearExplorer::Triangle> triangles;
            for ( int k = 0; k < this->M_primitiveInfo.size(); ++k )
            {
                int id = this->M_primitiveInfo[k].meshEntity().id();
                auto const& primInfo = this->M_primitiveInfo[k];
                auto const& meshEntity = primInfo.meshEntity();
                auto const& pt0 = meshEntity.point( 0 );
                auto const& pt1 = meshEntity.point( 1 );
                auto const& pt2 = meshEntity.point( 2 );
                bvhLinearExplorer::Triangle ltri;
                ltri.v1 = make_float4( pt0[0], pt0[1], pt0[2], 1.0f );
                ltri.v2 = make_float4( pt1[0], pt1[1], pt1[2], 1.0f );
                ltri.v3 = make_float4( pt2[0], pt2[1], pt2[2], 1.0f );
                ltri.id = id;
                triangles.push_back( ltri );
            }
            bvhl = lbvh::bvh<float, bvhLinearExplorer::Triangle, bvhLinearExplorer::aabb_getter>( triangles.begin(), triangles.end(), true );
            bvhl_dev = bvhl.get_device_repr();
            lbvh::aabb<float> global_aabb = thrust::reduce(
                thrust::device, bvhl_dev.aabbs, bvhl_dev.aabbs + bvhl_dev.num_objects,
                lbvh::aabb<float>(), bvhLinearExplorer::merge_aabb() );
            bvhLinearExplorer::CenterGlobalSpaceBox h_gBox;
            // bvhLinearExplorer::CenterGlobalSpaceBox *d_gBox;
            bvhLinearExplorer::buildInformationGlobalAABB( global_aabb, h_gBox );
            hipMalloc( (void**)&d_gBox, sizeof( bvhLinearExplorer::CenterGlobalSpaceBox ) );
            hipMemcpy( d_gBox, &h_gBox, sizeof( bvhLinearExplorer::CenterGlobalSpaceBox ), hipMemcpyHostToDevice );
        } // END modeGPU==4
    }

  private:
    // Value indicating whether we are in GPU mode
    bool isGPUHip() override
    {
        return ( true );
    }

    // Specifies the GPU device number
    void setNumDevice( int v )
    {
        int nbDevices = 0;
        hipGetDeviceCount( &nbDevices );
        if ( v > nbDevices )
        {
            v = 0;
        }
        numDevice = v;
    }

    // using rayintersection_result_type = BVHRayIntersectionResult;
    std::vector<std::vector<rayintersection_result_type>> intersectAllRaysWithGPU( std::vector<ray_type> const& rayons ) override
    {
        bool isModeDirectInDevice = false;
        isModeDirectInDevice = true;
        // static constexpr bool isAnyHit = false;
        size_t numRays = rayons.size();
        // isView = true;

        std::vector<std::vector<rayintersection_result_type>> resALL;
        // resALL.reserve(numRays);
        std::vector<rayintersection_result_type> res;
        // res.reserve(numRays);

        if ( modeGPU == 1 ) // mode hip
        {

            if ( isView ) std::cout << " [INFO GPU]: [BEGIN::LIST RAYs]"
                                    << " [" << rank << ":" << numDevice << "]"
                                    << "\n";

            bvhHip::Ray* deviceHipRays;
            // isUnifiedMemory=true;

            tic();
            hipEvent_t start1, stop1;
            hipEventCreate( &start1 );
            hipEventCreate( &stop1 );
            hipEventRecord( start1 );

            if ( !isUnifiedMemory )
            {
                std::vector<bvhHip::Ray> hostHipRays;
                // Load Ray
                for ( size_t k = 0; k < numRays; ++k )
                {
                    bvhHip::Ray ray;
                    ray.origin = bvhHip::Vec3( rayons[k].origin()[0], rayons[k].origin()[1], rayons[k].origin()[2] );
                    ray.direction = bvhHip::Vec3( rayons[k].dir()[0], rayons[k].dir()[1], rayons[k].dir()[2] );
                    ray.id = rayons[k].id;
                    hostHipRays.push_back( ray );
                }
                HIP_ASSERT( hipMalloc( &deviceHipRays, hostHipRays.size() * sizeof( bvhHip::Ray ) ) );
                HIP_ASSERT( hipMemcpy( deviceHipRays, hostHipRays.data(), hostHipRays.size() * sizeof( bvhHip::Ray ), hipMemcpyHostToDevice ) );
            }
            else
            {
                HIP_ASSERT( hipMallocManaged( &deviceHipRays, numRays * sizeof( bvhHip::Ray ) ) );
                for ( size_t k = 0; k < numRays; ++k )
                {
                    deviceHipRays[k].origin = bvhHip::Vec3( rayons[k].origin()[0], rayons[k].origin()[1], rayons[k].origin()[2] );
                    deviceHipRays[k].direction = bvhHip::Vec3( rayons[k].dir()[0], rayons[k].dir()[1], rayons[k].dir()[2] );
                }
            }

            if ( isView ) std::cout << "[INFO GPU]: [END::LIST RAYs]"
                                    << " [" << rank << ":" << numDevice << "]"
                                    << "\n";
            if ( isView ) std::cout << "[INFO GPU]: [BEGIN::RAYS TRACING]"
                                    << " [" << rank << ":" << numDevice << "]"
                                    << "\n";
            int* deviceHipHitTriangles;
            bvhHip::Vec3* deviceHipIntersectionPoint;
            float* deviceHipDistanceResults;
            int* deviceHipIdResults;

            HIP_ASSERT( hipMalloc( &deviceHipHitTriangles, numRays * sizeof( size_t ) ) );
            HIP_ASSERT( hipMalloc( &deviceHipIntersectionPoint, numRays * sizeof( bvhHip::Vec3 ) ) );
            HIP_ASSERT( hipMalloc( &deviceHipDistanceResults, numRays * sizeof( float ) ) );
            HIP_ASSERT( hipMalloc( &deviceHipIdResults, numRays * sizeof( int ) ) );

            hipEventRecord( stop1 );
            hipEventSynchronize( stop1 );
            float milliseconds1 = 0;
            hipEventElapsedTime( &milliseconds1, start1, stop1 );
            printf( "[INFO GPU]: dt1 duration data CPU to GPU =%f s\n", milliseconds1 / 1000.0 );
            hipEventDestroy( start1 );
            hipEventDestroy( stop1 );
            auto timeDataTransfertAllRaysCPUtoGPU = toc( "timeDataTransfertAllRaysCPUtoGPU" );

            tic();
            hipEvent_t start2, stop2;
            hipEventCreate( &start2 );
            hipEventCreate( &stop2 );
            hipEventRecord( start2 );

            // size_t blockSize = 512;
            size_t blockSize = 1024;
            size_t numBlocks = ( numRays + blockSize - 1 ) / blockSize;

            if (( numVersion == 0 ) || ( numVersion == 1 )  || ( numVersion == 4 ))
            {
                hipLaunchKernelGGL( bvhHip::raytraceKernel, dim3( numBlocks ), dim3( blockSize ), 0, 0,
                                    deviceHipRays,
                                    numRays,
                                    devicebvhHipNodes,
                                    deviceHipTriangles,
                                    deviceHipHitTriangles,
                                    deviceHipDistanceResults,
                                    deviceHipIntersectionPoint,
                                    deviceHipIdResults );
            }

            if ( numVersion == 2 )
            {
                hipLaunchKernelGGL( bvhHip::raytraceKernel_Parallel3, dim3( numBlocks ), dim3( blockSize ), 0, 0,
                                    deviceHipRays,
                                    numRays,
                                    devicebvhHipNodes,
                                    nbNodes,
                                    deviceHipTriangles,
                                    nbTriangles,
                                    deviceHipHitTriangles,
                                    deviceHipDistanceResults,
                                    deviceHipIntersectionPoint,
                                    deviceHipIdResults );
            }

            hipEventRecord( stop2 );
            hipEventSynchronize( stop2 );
            float milliseconds2 = 0;
            hipEventElapsedTime( &milliseconds2, start2, stop2 );
            printf( "[INFO GPU]: dt2 duration RT in Kernel=%f s\n", milliseconds2 / 1000.0 );
            hipEventDestroy( start2 );
            hipEventDestroy( stop2 );
            auto timeDurationInKernel = toc( "timeDurationInKernel" );

            tic();
            hipEvent_t start3, stop3;
            hipEventCreate( &start3 );
            hipEventCreate( &stop3 );
            hipEventRecord( start3 );

            std::vector<size_t> hostHipHitTriangles( numRays );
            hipMemcpy( hostHipHitTriangles.data(), deviceHipHitTriangles, numRays * sizeof( size_t ), hipMemcpyDeviceToHost );

            std::vector<bvhHip::Vec3> hostHipIntersectionPoint( numRays );
            hipMemcpy( hostHipIntersectionPoint.data(), deviceHipIntersectionPoint, numRays * sizeof( bvhHip::Vec3 ), hipMemcpyDeviceToHost );

            std::vector<float> hostHipDistanceResults( numRays );
            hipMemcpy( hostHipDistanceResults.data(), deviceHipDistanceResults, numRays * sizeof( float ), hipMemcpyDeviceToHost );

            std::vector<int> hostHipIdResults( numRays );
            hipMemcpy( hostHipIdResults.data(), deviceHipIdResults, numRays * sizeof( int ), hipMemcpyDeviceToHost );

            hipEventRecord( stop3 );
            hipEventSynchronize( stop3 );
            float milliseconds3 = 0;
            hipEventElapsedTime( &milliseconds3, start3, stop3 );
            printf( "[INFO GPU]: dt3 duration data GPU to CPU=%f s \n", milliseconds3 / 1000.0 );
            hipEventDestroy( start3 );
            hipEventDestroy( stop3 );
            auto timeDurationDataTransfertResultsGPUtoCPU = toc( "timeDurationDataTransfertResultsGPUtoCPU" );

            if ( isView ) std::cout << "[INFO GPU]: [END::RAYS TRACING]"
                                    << " [" << rank << ":" << numDevice << "]"
                                    << "\n";
            if ( isView ) std::cout << "[INFO GPU]: [BEGIN::DEBRIFING COLLISION]"
                                    << " [" << rank << ":" << numDevice << "]"
                                    << "\n";

            // Reading the results and transmitting the information that will be used later
            tic();

            for ( size_t i = 0; i < numRays; ++i )
            {
                double M_distance = std::numeric_limits<double>::max();
                size_t numId = -1;

                // if (hostHipHitResults[i]!=-1)
                if ( hostHipIdResults[i] != -1 )
                {
                    if ( isViewDataRT )
                    {
                        // std::cout<<"      Intersection found with Num Ray ["<<i<<"] ori= <"<<hostHipRays[i].origin.x<<","<<hostHipRays[i].origin.y<<","<<hostHipRays[i].origin.z<<"> ";
                        // std::cout<<" dir= <"<<hostHipRays[i].direction.x<<","<<hostHipRays[i].direction.y<<","<<hostHipRays[i].direction.z<<"> ";
                        std::cout << " [INFO GPU]: dist (min)=" << hostHipDistanceResults[i];
                        std::cout << " IntersectionPoint= <" << hostHipIntersectionPoint[i].x << "," << hostHipIntersectionPoint[i].y << "," << hostHipIntersectionPoint[i].z << "> ";
                        std::cout << " IdObject= " << hostHipIdResults[i] << " [" << rank << ":" << numDevice << "]"
                                  << "\n";
                    }

                    M_distance = fabs( double( hostHipDistanceResults[i] ) );

                    res.push_back( rayintersection_result_type( this->worldComm().rank(), hostHipIdResults[i], M_distance ) ); //(rank,idPrimitiv,distance)
                    res.back().setCoordinates( vector_realdim_type{ { hostHipIntersectionPoint[i].x, hostHipIntersectionPoint[i].y, hostHipIntersectionPoint[i].z } } );
                    res.back().M_id = rayons[i].id;
                    res.resize( 1 );
                    resALL.push_back( std::move( res ) );
                }
                else
                {
                    // TODO: Define what is returned, if there is no intersection point.
                    // res.push_back( rayintersection_result_type(this->worldComm().rank(),0, M_distance)); // No Collision
                    // res.back().setCoordinates(vector_realdim_type{{M_distance,M_distance,M_distance}});
                    // res.resize(1);
                    //
                    resALL.push_back( std::move( res ) );
                }
            }

            if ( isView ) std::cout << "[INFO GPU]: [END::DEBRIFING COLLISION]"
                                    << " [" << rank << ":" << numDevice << "]"
                                    << "\n";

            // Memory cleaning
            if ( isView ) std::cout << "[INFO GPU]: [END::MEMORY CLEANING]"
                                    << " [" << rank << ":" << numDevice << "]"
                                    << "\n";
            hipFree( deviceHipHitTriangles );
            hipFree( deviceHipDistanceResults );
            hipFree( deviceHipIntersectionPoint );
            hipFree( deviceHipHitTriangles );

            hostHipHitTriangles.clear();
            hostHipIntersectionPoint.clear();
            hostHipDistanceResults.clear();
            hostHipIdResults.clear();

            auto timeDataResultsDuration = toc( "timeDataResultsDuration" );

        } // END mode hip

        if ( modeGPU == 4 ) // mode LBVH Explorer
        {

            if ( isView ) std::cout << "[BEGIN::LIST RAYs]"
                                    << "\n";

            std::vector<bvhLinearExplorer::Ray> hostRays;
            // Load Ray
            for ( int k = 0; k < numRays; ++k )
            {
                bvhLinearExplorer::Ray r;
                r.origin = make_float4( rayons[k].origin()[0], rayons[k].origin()[1], rayons[k].origin()[2], 1.0f );
                r.direction = make_float4( rayons[k].dir()[0], rayons[k].dir()[1], rayons[k].dir()[2], 0.0f );
                r.id = rayons[k].id; // add
                normalizeRayDirection( r );

                hostRays.push_back( r );
            }

            bvhLinearExplorer::Ray* deviceRays;
            HIP_ASSERT( hipMalloc( &deviceRays, hostRays.size() * sizeof( bvhLinearExplorer::Ray ) ) );
            HIP_ASSERT( hipMemcpy( deviceRays, hostRays.data(), numRays * sizeof( bvhLinearExplorer::Ray ), hipMemcpyHostToDevice ) );

            hostRays.clear();
            if ( isView ) std::cout << "[END::LIST RAYs]"
                                    << "\n";

            if ( isView ) std::cout << "[BEGIN::RAYS TRACING]"
                                    << "\n";
            bvhLinearExplorer::HitRay* deviceHitRays;
            HIP_ASSERT( hipMalloc( &deviceHitRays, numRays * sizeof( bvhLinearExplorer::HitRay ) ) );

            int threadsPerBlock = 1024;
            int blocksPerGrid = ( numRays + threadsPerBlock - 1 ) / threadsPerBlock;

            ////
            const int numDirections = 14;
            float4* h_directions;
            float4* d_directions;
            h_directions = new float4[numDirections];
            bvhLinearExplorer::initializeDirections( h_directions );
            hipMalloc( &d_directions, numDirections * sizeof( float4 ) );
            bvhLinearExplorer::initializeDirectionsKernel<<<1, 1>>>( d_directions );

            bvhLinearExplorer::rayTracingKernelExplorationOptimized<float, bvhLinearExplorer::Triangle>
                <<<blocksPerGrid, threadsPerBlock>>>( bvhl_dev, deviceRays, deviceHitRays, numRays, d_directions, d_gBox );

            hipFree( d_directions );
            delete[] h_directions;
            ////

            hipDeviceSynchronize();
            std::vector<bvhLinearExplorer::HitRay> hostHitRays( numRays );
            HIP_ASSERT( hipMemcpy( hostHitRays.data(), deviceHitRays, numRays * sizeof( bvhLinearExplorer::HitRay ), hipMemcpyDeviceToHost ) );
            if ( isView ) std::cout << "[END::RAYS TRACING]"
                                    << "\n";

            if ( isView ) std::cout << "[BEGIN::DEBRIFING COLLISION]"
                                    << "\n";
            for ( int i = 0; i < numRays; ++i )
            {
                double M_distance = std::numeric_limits<double>::max();
                int numId = -1;

                if ( hostHitRays[i].idResults != -1 )
                {
                    if ( isView )
                    {
                        std::cout << " dist (min)=" << hostHitRays[i].distanceResults;
                        std::cout << " IntersectionPoint= <" << hostHitRays[i].intersectionPoint.x << "," << hostHitRays[i].intersectionPoint.y << "," << hostHitRays[i].intersectionPoint.z << "> ";
                        std::cout << " IdObject= " << hostHitRays[i].idResults << "\n";
                    }

                    M_distance = double( hostHitRays[i].distanceResults );

                    res.push_back( rayintersection_result_type( this->worldComm().rank(), hostHitRays[i].idResults, M_distance ) ); //(rank,idPrimitiv,distance)
                    res.back().setCoordinates( vector_realdim_type{ { hostHitRays[i].intersectionPoint.x, hostHitRays[i].intersectionPoint.y, hostHitRays[i].intersectionPoint.z } } );
                    res.back().M_id = rayons[i].id; // add
                    res.resize( 1 );
                    resALL.push_back( std::move( res ) );
                }
                else
                {
                    // TODO: Define what is returned, if there is no intersection point.
                    // res.push_back( rayintersection_result_type(this->worldComm().rank(),0, M_distance)); // No Collision
                    // res.back().setCoordinates(vector_realdim_type{{M_distance,M_distance,M_distance}});
                    // res.resize(1);
                    //
                    resALL.push_back( std::move( res ) );
                }
            }

            if ( isView ) std::cout << "[END::DEBRIFING COLLISION]"
                                    << "\n";

            // Memory cleaning
            if ( isView ) std::cout << "[END::MEMORY CLEANING]"
                                    << "\n";
            hostHitRays.clear();
            hipFree( deviceRays );
            hipFree( deviceHitRays );

        } ////END_LBVH

        return ( resALL );
    }

    vector_realdim_type CartesianCoordinates( float x, float y, float z ) const
    { // Deleted maybe later
        return vector_realdim_type{ {
            x,
            y,
            z,
        } };
    }

    std::vector<rayintersection_result_type> intersectSequential( ray_type const& ray, bool useRobustTraversal = true ) override
    {
        std::vector<rayintersection_result_type> res;
        return ( res );
    };
};
#endif

/***************************************************************************************************************************************************/
// NOTA : The objective of this function is to initiate the calculation from a CPU across multiple GPUs.
// The first step is to compute the BVH on one GPU. Retrieve the results and send them to all GPUs.
// The second step is ray tracing. We distribute the task of ray tracing across multiple GPUs, then we collect the result.

#if defined( FEELPP_HAS_HIP )

struct GpuData
{
    bvhHip::Triangle* triangles;
    bvhHip::BVHNode* nodes;
    size_t num_triangles;
    size_t num_nodes;
    size_t num_device;
};

template <typename MeshEntityType>
class BVH_HIP_CPU_GPUs_Party : public BVH<MeshEntityType>
{
    using super_type = BVH<MeshEntityType>;
    using mesh_entity_type = typename super_type::mesh_entity_type;
    using vector_realdim_type = typename super_type::vector_realdim_type;
    static constexpr uint16_type nRealDim = super_type::nRealDim;

  public:
    using ray_type = typename super_type::ray_type;
    using rayintersection_result_type = typename super_type::rayintersection_result_type;

    bvhHip::BVHNode* devicebvhHipNodes;
    bvhHip::Triangle* deviceHipTriangles;
    size_t nbTriangles;
    size_t nbNodes;

    int numDevice;
    int numVersion;
    int modeGPU; // 1- HIP 2-
    bool isUnifiedMemory;
    int rank;
    int num_gpus;
    int nbWorkDistributionGPUs;
    int deviceIdBegin;
    std::vector<GpuData> deviceInformationForGPUs;
    bool isView;
    bool isViewDataRT;

    BVH_HIP_CPU_GPUs_Party( BVHEnum::Quality quality, worldcomm_ptr_t worldComm )
        : super_type( quality, worldComm )
    {
        num_gpus = 0;
        hipGetDeviceCount( &num_gpus );
        // numDevice = 0;
        // numDevice = 3;
        rank = worldComm->rank();
        numDevice = rank % num_gpus;
        printf( "[INFO GPU]: constructor rank %i numDevice=%i  num_gpus=%i\n", rank, numDevice, num_gpus );
        numVersion = 2;
        modeGPU = 1;
        isUnifiedMemory = true; // isUnifiedMemory = false;

        nbWorkDistributionGPUs = num_gpus;
        deviceIdBegin = 0;
        isView = true;
        isView = false;
        isViewDataRT = false;
    }

    ~BVH_HIP_CPU_GPUs_Party()
    {
        // isView=true;
        //  Memory cleaning
        for ( int device_id = 0; device_id < nbWorkDistributionGPUs; ++device_id )
        {
            HIP_ASSERT( hipFree( deviceInformationForGPUs[device_id].triangles ) );
            HIP_ASSERT( hipFree( deviceInformationForGPUs[device_id].nodes ) );
            if ( isView ) std::cout << "[INFO GPU]: [MEMORY CLEANING MULTI-GPU]"
                                    << " [" << rank << "] Device: " << device_id << "\n";
        }
        if ( isView ) std::cout << "[INFO GPU]: [MULTI-GPU MEMORY CLEANING DONE]"
                                << "\n";
    }

    BVH_HIP_CPU_GPUs_Party( BVH_HIP_CPU_GPUs_Party&& ) = default;

    template <typename RangeType>
    void
    updateForUse( RangeType const& range )
    {

        // isView = true;
        //  up primitiveinfos
        super_type::updateForUse( range );
        // init bvh backend
        if ( isView ) std::cout << "[INFO GPU]: Size primitiveInfo=" << this->M_primitiveInfo.size() << " [" << rank << ":" << numDevice << "]"
                                << "\n";
        //...

        // hip device used
        hipSetDevice( numDevice );
        int numDeviceActivated;
        hipGetDevice( &numDeviceActivated );
        if ( isView ) std::cout << "[INFO GPU]: Num Device Activated=" << numDeviceActivated << " [" << rank << ":" << numDevice << "]"
                                << "\n";

        // Definition of operating modes
        bool isModeBox = true;
        isModeBox = false;
        bool isModeDirectInDevice = false;
        // isModeDirectInDevice = true; // Todo CTRL in infinity and size max of HIP GPU

        std::chrono::steady_clock::time_point t_begin_bvh_gpu, t_end_bvh_gpu;
        long int t_laps_bvh_gpu = 0;

        if ( modeGPU == 1 )
        {

            tic();
            size_t numTriangles = this->M_primitiveInfo.size();

#if 1
            HIP_ASSERT( hipMallocManaged( &deviceHipTriangles, numTriangles * sizeof( bvhHip::Triangle ) ) );
            HIP_ASSERT( hipMallocManaged( &devicebvhHipNodes, ( 2 * numTriangles - 1 ) * sizeof( bvhHip::BVHNode ) ) );
            // Load Mesh  in host-device
            for ( size_t k = 0; k < this->M_primitiveInfo.size(); ++k )
            {
                size_t id = this->M_primitiveInfo[k].meshEntity().id();
                auto const& primInfo = this->M_primitiveInfo[k];
                auto const& meshEntity = primInfo.meshEntity();
                auto const& pt0 = meshEntity.point( 0 );
                auto const& pt1 = meshEntity.point( 1 );
                auto const& pt2 = meshEntity.point( 2 );
                deviceHipTriangles[k].v0 = bvhHip::Vec3( pt0[0], pt0[1], pt0[2] );
                deviceHipTriangles[k].v1 = bvhHip::Vec3( pt1[0], pt1[1], pt1[2] );
                deviceHipTriangles[k].v2 = bvhHip::Vec3( pt2[0], pt2[1], pt2[2] );
                deviceHipTriangles[k].id = id;
            }
#endif

#if 0
            HIP_ASSERT( hipMallocManaged( &deviceHipTriangles, numTriangles * sizeof( bvhHip::Triangle ) ) );
            HIP_ASSERT( hipMallocManaged( &devicebvhHipNodes, ( 2 * numTriangles - 1 ) * sizeof( bvhHip::BVHNode ) ) );
            auto initializeTrianglesRange = [this](bvhHip::Triangle* deviceHipTriangles,
                                               size_t start, size_t end) {
            for (size_t k = start; k < end; ++k)
            {
                auto const& primInfo = this->M_primitiveInfo[k];
                auto const& meshEntity = primInfo.meshEntity();
                auto const& pt0 = meshEntity.point(0);
                auto const& pt1 = meshEntity.point(1);
                auto const& pt2 = meshEntity.point(2);
                deviceHipTriangles[k].v0 = bvhHip::Vec3(pt0[0], pt0[1], pt0[2]);
                deviceHipTriangles[k].v1 = bvhHip::Vec3(pt1[0], pt1[1], pt1[2]);
                deviceHipTriangles[k].v2 = bvhHip::Vec3(pt2[0], pt2[1], pt2[2]);
                deviceHipTriangles[k].id = meshEntity.id();
            }
            };

            size_t numThreads = std::thread::hardware_concurrency();
            //size_t numThreads = 3;
            std::vector<std::thread> threads;
            size_t chunkSize = numTriangles / numThreads;
            size_t remainder = numTriangles % numThreads;
            size_t start = 0;
            for (size_t i = 0; i < numThreads; ++i)
            {
                size_t end = start + chunkSize + (i < remainder ? 1 : 0);
                threads.emplace_back(initializeTrianglesRange, deviceHipTriangles, start, end);
                start = end;
            }

            for (auto& thread : threads)
            {
                thread.join();
            }
#endif

#if 0

            std::vector<bvhHip::Triangle> hostTriangles( numTriangles );
            auto initializeTrianglesRange = [this, &hostTriangles](
                                                size_t start, size_t end )
            {
                for ( size_t k = start; k < end; ++k )
                {
                    auto const& primInfo = this->M_primitiveInfo[k];
                    auto const& meshEntity = primInfo.meshEntity();
                    size_t id = meshEntity.id();
                    auto const& pt0 = meshEntity.point( 0 );
                    auto const& pt1 = meshEntity.point( 1 );
                    auto const& pt2 = meshEntity.point( 2 );
                    hostTriangles[k].v0 = bvhHip::Vec3( pt0[0], pt0[1], pt0[2] );
                    hostTriangles[k].v1 = bvhHip::Vec3( pt1[0], pt1[1], pt1[2] );
                    hostTriangles[k].v2 = bvhHip::Vec3( pt2[0], pt2[1], pt2[2] );
                    hostTriangles[k].id = id;
                }
            };

            size_t numThreads = std::thread::hardware_concurrency();
            std::vector<std::thread> threads;
            size_t chunkSize = numTriangles / numThreads;
            size_t remainder = numTriangles % numThreads;
            size_t start = 0;
            for ( size_t i = 0; i < numThreads; ++i )
            {
                size_t end = start + chunkSize + ( i < remainder ? 1 : 0 );
                threads.emplace_back( initializeTrianglesRange, start, end );
                start = end;
            }

            for ( auto& thread : threads )
            {
                thread.join();
            }

            size_t trianglesSize = numTriangles * sizeof( bvhHip::Triangle );
            size_t nodesSize = ( 2 * numTriangles - 1 ) * sizeof( bvhHip::BVHNode );
            HIP_CHECK( hipMalloc( &deviceHipTriangles, trianglesSize ) );
            HIP_CHECK( hipMalloc( &devicebvhHipNodes, nodesSize ) );
            HIP_CHECK( hipMemcpy( deviceHipTriangles, hostTriangles.data(), trianglesSize, hipMemcpyHostToDevice ) );
#endif

            auto timeDataTransfertDuration = toc( "timeDataTransfertTriangleDuration" );
            if ( isView ) std::cout << "[INFO GPU]: TransfertDataTriangle With Unified Memory=" << numTriangles << " [" << rank << ":" << numDevice << "]"
                                    << "\n";

            t_begin_bvh_gpu = std::chrono::steady_clock::now();
            if ( numVersion == 0 ) bvhHip::buildBVH_GPU_Version2( deviceHipTriangles, devicebvhHipNodes, numTriangles );
            if ( numVersion == 1 ) bvhHip::buildBVH_GPU_Version3( deviceHipTriangles, devicebvhHipNodes, numTriangles );
            if ( numVersion == 2 ) bvhHip::buildBVH_GPU_Parallel_Best_Axis( deviceHipTriangles, devicebvhHipNodes, numTriangles );
            t_end_bvh_gpu = std::chrono::steady_clock::now();

            if ( isView ) std::cout << "[INFO GPU]: BVH GPU FINISHED"
                                    << " [" << rank << ":" << numDevice << "]"
                                    << "\n";
            t_laps_bvh_gpu = std::chrono::duration_cast<std::chrono::microseconds>( t_end_bvh_gpu - t_begin_bvh_gpu ).count();
            if ( isView ) std::cout << "[INFO GPU]: Elapsed microseconds inside BVH  : " << t_laps_bvh_gpu << " us"
                                    << " [" << rank << ":" << numDevice << "]"
                                    << "\n";
            nbTriangles = numTriangles;
            nbNodes = 2 * nbTriangles - 1;

            // BEGIN::Dispach Informations
            std::cout << "[INFO]: Number of HIP devices found: " << num_gpus << "\n";
            std::cout << "[INFO]: Dispatch Number of HIP devices: " << nbWorkDistributionGPUs << "\n";

            std::chrono::steady_clock::time_point t_begin_stream_gpu, t_end_stream_gpu;
            t_begin_stream_gpu = std::chrono::steady_clock::now();
            long int t_laps_stream_gpu = 0;

            hipStream_t stream;
            HIP_ASSERT( hipStreamCreate( &stream ) );

            // Prefetch data to each GPU
            for ( int device_id = 0; device_id < nbWorkDistributionGPUs; ++device_id )
            {
                HIP_ASSERT( hipSetDevice( device_id ) );
                HIP_ASSERT( hipMemPrefetchAsync( deviceHipTriangles, numTriangles * sizeof( bvhHip::Triangle ), device_id, stream ) );
                HIP_ASSERT( hipMemPrefetchAsync( devicebvhHipNodes, ( 2 * numTriangles - 1 ) * sizeof( bvhHip::BVHNode ), device_id, stream ) );
            }

            // Wait for prefetching to complete
            HIP_ASSERT( hipStreamSynchronize( stream ) );

            for ( int device_id = 0; device_id < nbWorkDistributionGPUs; ++device_id )
            {
                // Configuration du GPU courant
                HIP_ASSERT( hipSetDevice( device_id ) );
                std::cout << "[INFO]: Sending data to device: " << device_id << "\n";
                GpuData tmp;
                // Memory allocation for each GPU
                tmp.num_device = device_id;
                tmp.triangles = deviceHipTriangles;
                tmp.nodes = devicebvhHipNodes;
                tmp.num_triangles = numTriangles;
                tmp.num_nodes = 2 * numTriangles - 1;
                deviceInformationForGPUs.push_back( tmp );
                std::cout << "[INFO]: Data is available on device " << device_id << " (Unified Memory).\n";
            }
            HIP_ASSERT( hipStreamDestroy( stream ) );

            t_end_stream_gpu = std::chrono::steady_clock::now();

            t_laps_stream_gpu = std::chrono::duration_cast<std::chrono::microseconds>( t_end_stream_gpu - t_begin_stream_gpu ).count();
            if ( isView ) std::cout << "[INFO GPU]: Elapsed microseconds inside Sending all data to device  : " << t_laps_stream_gpu << " us"
                                    << " [" << rank << ":" << numDevice << "]"
                                    << "\n";

            // END::Dispach Informations

            // Back to numDevice
            HIP_ASSERT( hipSetDevice( numDevice ) );

        } // END modeGPU==1
    }

  private:
    // Value indicating whether we are in GPU mode
    bool isGPUHip() override
    {
        return ( true );
    }

    // Specifies the GPU device number
    void setNumDevice( int v )
    {
        int nbDevices = 0;
        hipGetDeviceCount( &nbDevices );
        if ( v > nbDevices )
        {
            v = 0;
        }
        numDevice = v;
    }

    // Specifies the number of desired working GPUs
    void setNbDesiredWorkingGPU( int v )
    {
        if ( v > num_gpus ) v = num_gpus;
        if ( v < 1 ) v = num_gpus;
        nbWorkDistributionGPUs = v;
    }

    std::vector<std::vector<rayintersection_result_type>> processRaysOnGPU( const std::vector<ray_type>& rayons, int deviceId )
    {
        bool isModeDirectInDevice = true;
        size_t numRays = rayons.size();
        // isView = true;

        std::vector<std::vector<rayintersection_result_type>> resALL;
        // resALL.reserve(numRays);

        if ( modeGPU == 1 ) // mode hip
        {
            if ( isView ) std::cout << " [INFO GPU]: [BEGIN::LIST RAYs]"
                                    << " [" << rank << ":" << numDevice << "] Device: " << deviceId << "\n";

            bvhHip::Ray* deviceHipRays;
            tic();
            hipEvent_t start1, stop1;
            HIP_ASSERT( hipEventCreate( &start1 ) );
            HIP_ASSERT( hipEventCreate( &stop1 ) );
            HIP_ASSERT( hipEventRecord( start1 ) );

            if ( !isUnifiedMemory )
            {
                std::vector<bvhHip::Ray> hostHipRays;
                // Load Ray
                for ( size_t k = 0; k < numRays; ++k )
                {
                    bvhHip::Ray ray;
                    ray.origin = bvhHip::Vec3( rayons[k].origin()[0], rayons[k].origin()[1], rayons[k].origin()[2] );
                    ray.direction = bvhHip::Vec3( rayons[k].dir()[0], rayons[k].dir()[1], rayons[k].dir()[2] );
                    ray.id = rayons[k].id;
                    hostHipRays.push_back( ray );
                }
                HIP_ASSERT( hipMalloc( &deviceHipRays, hostHipRays.size() * sizeof( bvhHip::Ray ) ) );
                HIP_ASSERT( hipMemcpy( deviceHipRays, hostHipRays.data(), hostHipRays.size() * sizeof( bvhHip::Ray ), hipMemcpyHostToDevice ) );
            }
            else
            {
                HIP_ASSERT( hipMallocManaged( &deviceHipRays, numRays * sizeof( bvhHip::Ray ) ) );
                for ( size_t k = 0; k < numRays; ++k )
                {
                    deviceHipRays[k].origin = bvhHip::Vec3( rayons[k].origin()[0], rayons[k].origin()[1], rayons[k].origin()[2] );
                    deviceHipRays[k].direction = bvhHip::Vec3( rayons[k].dir()[0], rayons[k].dir()[1], rayons[k].dir()[2] );
                }
            }

            if ( isView ) std::cout << "[INFO GPU]: [END::LIST RAYs]"
                                    << " [" << rank << ":" << numDevice << "] Device: " << deviceId << "\n";
            if ( isView ) std::cout << "[INFO GPU]: [BEGIN::RAYS TRACING]"
                                    << " [" << rank << ":" << numDevice << "] Device: " << deviceId << "\n";

            int* deviceHipHitTriangles;
            bvhHip::Vec3* deviceHipIntersectionPoint;
            float* deviceHipDistanceResults;
            int* deviceHipIdResults;

            HIP_ASSERT( hipMalloc( &deviceHipHitTriangles, numRays * sizeof( size_t ) ) );
            HIP_ASSERT( hipMalloc( &deviceHipIntersectionPoint, numRays * sizeof( bvhHip::Vec3 ) ) );
            HIP_ASSERT( hipMalloc( &deviceHipDistanceResults, numRays * sizeof( float ) ) );
            HIP_ASSERT( hipMalloc( &deviceHipIdResults, numRays * sizeof( int ) ) );

            HIP_ASSERT( hipEventRecord( stop1 ) );
            HIP_ASSERT( hipEventSynchronize( stop1 ) );
            float milliseconds1 = 0;
            HIP_ASSERT( hipEventElapsedTime( &milliseconds1, start1, stop1 ) );
            if ( isView ) printf( "[INFO GPU]: dt1 duration data CPU to GPU =%f s, Device %d\n", milliseconds1 / 1000.0, deviceId );
            HIP_ASSERT( hipEventDestroy( start1 ) );
            HIP_ASSERT( hipEventDestroy( stop1 ) );
            auto timeDataTransfertAllRaysCPUtoGPU = toc( "timeDataTransfertAllRaysCPUtoGPU" );

            tic();
            hipEvent_t start2, stop2;
            HIP_ASSERT( hipEventCreate( &start2 ) );
            HIP_ASSERT( hipEventCreate( &stop2 ) );
            HIP_ASSERT( hipEventRecord( start2 ) );

            size_t blockSize = 1024;
            size_t numBlocks = ( numRays + blockSize - 1 ) / blockSize;

            if (( numVersion == 0 ) || ( numVersion == 1 )  || ( numVersion == 4 ))
            {
                hipLaunchKernelGGL( bvhHip::raytraceKernel, dim3( numBlocks ), dim3( blockSize ), 0, 0,
                                    deviceHipRays, numRays, devicebvhHipNodes, deviceHipTriangles,
                                    deviceHipHitTriangles, deviceHipDistanceResults,
                                    deviceHipIntersectionPoint, deviceHipIdResults );
            }
            
            if ( numVersion == 2 )
            {
                hipLaunchKernelGGL( bvhHip::raytraceKernel_Parallel3, dim3( numBlocks ), dim3( blockSize ), 0, 0,
                                    deviceHipRays, numRays, devicebvhHipNodes, nbNodes,
                                    deviceHipTriangles, nbTriangles, deviceHipHitTriangles,
                                    deviceHipDistanceResults, deviceHipIntersectionPoint, deviceHipIdResults );
            }

            HIP_ASSERT( hipEventRecord( stop2 ) );
            HIP_ASSERT( hipEventSynchronize( stop2 ) );
            float milliseconds2 = 0;
            HIP_ASSERT( hipEventElapsedTime( &milliseconds2, start2, stop2 ) );
            if ( isView ) printf( "[INFO GPU]: dt2 duration RT in Kernel=%f s, Device %d\n", milliseconds2 / 1000.0, deviceId );
            HIP_ASSERT( hipEventDestroy( start2 ) );
            HIP_ASSERT( hipEventDestroy( stop2 ) );
            auto timeDurationInKernel = toc( "timeDurationInKernel" );

            tic();
            hipEvent_t start3, stop3;
            HIP_ASSERT( hipEventCreate( &start3 ) );
            HIP_ASSERT( hipEventCreate( &stop3 ) );
            HIP_ASSERT( hipEventRecord( start3 ) );

            std::vector<size_t> hostHipHitTriangles( numRays );
            std::vector<bvhHip::Vec3> hostHipIntersectionPoint( numRays );
            std::vector<float> hostHipDistanceResults( numRays );
            std::vector<int> hostHipIdResults( numRays );

            HIP_ASSERT( hipMemcpy( hostHipHitTriangles.data(), deviceHipHitTriangles, numRays * sizeof( size_t ), hipMemcpyDeviceToHost ) );
            HIP_ASSERT( hipMemcpy( hostHipIntersectionPoint.data(), deviceHipIntersectionPoint, numRays * sizeof( bvhHip::Vec3 ), hipMemcpyDeviceToHost ) );
            HIP_ASSERT( hipMemcpy( hostHipDistanceResults.data(), deviceHipDistanceResults, numRays * sizeof( float ), hipMemcpyDeviceToHost ) );
            HIP_ASSERT( hipMemcpy( hostHipIdResults.data(), deviceHipIdResults, numRays * sizeof( int ), hipMemcpyDeviceToHost ) );

            HIP_ASSERT( hipEventRecord( stop3 ) );
            HIP_ASSERT( hipEventSynchronize( stop3 ) );
            float milliseconds3 = 0;
            HIP_ASSERT( hipEventElapsedTime( &milliseconds3, start3, stop3 ) );
            if ( isView ) printf( "[INFO GPU]: dt3 duration data GPU to CPU=%f s , Device %d\n", milliseconds3 / 1000.0, deviceId );
            HIP_ASSERT( hipEventDestroy( start3 ) );
            HIP_ASSERT( hipEventDestroy( stop3 ) );
            auto timeDurationDataTransfertResultsGPUtoCPU = toc( "timeDurationDataTransfertResultsGPUtoCPU" );

            if ( isView ) std::cout << "[INFO GPU]: [END::RAYS TRACING]"
                                    << " [" << rank << ":" << numDevice << "] Device: " << deviceId << "\n";
            if ( isView ) std::cout << "[INFO GPU]: [BEGIN::DEBRIFING COLLISION]"
                                    << " [" << rank << ":" << numDevice << "] Device: " << deviceId << "\n";
            tic();
            for ( size_t i = 0; i < numRays; ++i )
            {
                double M_distance = std::numeric_limits<double>::max();
                size_t numId = -1;

                if ( hostHipIdResults[i] != -1 )
                {
                    if ( isViewDataRT )
                    {
                        std::cout << " [INFO GPU]: dist (min)=" << hostHipDistanceResults[i];
                        std::cout << " IntersectionPoint= <" << hostHipIntersectionPoint[i].x << "," << hostHipIntersectionPoint[i].y << "," << hostHipIntersectionPoint[i].z << "> ";
                        std::cout << " IdObject= " << hostHipIdResults[i] << " [" << rank << ":" << numDevice << "] Device: " << deviceId << "\n";

                        // deviceHipRays[k].origin.x
                    }

                    M_distance = fabs( double( hostHipDistanceResults[i] ) );

                    std::vector<rayintersection_result_type> res;
                    res.push_back( rayintersection_result_type( this->worldComm().rank(), hostHipIdResults[i], M_distance ) );
                    res.back().setCoordinates( vector_realdim_type{ { hostHipIntersectionPoint[i].x, hostHipIntersectionPoint[i].y, hostHipIntersectionPoint[i].z } } );
                    res.back().M_id = rayons[i].id;
                    resALL.push_back( std::move( res ) );
                }
                else
                {
                    resALL.push_back( std::vector<rayintersection_result_type>() );
                }
            }

            if ( isView ) std::cout << "[INFO GPU]: Nb resALL=" << resALL.size() << " [" << rank << ":" << numDevice << "] Device: " << deviceId << "\n";

            if ( isView ) std::cout << "[INFO GPU]: [END::DEBRIFING COLLISION]"
                                    << " [" << rank << ":" << numDevice << "] Device: " << deviceId << "\n";

            // Memory cleaning
            if ( isView ) std::cout << "[INFO GPU]: [BEGIN::MEMORY CLEANING]"
                                    << " [" << rank << ":" << numDevice << "] Device: " << deviceId << "\n";
            HIP_ASSERT( hipFree( deviceHipRays ) );
            HIP_ASSERT( hipFree( deviceHipHitTriangles ) );
            HIP_ASSERT( hipFree( deviceHipDistanceResults ) );
            HIP_ASSERT( hipFree( deviceHipIntersectionPoint ) );
            HIP_ASSERT( hipFree( deviceHipIdResults ) );

            if ( isView ) std::cout << "[INFO GPU]: [END::MEMORY CLEANING]"
                                    << " [" << rank << ":" << numDevice << "] Device: " << deviceId << "\n";

            auto timeDataResultsDuration = toc( "timeDataResultsDuration" );
        }
        else
        {
            // mode lbvh explorer....
        }
        return resALL;
    }

    std::vector<std::vector<rayintersection_result_type>> intersectAllRaysWithGPU( std::vector<ray_type> const& rayons ) override
    {
        size_t numRays = rayons.size();
        // isView = true;

        // Calculate the number of rays per GPU
        if ( nbWorkDistributionGPUs > numRays ) nbWorkDistributionGPUs = 1;

        size_t raysPerGPU = numRays / nbWorkDistributionGPUs;
        size_t remainingRays = numRays % nbWorkDistributionGPUs;

        std::vector<std::thread> threads;
        std::vector<std::vector<std::vector<rayintersection_result_type>>> partialResults( nbWorkDistributionGPUs );

        for ( size_t deviceId = 0; deviceId < nbWorkDistributionGPUs; ++deviceId )
        {
            size_t startIdx = deviceId * raysPerGPU + std::min( deviceId, remainingRays );
            size_t endIdx = startIdx + raysPerGPU + ( deviceId < remainingRays ? 1 : 0 );

            threads.emplace_back( [this, deviceId, &rayons, startIdx, endIdx, &partialResults]()
                                  {
            HIP_ASSERT(hipSetDevice(deviceId));

            std::vector<ray_type> deviceRayons(rayons.begin() + startIdx, rayons.begin() + endIdx);
            partialResults[deviceId] = this->processRaysOnGPU(deviceRayons, deviceId); } );
        }

        // Wait for all threads to finish
        for ( auto& thread : threads )
        {
            thread.join();
        }

        // Merge all the results
        std::vector<std::vector<rayintersection_result_type>> resALL;
        // resALL.reserve(numRays);
        for ( const auto& partialResult : partialResults )
        {
            resALL.insert( resALL.end(), partialResult.begin(), partialResult.end() );
        }

        if ( isView ) std::cout << "[INFO GPU]: Nb resALL Total=" << resALL.size() << " [" << rank << ":" << numDevice << "]"
                                << "\n";

        return resALL;
    }

    vector_realdim_type CartesianCoordinates( float x, float y, float z ) const
    { // Deleted maybe later
        return vector_realdim_type{ {
            x,
            y,
            z,
        } };
    }

    std::vector<rayintersection_result_type> intersectSequential( ray_type const& ray, bool useRobustTraversal = true ) override
    {
        std::vector<rayintersection_result_type> res;
        return ( res );
    };
};
#endif

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
        BVHNode* parent() const { return M_parent; }

        vector_realdim_type const& boundMin() const noexcept { return M_bounds_min; }
        vector_realdim_type const& boundMax() const noexcept { return M_bounds_max; }
        vector_realdim_type centroid() const { return 0.5 * ( M_bounds_min + M_bounds_max ); }
        int splitAxis() const noexcept { return M_splitAxis; }
        int nPrimitives() const noexcept { return M_nPrimitives; }
        int firstPrimOffset() const noexcept { return M_firstPrimOffset; }

        BVHNode* child( int k ) const { return M_children[k].get(); }

        bool isLeaf() const { return !M_children[0] && !M_children[1]; }

        BVHNode* nearChild( ray_type const& ray ) const
        {
            if ( ray.dir()( this->splitAxis() ) > 0 )
                return this->child( 0 );
            else
                return this->child( 1 );
        }

        BVH_InHouse::BVHNode* siblingNode() const
        {
            if ( !M_parent )
                return nullptr;
            return M_parent->child( this == M_parent->child( 0 ) ? 1 : 0 );
        }

        bool checkIntersection( ray_type const& rayon )
        {
            double tmin = 0.0;
            double tmax = FLT_MAX;

            for ( int i = 0; i < nRealDim; i++ )
            {
                double ratio = 1.0 / ( rayon.dir()[i] + 2 * FLT_MIN );
                double t1 = ( M_bounds_min[i] - rayon.origin()[i] ) * ratio;
                double t2 = ( M_bounds_max[i] - rayon.origin()[i] ) * ratio;
                if ( t1 > t2 )
                {
                    double tTemp = t1;
                    t1 = t2;
                    t2 = tTemp;
                }
                if ( t1 > tmin )
                    tmin = t1;
                if ( t2 > tmax )
                    tmax = t2;
                if ( tmin > tmax )
                    return false;
            }

            return true;
        }

        std::pair<bool, double> checkIntersectionWithSegment( ray_type const& ray, std::vector<primitiveinfo_type> const& primitiveInfo ) const
        {
            auto const& meshElt = primitiveInfo[this->firstPrimOffset()].meshEntity();
            auto p1 = Eigen::Map<const Eigen::Matrix<double, nRealDim, 1>>( meshElt.point( 0 ).node().data().begin() );
            auto p2 = Eigen::Map<const Eigen::Matrix<double, nRealDim, 1>>( meshElt.point( 1 ).node().data().begin() );

            auto const& origin = ray.origin();
            auto const& direction = ray.dir();

            vector_realdim_type v1 = origin - p1;
            vector_realdim_type v2 = p2 - p1;
            vector_realdim_type v3{ -direction[1], direction[0] };

            double dot = v2.dot( v3 );
            if ( math::abs( dot ) < 1e-6 )
                return std::make_pair( false, 0 );

            double t1 = ( v2[0] * v1[1] - v2[1] * v1[0] ) / dot;
            double t2 = v1.dot( v3 ) / dot;

            if ( t1 > 2 * FLT_MIN && ( t2 >= 0.0 && t2 <= 1.0 ) )
            {
#if 0
                    vector_realdim_type w_{
                        origin[0] + direction[0]*t1,
                        origin[1] + direction[1]*t1; };
#endif
                return std::make_pair( true, t1 );
            }
            return std::make_pair( false, t1 );
        }

        // Verify if the ray intersects the element
        std::pair<bool, double> checkIntersectionWithTriangle( ray_type const& ray, std::vector<primitiveinfo_type> const& primitiveInfo ) const
        {
            DCHECK( this->isLeaf() ) << "should be a leaf: ";

            auto const& meshElt = primitiveInfo[this->firstPrimOffset()].meshEntity();
            auto p1 = Eigen::Map<const Eigen::Matrix<double, nRealDim, 1>>( meshElt.point( 0 ).node().data().begin() );
            auto p2 = Eigen::Map<const Eigen::Matrix<double, nRealDim, 1>>( meshElt.point( 1 ).node().data().begin() );
            auto p3 = Eigen::Map<const Eigen::Matrix<double, nRealDim, 1>>( meshElt.point( 2 ).node().data().begin() );

            auto const& origin = ray.origin();
            auto const& direction = ray.dir();

            // // normal vector
            auto n1 = ( p2 - p1 ).cross( p3 - p1 );
            n1 = n1 / n1.norm();
            double n_dot_dir = direction.dot( n1 );
            // Ray is parallel to the triangle's plane
            if ( math::abs( n_dot_dir ) < 1e-6 )
            {
                return std::make_pair( false, 0 );
            }
            double d = -p1.dot( n1 );
            double t_line = -( origin.dot( n1 ) + d ) / n_dot_dir;
            if ( t_line <= 1e-10 ) // intersection not in the same direction as the ray
                return std::make_pair( false, 0 );
            // intersection point
            auto w = origin + direction * t_line;

            Eigen::Matrix<double, 3, 3> m;
            m.col( 0 ) = p2 - p1;
            m.col( 1 ) = p3 - p1;
            m.col( 2 ) = n1;
            auto w_ = m.inverse() * ( w - p1 );

            return std::make_pair( ( w_( 0 ) > 2 * FLT_MIN ) && ( w_( 1 ) > 0 ) && ( w_( 0 ) + w_( 1 ) < 1 ), t_line );
        }

        std::pair<bool, double> checkLeafIntersection( ray_type const& rayon, std::vector<primitiveinfo_type> const& primitiveInfo )
        {
            if constexpr ( nRealDim == 2 )
                return checkIntersectionWithSegment( rayon, primitiveInfo );
            else if constexpr ( nRealDim == 3 )
                return checkIntersectionWithTriangle( rayon, primitiveInfo );
        }

      private:
        BVHNode* setChild( uint16_type k, std::unique_ptr<BVHNode>&& childNode )
        {
            if ( childNode->M_parent )
            { /*TODO remove child in this parent*/
            }

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
        std::array<std::unique_ptr<BVHNode>, 2> M_children;
        BVHNode* M_parent = nullptr;
        int M_splitAxis = 0, M_nPrimitives = 0, M_firstPrimOffset = 0;
        vector_realdim_type M_bounds_min, M_bounds_max;
    };

    BVH_InHouse( worldcomm_ptr_t worldComm )
        : super_type( BVHEnum::Quality::High, worldComm ) {}

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
    // Value indicating whether we are in GPU mode
    bool isGPUHip() override
    {
        return ( false );
    }

    std::vector<std::vector<rayintersection_result_type>> intersectAllRaysWithGPU( std::vector<ray_type> const& rayons ) override
    {
        CHECK( false ) << "no implementation";
        std::vector<std::vector<rayintersection_result_type>> res;
        return ( res );
    }

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
        if ( M_rootNode->checkIntersection( rayon ) )
        {
            traverse_stackless( M_rootNode.get(), rayon );
        }
        if ( !M_intersected_leaf.empty() )
        {
            int argmin_lengths = std::distance( M_lengths.begin(), std::min_element( M_lengths.begin(), M_lengths.end() ) );
            res.push_back( rayintersection_result_type( this->worldComm().rank(), M_intersected_leaf[argmin_lengths], M_lengths[argmin_lengths] ) );
        }
        return res;
    }

    void buildTree()
    {
        if ( M_rootNode )
            return;

        M_rootNode = std::make_unique<BVHNode>();

        std::stack<std::tuple<BVHNode*, int, int, int>> stack;
        stack.push( std::make_tuple( M_rootNode.get(), 0, 0, this->M_primitiveInfo.size() ) );
        // TODO case only one 1 element
        while ( !stack.empty() )
        {
            auto [currentNode, cut_dimension, start_index_primitive, end_index_primitive] = stack.top();
            stack.pop();

            int nPrimitives = end_index_primitive - start_index_primitive;
            auto [bound_min_node, bound_max_node] = nPrimitives > 0 ? this->bounds( start_index_primitive, end_index_primitive ) : std::make_tuple( vector_realdim_type{}, vector_realdim_type{} );

            if ( nPrimitives <= 1 )
            {
                // Create a leaf, since there is only one primitive in the list
                int firstPrimOffset = M_orderedPrims.size();
                for ( int i = start_index_primitive; i < end_index_primitive; ++i )
                {
                    int primNum = this->M_primitiveInfo[i].meshEntity().id();
                    M_orderedPrims.push_back( primNum );
                }
                currentNode->updateForUse( firstPrimOffset, nPrimitives, -1, bound_min_node, bound_max_node );
            }
            else
            {
                CHECK( start_index_primitive >= 0 && end_index_primitive <= this->M_primitiveInfo.size() ) << start_index_primitive << " " << end_index_primitive;
                auto mid = ( start_index_primitive + end_index_primitive ) / 2;
                std::nth_element( &this->M_primitiveInfo[start_index_primitive], &this->M_primitiveInfo[mid],
                                  &this->M_primitiveInfo[end_index_primitive - 1] + 1,
                                  [cut_dimension = cut_dimension]( primitiveinfo_type const& a, primitiveinfo_type const& b )
                                  {
                                      return a.centroid()[cut_dimension] < b.centroid()[cut_dimension];
                                  } );

                int next_cut_dimension = ( cut_dimension + 1 ) % nRealDim;
                auto childNode0 = currentNode->setChild( 0, std::make_unique<BVHNode>() );
                stack.push( std::make_tuple( childNode0, next_cut_dimension, start_index_primitive, mid ) );
                auto childNode1 = currentNode->setChild( 1, std::make_unique<BVHNode>() );
                stack.push( std::make_tuple( childNode1, next_cut_dimension, mid, end_index_primitive ) );

                currentNode->updateForUse( -1, nPrimitives, next_cut_dimension, bound_min_node, bound_max_node );
            }
        }
    }

    std::tuple<vector_realdim_type, vector_realdim_type> bounds( int start_index_primitive, int end_index_primitive ) const
    {
        if ( start_index_primitive >= end_index_primitive )
            throw std::logic_error( "Error in BVHNode : compute bounds with no elemnent" );

        // vector_realdim_type newBoundsMin, newBoundsMax;
        vector_realdim_type newBoundsMin = this->M_primitiveInfo[start_index_primitive].boundMin();
        vector_realdim_type newBoundsMax = this->M_primitiveInfo[start_index_primitive].boundMax();
        for ( int i = start_index_primitive + 1; i < end_index_primitive; ++i )
        {
            auto const& primitiveInfo = this->M_primitiveInfo[i];
            for ( uint8_type d = 0; d < vector_realdim_type::SizeAtCompileTime; ++d )
            {
                newBoundsMin[d] = std::min( newBoundsMin[d], primitiveInfo.boundMin()[d] );
                newBoundsMax[d] = std::max( newBoundsMax[d], primitiveInfo.boundMax()[d] );
            }
        }
        return std::make_tuple( std::move( newBoundsMin ), std::move( newBoundsMax ) );
    }

    void traverse_stackless( BVH_InHouse::BVHNode* tree, ray_type const& rayon )
    {
        auto current_node = M_rootNode->nearChild( rayon );
        if ( !current_node ) // case where root is leaf
        {
            auto [has_intersected_leaf, distance] = M_rootNode->checkLeafIntersection( rayon, this->M_primitiveInfo );
            if ( has_intersected_leaf )
            {
                M_intersected_leaf.push_back( M_rootNode->firstPrimOffset() );
                M_lengths.push_back( distance );
            }
            return;
        }
        char state = 'P'; // the current node is being traversed from its Parent ('P')

        while ( true )
        {
            switch ( state )
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
                    auto [has_intersected_leaf, distance] = current_node->checkLeafIntersection( rayon, this->M_primitiveInfo );
                    if ( has_intersected_leaf )
                    {
                        // if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), this->M_primitiveInfo[current_node->firstPrimOffset()].meshEntity().id()) == M_intersected_leaf.end() )
                        if ( std::find( M_intersected_leaf.begin(), M_intersected_leaf.end(), current_node->firstPrimOffset() ) == M_intersected_leaf.end() )
                        {
                            // M_intersected_leaf.push_back(this->M_primitiveInfo[current_node->firstPrimOffset()].meshEntity().id());
                            M_intersected_leaf.push_back( current_node->firstPrimOffset() );
                            M_lengths.push_back( distance );
                        }
                    }
                    current_node = current_node->parent();
                    state = 'C'; // the current node is being accessed from its child
                }
                else
                {
                    current_node = current_node->nearChild( rayon );
                    state = 'P'; // the current node has been accessed from its parent
                }
                break;

            case 'P':
                if ( current_node->checkIntersection( rayon ) == false )
                {
                    current_node = current_node->siblingNode();
                    state = 'S'; // the current node has been accessed from its sibling
                }
                else if ( current_node->isLeaf() )
                {
                    auto [has_intersected_leaf, distance] = current_node->checkLeafIntersection( rayon, this->M_primitiveInfo );
                    if ( has_intersected_leaf )
                    {
                        // if ( std::find(M_intersected_leaf.begin(), M_intersected_leaf.end(), this->M_primitiveInfo[current_node->firstPrimOffset()].meshEntity().id()) == M_intersected_leaf.end() )
                        if ( std::find( M_intersected_leaf.begin(), M_intersected_leaf.end(), current_node->firstPrimOffset() ) == M_intersected_leaf.end() )
                        {
                            // M_intersected_leaf.push_back(this->M_primitiveInfo[current_node->firstPrimOffset()].meshEntity().id());
                            M_intersected_leaf.push_back( current_node->firstPrimOffset() );
                            M_lengths.push_back( distance );
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

                LOG( ERROR ) << "ERROR: None of the previous cases has been traversed";

                throw std::logic_error( "Error in BVH traversal: none of the previous cases has been traversed." );

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

template <typename... Ts>
auto boundingVolumeHierarchy( Ts&&... v )
{

    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto&& range = args.get( _range );
    using mesh_entity_type = std::remove_const_t<entity_range_t<std::decay_t<decltype( range )>>>;
    std::string const& kind = args.get_else( _kind, mesh_entity_type::nRealDim == 3 ? "third-party" : "in-house" );
    BVHEnum::Quality quality = args.get_else( _quality, BVHEnum::Quality::High );
    worldcomm_ptr_t worldcomm = args.get_else( _worldcomm, Environment::worldCommPtr() ); // TODO : use default worldcomm from range

    using bvh_type = BVH<mesh_entity_type>;
    std::unique_ptr<bvh_type> bvh;

    if ( kind == "in-house" )
    {
        using bvh_inhouse_type = BVH_InHouse<mesh_entity_type>;
        auto bvhInHouse = std::make_unique<bvh_inhouse_type>( worldcomm );
        bvhInHouse->updateForUse( range );
        bvh = std::move( bvhInHouse );
    }

    else

        if ( kind == "third-party" )
    {
        if constexpr ( mesh_entity_type::nRealDim != 3 )
            throw std::invalid_argument( "third-party only implement with triangle in 3D" );
        auto bvhThirdParty = std::make_unique<BVH_ThirdParty<mesh_entity_type>>( quality, worldcomm );
        bvhThirdParty->updateForUse( range );
        bvh = std::move( bvhThirdParty );
    }

#ifdef FEELPP_HAS_HIP
    else if ( kind == "hip-party" )
    {
        if constexpr ( mesh_entity_type::nRealDim != 3 )
            throw std::invalid_argument( "hip-party only implement with triangle in 3D" );
        auto bvhHIPParty = std::make_unique<BVH_HIP_Party<mesh_entity_type>>( quality, worldcomm );
        bvhHIPParty->updateForUse( range );
        bvh = std::move( bvhHIPParty );
    }
    else if ( kind == "hip-multi-gpu-party" )
    {
        printf( "<<<<<<<<<<> hip-multi-gpu-party <>>>>>>>>>>>>>>\n" );
        if constexpr ( mesh_entity_type::nRealDim != 3 )
            throw std::invalid_argument( "hip-multi-gpu-party only implement with triangle in 3D" );
        auto bvhHIPMultiGPUParty = std::make_unique<BVH_HIP_CPU_GPUs_Party<mesh_entity_type>>( quality, worldcomm );
        bvhHIPMultiGPUParty->updateForUse( range );
        bvh = std::move( bvhHIPMultiGPUParty );
    }
#endif

    else
        throw std::invalid_argument( fmt::format( "invalid bvh arg kind {} (should be third-party or in-house)", kind ) );

    return bvh;
}

} // namespace Feel
