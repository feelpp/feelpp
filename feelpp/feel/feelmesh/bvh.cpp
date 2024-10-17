
#include <feel/feelmesh/bvh.hpp>

/*
__global__ void kernel(float* x,float* y,int n){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) y[idx] += 1;
}
*/


namespace bvhhip
{
// Definition of manipulation tools, modifications for the creation of the BVH and doing Ray Tracing on the GPU

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
    int id;
};



struct F3TriangleInit {
    float3 _v0, _v1, _v2;
    int _id;
    
    F3TriangleInit(float3 v0, float3 v1,float3 v2,int id) 
        : _v0(v0), _v1(v1), _v2(v2),_id(id) {}

    __host__ __device__
    F3Triangle operator()(int idx) {
        F3Triangle t;
        t.v0 = _v0;
        t.v1 = _v1;
        t.v2 = _v2;
        t.id = _id;
        return t;
    }
};



struct Triangle {
    Vec3 v0, v1, v2;
    int id;
};



struct Box {
    Vec3 min;
    Vec3 max;
    int id;
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
    float tMin;
    float tMax;
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



// Function to calculate the AABB of a triangle
struct AABB {
    float3 min, max;
};



// Structure to represent a BVH node AABB
struct BVHNodeAABB {
    AABB bounds;
    int leftChild;
    int rightChild;
    int firstTriangle;
    int triangleCount;
};



__device__ float3 make_float3_device(float x, float y, float z)
{
    return make_float3(x, y, z);
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


__device__ __forceinline__ 
float3 make_float3_fast(float x, float y, float z) {
    //return make_float3(__float2int_rn(x), __float2int_rn(y), __float2int_rn(z));
    return make_float3(x,y,z); // A VOIR
}


float3 toFloat3(const Vec3& v) { return {v.x, v.y, v.z}; }

std::vector<F3Triangle> boxToTriangles(const Box& box,const int& id) {
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
        tri.id = id;
        triangles.push_back(tri);
        triangles.back().id=id;
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
// Nota : The problem is that it takes a lot of time and this version cannot be used for HPC.
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

// Function to build a BVH
// Nota : It is faster, I use GPU parallelism
void buildBVHWithTriangleVersion2(thrust::device_vector<F3Triangle>& triangles, thrust::device_vector<BVHNode>& nodes)
{
    int numTriangles = triangles.size();
    nodes.resize(2 * numTriangles - 1);

    // Initialize the leaves in parallel
    thrust::for_each(thrust::device, 
                     thrust::make_counting_iterator(0),
                     thrust::make_counting_iterator(numTriangles),
                     [triangles = thrust::raw_pointer_cast(triangles.data()),
                      nodes = thrust::raw_pointer_cast(nodes.data()),
                      numTriangles] __device__ (int i) {
        BVHNode& node = nodes[numTriangles - 1 + i];
        calculateBoundingBox(triangles[i], node.min, node.max);
        node.triangleIndex = i;
        node.leftChild = node.rightChild = -1;
    });

    for (int i = numTriangles - 2; i >= 0; --i) {
        BVHNode* raw_ptr = thrust::raw_pointer_cast(nodes.data());
        BVHNode& node = raw_ptr[i];
        int leftChild = 2 * i + 1;
        int rightChild = 2 * i + 2;

        node.leftChild = leftChild;
        node.rightChild = rightChild;
        node.triangleIndex = -1;

        const BVHNode& leftNode = raw_ptr[leftChild];
        const BVHNode& rightNode = raw_ptr[rightChild];

        node.min = make_float3(fminf(leftNode.min.x, rightNode.min.x),
                               fminf(leftNode.min.y, rightNode.min.y),
                               fminf(leftNode.min.z, rightNode.min.z));

        node.max = make_float3(fmaxf(leftNode.max.x, rightNode.max.x),
                               fmaxf(leftNode.max.y, rightNode.max.y),
                               fmaxf(leftNode.max.z, rightNode.max.z));
    }
}



__global__ void initializeLeaves3(F3Triangle* triangles, BVHNode* nodes, int numTriangles) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < numTriangles) {
        BVHNode& node = nodes[numTriangles - 1 + i];
        node.min = make_float3(INFINITY, INFINITY, INFINITY);
        node.max = make_float3(-INFINITY, -INFINITY, -INFINITY);
        calculateBoundingBox(triangles[i], node.min, node.max);
        
        node.triangleIndex = i;
        node.leftChild = node.rightChild = -1;
        
    }
    //__syncthreads();
}

__global__ void buildInternalNodes3(BVHNode* nodes, int numTriangles) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx <= numTriangles - 2) 
    {
        int nodeIdx = numTriangles - 2 - idx;
        //int nodeIdx = idx;
        BVHNode& node = nodes[nodeIdx];
        node.min = make_float3(INFINITY, INFINITY, INFINITY);
        node.max = make_float3(-INFINITY, -INFINITY, -INFINITY);

        int leftChild = 2 * nodeIdx + 1;
        int rightChild = 2 * nodeIdx + 2;

        node.leftChild = leftChild;
        node.rightChild = rightChild;
        node.triangleIndex = -1;

        const BVHNode& leftNode  = nodes[leftChild];
        const BVHNode& rightNode = nodes[rightChild];

        node.min = make_float3(fminf(leftNode.min.x, rightNode.min.x),
                               fminf(leftNode.min.y, rightNode.min.y),
                               fminf(leftNode.min.z, rightNode.min.z));

        node.max = make_float3(fmaxf(leftNode.max.x, rightNode.max.x),
                               fmaxf(leftNode.max.y, rightNode.max.y),
                               fmaxf(leftNode.max.z, rightNode.max.z)); 
    }
    __syncthreads();
}


void buildBVHWithTriangleVersion3(thrust::device_vector<F3Triangle>& triangles, thrust::device_vector<BVHNode>& nodes)
{
    int numTriangles = triangles.size();
    nodes.resize(2 * numTriangles - 1);
    // Initialize leaves
    int blockSize = 512;
    int numBlocks = (numTriangles + blockSize - 1) / blockSize;
    initializeLeaves3<<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(triangles.data()),
                                               thrust::raw_pointer_cast(nodes.data()),
                                               numTriangles);
    buildInternalNodes3<<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(nodes.data()),numTriangles);
    hipDeviceSynchronize();
     
}

// Function to calculate the height of the tree
int calculateTreeHeight(int numTriangles) {
    return static_cast<int>(std::ceil(std::log2(numTriangles)));
}

// Function to get the start index for a given level
int getStartIndexForLevel(int level, int numTriangles) {
    return std::max(0, (1 << level) - 1);
}

// Function to get the end index for a given level
int getEndIndexForLevel(int level, int numTriangles) {
    return std::min((1 << (level + 1)) - 1, numTriangles - 1);
}

// Construction of internal nodes
struct BuildInternalNodesFunctor {
    BVHNode* nodes;
    int numTriangles;

    __host__ __device__
    BuildInternalNodesFunctor(BVHNode* _nodes, int _numTriangles) 
        : nodes(_nodes), numTriangles(_numTriangles) {}

    __host__ __device__
    void operator()(int i) const {
        if (i >= numTriangles - 1) return;  

        BVHNode& node = nodes[i];
        int leftChild = 2 * i + 1;
        int rightChild = 2 * i + 2;

        node.leftChild = leftChild;
        node.rightChild = rightChild;
        node.triangleIndex = -1;

        const BVHNode& leftNode = nodes[leftChild];
        const BVHNode& rightNode = nodes[rightChild];

        node.min = make_float3(fminf(leftNode.min.x, rightNode.min.x),
                               fminf(leftNode.min.y, rightNode.min.y),
                               fminf(leftNode.min.z, rightNode.min.z));

        node.max = make_float3(fmaxf(leftNode.max.x, rightNode.max.x),
                               fmaxf(leftNode.max.y, rightNode.max.y),
                               fmaxf(leftNode.max.z, rightNode.max.z));
    }
};


// Function to build a BVH
// Nota : Even faster and more optimized
void buildBVHWithTriangleVersion5(thrust::device_vector<F3Triangle>& triangles, thrust::device_vector<BVHNode>& nodes)
{
    int numTriangles = triangles.size();
    nodes.resize(2 * numTriangles - 1);

    // Initialize the leaves in parallel
    thrust::for_each(thrust::device, 
                     thrust::make_counting_iterator(0),
                     thrust::make_counting_iterator(numTriangles),
                     [triangles = thrust::raw_pointer_cast(triangles.data()),
                      nodes = thrust::raw_pointer_cast(nodes.data()),
                      numTriangles] __device__ (int i) {
        BVHNode& node = nodes[numTriangles - 1 + i];
        calculateBoundingBox(triangles[i], node.min, node.max);
        node.triangleIndex = i;
        node.leftChild = node.rightChild = -1;
    });

    // Synchronize to ensure all leaves are initialized. Very important !!!
    hipDeviceSynchronize();

    // Build the internal nodes level by level
    int treeHeight = calculateTreeHeight(numTriangles);
    BVHNode* raw_ptr = thrust::raw_pointer_cast(nodes.data());

    for (int level = treeHeight - 1; level >= 0; --level) {
        int startIdx = getStartIndexForLevel(level, numTriangles);
        int endIdx = getEndIndexForLevel(level, numTriangles);

        thrust::for_each(
            thrust::device,
            thrust::make_counting_iterator(startIdx),
            thrust::make_counting_iterator(endIdx + 1),
            BuildInternalNodesFunctor(raw_ptr, numTriangles)
        );

        // Synchronize after each level. Very important !!!
        hipDeviceSynchronize();
    }
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

__device__ __forceinline__ bool rayTriangleIntersectVersion2(const F3Ray& ray, const F3Triangle& triangle, float& t, float3& intersectionPoint) {
 
    const float3 v0 = triangle.v0;
    const float3 edge1 = triangle.v1 - v0;
    const float3 edge2 = triangle.v2 - v0;
    const float3 pvec = cross(ray.direction, edge2);
    const float det = dot(edge1, pvec);
    if (fabsf(det) < 1e-6f) return false;
    const float invDet = __frcp_rn(det);
    const float3 tvec = ray.origin - v0;
    const float u = dot(tvec, pvec) * invDet;
    if (u < 0.0f || u > 1.0f) return false;
    const float3 qvec = cross(tvec, edge1);
    const float v = dot(ray.direction, qvec) * invDet;
    if (v < 0.0f || u + v > 1.0f) return false;
    t = dot(edge2, qvec) * invDet;
    if (t > 1e-6f) {
        intersectionPoint = ray.origin + t * ray.direction;
        return true;
    }
    intersectionPoint =make_float3(INFINITY, INFINITY, INFINITY);
    return false;
}


__global__ void rayTracingKernelVersion1(
    BVHNode* nodes, 
    F3Triangle* triangles, 
    F3Ray* rays, 
    int* hitResults, 
    float* distance,
    float3* intersectionPoint,
    int* hitId, 
    int numRays
) 

{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numRays) return;

    F3Ray ray = rays[idx];
    int stack[64];
    int stackPtr = 0;
    stack[stackPtr++] = 0;

    float closestT = INFINITY;
    int closestTriangle = -1;
    int closesIntersectionId = -1;
    float3 closestIntersectionPoint=make_float3(INFINITY, INFINITY, INFINITY);
    bool isView=false; //isView=true;

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

        int numIdNodeTrianggleIndex=node.triangleIndex;

        if (node.triangleIndex != -1) {
            // Sheet: test the intersection with the triangle
            float t;
            float3 intersectionPointT;
            if (rayTriangleIntersect(ray, triangles[node.triangleIndex], t,intersectionPointT)) {

                //To view all intersections
                if (isView) printf("      Node Idx [%i] Num Ray[%i] <%f %f %f>\n",nodeIdx,idx,intersectionPointT.x,intersectionPointT.y,intersectionPointT.z);

                if (t < closestT) {
                    closestT = t;
                    closestTriangle = node.triangleIndex;
                    closestIntersectionPoint=intersectionPointT;
                    closesIntersectionId=triangles[numIdNodeTrianggleIndex].id;
                    //printf("      NodeTriangleIndex=%i %i\n",numIdNodeTrianggleIndex,triangles[numIdNodeTrianggleIndex].id);
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
    hitId[idx]             = closesIntersectionId;
    //if (closestTriangle!=-1) { printf("t=%f\n",closestT); } OK
}



__global__ void rayTracingKernelVersion3(
    BVHNode* nodes, 
    F3Triangle* __restrict__ triangles, 
    F3Ray* rays, 
    int* __restrict__ hitResults, 
    float* __restrict__ distance,
    float3* __restrict__ intersectionPoint,
    int* __restrict__ hitId, 
    int numRays
) 

{
    __shared__ BVHNode sharedNodes[64]; 
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numRays) return;
    // On vharge les premiers nœuds BVH dans la mémoire partagée
    if (threadIdx.x < 64) {
        sharedNodes[threadIdx.x] = nodes[threadIdx.x];
    }
    __syncthreads();

    F3Ray ray = rays[idx];
    constexpr int MAX_STACK_SIZE = 64;
    int stack[MAX_STACK_SIZE];
    int stackPtr = 0;
    stack[stackPtr++] = 0;

    float closestT = INFINITY;
    int closestTriangle = -1;
    int closesIntersectionId = -1;
    float3 closestIntersectionPoint=make_float3_fast(INFINITY, INFINITY, INFINITY);
    const float3 invDir = make_float3_fast(__frcp_rn(ray.direction.x), __frcp_rn(ray.direction.y), __frcp_rn(ray.direction.z));
    bool isView=false; //isView=true;


    while (stackPtr > 0) {
        const int nodeIdx = stack[--stackPtr];
        const BVHNode& node = nodeIdx < 64 ? sharedNodes[nodeIdx] : nodes[nodeIdx];

        float tmin = fmaxf(fminf((node.min.x - ray.origin.x) * invDir.x, 
                                 (node.max.x - ray.origin.x) * invDir.x),
                           fmaxf((node.min.y - ray.origin.y) * invDir.y, 
                                 (node.max.y - ray.origin.y) * invDir.y));
        float tmax = fminf(fmaxf((node.min.x - ray.origin.x) * invDir.x, 
                                 (node.max.x - ray.origin.x) * invDir.x),
                           fminf((node.min.y - ray.origin.y) * invDir.y, 
                                 (node.max.y - ray.origin.y) * invDir.y));
        
        tmin = fmaxf(tmin, fminf((node.min.z - ray.origin.z) * invDir.z, 
                                 (node.max.z - ray.origin.z) * invDir.z));
        tmax = fminf(tmax, fmaxf((node.min.z - ray.origin.z) * invDir.z, 
                                 (node.max.z - ray.origin.z) * invDir.z));

        if (tmax < 0 || tmin > tmax || tmin > closestT) continue;

        if (node.triangleIndex != -1) {
            float t;
            float3 intersectionPointT;
            if (rayTriangleIntersectVersion2(ray, triangles[node.triangleIndex], t, intersectionPointT)) 
            {
                //To view all intersections
                if (isView) printf("      Node Idx [%i] Num Ray[%i] <%f %f %f>\n",nodeIdx,idx,intersectionPointT.x,intersectionPointT.y,intersectionPointT.z);
                if (t < closestT) {
                    closestT = t;
                    closestTriangle = node.triangleIndex;
                    closestIntersectionPoint = intersectionPointT;
                    closestIntersectionId = triangles[node.triangleIndex].id;
                }
            }
        } else if (stackPtr < MAX_STACK_SIZE - 1) {
            stack[stackPtr++] = node.leftChild;
            stack[stackPtr++] = node.rightChild;
        }
    }


 

    hitResults[idx]        = closestTriangle;
    distance[idx]          = closestT;
    intersectionPoint[idx] = closestIntersectionPoint;
    hitId[idx]             = closesIntersectionId;
    //if (closestTriangle!=-1) { printf("t=%f\n",closestT); } OK
}


void writeBVHNodes(const thrust::device_vector<BVHNode>& nodes) {
    std::vector<BVHNode> hostNodes(nodes.size());
    thrust::copy(nodes.begin(), nodes.end(), hostNodes.begin());

    std::cout << "BVH Nodes:" << std::endl;
    for (size_t i = 0; i < hostNodes.size(); ++i) {
        const BVHNode& node = hostNodes[i];
        std::cout << "Node " << i << ":" << std::endl;
        std::cout << "  Min: (" << node.min.x << ", " << node.min.y << ", " << node.min.z << ")" << std::endl;
        std::cout << "  Max: (" << node.max.x << ", " << node.max.y << ", " << node.max.z << ")" << std::endl;
        std::cout << "  Left Child: " << node.leftChild << std::endl;
        std::cout << "  Right Child: " << node.rightChild << std::endl;
        std::cout << "  Triangle Index: " << node.triangleIndex << std::endl;
        std::cout << std::endl;
    }
}

// Function that loads the BVH node
void loadBVH(const std::string& filename, thrust::device_vector<BVHNode>& nodes) {
    std::ifstream inFile(filename, std::ios::binary);
    
    // Returns an error message if the file does not exist.
    if (!inFile) {
        throw std::runtime_error("Could not open file for reading: " + filename);
    }

    int nodeCount;
    // Read the number of nodes
    inFile.read(reinterpret_cast<char*>(&nodeCount), sizeof(int));

    // Resize the device vector to hold the nodes
    nodes.resize(nodeCount);

    // Read each node
    for (int i = 0; i < nodeCount; ++i) {
        BVHNode* raw_ptr = thrust::raw_pointer_cast(nodes.data());
        BVHNode& node = raw_ptr[i];
        inFile.read(reinterpret_cast<char*>(&node.min), sizeof(float3));
        inFile.read(reinterpret_cast<char*>(&node.max), sizeof(float3));
        inFile.read(reinterpret_cast<char*>(&node.leftChild), sizeof(int));
        inFile.read(reinterpret_cast<char*>(&node.rightChild), sizeof(int));
        inFile.read(reinterpret_cast<char*>(&node.triangleIndex), sizeof(int));
        inFile.read(reinterpret_cast<char*>(&node.boxIndex), sizeof(int)); 
    }

    inFile.close();
}

// Function that saves the BVH node
void saveBVH(const std::string& filename, const thrust::device_vector<BVHNode>& nodes) {
    std::ofstream outFile(filename, std::ios::binary);
    
    // Returns an error message if the file does not exist.
    if (!outFile) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    int nodeCount = nodes.size();
    // Write the number of nodes
    outFile.write(reinterpret_cast<const char*>(&nodeCount), sizeof(int));

    // Write each node
    for (int i = 0; i < nodeCount; ++i) {
        const BVHNode& node = nodes[i];
        outFile.write(reinterpret_cast<const char*>(&node.min), sizeof(float3));
        outFile.write(reinterpret_cast<const char*>(&node.max), sizeof(float3));
        outFile.write(reinterpret_cast<const char*>(&node.leftChild), sizeof(int));
        outFile.write(reinterpret_cast<const char*>(&node.rightChild), sizeof(int));
        outFile.write(reinterpret_cast<const char*>(&node.triangleIndex), sizeof(int));
        outFile.write(reinterpret_cast<const char*>(&node.boxIndex), sizeof(int)); 
    }

    outFile.close();
}


struct Camera {
    float3 position;  // Camera position
    float3 target;    // Point the camera is looking at
    float3 up;        // Up vector for the camera
    float fov;        // Field of view in radians
    float aspect;     // Aspect ratio (width/height)

    __device__ void getRay(int x, int y, int width, int height, float3* rayOrigin, float3* rayDirection) const {
        // Calculate normalized device coordinates (NDC)
        float ndcX = (2.0f * (x + 0.5f) / width - 1.0f) * aspect;
        float ndcY = 1.0f - 2.0f * (y + 0.5f) / height;

        // Calculate the direction of the ray in world space
        float3 forward = make_float3(target.x - position.x, target.y - position.y, target.z - position.z);
        forward = normalize(forward);
        
        float3 right = cross(forward, up);
        right = normalize(right);
        
        float3 cameraUp = cross(right, forward);

        // Calculate the ray direction
        float3 horizontal = right * tan(fov / 2.0f);
        float3 vertical = cameraUp * tan(fov / 2.0f);

        *rayDirection = forward + horizontal * ndcX + vertical * ndcY;
        *rayDirection = normalize(*rayDirection);
        
        *rayOrigin = position;
    }
};


__global__ void rayTracingImgKernel(
    unsigned char* image,
    int width,
    int height,
    Camera camera,
    BVHNode* nodes, 
    F3Triangle* triangles, 
    int* hitResults, 
    float* distance,
    float3* intersectionPoint,
    int* hitId 
)
 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= width * height) return;


    int x = idx % width;
    int y = idx / width;

    F3Ray ray;
    camera.getRay(x,y,width,height,&ray.origin,&ray.direction);

    int stack[64];
    int stackPtr = 0;
    stack[stackPtr++] = 0;

    float closestT = INFINITY;
    int closestTriangle = -1;
    int closesIntersectionId = -1;
    float3 closestIntersectionPoint=make_float3(INFINITY, INFINITY, INFINITY);
    bool isView=false; //isView=true;

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

        int numIdNodeTrianggleIndex=node.triangleIndex;

        if (node.triangleIndex != -1) {
            // Sheet: test the intersection with the triangle
            float t;
            float3 intersectionPointT;
            if (rayTriangleIntersect(ray, triangles[node.triangleIndex], t,intersectionPointT)) {

                //To view all intersections
                //if (isView) printf("      Node Idx [%i] Num Ray[%i] <%f %f %f>\n",nodeIdx,idx,intersectionPointT.x,intersectionPointT.y,intersectionPointT.z);
                if (isView) printf("      Num Ray[%i] <%f %f %f>\n",idx,intersectionPointT.x,intersectionPointT.y,intersectionPointT.z);

                if (t < closestT) {
                    closestT = t;
                    closestTriangle = node.triangleIndex;
                    closestIntersectionPoint=intersectionPointT;
                    closesIntersectionId=triangles[numIdNodeTrianggleIndex].id;
                    //printf("      NodeTriangleIndex=%i %i\n",numIdNodeTrianggleIndex,triangles[numIdNodeTrianggleIndex].id);
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
    hitId[idx]             = closesIntersectionId;

    //if (closestTriangle!=-1) { printf("t=%f\n",closestT); } OK

    if (hitId[idx]!=-1)
    {
           int value= 255;
           image[(y * width + x) * 3]     = value;
           image[(y * width + x) * 3 + 1] = value;
           image[(y * width + x) * 3 + 2] = value;
    }
    else
    {
           image[(y * width + x) * 3]     = 0; // Red chanel
           image[(y * width + x) * 3 + 1] = 0; // Green chanel
           image[(y * width + x) * 3 + 2] = 0; // Blue chanel
     
    }

}

void savePPM(const std::string& filename, unsigned char* data, int width, int height) {
    std::ofstream file(filename, std::ios::binary);
    file << "P6\n" << width << " " << height << "\n255\n";
    file.write(reinterpret_cast<char*>(data), width * height * 3);
}

void buildPicturRayTracingPPM(
	thrust::device_vector<F3Triangle>& triangles, 
	thrust::device_vector<BVHNode>& nodes,
	Camera camera,
    int width,
    int height,
	const std::string& filename)
{
    // Before using this function, the BVH must already be calculated and the triangles must be in the device.
	// Ray Tracing
    const int threadsPerBlock = 256;
    const int numRays = width * height; // Total number of rays based on image dimensions
    int blocksPerGrid = (numRays + threadsPerBlock - 1) / threadsPerBlock;
    
    //...
    thrust::device_vector<unsigned char> deviceImage(width * height * 3);
    thrust::device_vector<int>           deviceHitResults(numRays);
    thrust::device_vector<float>         deviceDistanceResults(numRays);
    thrust::device_vector<float3>        deviceIntersectionPoint(numRays);
    thrust::device_vector<int>           deviceIdResults(numRays);

    //...
    rayTracingImgKernel<<<blocksPerGrid, threadsPerBlock>>>(
        thrust::raw_pointer_cast(deviceImage.data()),
        width,
        height,
        camera,
        thrust::raw_pointer_cast(nodes.data()),
        thrust::raw_pointer_cast(triangles.data()),
        thrust::raw_pointer_cast(deviceHitResults.data()),
        thrust::raw_pointer_cast(deviceDistanceResults.data()),
        thrust::raw_pointer_cast(deviceIntersectionPoint.data()),
        thrust::raw_pointer_cast(deviceIdResults.data())
    );

    
    //...
    thrust::host_vector<unsigned char> hostImage = deviceImage;
    savePPM(filename, hostImage.data(), width, height);

	// Memory cleaning
    deviceHitResults.clear(); 
    deviceDistanceResults.clear(); 
    deviceIntersectionPoint.clear(); 
    deviceIdResults.clear(); 
    deviceImage.clear(); 
}


// Function to calculate the AABB of a triangle
__device__ AABB calculateAABB(const F3Triangle& triangle) {
    AABB aabb;
    aabb.min = make_float3(
        fminf(fminf(triangle.v0.x, triangle.v1.x), triangle.v2.x),
        fminf(fminf(triangle.v0.y, triangle.v1.y), triangle.v2.y),
        fminf(fminf(triangle.v0.z, triangle.v1.z), triangle.v2.z)
    );
    aabb.max = make_float3(
        fmaxf(fmaxf(triangle.v0.x, triangle.v1.x), triangle.v2.x),
        fmaxf(fmaxf(triangle.v0.y, triangle.v1.y), triangle.v2.y),
        fmaxf(fmaxf(triangle.v0.z, triangle.v1.z), triangle.v2.z)
    );
    return aabb;
}

struct CalculateAABB {
    __host__ __device__
    AABB operator()(const F3Triangle& triangle) const {
        AABB aabb;
        aabb.min = make_float3(
            fminf(fminf(triangle.v0.x, triangle.v1.x), triangle.v2.x),
            fminf(fminf(triangle.v0.y, triangle.v1.y), triangle.v2.y),
            fminf(fminf(triangle.v0.z, triangle.v1.z), triangle.v2.z)
        );
        aabb.max = make_float3(
            fmaxf(fmaxf(triangle.v0.x, triangle.v1.x), triangle.v2.x),
            fmaxf(fmaxf(triangle.v0.y, triangle.v1.y), triangle.v2.y),
            fmaxf(fmaxf(triangle.v0.z, triangle.v1.z), triangle.v2.z)
        );
    return aabb;
    }
};

// Function to merge two AABBs
__device__ AABB mergeAABB(const AABB& a, const AABB& b) {
    AABB result;
    result.min = make_float3(
        fminf(a.min.x, b.min.x),
        fminf(a.min.y, b.min.y),
        fminf(a.min.z, b.min.z)
    );
    result.max = make_float3(
        fmaxf(a.max.x, b.max.x),
        fmaxf(a.max.y, b.max.y),
        fmaxf(a.max.z, b.max.z)
    );
    return result;
}

__device__ bool intersectAABB(const F3Ray& ray, const AABB& aabb) {
    float3 invDir = make_float3(1.0f / ray.direction.x, 1.0f / ray.direction.y, 1.0f / ray.direction.z);
    float3 t0 = make_float3((aabb.min.x - ray.origin.x) * invDir.x,
                            (aabb.min.y - ray.origin.y) * invDir.y,
                            (aabb.min.z - ray.origin.z) * invDir.z);
    float3 t1 = make_float3((aabb.max.x - ray.origin.x) * invDir.x,
                            (aabb.max.y - ray.origin.y) * invDir.y,
                            (aabb.max.z - ray.origin.z) * invDir.z);
    float tmin = fmaxf(fmaxf(fminf(t0.x, t1.x), fminf(t0.y, t1.y)), fminf(t0.z, t1.z));
    float tmax = fminf(fminf(fmaxf(t0.x, t1.x), fmaxf(t0.y, t1.y)), fmaxf(t0.z, t1.z));
    return tmax >= tmin && tmin < ray.tMax && tmax > ray.tMin;
}


__device__ bool intersectTriangleVersion2(const F3Ray& ray, const F3Triangle& triangle, float& t, float3& intersectionPoint) {
        float3 edge1 = triangle.v1 - triangle.v0;
        float3 edge2 = triangle.v2 - triangle.v0;
        float3 h = cross(ray.direction, edge2);
        float a = dot(edge1, h);
        
        if (a > -1e-6f && a < 1e-6f) return false;
        float f = 1.0f / a;
        float3 s = ray.origin - triangle.v0;
        float u = f * dot(s, h);
        
        if (u < 0.0f || u > 1.0f) return false;
        
        float3 q = cross(s, edge1);
        float v = f * dot(ray.direction, q);
        
        if (v < 0.0f || u + v > 1.0f) return false;
        
        t = f * dot(edge2, q);
        
        if (t > 1e-6f) {
            intersectionPoint = ray.origin + ray.direction * t;
        }
    return t > ray.tMin && t < ray.tMax;
}


__device__ float3 calculateCentroid(const F3Triangle& triangle) {
    return make_float3(
        (triangle.v0.x + triangle.v1.x + triangle.v2.x) / 3.0f,
        (triangle.v0.y + triangle.v1.y + triangle.v2.y) / 3.0f,
        (triangle.v0.z + triangle.v1.z + triangle.v2.z) / 3.0f
    );
}


struct MergeAABB {
    __host__ __device__
    AABB operator()(const AABB& a, const AABB& b) const {
        AABB result;
        result.min.x = fminf(a.min.x, b.min.x);
        result.min.y = fminf(a.min.y, b.min.y);
        result.min.z = fminf(a.min.z, b.min.z);
        result.max.x = fmaxf(a.max.x, b.max.x);
        result.max.y = fmaxf(a.max.y, b.max.y);
        result.max.z = fmaxf(a.max.z, b.max.z);
        return result;
    }

    __host__ __device__
    AABB identity() const {
        AABB result;
        result.min = make_float3( FLT_MAX,  FLT_MAX,  FLT_MAX);
        result.max = make_float3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
        return result;
    }
};

struct CalculateCentroid {
    __host__ __device__
    float3 operator()(const F3Triangle& triangle) const {
        float3 centroid;
        centroid.x = (triangle.v0.x + triangle.v1.x + triangle.v2.x) / 3.0f;
        centroid.y = (triangle.v0.y + triangle.v1.y + triangle.v2.y) / 3.0f;
        centroid.z = (triangle.v0.z + triangle.v1.z + triangle.v2.z) / 3.0f;
        return centroid;
    }
};

__device__ __host__ inline float getComponent(const float3& vec, int axis) {
    switch(axis) {
        case 0: return vec.x;
        case 1: return vec.y;
        case 2: return vec.z;
        default: return 0.0f; // or handle error
    }
}


// Recursive function to build the BVH
void buildBVH_AABB_Recursive(thrust::device_vector<F3Triangle>& triangles,
                       thrust::device_vector<AABB>& aabbs,
                       thrust::device_vector<float3>& centroids,
                       thrust::device_vector<BVHNodeAABB>& nodes,
                       int& nodeIndex,
                       int start,
                       int end)
{
    BVHNodeAABB* raw_ptr = thrust::raw_pointer_cast(nodes.data());
    BVHNodeAABB& node = raw_ptr[nodeIndex];
    node.firstTriangle = start;
    node.triangleCount = end - start;
    node.bounds = thrust::reduce(thrust::device, aabbs.begin() + start, aabbs.begin() + end, AABB(), MergeAABB());

    if (node.triangleCount <= 2) {
        // Leaf node
        node.leftChild = -1;
        node.rightChild = -1;
    } else {
        // Internal node
        int axis = 0;
        //float splitPos = (node.bounds.min[axis] + node.bounds.max[axis]) * 0.5f;
        float splitPos = 0.5f * (getComponent(node.bounds.min, axis) + getComponent(node.bounds.max, axis));

        // Partition the triangles
        auto splitIter = thrust::partition(thrust::device,
            thrust::make_zip_iterator(thrust::make_tuple(triangles.begin() + start, aabbs.begin() + start, centroids.begin() + start)),
            thrust::make_zip_iterator(thrust::make_tuple(triangles.begin() + end, aabbs.begin() + end, centroids.begin() + end)),
                [=] __device__ (const thrust::tuple<F3Triangle, AABB, float3>& t) {
                    return getComponent(thrust::get<2>(t), axis) < splitPos;
                }
        );
        
        int mid = start + thrust::distance(
                thrust::make_zip_iterator(thrust::make_tuple(triangles.begin() + start, aabbs.begin() + start, centroids.begin() + start)),
                splitIter
            );
            
        // Check if the partition actually divided the triangles
        if (mid == start || mid == end) {
                // If the partition did not divide the triangles, force a division in the middle
                mid = start + (end - start) / 2;
        }
        
        //std::cout<<"mid="<<mid<<"\n"; CTRL OK

        // Create the child nodes
        node.leftChild = ++nodeIndex;
        buildBVH_AABB_Recursive(triangles, aabbs, centroids, nodes, nodeIndex, start, mid);
        node.rightChild = ++nodeIndex;
        buildBVH_AABB_Recursive(triangles, aabbs, centroids, nodes, nodeIndex, mid, end);
    }
}


struct BVHBuildTask {
    int nodeIndex;
    int start;
    int end;
    BVHBuildTask(int ni, int s, int e) : nodeIndex(ni), start(s), end(e) {}
};


void buildBVH_AABB_Iterative(thrust::device_vector<F3Triangle>& triangles,
                       thrust::device_vector<AABB>& aabbs,
                       thrust::device_vector<float3>& centroids,
                       thrust::device_vector<BVHNodeAABB>& nodes,
                       int& nodeIndex,
                       int start,
                       int end) {
                        
    std::stack<BVHBuildTask> taskStack;
    taskStack.push(BVHBuildTask(nodeIndex, start, end));

    while (!taskStack.empty()) {
        BVHBuildTask task = taskStack.top();
        taskStack.pop();

        int currentNodeIndex = task.nodeIndex;
        int currentStart = task.start;
        int currentEnd = task.end;

    
        BVHNodeAABB* raw_ptr = thrust::raw_pointer_cast(nodes.data());
        BVHNodeAABB& node = raw_ptr[currentNodeIndex];

        node.firstTriangle = currentStart;
        node.triangleCount = currentEnd - currentStart;

        //std::cout<<"current nodeIndex="<<currentNodeIndex<<" Start= "<<currentStart<<" End="<<currentEnd <<"\n";
        //std::cout<<"node.triangleCount="<<node.triangleCount<<"\n";
        //getchar();

        // Calculates the AABB of the node
        node.bounds = thrust::reduce(thrust::device, aabbs.begin() + currentStart, aabbs.begin() + currentEnd, AABB(), MergeAABB());

        if (node.triangleCount <= 2) {
            // Leaf node
            node.leftChild = -1;
            node.rightChild = -1;
        } else {
            // Internal node
            int axis = 0; // // separation axis (can be optimized). To be seen later, if ...
            float splitPos = 0.5f * (getComponent(node.bounds.min, axis) + getComponent(node.bounds.max, axis));
            // Partition the triangles
            auto splitIter = thrust::partition(thrust::device,
                thrust::make_zip_iterator(thrust::make_tuple(triangles.begin() + currentStart, aabbs.begin() + currentStart, centroids.begin() + currentStart)),
                thrust::make_zip_iterator(thrust::make_tuple(triangles.begin() + currentEnd, aabbs.begin() + currentEnd, centroids.begin() + currentEnd)),
                [=] __device__ (const thrust::tuple<F3Triangle, AABB, float3>& t) {
                    return getComponent(thrust::get<2>(t), axis) < splitPos;
                }
            );

            int mid = currentStart + thrust::distance(
                thrust::make_zip_iterator(thrust::make_tuple(triangles.begin() + currentStart, aabbs.begin() + currentStart, centroids.begin() + currentStart)),
                splitIter
            );

            // Check if the partition actually divided the triangles
            if (mid == currentStart || mid == currentEnd) {
                // If the partition did not divide the triangles, force a division in the middle
                mid = currentStart + (currentEnd - currentStart) / 2;
            }

            //std::cout<<"mid="<<mid<<"\n";
            node.leftChild = ++nodeIndex;
            node.rightChild = ++nodeIndex;

            taskStack.push(BVHBuildTask(node.rightChild, mid, currentEnd));
            taskStack.push(BVHBuildTask(node.leftChild, currentStart, mid));
        }
    }
}




void buildBVH_AABB(thrust::device_vector<F3Triangle>& triangles, thrust::device_vector<BVHNodeAABB>& nodes) {
    int numTriangles = triangles.size();
    nodes.resize(2 * numTriangles - 1);
    //Calculate AABBs and centroids for all triangles
    thrust::device_vector<AABB>   aabbs(numTriangles);
    thrust::device_vector<float3> centroids(numTriangles);

    thrust::transform(thrust::device, triangles.begin(), triangles.end(), aabbs.begin(), CalculateAABB());
    thrust::transform(thrust::device, triangles.begin(), triangles.end(), centroids.begin(), CalculateCentroid());

    // Build the BVH recursively or iteratively
    int rootNodeIndex = 0;
    buildBVH_AABB_Recursive(triangles, aabbs, centroids, nodes, rootNodeIndex, 0, numTriangles); //Nota: it is a bit faster than compared to the iterative
    //buildBVH_AABB_Iterative(triangles, aabbs, centroids, nodes, rootNodeIndex, 0, numTriangles);
}


__global__ void rayTracingImgKernel_AABB(
    const BVHNodeAABB* nodes,
    const F3Triangle* triangles,
    const F3Ray* rays,
    int* hitResults,
    float* hitDistances,
    float3* intersectionPoint,
    int* hitId,
    int numRays
)

{
    int rayIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (rayIdx >= numRays) return;

    F3Ray ray = rays[rayIdx];
    int stack[64];
    int stackPtr = 0;
    stack[stackPtr++] = 0;

    //float closestHit = ray.tMax;
    float closestHit = INFINITY;
    int closestTriangle = -1;
    int closesIntersectionId = -1;
    float3 closestIntersectionPoint = make_float3(INFINITY, INFINITY, INFINITY);
    bool isView=false; //isView=true;

    while (stackPtr > 0) {
        int nodeIdx = stack[--stackPtr];
        const BVHNodeAABB& node = nodes[nodeIdx];

        if (intersectAABB(ray, node.bounds)) {
            //printf("nodeIdx=%i\n",nodeIdx);
            if (node.leftChild == -1 && node.rightChild == -1) {
                // Leaf node
                for (int i = 0; i < node.triangleCount; ++i) {
                    const F3Triangle& tri = triangles[node.firstTriangle + i];
                    float t;
                    float3 intersectionPointT;

                    if (intersectTriangleVersion2(ray, tri, t,intersectionPointT)) {
                        if (isView) printf("[%i] %f \n",rayIdx,t);
                        if (t < closestHit) {
                            closestHit = t;
                            closestTriangle = node.firstTriangle + i;
                            closestIntersectionPoint = intersectionPointT;
                            closesIntersectionId=triangles[closestTriangle].id;
                        }
                    }
                }
            } else {
                if (node.rightChild != -1) stack[stackPtr++] = node.rightChild;
                if (node.leftChild != -1) stack[stackPtr++] = node.leftChild;
            }
        }
    }

    hitResults[rayIdx]        = closestTriangle;
    hitDistances[rayIdx]      = fabs(closestHit);
    intersectionPoint[rayIdx] = closestIntersectionPoint;
    hitId[rayIdx]             = closesIntersectionId;
}


}





