//#include <future>

//#include <feel/feelcore/testsuite.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/bvh.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <unordered_map>
#include <functional>


#include "tools.hpp"
#include <feel/feeltask/taskpu.hpp>

#include <boost/filesystem.hpp>


namespace Feel {


template <typename MeshType>
class ShadingMask
{
    using mesh_type = MeshType;
    typedef typename MeshType::ptrtype mesh_ptrtype;

    using tr_mesh_type = typename std::conditional<MeshType::nDim==MeshType::nRealDim,
                                                //typename MeshType::trace_mesh_type,
                                                typename MeshType::template trace_mesh_type<0>,
                                                typename MeshType::type >::type ;
    using tr_mesh_ptrtype = typename std::conditional<MeshType::nDim==MeshType::nRealDim,
                                                //trace_mesh_t<MeshType>,//typename MeshType::trace_mesh_ptrtype,
                                                typename MeshType::template trace_mesh_ptrtype<0>,
                                                typename MeshType::ptrtype >::type;
/*
    using mesh_entity_type = typename std::conditional<MeshType::nDim==MeshType::nRealDim,
                                            entity_range_t<elements_reference_wrapper_t<MeshType>>,
                                            entity_range_t<faces_reference_wrapper_t<MeshType>> >::type;
                                            //entity_range_t<elements_reference_wrapper_t<MeshType>>,
                                            //entity_range_t<faces_reference_wrapper_t<MeshType>> >::type;
*/


    using wrapper_type = elements_reference_wrapper_t<MeshType>;
                                        
    typedef typename matrix_node<double>::type matrix_node_type;

public:
    using value_type = double;
    bool QSaveTypeThreadDotON;
    bool QSaveControlFiles;
    int  numTypeThread;
    int  numModeTaskUsed; 


    ShadingMask(mesh_ptrtype mesh, nl::json const& specs, int intervalsAzimuth=72, int intervalsAltitude=10 );

    auto commonComputePartCTRL(int NumOption,matrix_node_type const& element_points,int n_rays_thread, int id_thread);

    bool computePartMarker(std::vector<std::string> marker_list_thread, int id_thread, int start_index);



    // Create the random number generators
    void makeRandomNumberGeneratorsSeed(bool QCTRL_SAVE_SEED,bool QCTRL_LOAD_SEED);

     // Create and store the directions of the M_Nrays
    void makeCreateM_NraysMatrix(int intervalsAzimuth, int intervalsAltitude);

    // For each building, save the surface mesh and build the corresponding BVH tree for ray search
    void loadMeshData(mesh_ptrtype mesh, nl::json const& specs);
        void loadMeshDataSubPartList          (mesh_ptrtype mesh);
        void loadMeshDataSubPartVolumes       (mesh_ptrtype mesh,int numOp);
        void loadMeshDataSubPartSurfacesFaces (mesh_ptrtype mesh,int numOp);
        void loadMeshDataSubPartMarkers       (mesh_ptrtype mesh);

    // Subdivide the azimuth angles [0,360]° and altitude angles [0,90]° in subsets for easier computation of the shading masks
    void fixAzimuthAltitudeDiscretization(int intervalsAzimuth=72, int intervalsAltitude=10);

    // Choose a random pair of indices in the discretized azimuth and altitude vectors
    void getRandomDirectionSM(std::vector<double> &random_direction, std::mt19937 & M_gen, std::mt19937 & M_gen2, int& index_azimuth, int& index_altitude);

    // Get a random point from the surface of the triangle
    Eigen::VectorXd get_random_point(matrix_node_type const& element_points);

    Eigen::VectorXd get_element_normal (Eigen::VectorXd p1,Eigen::VectorXd p2,Eigen::VectorXd p3);
   

    // 3D case
    // Compute the sum of the areas of three subtriangles
    double elementArea(Eigen::VectorXd const& point,Eigen::VectorXd const& el_p1,Eigen::VectorXd const& el_p2,Eigen::VectorXd const& el_p3);

    // 3D case
    // Compare the area of the 2d simplex V1V2V3 (as sum of 3 subtriangles V_iV_jB) and the 2d triangle
    // created by the intersection P of the ray with the plane of  V1V2V3 (as sum of V_iV_jP)
    bool isOnSurface(Eigen::VectorXd const &point,Eigen::VectorXd const &el_p1,Eigen::VectorXd const &el_p2,Eigen::VectorXd const &el_p3);


    // Compute shading masks for the buildings in the json file
    void computeMasksMaster();
        void computeMasksSubPartList();
        void computeMasksSubPartSurfaceVolumes(int numOp);
        void computeMasksSubPartMarkersCTRL();

    void computeThread(Eigen::MatrixXd SM_table_marker,Eigen::MatrixXd Angle_table_marker,matrix_node_type const& element_points);

    // Save Compute shading masks
    void computeSaveMasks(std::vector<double> SM_tables);

    // Compute shading masks for one building only
    void computeMasksOneBuildingCTRL(std::string building_name);

    //void computeMasksOneBuildingOld(std::string building_name);

    // Save the shading mask table to a CSV file
    void saveShadingMask(std::string prefix_name,std::string building_name, std::string marker_name, const Eigen::Ref<const Eigen::MatrixXd>& M);

    // Save the shading mask table metadata to json
    void saveMetadata(std::string name);

    void saveMetadataInfoPart();

    void testComparisonAllMasksValidation();
    bool testComparisonMaskValidationLevel1(std::string matrix_filename_NEW,std::string matrix_filename_CTRL);
    bool testComparisonMaskValidationLevel2(std::string matrix_filename_NEW,std::string matrix_filename_CTRL);


private:
    
    std::map<std::string,std::unique_ptr<BVH<typename tr_mesh_type::element_type>>> M_bvh_tree_vector;
    std::map<std::string,tr_mesh_ptrtype> M_submeshes;
    std::map<int,node_type> M_faces_to_normals;

    //std::unique_ptr<BVH<typename tr_mesh_type::element_type>> M_bvh;

    using bvh_type = BVH<typename tr_mesh_type::element_type>;
    std::unique_ptr<bvh_type> M_bvh;

    std::unordered_map<int,std::string> M_mapEntityToBuildingFace;
    wrapper_type M_rangeFaces;
    std::vector<std::string> M_listFaceMarkers;

    Eigen::VectorXd M_azimuthAngles, M_altitudeAngles;

    std::vector< std::tuple< std::vector<double> /* random direction */,int /* index azimuth */, int /* index altitude */ > > M_raysdirections;

    std::map<std::string /* marker name */,std::vector<typename tr_mesh_type::element_type> /* list of references to face entities */> M_listMarkerFaceEntity;

    int M_azimuthSize;
    int M_altitudeSize;
    int M_Nrays;
    int M_Nthreads;
    std::string M_mthreadtype;
    bool M_saveMasks;

    int matrixSize;
    int dim=3;

    nl::json j_,M_metadataJson;

    std::random_device M_rd,M_rd2;
    std::mt19937 M_gen,M_gen2;

    std::vector<double> SM_tables_Alpha;
    std::vector<double> Angle_tables_Alpha;

    std::time_t beginning_time;

    
};
} // namespace Feel

#include "shading_mask_init.hpp"
#include "shading_mask_geometry.hpp"
#include "shading_mask_compute.hpp"
#include "shading_mask_save.hpp"
#include "shading_mask_control.hpp"






