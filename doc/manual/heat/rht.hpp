#include <iostream>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelviewfactor/unobstructedplanarviewfactor.hpp>
#include <assert.h> 

namespace Feel
{      
    template<int Dim,int Order=1>
    class RHT  
    {
    public:
        using mesh_t = Mesh<Simplex<Dim>>;
        using mesh_ptr_t = std::shared_ptr<mesh_t>;
        using mesh_trace_t = typename mesh_t::trace_mesh_type;
        using mesh_trace_ptr_t = typename mesh_t::trace_mesh_ptrtype;
        using spacevect_t = Pchv_type<mesh_t, Order>;
        using spacevect_ptr_t = Pchv_ptrtype<mesh_t, Order>;
        using space_t = Pch_type<mesh_t, Order>;
        using space_ptr_t = Pch_ptrtype<mesh_t, Order>;
        using spacedisc_t = Pdh_type<mesh_trace_t, Order-1>;
        using spacedisc_ptr_t = Pdh_ptrtype<mesh_trace_t, Order-1>;
        using spacedisc_interm_ptr_t = Pdh_ptrtype<mesh_trace_t, Order>;
        static constexpr int nDim = Dim;
        using scalar_t =  double;
        using coord_t = Eigen::Matrix<double, Dim, 1>;
        using tensor_t = Eigen::Matrix<double, Dim, Dim>;
        using elementvect_t = typename spacevect_t::element_type;
        using elementvect_ptr_t = typename spacevect_t::element_ptrtype;
        using element_t = typename space_t::element_type;
        using element_ptr_t = typename space_t::element_ptrtype;
        using elementdisc_ptr_t = typename spacedisc_t::element_ptrtype;

        RHT(nl::json specs)
        {
            // Assign the json structures to the members of the class
            this->specs=specs;            
        }    

        void init();
        void solveHeatEquation(element_ptr_t T, element_ptr_t q );
        void solveHeatEquationNonLinear(element_ptr_t T );

        typedef Backend<double> backend_type;
        typedef std::shared_ptr<backend_type> backend_ptrtype;

        /*matrix*/
        typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
        typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
        typedef typename backend_type::vector_type vector_type;
        typedef typename backend_type::vector_ptrtype vector_ptrtype;

        typedef Bdf<space_t>  bdf_type;
        typedef std::shared_ptr<bdf_type> bdf_ptrtype;

        typedef Exporter<mesh_t,1> exporter_type;
        typedef std::shared_ptr <exporter_type> exporter_ptrtype;

        void rebuild_D(element_ptr_t T);
        void rebuild_N(element_ptr_t q);
        void build_M();
        void execute();
        void executeNonLinear();
        void computeVF_and_save();
        void saveVF(std::string cavity_name,const Eigen::Ref<const Eigen::MatrixXd>&  M);
        void loadVF(std::string cavity_name,std::string filename);
        void computeVF(std::string cavity_name,std::string filename);
        void checkResults();
        void read_Tsky();
        void readTskyFromCSV(std::string const& meteo_station_name, std::string const& filename);
        double accessTsky(std::string const& meteo_station_name,double const& time, double const& dt);
        void read_Text();
        void readTextFromCSV(std::string const& meteo_station_name, std::string const& filename);
        double accessText(std::string const& meteo_station_name,double const& time, double const& dt);
        void read_solarRadiation();
        void readSolarRadiationFromCSV(std::string const& meteo_station_name, std::string const& building_name, std::string const& surface_name, std::string const& directory);
        double accessSolarRadiation(std::string const& meteo_station_name, std::string const& building_name, std::string const& surface_name, double const& time, double const& dt);

        struct TandQ
        {
            element_ptr_t T_;
            element_ptr_t q_;

            element_ptr_t T(){return T_;}
            element_ptr_t q(){return q_;}

            // void setT(element_ptr_t& T){T_=boost::unwrap_ref(T);}
            // void setq(element_ptr_t& q){q_=boost::unwrap_ref(q);}
            void setT(element_ptr_t& T){T_=T;}
            void setq(element_ptr_t& q){q_=q;}
        };

        void solvePicardIteration();
        void exportHeat();
        void initHeatEquation();        
        nl::json specs;
        nl::json j_viewfactor;
        space_ptr_t M_Xh;
        spacevect_ptr_t M_Xhvec;
        spacedisc_ptr_t M_Xhd0;
        spacedisc_interm_ptr_t M_Xhd1;
        std::map< std::string, spacedisc_ptr_t > M_Xhd0_map;
        std::map< std::string, spacedisc_interm_ptr_t > M_Xhd1_map;
        mesh_ptr_t M_mesh;
        std::map< std::string, mesh_trace_ptr_t> M_surface_submesh_map;
        std::map< std::string, std::vector<std::string> > M_markers_map;   
        std::map< std::string, std::vector<double> > M_solarRadiation_structure;
        std::map< std::string, std::vector<double> > M_station_Tsky;
        std::map< std::string, std::vector<double> > M_station_Text;

        bdf_ptrtype M_bdf;
        TandQ M_currentTempAndFlux;        

        sparse_matrix_ptrtype M_M, M_M_block; // matrix for problem on multiple irradiating surfaces
        std::map< std::string, Eigen::MatrixXd > M_matrix_vf_map; // matrix storing view factors
        vector_ptrtype M_D, M_N, M_D_block,M_N_block; // vectors for problem on multiple irradiating surfaces

        sparse_matrix_ptrtype M_a,M_at; // matrices for heat transfer PDE
        vector_ptrtype M_l,M_lt; // right-hand sides for heat transfer PDE

        exporter_ptrtype M_e; // BDF exporter
    };

}