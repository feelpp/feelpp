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
#include <cassert> 

namespace Feel
{      
    // Radiative heat transfer in a cavity
    template<int Dim,int Order=1>
    class RHT  
    {
    public:
        using mesh_t = Mesh<Simplex<Dim>>;
        using mesh_ptr_t = std::shared_ptr<mesh_t>;
        using mesh_trace_t = trace_mesh_t<mesh_t>;
        using mesh_trace_ptr_t = trace_mesh_ptr_t<mesh_t>;
        using space_t = Pch_type<mesh_t, Order>;
        using space_ptr_t = Pch_ptrtype<mesh_t, Order>;
        using spacedisc_t = Pdh_type<mesh_trace_t, Order-1>;
        using spacedisc_ptr_t = Pdh_ptrtype<mesh_trace_t, Order-1>;
        using spacedisc_surf_t = Pdh_type<mesh_t, 0>;
        using spacedisc_surf_ptr_t = Pdh_ptrtype<mesh_t, 0>;
        static constexpr int nDim = Dim;
        using scalar_t =  double;
        using coord_t = Eigen::Matrix<double, Dim, 1>;
        using tensor_t = Eigen::Matrix<double, Dim, Dim>;
        using element_t = typename space_t::element_type;
        using element_ptr_t = typename space_t::element_ptrtype;
        using elementdisc_ptr_t = typename spacedisc_t::element_ptrtype;
        using elementdisc_surf_ptr_t = typename spacedisc_surf_t::element_ptrtype;

        // Physics constants for radiative heat transfer
        static constexpr double STEFAN_BOLTZMANN_DERIVATIVE_COEFF = 4.0;  // d(T⁴)/dT = 4T³
        static constexpr int RADIATIVE_POWER = 4;                          // T^4 in Stefan-Boltzmann law

        RHT(nl::json specs)
        {
            // Assign the json structures to the members of the class
            this->specs=specs;            
        }    

        void init();
        void solveHeatEquationNonLinear(element_ptr_t T );

        using backend_type = Backend<double>;
        using backend_ptrtype = std::shared_ptr<backend_type>;

        /*matrix*/
        using sparse_matrix_type = typename backend_type::sparse_matrix_type;
        using sparse_matrix_ptrtype = typename backend_type::sparse_matrix_ptrtype;
        using vector_type = typename backend_type::vector_type;
        using vector_ptrtype = typename backend_type::vector_ptrtype;

        using bdf_type = Bdf<space_t>;
        using bdf_ptrtype = std::shared_ptr<bdf_type>;

        using exporter_type = Exporter<mesh_t,1>;
        using exporter_ptrtype = std::shared_ptr<exporter_type>;

        void executeNonLinear();
        void computeVF_and_save();
        void saveVF(const std::string& cavity_name, const Eigen::Ref<const Eigen::MatrixXd>& M);
        void loadVF(const std::string& cavity_name, const std::string& filename);
        void computeVF(const std::string& cavity_name, const std::string& filename);
        void checkResults();
        
        // Helper method to find coating emissivity for a given marker
        std::optional<std::string> getCoatingEpsilon(const std::string& marker) const;
        
        struct Tstruct
        {
            element_ptr_t T_;

            const element_ptr_t& T() const { return T_; }
            
            void setT(const element_ptr_t& T) { T_ = T; }

        };

        // void solvePicardIteration();
        void exportHeat();
        void initHeatEquation();        
        nl::json specs;
        nl::json j_viewfactor;
        space_ptr_t M_Xh;
        spacedisc_ptr_t M_Xhd0;
        std::map< std::string, spacedisc_ptr_t > M_Xhd0_map;
        spacedisc_surf_ptr_t M_Xhds0;
        elementdisc_surf_ptr_t M_conductivity;        
        mesh_ptr_t M_mesh;
        std::map< std::string, std::vector<std::string> > M_markers_map;   

        bdf_ptrtype M_bdf;
        Tstruct M_currentTemp;        

        std::map< std::string, Eigen::MatrixXd > M_matrix_vf_map; // matrix storing view factors

        sparse_matrix_ptrtype M_a,M_at; // matrices for heat transfer PDE
        vector_ptrtype M_l,M_lt; // right-hand sides for heat transfer PDE

        exporter_ptrtype M_e; // BDF exporter
    };

}