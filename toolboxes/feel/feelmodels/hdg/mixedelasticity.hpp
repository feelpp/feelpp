#ifndef _MIXEDELASTICITY_HPP
#define _MIXEDELASTICITY_HPP

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/complement.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <feel/feelts/newmark.hpp>
#include <boost/algorithm/string.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/projector.hpp>

// #define USE_SAME_MATH 1

namespace Feel {


template <typename SpaceType>
NullSpace<double> hdgNullSpace( SpaceType const& space, mpl::int_<2> /**/ )
{
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
    auto mode3 = space->element( vec(Py(),-Px()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3 } );
    return userNullSpace;
}

template <typename SpaceType>
NullSpace<double> hdgNullSpace( SpaceType const& space, mpl::int_<3> /**/ )
{
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
    auto mode3 = space->element( oneZ() );
    auto mode4 = space->element( vec(Py(),-Px(),cst(0.)) );
    auto mode5 = space->element( vec(-Pz(),cst(0.),Px()) );
    auto mode6 = space->element( vec(cst(0.),Pz(),-Py()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3,mode4,mode5,mode6 } );
    return userNullSpace;
}



namespace FeelModels {

inline
po::options_description
makeMixedElasticityOptions( std::string prefix = "mixedelasticity" )
{
    po::options_description mpOptions( "Mixed Elasticity HDG options");
    mpOptions.add_options()
        ( prefixvm( prefix, "gmsh.submesh").c_str(), po::value<std::string>()->default_value( "" ), "submesh extraction" )
        // ( "gmsh.submesh2", po::value<std::string>()->default_value( "" ), "submesh extraction" )
        ( prefixvm( prefix, "hface").c_str(), po::value<int>()->default_value( 0 ), "hface" )
        ( "lambda", po::value<std::string>()->default_value( "1" ), "lambda" )
        ( "mu", po::value<std::string>()->default_value( "1" ), "mu" )
        ( prefixvm( prefix, "model_json").c_str(), po::value<std::string>()->default_value("model.json"), "json file for the model")
        ( prefixvm( prefix,"tau_constant").c_str(), po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( prefixvm( prefix,"tau_order").c_str(), po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ( prefixvm( prefix, "use-sc").c_str(), po::value<bool>()->default_value(true), "use static condensation")
        ( prefixvm( prefix, "nullspace").c_str(), po::value<bool>()->default_value( false ), "add null space" )
        ;
    mpOptions.add( modelnumerical_options( prefix ) );
	mpOptions.add ( backend_options( prefix+".sc" ) );
    return mpOptions;
}

inline po::options_description
makeMixedElasticityLibOptions( std::string prefix = "mixedelasticity" )
{
    po::options_description mpLibOptions( "Mixed Elasticity HDG Lib options");
    // if ( !prefix.empty() )
    //     mpLibOptions.add( backend_options( prefix ) );
    return mpLibOptions;
}

template<int Dim, int Order, int G_Order = 1, int E_Order = 4>
class MixedElasticity    :	public ModelNumerical
{
public:
    typedef ModelNumerical super_type;


    static const uint16_type expr_order = Order + E_Order;
    //! numerical type is double
    typedef double value_type;
    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename std::shared_ptr<backend_type> backend_ptrtype ;

    typedef MixedElasticity<Dim,Order,G_Order,E_Order> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    using sparse_matrix_type = backend_type::sparse_matrix_type;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;
    using vector_type = backend_type::vector_type;
    using vector_ptrtype = backend_type::vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim,G_Order> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim-1,G_Order,Dim> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef std::shared_ptr<face_mesh_type> face_mesh_ptrtype;


    // ---- //
    using Vh_t =  Pdhms_type<mesh_type,Order>;
    using Vh_ptr_t =  Pdhms_ptrtype<mesh_type,Order>;
    using Vh_element_t = typename Vh_t::element_type;
    using Vh_element_ptr_t = typename Vh_t::element_ptrtype;
    // ---- //
    using Wh_t =  Pdhv_type<mesh_type,Order>;
    using Wh_ptr_t =  Pdhv_ptrtype<mesh_type,Order>;
    using Wh_element_t = typename Wh_t::element_type;
    using Wh_element_ptr_t = typename Wh_t::element_ptrtype;
    // ---- //
    using Mh_t =  Pdhv_type<face_mesh_type,Order>;
    using Mh_ptr_t =  Pdhv_ptrtype<face_mesh_type,Order>;
    using Mh_element_t = typename Mh_t::element_type;
    using Mh_element_ptr_t = typename Mh_t::element_ptrtype;
    // ---- //
    using M0h_t =  Pdh_type<face_mesh_type,0>;
    using M0h_ptr_t =  Pdh_ptrtype<face_mesh_type,0>;
    using M0h_element_t = typename M0h_t::element_type;
    using M0h_element_ptr_t = typename M0h_t::element_ptrtype;
    // ---- //
    using Ch_t = Pchv_type<face_mesh_type,0>;
    using Ch_ptr_t = Pchv_ptrtype<face_mesh_type,0>;
    using Ch_element_t = typename Ch_t::element_type;
    using Ch_element_ptr_t = typename Ch_t::element_ptrtype;
    using Ch_element_vector_type = std::vector<Ch_element_t>;



	using op_interp_ptrtype = std::shared_ptr<OperatorInterpolation<Wh_t, Pdhv_type<mesh_type,Order>>>;
    using opv_interp_ptrtype = std::shared_ptr<OperatorInterpolation<Vh_t, Pdhms_type<mesh_type,Order>>>;

    //! Model properties type
    using model_prop_type = ModelProperties;
    using model_prop_ptrtype = std::shared_ptr<model_prop_type>;

    /*
     using product_space_std = ProductSpaces<Vh_ptr_t,Wh_ptr_t,Mh_ptr_t>;
     using product_space_ptrtype = std::shared_ptr<product_space_std>;
     using bilinear_block_std = BlockBilinearForm<product_space_std>;
     */

    using product2_space_type = ProductSpaces2<Ch_ptr_t,Vh_ptr_t,Wh_ptr_t,Mh_ptr_t>;
    using product2_space_ptrtype = std::shared_ptr<product2_space_type>;
    using integral_boundary_list_type = std::vector<ExpressionStringAtMarker>;

    typedef Exporter<mesh_type,G_Order> exporter_type;
    typedef std::shared_ptr <exporter_type> exporter_ptrtype;

    // typedef Newmark<space_mixedelasticity_type>  newmark_type;
    typedef Newmark <Wh_t> newmark_type;
    typedef std::shared_ptr<newmark_type> newmark_ptrtype;

    //private:
protected:
    model_prop_ptrtype M_modelProperties;
    std::string M_prefix;

    mesh_ptrtype M_mesh;

    Vh_ptr_t M_Vh; // stress
    Wh_ptr_t M_Wh; // displacement
    Mh_ptr_t M_Mh; // displacement trace
    M0h_ptr_t M_M0h;
    Ch_ptr_t M_Ch; // Lagrange multiplier for IBC

	product2_space_ptrtype M_ps;

    backend_ptrtype M_backend;
    condensed_matrix_ptr_t<value_type> M_A_cst;
    condensed_vector_ptr_t<value_type> M_F;

    Vh_element_t M_up; // stress solution
    Wh_element_t M_pp; // displacement solution
    Ch_element_vector_type M_mup; // displacement solution on the IBC

    double M_tauCst;
    int M_tauOrder;
    int M_hFace;
    bool M_useSC;
    bool M_nullspace;

    int M_integralCondition;
    int M_useUserIBC;
    integral_boundary_list_type M_IBCList;

    // time discretization
    newmark_ptrtype M_nm_mixedelasticity;

    // save time
    std::map<std::string, std::vector<double> > M_timers;

public:

    // constructor
    MixedElasticity( std::string const& prefix = "mixedelasticity",
                     worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
                     std::string const& subPrefix = "",
                     ModelBaseRepository const& modelRep = ModelBaseRepository() );

    MixedElasticity( self_type const& ME ) = default;
    static self_ptrtype New( std::string const& prefix = "mixedelasticity",
                             worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                             std::string const& subPrefix = "",
                             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    // Get Methods
    mesh_ptrtype mesh() const { return M_mesh; }
    Vh_ptr_t fluxSpace() const { return M_Vh; }
    Wh_ptr_t potentialSpace() const { return M_Wh; }
    Mh_ptr_t traceSpace() const { return M_Mh; }
    M0h_ptr_t traceSpaceOrder0() const { return M_M0h; }
    Ch_ptr_t constantSpace() const {return M_Ch;}

    Vh_element_t fluxField() const { return M_up; }
    Wh_element_t potentialField() const { return M_pp; }
    model_prop_ptrtype modelProperties() { return M_modelProperties; }
    model_prop_ptrtype modelProperties() const { return M_modelProperties; }

    integral_boundary_list_type integralBoundaryList() const { return M_IBCList; }
    int integralCondition() const { return M_integralCondition; }
    void setIBCList(std::vector<std::string> markersIbc);
    product2_space_ptrtype getPS() const { return M_ps; }

    int tauOrder() const { return M_tauOrder; }
    void setTauOrder(int order) { M_tauOrder = order; }
    double tauCst() const { return M_tauCst; }
    void setTauCst(double cst) { M_tauCst = cst; }
    int hFace() const { return M_hFace; }
    void setHFace(int h) { M_hFace = h; }
    bool useSC() const { return M_useSC; }
    void setUseSC(bool sc) { M_useSC = sc; }
    bool nullspace() const { return M_nullspace; }
    void setNullSpace(bool nullspace) { M_nullspace = nullspace; }

    backend_ptrtype get_backend() { return M_backend; }
    condensed_vector_ptr_t<value_type> getF() {return M_F; }
    std::map<std::string, std::vector<double> > getTimers() {return M_timers; }

    // Exporter
    virtual void exportResults( mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr  )
        {
            this->exportResults (this->currentTime(), mesh , Idh, Idhv);
            M_exporter -> save();
        }

    void exportResults ( double Time, mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr  ) ;
    exporter_ptrtype M_exporter;
    exporter_ptrtype exporterME() { return M_exporter; }

    void init( mesh_ptrtype mesh = nullptr, mesh_ptrtype meshVisu = nullptr);

    virtual void initModel();
    virtual void initSpaces();
    virtual void initExporter( mesh_ptrtype meshVisu = nullptr );
    virtual void exportTimers();

    virtual void assemble();
    void assembleSTD();
    void assembleF();
    void assembleMatrixIBC(int i, std::string markerOpt = "" );
    void assembleRhsIBC(int i, std::string marker = "", double intjn = 0);

	void assembleCst();
	void assembleNonCst();

    void geometricTest();

    void solve();


    // time step scheme
    virtual void createTimeDiscretization() ;
    newmark_ptrtype timeStepNM() { return M_nm_mixedelasticity; }
    newmark_ptrtype const& timeStepNM() const { return M_nm_mixedelasticity; }
    std::shared_ptr<TSBase> timeStepBase() { return this->timeStepNM(); }
    std::shared_ptr<TSBase> timeStepBase() const { return this->timeStepNM(); }
    virtual void updateTimeStepNM();
    virtual void initTimeStep();
    void updateTimeStep() { this->updateTimeStepNM(); }

};

} // Namespace FeelModels

} // Namespace Feel

#endif
