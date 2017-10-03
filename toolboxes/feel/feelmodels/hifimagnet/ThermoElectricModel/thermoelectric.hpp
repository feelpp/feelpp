#ifndef __THERMOELECTRIC_3D_HPP
#define __THERMOELECTRIC_3D_HPP 1

#include <fstream>
#include <iostream>
#include <string>

/** include predefined feel command line options */
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>

/** include linear algebra backend */
#include <feel/feelalg/backend.hpp>

/** include function space class */
#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/product.hpp>

#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feeldiscr/dh.hpp>

/** include gmsh mesh importer */
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>

/** include exporter factory class */
#include <feel/feelfilters/exporter.hpp>

/** for high order visualisation **/
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>

/** include general hifimagnet applications options **/
#include <feel/feelmodels/hifimagnet/options.hpp>

/** include modelproperties **/
#include <feel/feelmodels/modelproperties.hpp>

/** include mesh_initializer for Hifimagnet **/
#include <feel/feelmodels/hifimagnet/Tools/mesh_initializer.hpp>
#include <feel/feelmodels/hifimagnet/Tools/stats.hpp>
#include <feel/feelmodels/hifimagnet/Tools/init_field.hpp>

/** include the header for mesh adaptation **/
#include <feel/feelmodels/hifimagnet/ThermoElectricModel/estimators.hpp>
#if defined( FEELPP_HAS_GMSH_ADAPT_H )
#include <feel/feelmesh/meshadaptation.hpp>
#endif

/** include icontains **/
#include <boost/algorithm/string.hpp>
//#include <boost/any.hpp>

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
Feel::po::options_description
ThermoElectricOptions()
{
    Feel::po::options_description thermoelectric_options("thermoelectric options");
    thermoelectric_options.add_options()
        ( "thermoelectric.verbosity" , Feel::po::value<int>()->default_value( 0 ), "verbosity level" )
        ( "thermoelectric.resolution", Feel::po::value<std::string>()->default_value( "picard"), "newton or picard resolution method for nonlinear model" )
        ( "thermoelectric.weakdir", Feel::po::value<bool>()->default_value( true ), "use weak Dirichlet condition" )
        ( "thermoelectric.penaldir", Feel::po::value<double>()->default_value( 50 ), "penalisation parameter for the weak boundary Dirichlet" )
        ( "thermoelectric.model_json", Feel::po::value<std::string>()->default_value( "model.json" ), "use json file for model properties" )
        ( "thermoelectric.print_info", Feel::po::value<bool>()->default_value( true ), "print some info on the model" )
        // intensity
        ( "thermoelectric.update_intensity", Feel::po::value<bool>()->default_value( false ), "update boundary condition to have target intensity" )
        ( "thermoelectric.target_intensity", Feel::po::value< std::vector<double> >()->default_value(std::vector<double>()), "vector of target intensities" )
        ( "thermoelectric.marker_intensity", Feel::po::value< std::vector<std::string> >()->default_value(std::vector<std::string>()), "vecotor of markers for which intensity is given" )
        ( "thermoelectric.eps_intensity", Feel::po::value<double>()->default_value( 1e-10 ), "error criterion for intensity" )

        ( "thermoelectric.itmax_picard", Feel::po::value<int>()->default_value( 10 ), "max iteration for Picard, Newton and current" )
        ( "thermoelectric.eps_potential", Feel::po::value<double>()->default_value( 1e-10 ), "error criterion for potential" )
        ( "thermoelectric.eps_temperature", Feel::po::value<double>()->default_value( 1e-10 ), "error criterion for temperature" )
        // ( "thermoelectric.V_guess", Feel::po::value<double>()->default_value( 0 ), "initial guess for potential (Volt)" )
        // ( "thermoelectric.T_guess", Feel::po::value<double>()->default_value( 293 ), "initial guess for temperature (Kelvin)" )
        // ( "thermoelectric.V_guess_filename", Feel::po::value<std::string>()->default_value( "" ), "initial guess for potential (h5 file)" )
        // ( "thermoelectric.T_guess_filename", Feel::po::value<std::string>()->default_value( "" ), "initial guess for temperature (h5 file)" )
        ( "thermoelectric.V_guess", Feel::po::value<std::string>()->default_value( "0" ), "initial guess for potential (Volt)" )
        ( "thermoelectric.T_guess", Feel::po::value<std::string>()->default_value( "293" ), "initial guess for temperature (Kelvin)" )
        ( "thermoelectric.restart", Feel::po::value<bool>()->default_value( true ), "restart from initial solutions (id provided by V_guess and T_guess)" )

        ( "thermoelectric.estimator_type",Feel::po::value<std::string>()->default_value( "zz" ),"residual estimator=res, ZZ estimator=zz, anisotropic")
        ;
    return thermoelectric_options;
}

inline
Feel::po::options_description
ThermoElectricOptionsLib()
{
    Feel::po::options_description thermoelectric_optionsLib("thermoelectric options");
    thermoelectric_optionsLib.add(Feel::backend_options("electro"));
    thermoelectric_optionsLib.add(Feel::backend_options("thermal"));
    thermoelectric_optionsLib.add(Feel::backend_options("newton"));
    return thermoelectric_optionsLib;
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline
Feel::AboutData
makeAboutThermoElectric()
{
    Feel::AboutData about( "thermoelectric" ,
                           "thermoelectric" ,
                           "0.1",
                           "nD(n=1,2,3) Potential and Temperature on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008-2009 Universite Joseph Fourier"
                           "Copyright (c) CNRS");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "");
    about.addAuthor("Cecile Daversin", "developer", "cecile.daversin@lncmi.cnrs.fr", "");

    return about;

}

// Class thermoelectric
namespace Feel
{
    template<int Dim, int OrderV, int OrderT, int G_order>
    class thermoelectric
    {
    public:
        using self_type = thermoelectric<Dim, OrderV, OrderT, G_order>;
        using self_ptrtype = boost::shared_ptr<self_type>;

        //! numerical type is double
        typedef double value_type;

        //! linear algebra backend factory
        typedef Backend<value_type> backend_type;
        //! linear algebra backend factory shared_ptr<> type
        typedef boost::shared_ptr<backend_type> backend_ptrtype;

        //! sparse matrix type associated with backend
        typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
        //! sparse matrix type associated with backend (shared_ptr<> type)
        typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
        //! vector type associated with backend
        typedef typename backend_type::vector_type vector_type;
        //! vector type associated with backend (shared_ptr<> type)
        typedef typename backend_type::vector_ptrtype vector_ptrtype;

        //! geometry entities type composing the mesh, here Simplex in Dimension Dim (Order G_order)
        typedef Simplex<Dim,G_order> convex_type;
        typedef Mesh<convex_type> mesh_type;
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

        //POTENTIAL
        using V_space_type = Pch_type<mesh_type,OrderV>;
        using V_space_ptrtype =  boost::shared_ptr<V_space_type>;
        typedef typename V_space_type::element_type V_element_type;

        //CURRENT
        using C_space_type = Pdhv_type<mesh_type,OrderV-1>;
        using C_space_ptrtype =  boost::shared_ptr<C_space_type>;
        typedef typename C_space_type::element_type C_element_type;

        //TEMPERATURE
        using  T_space_type = Pch_type<mesh_type,OrderT>;
        typedef boost::shared_ptr<T_space_type> T_space_ptrtype;
        typedef typename T_space_type::element_type T_element_type;

        // //! Product space types
        // typedef Lagrange<OrderV,Scalar> single_basis_typeV;
        // typedef Lagrange<OrderT,Scalar> single_basis_typeT;
        // typedef bases<single_basis_typeV, single_basis_typeT> prod_basis_type;
        // typedef FunctionSpace< mesh_type, prod_basis_type, value_type > prod_space_type;
        // typedef boost::shared_ptr<prod_space_type> prod_space_ptrtype;

        //! Scalar P0 space
        using p0_space_type = Pdh_type<mesh_type,0>;
        typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
        typedef typename p0_space_type::element_type p0_element_type;

        //! the exporter factory type
        typedef Exporter<mesh_type> export_type;
        //! the exporter factory (shared_ptr<> type)
        typedef boost::shared_ptr<export_type> export_ptrtype;

        //! Model properties type
        typedef ModelProperties model_prop_type;
        typedef std::shared_ptr<model_prop_type> model_prop_ptrtype;
        //! Material properties type
        typedef ModelMaterials material_map_type;
        //! Boundary conditions type
        typedef BoundaryConditions condition_map_type;

        //For estimator class
        typedef Estimators<Dim, OrderV, OrderT, G_order> estimator;

        // For mesh adaptation
        #if defined( FEELPP_HAS_GMSH_ADAPT_H )
        typedef Eigen::Matrix<double, Dim, 1> vectorN_type;
        typedef Eigen::Matrix<double, Dim, Dim> matrixN_type;

        typedef MeshAdaptation<Dim, OrderT, G_order> MeshAdapt;
        #endif

    private:
        //! mesh characteristic size
        double M_meshSize; // Coil

        //! shape of the domain
        std::string M_shape;

        //! Geofiles input (conductor)
        std::string M_geofile;
        std::string M_geoDepends;
        std::string M_geofilePath;

        std::vector<std::string> M_domains;
        std::vector<std::string> M_conductors;

        bool M_restart;
        
        //! Mesh adaptation method / type
        bool M_meshadapt;
        int M_meshadapt_itMax;
        double M_meshadapt_tol;
        std::string M_meshadaptMethod;
        std::string M_meshadaptType;

        /* mesh, spaces and elements */
        mesh_ptrtype M_mesh;
        mesh_ptrtype M_mesh_V; //mesh for conductors
        V_space_ptrtype M_XhV;
        C_space_ptrtype M_XhC;
        T_space_ptrtype M_XhT;
        V_element_type M_V;
        C_element_type M_C;
        T_element_type M_T;
        V_element_type M_VOld;
        T_element_type M_TOld;

        // options
        int M_verbosity;
        bool M_weakdir;
        double M_penaldir;
        std::string M_resolution;
        double M_epsilonV;
        double M_epsilonT;
        int M_itMax;

        // intensity
        double M_sigmaMax;
        bool M_updateIntensity;
        std::vector<std::string> M_markerIntensity;
        std::vector<double> M_vIntensity;
        std::vector<double> M_targetIntensity;
        double M_epsilonIntensity;


        // properties
        model_prop_ptrtype M_modelProps;
        material_map_type M_matProps;
        bool M_isNonLinear;

        //backends
        backend_ptrtype M_V_backend;
        backend_ptrtype M_T_backend;

        sparse_matrix_ptrtype M_aV;
        vector_ptrtype M_rV;
        sparse_matrix_ptrtype M_aT;
        vector_ptrtype M_rT;

        // ho exporter
        bool M_isHOVisu;

    public:
        thermoelectric();
        static self_ptrtype New();
        void changeRepository();

        void init( mesh_ptrtype mesh = NULL, std::vector<std::string> = {} );
        void solve();
        void solveLinear();
        void solvePicard();
        void solveNewton();
        void adapt();
        void exportResults(double time = 0 );
        void V_exportResults(double time = 0 );
        void T_exportResults(double time = 0 );
        void exportHOResults(double time = 0 );
        bool isValidModel();

        bool is_P1V() const {return OrderV == 1;};
        bool is_P1T() const {return OrderT == 1;};
        bool is_adaptmesh() const {return M_meshadapt;};

        V_element_type electricField() const { return M_V; }
        T_element_type temperature() const { return M_T; }
        C_element_type current() const { return M_C; }

        mesh_ptrtype potentialFieldMesh() const { return M_mesh_V; }

        model_prop_ptrtype modelProperties() const { return M_modelProps; }

        bool isConductor(std::string );
        bool isDomain(std::string);

    private:
        bool initOldVT();
        void initParametersV(bool isRestart);

        void computeV( int iter = -1 );
        std::vector<double> updateIntensityBC();
        void computeVBilinear( int iter = -1 );
        void computeVBoundaryCond( int iter = -1 );
        void computeVStrongDirichlet();

        void computeVFlux(); // current density
        void computeCurrentIntensities(); // current density

        void computeT( int iter = -1 );
        void computeTBilinear( int iter = -1 );
        void computeTLinear( int iter = -1 );
        void computeTBoundaryCond( int iter = -1 );

        void computeTFlux(); // heat flux

    }; // class thermoelectric

} // namespace

template<int Dim, int OrderV, int OrderT, int G_order>
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::thermoelectric()
{
    M_meshSize = doption( "hsize" );
    M_shape = soption( "shape" );
    M_geofile = soption( "geofile" );
    M_geoDepends = soption( "geo_depends" );
    M_geofilePath = soption( "geofile-path" );
    M_meshadapt = boption( "meshadapt" );
    M_meshadapt_itMax = ioption("meshadapt_maxiter");
    M_meshadapt_tol = doption("meshadapt_tol");
    M_meshadaptMethod = soption( "meshadapt_method" );
    M_meshadaptType = soption( "meshadapt_type" );
    M_verbosity = ioption( "thermoelectric.verbosity" );
    M_weakdir = boption( "thermoelectric.weakdir" );
    M_penaldir = doption( "thermoelectric.penaldir" );
    M_resolution = soption( "thermoelectric.resolution");
    M_epsilonV = doption( "thermoelectric.eps_potential");
    M_epsilonT = doption( "thermoelectric.eps_temperature");
    M_itMax = ioption( "thermoelectric.itmax_picard");
    M_updateIntensity = boption("thermoelectric.update_intensity");
    M_modelProps = std::make_shared<model_prop_type>( Environment::expand( soption("thermoelectric.model_json") ) );
    if ( !isValidModel() )
    {
        Feel::cout << "ThermoElectric Model " << M_modelProps->model() << " not supported! Aborting" << std::endl;
        LOG(FATAL) << "ThermoElectric Model " << M_modelProps->model() << " not supported! Aborting";
    }
    M_matProps = M_modelProps->materials();
    M_isNonLinear = boost::icontains(M_modelProps->model(), "nonlinear");

    M_isHOVisu = (OrderV > 1) || (OrderT > 1);
    if ( M_isHOVisu )
        M_isHOVisu = boption("export-highorders");
}

template<int Dim, int OrderV, int OrderT, int G_order>
typename Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::self_ptrtype
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::New()
{
    return boost::make_shared<self_type>();
}

template<int Dim, int OrderV, int OrderT, int G_order>
bool
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::isValidModel()
{
    bool b = false;
    b = b || boost::icontains(M_modelProps->model(),"thermoelectric-linear")
        || boost::icontains(M_modelProps->model(),"thermoelectric-nonlinear")
        || boost::icontains(M_modelProps->model(),"coupled")
        || boost::icontains(M_modelProps->model(),"coupled-nonlinear");

    // additionnal check when M_meshadapt
    return b;

}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::changeRepository( )
{
    /*** eventually change output dir***/
    std::string appname = "ThermoElectricModel";
    Environment::changeRepository( boost::format( "hifimagnet/%1%/%2%/%3%_V%4%T%5%_N%6%/" )
                                   % appname
                                   % M_geofile
                                   % M_modelProps->model()
                                   % OrderV
                                   % OrderT
                                   % G_order );

}
template<int Dim, int OrderV, int OrderT, int G_order>
bool Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::isConductor(std::string name)
{
    bool check = false;
    if ( !M_conductors.empty() )
        check = check || std::find(std::begin(M_conductors), std::end(M_conductors), name) != std::end(M_conductors);
    return check;
}

template<int Dim, int OrderV, int OrderT, int G_order>
bool Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::isDomain(std::string name)
{
    bool check = false;
    if ( !M_domains.empty() )
        check = check || std::find(std::begin(M_domains), std::end(M_domains), name) != std::end(M_domains);
    return check;
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::init( mesh_ptrtype mesh, std::vector<std::string> th_domains  )
{
    tic();
    // Feel::cout << "M_meshadapt=" << M_meshadapt << " (" << boption("meshadapt") << ")" << std::endl;
    // Feel::cout << "M_verbosity=" << M_verbosity << "(" <<  ioption( "thermoelectric.verbosity" ) << ")" << std::endl;

    if ( !mesh )
    {
        mesh_initializer<Dim, G_order> mesh_init;
        M_mesh = mesh_init.initializeMesh(Environment::vm());
    } else
    {
        M_mesh = mesh;
    }

    // create submesh for conductor only
    // Feel::cout << "loading conductors : " << std::flush;
    M_conductors = vsoption(_name="conductor_volume");

    // check if conductor are properly declared as markers
    // Feel::cout << "\n" << "Check conductor_list: ";
    bool check_conductor_list = false;
    for( auto marker: M_mesh->markerNames() )
        {
            auto name = marker.first;
            auto data = marker.second;
            // Feel::cout << name << ", dim=" << data[1] << ", Dim=" << Dim << ",  M_mesh->dimension()=" <<  M_mesh->dimension() << std::endl;
            if ( data[1] == M_mesh->dimension() )
                {
                    check_conductor_list = check_conductor_list || isConductor(name);
                    // Feel::cout << name << "[" ;
                    // Feel::cout << (std::find(std::begin(conductor_list), std::end(conductor_list), name) != std::end(conductor_list)) << ",";
                    // Feel::cout << (std::find(std::begin(conductor), std::end(conductor), name) != std::end(conductor)) << "]\n";
                }
        }
    // Feel::cout << std::endl;

    if (!check_conductor_list)
        {
            Feel::cout << "ThermoElectric Model : invalid conductor_list [";
            for (std::string item : M_conductors)
                Feel::cout << item << ",";
            Feel::cout << "]\n" << std::flush;
            LOG(FATAL) << "ThermoElectric Model : invalid conductor_list! Aborting";
        }

    LOG(INFO) << "Vh_domains:" << std::endl;
    for (auto domain: M_conductors)
        LOG(INFO) << domain << ", ";
    LOG(INFO) << std::endl;

    if ( M_conductors.empty() )
        M_mesh_V = M_mesh;
    else
        M_mesh_V = createSubmesh(M_mesh, markedelements(M_mesh, M_conductors));

    // Th Domains
    if ( !th_domains.empty() )
        std::copy( th_domains.begin(), th_domains.end(), std::back_inserter( M_domains ) );

    if ( M_domains.empty() )
        for( auto marker: M_mesh->markerNames() )
            {
                auto name = marker.first;
                M_domains.push_back(name);
            }
    LOG(INFO) << "M_domains:" << std::endl;
    for (auto domain: M_domains)
        LOG(INFO) << domain << ", ";
    LOG(INFO) << std::endl;

    // V for Potential, C for Current density and T for Temperature
    M_XhV = V_space_type::New( M_mesh_V );
    M_XhC = C_space_type::New( M_mesh_V );
    M_XhT = T_space_type::New( M_mesh );
    M_V = M_XhV->element();
    M_C = M_XhC->element();
    M_T = M_XhT->element();

    M_VOld = M_XhV->element();
    M_TOld = M_XhT->element();

    M_V_backend = backend(_name="electro", _rebuild=true);
    M_aV = M_V_backend->newMatrix( _test=M_XhV, _trial=M_XhV );
    M_rV = M_V_backend->newVector( _test=M_XhV );
    M_T_backend = backend(_name="thermal", _rebuild=true);
    M_aT = M_T_backend->newMatrix( _test=M_XhT, _trial=M_XhT );
    M_rT = M_T_backend->newVector( _test=M_XhT );

    M_restart = boption("thermoelectric.restart");
    bool M_init = this->initOldVT(); // should only be done once if ( M_meshadaptMethod != "no" )
    this->initParametersV(M_restart && M_init);

    int proc_rank = Environment::worldComm().globalRank();
    if( boption( "thermoelectric.print_info" ) && proc_rank==0 )
    {
        cout << std::endl;
        cout << "model : " << M_modelProps->model() << std::endl;
        cout << "geofile name : " << M_geofile << std::endl;
        cout << "Number of dofs for potential : " << M_XhV->nDof() << std::endl;
        cout << "Number of dofs for temperature : " << M_XhT->nDof() << std::endl;
        cout << "Number of vertices in mesh : " <<  M_mesh->numGlobalVertices() << std::endl;
        cout << "Number of faces in mesh : " <<  ( (Dim==3) ? M_mesh->numGlobalFaces() : M_mesh->numGlobalElements() ) << std::endl;
        cout << "Number of edges in mesh : " << ( (Dim==3) ? M_mesh->numGlobalEdges() :M_mesh->numGlobalFaces() ) << std::endl;
        if ( Dim==3 )
            cout << "Number of elts in mesh : " << M_mesh->numGlobalElements() << std::endl;

        if( M_weakdir )
            {
                cout << "Use weak Dirichlet conditions " << std::endl;
                cout << "penalisation coefficient : " << M_penaldir << std::endl;
            }
        else
            cout << "Use strong Dirichlet conditions " << std::endl;

        if( M_updateIntensity )
            for (unsigned int i=0; i<M_markerIntensity.size(); i++)
                {
                    cout << "[" <<  M_markerIntensity[i] << "] Update intensity to target : " << M_targetIntensity[i] << std::endl;
                }

        if ( M_meshadapt  )
            {
                cout << "Mesh adaptation method : " << M_meshadaptMethod  << std::endl;
                cout << "Mesh adaptation : " << M_meshadaptType  << std::endl;
            }

        if ( M_isHOVisu )
            cout << "High Order Visu : " << M_isHOVisu << std::endl;

        cout << "High Order Visu : " << M_isHOVisu << std::endl;
        cout << "export-highorders =" << boption("export-highorders") << std::endl;
        cout << "OrderV = " << OrderV << ", OrderT = " << OrderT << std::endl;
        //cout << "Environment::vm().count(export-highorders) = " << Environment::vm().count("export-highorders") << std::endl;
    }
    toc("init", M_verbosity > 0);
}

/** Initialize old temperature and potential with initial guess or from a file */
template<int Dim, int OrderV, int OrderT, int G_order>
bool
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::initOldVT()
{
    tic();
    bool result = false;

    bool isV_double = init_field(M_VOld, "V", soption("thermoelectric.V_guess"));
    bool isT_double = init_field(M_TOld, "T", soption("thermoelectric.T_guess"));

    M_T = M_TOld;
    M_V = M_VOld;

    if ( !isV_double && !isT_double )
        {
            Feel::cout << "**** Checking Inititial States ****" << std::endl;
            stats(M_V, "Init Potential", "V");
            // auto V_min = M_VOld.min();
            // auto V_max = M_VOld.max();
            // Feel::cout << "Vmin=" << V_min << " V, ";
            // Feel::cout << "Vmax=" << V_max << " V" << std::endl;

            stats(M_T, "Init Temperature", "K");
            // auto T_min = M_TOld.min();
            // auto T_max = M_TOld.max();
            // Feel::cout << "Tmin=" << T_min << " K, ";
            // Feel::cout << "Tmax=" << T_max << " K" << std::endl;

            computeCurrentIntensities();
            Feel::cout << "*********************************" << std::endl;

            result = true;
        }

    toc("init V and T", M_verbosity > 0);
    return result;
}

/** Initialize max of relative electric conductivity and BC */
template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::initParametersV(bool isRestart)
{
    tic();
    M_sigmaMax = 0;
    for( auto const& pairMat : M_matProps )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        if ( isConductor(marker) )
            {
                auto sigmaProj = vf::project( _range=markedelements(M_mesh_V, marker),
                                              _space=M_XhV,
                                              _expr=material.getScalar("sigma0") );
                auto norm = sigmaProj.linftyNorm(); // M_sigmaMax i no longer a sigma?? should divide by volume??
                if ( norm > M_sigmaMax )
                    M_sigmaMax = norm;
            }
    }

    if ( M_updateIntensity)
    {
        M_targetIntensity = vdoption("thermoelectric.target_intensity");
        M_markerIntensity = vsoption("thermoelectric.marker_intensity");
        M_epsilonIntensity = doption("thermoelectric.eps_intensity");

        CHECK( M_targetIntensity.size() == M_markerIntensity.size() ) << "incoherent set of imposed current";
        M_vIntensity.resize(M_targetIntensity.size(), 0);
        // int proc_rank = Environment::worldComm().globalRank();
        // std::cout << "np=" << proc_rank << ":" <<  " M_vIntensity size=" << M_vIntensity.size() << std::endl;

        if (isRestart)
            {
                // should pick a point on marker_intensity and retreive associated V when restarting from previous calculations
                // should make it optionnal
                for (unsigned int i=0; i<M_markerIntensity.size(); i++)
                    {
                        double Mean_V = integrate(markedfaces(M_mesh_V, M_markerIntensity[i]), idv(M_V) ).evaluate()(0,0);
                        double Area =  integrate(markedfaces(M_mesh_V, M_markerIntensity[i]), cst(1.0) ).evaluate()(0,0);
                        M_vIntensity[i] = Mean_V/Area;
                        Feel::cout << "restart using V[" << M_markerIntensity[i] << "] : "  << M_vIntensity[i] << std::endl;
                    }
            }
        else
            {
                auto itField = M_modelProps->boundaryConditions().find( "potential");
                if ( itField != M_modelProps->boundaryConditions().end() )
                    {
                        auto mapField = (*itField).second;
                        auto itType = mapField.find( "Dirichlet" );
                        if ( itType != mapField.end() )
                            {
                                for ( auto const& exAtMarker : (*itType).second )
                                    {
                                        std::string marker = exAtMarker.marker();
                                        auto g = expr(exAtMarker.expression1());
                                        std::string domain = exAtMarker.expression2();

                                        for (unsigned int i=0; i<M_markerIntensity.size(); i++)
                                            {
                                                if ( marker == M_markerIntensity[i] )
                                                    {
                                                        // std::cout << "np=" << proc_rank << ": " << std::flush;
                                                        // std::cout << "M_markerIntensity[" << i << "]=" << M_markerIntensity[i] << " " << std::flush;
                                                        // std::cout << "marker=" << marker << std::flush;
                                                        // std::cout << "V=" << exAtMarker.expression() << std::flush; // if V is null pb!!!
                                                        M_vIntensity[i] = std::stod(exAtMarker.expression1());
                                                        Feel::cout << "V[" << M_markerIntensity[i] << "] : "  << M_vIntensity[i] << std::endl;
                                                        // std::cout << "M_vIntensity[" << i << "]=" << M_vIntensity[i] << " " << std::flush;
                                                        // std::cout << std::endl;
                                                    }
                                            }
                                    }
                            }
                    }
            }
    }
    toc("initParametersV", M_verbosity > 0);
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::solve()
{
    tic();
    Feel::cout << "thermoelectric::solve(...)" << std::endl;
    if ( !M_isNonLinear )
        this->solveLinear();
    else if ( boost::iequals(M_resolution, "picard") )
        this->solvePicard();
    else if ( boost::iequals(M_resolution, "newton") )
        this->solveNewton();
    else
        Feel::cout << "WARNING! resolution strategy unknwon ( Newton or Picard in non-linear case )";

    computeVFlux();
    computeCurrentIntensities();

    auto V_min = M_V.min();
    auto V_max = M_V.max();
    Feel::cout << "Vmin=" << V_min << " V, ";
    Feel::cout << "Vmax=" << V_max << " V" << std::endl;

    auto T_min = M_T.min();
    auto T_max = M_T.max();
    Feel::cout << "Tmin=" << T_min << " K, ";
    Feel::cout << "Tmax=" << T_max << " K" << std::endl;

    toc("solve",  M_verbosity > 0);
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::adapt()
{
    tic();
    Feel::cout << "Mesh adaptation [method=" << M_meshadaptMethod
                           << ", type=" << M_meshadaptType
                           << ", max_iter=" << M_meshadapt_itMax << "]...\n" << std::endl;
    estimator estimator;
    double mesh_eps = 1.0;// _tol = geometric tolerance (for adaptation from hessian) - Alauzet formula - p35 Metric-Based Anisotropic Mesh Adaptation

    std::string access_geofile = M_geofile;
    if ( !M_geofilePath.empty() )
        access_geofile = ( boost::format( "%1%/%2%" ) % Environment::expand(M_geofilePath) % M_geofile ).str();
    Feel::cout << "access_geofile =" << access_geofile << std::endl;

    std::string geo_depends = soption("geo_depends");
    Feel::cout << "geo_depends =" << geo_depends << std::endl;

#if defined( FEELPP_HAS_GMSH_ADAPT_H )

    double criterion_T; // Relative error
    bool criterion = false;
    int iter = 0;
    MeshAdapt mesh_adaptation;
    do {
        solve();

        auto M_V0 = M_V;
        auto M_T0 = M_T;

        // error estimation: only ZZ, should consider M_meshapdatMethod
        boost::tuple<double, double, p0_element_type> estimator_T;
        estimator_T = estimator.zz_estimator(M_T, M_mesh);

        // std::vector<vectorN_type> measures = mesh_adaptation.measures();
        // std::vector<matrixN_type> directions = mesh_adaptation.directions();
        // boost::tuple<double, p0_element_type> estimator_T = estimator.anisotropic_estimator_T(M_T, M_V, M_mesh, M_matProps, M_modelProps, M_weakdir, M_penaldir, measures, directions);

        auto norm_T = normL2( elements(M_mesh), idv(M_T) );
        criterion_T = math::abs( estimator_T.template get<0>() - norm_T  )/norm_T;
        criterion = (criterion_T >= (1-M_meshadapt_tol) && criterion_T <= (1+M_meshadapt_tol) );

        Feel::cout << M_meshadaptType << "_Adapt#" << iter << " ; ";
        Feel::cout << "Potential=" << estimator.zz_estimator(M_V, M_mesh_V).template get<0>() << ", ";
        Feel::cout << "Temperature=" << estimator_T.template get<0>() << ", ";
        Feel::cout << "criterion =" << criterion_T << "(";
        Feel::cout << "tol=[" << (1-M_meshadapt_tol) << "," << (1+M_meshadapt_tol) << "])";
        Feel::cout << "criterion=" << criterion;
        Feel::cout <<std::endl;
        if ( criterion )
            break;

        std::list<std::pair<T_element_type, std::string> > var_list;
        std::pair<T_element_type, std::string> T_pair = std::make_pair( M_T, "temperature");
        var_list.push_back(T_pair); //adapt from temperature values

        Feel::cout << "adaptmesh";
        M_mesh = mesh_adaptation.adaptMesh(_initMesh=M_mesh, _geofile=access_geofile, _adaptType=M_meshadaptType, _var=var_list, _tol=mesh_eps);
        Feel::cout << "Number of faces in mesh : " << M_mesh->numFaces() << std::endl;
        Feel::cout << "Number of edges in mesh : " << M_mesh->numEdges() << std::endl;
        Feel::cout << "Number of elts in mesh : " << M_mesh->numElements() << std::endl;

        Feel::cout << "init" << std::endl;
        init(M_mesh);

        // how to interpolate M_V0 on new mesh M_mesh_V and get a new M_VOld??
    } while ( !criterion && iter++  < M_meshadapt_itMax );

#else
    Feel::cout << "ThermoElectricModel: adapt not available - requires gmsh to be patched \n";
#endif

    toc("adapt", M_verbosity > 0);
}


template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::computeVFlux()
{
    // Feel::cout << "ComputeVFlux" << std::flush;
    tic();
    for( auto const& pairMat : M_matProps )
        {
            auto marker = pairMat.first;
            auto material = pairMat.second;
            // check if marker is also defined in M_mesh_V
            if ( isConductor(marker) ) //nelements( markedelements(M_mesh_V, marker) ) > 0 )
                {
                    auto sigma0 = material.getDouble("sigma0");
                    if ( M_isNonLinear )
                        {
                            auto alpha = material.getDouble("alpha");
                            auto T0 = material.getDouble("T0");
                            auto sigma = material.getScalar("sigma", {"T"}, {idv(M_T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
                            M_C += vf::project(_range=markedelements(M_mesh_V, marker), _space=M_XhC, _expr=-sigma*trans(gradv(M_V)) );
                        }
                    else
                        M_C += vf::project(_range=markedelements(M_mesh_V, marker), _space=M_XhC, _expr=-cst(sigma0)*trans(gradv(M_V)) );
                }
        }
    // Feel::cout << " done" << std::endl;
    toc("ComputeVFlux", M_verbosity > 0);
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::computeTFlux()
{
    tic();
    cout << "thermoelectric<Dim, OrderV, OrderT, G_order>::computeTFlux : not implemented yet\n";
    // for( auto marker: M_mesh->markerNames() )
    //     {
    //         auto name = marker.first;
    //         auto data = marker.second;
    //         if ( data[1] == M_mesh_T->dimension()-1 and isDomain(name) )
    //             {
    //                 auto material = M_matProps[name];
    //                 auto alpha = material.getDouble("alpha");
    //                 auto T0 = material.getDouble("T0");
    //                 auto k0 = material.getDouble("k0");

    //                 if ( M_isNonLinear )
    //                     {
    //                          auto k = material.getScalar("k", {"T"}, {idv(M_T)}, {{"k0",k0},{"T0",T0},{"alpha",alpha}});
    //                          auto flux =  integrate( markedfaces(M_mesh, marker), -k*gradv(M_T)*N() ).evaluate()(0,0);
    //                     }
    //                 else
    //                     auto flux =   integrate( markedfaces(M_mesh, marker), -cst(k0)*gradv(M_T)*N() ).evaluate()(0,0);
    //             }
    //     }
    toc("ComputeTFlux", M_verbosity > 0);
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::solveLinear()
{
    // Feel::cout << "linear solver" << std::endl;
    double errorIntensity = 0;
    int i = 0;
    do {
        tic();
        this->computeV();
        toc("computeV", M_verbosity > 0);

        if ( M_updateIntensity )
        {
            std::vector<double> errorI = this->updateIntensityBC();
            std::vector<double>::iterator res = std::max_element(errorI.begin(), errorI.end()); // does it work in parallel??
            if( M_verbosity > 2 )
                {
                    for (unsigned int j=0; j<M_markerIntensity.size(); j++)
                        {
                            Feel::cout << "ErrorI[" << M_markerIntensity[j] << "] = " << errorI[j] << std::endl;
                        }
                    Feel::cout  << std::endl;
                }
            errorIntensity = *res;
        }
    } while ( errorIntensity > M_epsilonIntensity && i++ < M_itMax );

    tic();
    this->computeT();
    toc("computeT", M_verbosity > 0);
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::solvePicard()
{
    // Feel::cout << "Picard solver" << std::endl;

    int i=0, j=0;
    double errorV = std::numeric_limits<double>::max();
    double errorT = std::numeric_limits<double>::max();
    double errorI;
    if ( M_updateIntensity )
        errorI = std::numeric_limits<double>::max();
    else
        errorI = 0;

    std::string msg;

    do {
        tic();
        this->computeV(i);
        toc("compute potential", M_verbosity > 0);
        errorV = normL2(elements(M_mesh_V), idv(M_V)-idv(M_VOld)) / normL2(elements(M_mesh_V), idv(M_V));
        M_VOld = M_V;
        msg = (boost::format("Picard Iter=%1% ErrorV=%2% (%3%)") % i % errorV % M_epsilonV).str();

        // j = 0;
        // do {
        //     tic();
        //     this->computeT(j);
        //     toc("compute temperature", M_verbosity > 0 );
        //     errorT = normL2(elements(M_mesh), idv(M_T)-idv(M_TOld)) / normL2(elements(M_mesh), idv(M_T));
        //     M_TOld = M_T;
        // } while ( (errorT > M_epsilonT) && (j++ < M_itMax) );
        // msg += (boost::format(" SubIter(T)=%1% ErrorT=%2% (%3%)") % j % errorT % M_epsilonT).str();

        tic();
        this->computeT(i);
        toc("compute temperature", M_verbosity > 0 );
        errorT = normL2(elements(M_mesh), idv(M_T)-idv(M_TOld)) / normL2(elements(M_mesh), idv(M_T));
        M_TOld = M_T;
        msg += (boost::format(" ErrorT=%1% (%2%)") % errorT % M_epsilonT).str();

        if (j > M_itMax )
            {
                std::string errmsg = msg + " Failed to converge";
                throw std::logic_error( errmsg );
            }

        if ( M_updateIntensity )
            {
                std::vector<double> errI = this->updateIntensityBC();
                std::vector<double>::iterator res = std::max_element(errI.begin(), errI.end()); // does it work properly in parallel??
                if ( M_verbosity > 2 )
                    {
                        for (unsigned int l=0;  l<M_targetIntensity.size(); l++)
                             Feel::cout << "ErrorI[" << M_markerIntensity[j] << "] = " << errI[j] << std::endl;
                    }
                errorI = *res;
                if ( M_verbosity > 2 )
                    Feel::cout  << std::endl;

                msg += (boost::format(" ErrorI=%1% (%2%)") % errorI % M_epsilonIntensity).str();
            }

        Feel::cout << tc::green << msg << tc::reset << std::endl;


    } while ( ( (errorV > M_epsilonV)
                || (errorT > M_epsilonT)
                || (errorI > M_epsilonIntensity) )
              && (i++ < M_itMax) );

    if (i > M_itMax )
        {
            std::string errmsg = msg + " Failed to converge";
            throw std::logic_error( errmsg );
        }
    Feel::cout << "\n[Picard] " << tc::green << msg << tc::reset << std::endl;
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::solveNewton()
{
    if ( M_verbosity > 0 )
        Feel::cout << "Newton solver not yet implemented" << std::endl;
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::computeV( int iter )
{
    tic();
    form1( _test=M_XhV, _vector=M_rV) = integrate( elements(M_mesh_V), cst(0.));
    form2( _test=M_XhV, _trial=M_XhV, _matrix=M_aV) = integrate( elements(M_mesh_V), cst(0.));

    this->computeVBilinear(iter);
    this->computeVBoundaryCond(iter);
    if ( !M_weakdir )
        this->computeVStrongDirichlet();

    auto near_nullspace = NullSpace<value_type>( {M_XhV->element(cst(1.))} );
    auto result = M_V_backend->solve(_matrix=M_aV, _rhs=M_rV, _solution=M_V, _near_null_space=near_nullspace);
    std::string msg = (boost::format("\n[Potential %1%] NbIter=%2% Residual=%3%") % soption("electro.pc-type") % result.nIterations() % result.residual() ).str();
    if (result.isConverged())
        {
            Feel::cout << tc::green << msg << tc::reset << std::endl;
        }
    else
        {
            std::string errmsg = msg + " Failed to converge";
            throw std::logic_error( errmsg );
        }
    toc("potential [compute V]", M_verbosity > 0);

#if defined(FEELPP_HAS_HDF5)
    tic();
    LOG(INFO) << " Saving results in HDF5 format (" << M_geofile << ")";
    std::string M_mshfile;
    if ( !M_geofilePath.empty() )
        M_mshfile = ( boost::format( "%1%/%2%" ) % Environment::expand( M_geofilePath ) % M_geofile ).str();
    fs::path M_mshfile_path(M_mshfile);
    M_V.saveHDF5(M_mshfile_path.stem().string()+"_V.h5");
    LOG(INFO) << " HDF5 data saved";
    toc("potential [save V to hdf5]", M_verbosity > 0);
#endif

    if ( M_verbosity > 1 )
        {
            stats(M_V, "Potential", "V");
            auto V_min = M_V.min();
            auto V_max = M_V.max();
            Feel::cout << "Vmin=" << V_min << " V, ";
            Feel::cout << "Vmax=" << V_max << " V" << std::endl;

            computeCurrentIntensities();
        }
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::computeCurrentIntensities( )
{
    /* compute current intensities */
    auto itField = M_modelProps->boundaryConditions().find( "potential");
    if ( itField != M_modelProps->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression1());
                std::string domain = exAtMarker.expression2();

                auto material = M_matProps[domain];
                // cout << marker << " " << material << "\n";
                auto sigma0 = material.getDouble("sigma0");

                double I, tmp;
                if ( M_isNonLinear )
                    {
                        auto alpha = material.getDouble("alpha");
                        auto T0 = material.getDouble("T0");
                        auto sigma = material.getScalar("sigma", {"T"}, {idv(M_T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
                        I = integrate( markedfaces(M_mesh_V, marker),
                                       -sigma*gradv(M_V)*N() ).evaluate()(0,0);
                    }
                else
                    I = integrate( markedfaces(M_mesh_V, marker),
                                   -sigma0*gradv(M_V)*N() ).evaluate()(0,0);

                cout << "I[" << marker << "," << domain << "] : " <<  I << " A\n" << std::flush;
            }
        }
    }

}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::computeVBilinear( int iter )
{
    tic();
    // Feel::cout << "computeVBilinear:\n" <<std::endl;
    auto V = M_XhV->element();

    for( auto const& pairMat : M_matProps )
    {
        auto marker = pairMat.first;
        // Feel::cout << marker << ":" << std::flush;
        auto material = pairMat.second;
        // Feel::cout << material << "," << std::flush;
        // Feel::cout << "isConductor=" << isConductor(marker) << std::flush;
        if ( isConductor(marker) )
            {
                // cout << material << ":" <<std::flush;
                auto sigma0 = material.getDouble("sigma0")/M_sigmaMax;  // cout << "sigma0\n" << std::flush;
                auto alpha = material.getDouble("alpha"); // cout << "alpha\n" << std::flush;
                auto T0 = material.getDouble("T0"); // cout << "T0\n" << std::flush;
                auto sigma = material.getScalar("sigma", {"T"}, {idv(M_T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}}); // cout << "sigma(T, T0)\n" << std::flush;

                if ( M_resolution == "picard" && iter > 0 )
                    {
                        form2( _test=M_XhV, _trial=M_XhV, _matrix=M_aV ) +=
                            integrate( markedelements(M_mesh_V, marker), sigma*inner(gradt(V),grad(V)) );
                    }
                else
                    {
                        form2( _test=M_XhV, _trial=M_XhV, _matrix=M_aV ) +=
                            integrate( markedelements(M_mesh_V, marker), sigma0*inner(gradt(V),grad(V)) );
                    }
            }
        // Feel::cout << "... Done" << std::endl;
    }
    toc("computeVBilinear", M_verbosity > 0);
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::computeVBoundaryCond( int iter )
{
    tic();
    auto a = form2( _test=M_XhV, _trial=M_XhV, _matrix=M_aV);
    auto l = form1( _test=M_XhV, _vector=M_rV);
    auto V = M_XhV->element();

    auto itField = M_modelProps->boundaryConditions().find( "potential");
    if ( itField != M_modelProps->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() && M_weakdir )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression1());
                std::string domain = exAtMarker.expression2();

                auto material = M_matProps[domain];
                cout << marker << " " << material << "\n";
                auto alpha = material.getDouble("alpha");
                auto T0 = material.getDouble("T0");
                auto sigma0 = material.getDouble("sigma0");
                auto sigma = material.getScalar("sigma", {"T"}, {idv(M_T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});

                if( !M_updateIntensity || std::find(M_markerIntensity.begin(), M_markerIntensity.end(), marker) == M_markerIntensity.end() )
                    {
                        if ( M_resolution == "picard" && iter > 0 )
                            {
                                a += integrate( markedfaces(M_mesh_V, marker),
                                                -sigma/M_sigmaMax*(gradt(V)*N())*id(V)
                                                -sigma/M_sigmaMax*(grad(V)*N())*idt(V)
                                                +sigma/M_sigmaMax*M_penaldir*id(V)*idt(V)/hFace() );
                                l += integrate( markedfaces(M_mesh_V, marker),
                                                - sigma/M_sigmaMax*grad(V)*N()*g
                                                + sigma/M_sigmaMax*M_penaldir*id(V)*g/hFace() );
                            }
                        else
                            {
                                a += integrate( markedfaces(M_mesh_V, marker),
                                                -sigma0/M_sigmaMax*(gradt(V)*N())*id(V)
                                                -sigma0/M_sigmaMax*(grad(V)*N())*idt(V)
                                                +sigma0/M_sigmaMax*M_penaldir*id(V)*idt(V)/hFace() );
                                l += integrate( markedfaces(M_mesh_V, marker),
                                                - sigma0/M_sigmaMax*grad(V)*N()*g
                                                + sigma0/M_sigmaMax*M_penaldir*id(V)*g/hFace() );
                            }
                    }
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());
                l += integrate( markedfaces(M_mesh_V, marker), -g*id(V)/cst(M_sigmaMax) );
            }
        }
        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g1 = expr(exAtMarker.expression1());
                auto g2 = expr(exAtMarker.expression2());
                a += integrate( markedfaces(M_mesh_V, marker), g1*idt(V)*id(V)/cst(M_sigmaMax) );
                l += integrate( markedfaces(M_mesh_V, marker), g2*id(V)/cst(M_sigmaMax) );
            }
        }
    }
    if ( M_updateIntensity )
    {
        for (unsigned int i=0; i<M_markerIntensity.size(); i++)
            {
                auto itField = M_modelProps->boundaryConditions().find( "potential");
                if ( itField != M_modelProps->boundaryConditions().end() )
                    {
                        auto mapField = (*itField).second;
                        auto itType = mapField.find( "Dirichlet" );
                        if ( itType != mapField.end() && M_weakdir )
                            {
                                for ( auto const& exAtMarker : (*itType).second )
                                    {
                                        std::string marker = exAtMarker.marker();
                                        auto g = expr(exAtMarker.expression1());
                                        std::string domain = exAtMarker.expression2();

                                        if ( marker == M_markerIntensity[i] )
                                            {
                                                auto material = M_matProps[ domain ];
                                                auto alpha = material.getDouble("alpha");
                                                auto T0 = material.getDouble("T0");
                                                auto sigma0 = material.getDouble("sigma0");
                                                auto sigma = material.getScalar("sigma", {"T"}, {idv(M_T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});

                                                if ( M_resolution == "picard" && iter > 0 )
                                                    {
                                                        a += integrate( markedfaces(M_mesh_V, M_markerIntensity[i]),
                                                                        -sigma/M_sigmaMax*(gradt(V)*N())*id(V)
                                                                        -sigma/M_sigmaMax*(grad(V)*N())*idt(V)
                                                                        +sigma/M_sigmaMax*M_penaldir*id(V)*idt(V)/hFace() );
                                                        l += integrate( markedfaces(M_mesh_V, M_markerIntensity[i]),
                                                                        - sigma/M_sigmaMax*grad(V)*N()*cst(M_vIntensity[i])
                                                                        + sigma/M_sigmaMax*M_penaldir*id(V)*cst(M_vIntensity[i])/hFace() );
                                                    }
                                                else
                                                    {
                                                        a += integrate( markedfaces(M_mesh_V, M_markerIntensity[i]),
                                                                        -sigma0/M_sigmaMax*(gradt(V)*N())*id(V)
                                                                        -sigma0/M_sigmaMax*(grad(V)*N())*idt(V)
                                                                        +sigma0/M_sigmaMax*M_penaldir*id(V)*idt(V)/hFace() );
                                                        l += integrate( markedfaces(M_mesh_V, M_markerIntensity[i]),
                                                                        - sigma0/M_sigmaMax*grad(V)*N()*cst(M_vIntensity[i])
                                                                        + sigma0/M_sigmaMax*M_penaldir*id(V)*cst(M_vIntensity[i])/hFace() );
                                                    }
                                            }
                                    }
                            }
                    }
            }
    }
    toc("computeVBoundaryCond", M_verbosity > 0);
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::computeVStrongDirichlet()
{
    auto a = form2( _test=M_XhV, _trial=M_XhV, _matrix=M_aV);
    auto l = form1( _test=M_XhV, _vector=M_rV);

    auto itField = M_modelProps->boundaryConditions().find( "potential");
    if ( itField != M_modelProps->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if( !M_updateIntensity || std::find(M_markerIntensity.begin(), M_markerIntensity.end(), marker) == M_markerIntensity.end() )
                {
                    auto g = expr(exAtMarker.expression());
                    a += on( markedfaces(M_mesh_V, marker), M_V, M_rV, g);
                }
            }
        }
    }
    if ( M_updateIntensity )
        for (unsigned int i=0; i<M_markerIntensity.size(); i++)
            a += on( markedfaces(M_mesh_V, M_markerIntensity[i]), M_V, M_rV, cst(M_vIntensity[i]));
}

template<int Dim, int OrderV, int OrderT, int G_order>
std::vector<double>
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::updateIntensityBC()
{
    int proc_rank = Environment::worldComm().globalRank();
    // std::cout << "\nnp=" << proc_rank << ": updateIntensityBC()" << std::endl;

    std::vector<double> errorI(M_markerIntensity.size());
    // std::cout << "np" << proc_rank <<  ": M_markerIntensity.size()=" << M_markerIntensity.size() << std::endl;
    for (unsigned int i=0; i<M_markerIntensity.size(); i++)
        {
            // std::cout << "np" << proc_rank <<  ": M_markerIntensity[" << i << "]=" << M_markerIntensity[i] << std::endl;
            auto itField = M_modelProps->boundaryConditions().find( "potential");
            if ( itField != M_modelProps->boundaryConditions().end() )
                {
                    auto mapField = (*itField).second;
                    auto itType = mapField.find( "Dirichlet" );
                    if ( itType != mapField.end() )
                        {
                            for ( auto const& exAtMarker : (*itType).second )
                                {
                                    std::string marker = exAtMarker.marker();
                                    auto g = expr(exAtMarker.expression1());
                                    std::string domain = exAtMarker.expression2();
                                    // std::cout << "np" << proc_rank << ": marker=" << marker << " domain=" << domain <<  std::endl;

                                    if ( marker == M_markerIntensity[i] )
                                        {
                                            auto material = M_matProps[ domain ];
                                            auto sigma0 = material.getDouble("sigma0");

                                            double currentIntensity;
                                            if ( M_isNonLinear )
                                                {
                                                    auto alpha = material.getDouble("alpha");
                                                    auto T0 = material.getDouble("T0");
                                                    auto sigma = material.getScalar("sigma", {"T"}, {idv(M_T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
                                                    currentIntensity = integrate( markedfaces(M_mesh_V, M_markerIntensity[i]),
                                                                                  -sigma*gradv(M_V)*N() ).evaluate()(0,0);
                                                }
                                            else
                                                currentIntensity = integrate( markedfaces(M_mesh_V, M_markerIntensity[i]),
                                                                              -sigma0*gradv(M_V)*N() ).evaluate()(0,0);

                                            errorI[i] = std::abs(currentIntensity/M_targetIntensity[i]-1);
                                            M_vIntensity[i] *= M_targetIntensity[i]/currentIntensity;
                                            // std::cout << "np=" << proc_rank << ": I=" << currentIntensity << " errorI[" << i << "]=" << errorI[i] << " M_vIntensity[" << i << "]=" << M_vIntensity[i] << std::endl;

                                            if ( M_verbosity > 0 )
                                                {
                                                    Feel::cout << "I[" << M_markerIntensity[i] << "] : " << currentIntensity;
                                                    Feel::cout<< " target I=" << M_targetIntensity[i];
                                                    Feel::cout<< " new V=" << M_vIntensity[i] << std::endl;
                                                    // std::cout << "I[" << M_markerIntensity[i] << ",proc=" << proc_rank << "] : " << currentIntensity;
                                                    // std::cout<< " target I=" << M_targetIntensity[i];
                                                    // std::cout<< " new V=" << M_vIntensity[i] << std::flush << std::endl;
                                                }
                                        }
                                }
                        }
                }
        }
    return errorI;
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::computeT( int iter )
{
    form1( _test=M_XhT, _vector=M_rT) = integrate( elements(M_mesh), cst(0.));
    form2( _test=M_XhT, _trial=M_XhT, _matrix=M_aT) = integrate( elements(M_mesh), cst(0.));

    this->computeTBilinear(iter);
    this->computeTLinear(iter);
    this->computeTBoundaryCond(iter);

    auto near_nullspace = NullSpace<value_type>( {M_XhT->element(cst(1.))} );
    auto result = M_T_backend->solve(_matrix=M_aT, _rhs=M_rT, _solution=M_T, _near_null_space=near_nullspace);
    std::string msg = (boost::format("\n[Temperature %1%] NbIter=%2% Residual=%3%") % soption("thermal.pc-type") % result.nIterations() % result.residual() ).str();
    if (result.isConverged())
        {
            Feel::cout << tc::green << msg << tc::reset << std::endl;
        }
    else
        {
            std::string errmsg = msg + " Failed to converge";
            throw std::logic_error( errmsg );
        }

#if defined(FEELPP_HAS_HDF5)
    tic();
    LOG(INFO) << " Saving results in HDF5 format (" << M_geofile << ")";
    std::string M_mshfile;
    if ( !M_geofilePath.empty() )
        M_mshfile = ( boost::format( "%1%/%2%" ) % Environment::expand( M_geofilePath ) % M_geofile ).str();
    fs::path M_mshfile_path(M_mshfile);
    M_T.saveHDF5(M_mshfile_path.stem().string()+"_T.h5");
    LOG(INFO) << " HDF5 data saved";
    toc("temperature [save T to hdf5]", M_verbosity > 0);
#endif

    if ( M_verbosity > 1 )
        {
            stats(M_T, "Temperature", "K", (M_matProps.size() > 1) ? true : false);
            auto T_min = M_T.min();
            auto T_max = M_T.max();
            Feel::cout << "Tmin=" << T_min << " K, ";
            Feel::cout << "Tmax=" << T_max << " K" << std::endl;
        }
    if ( M_verbosity > 2 )
        stats_bcs(M_T, M_modelProps, "temperature", "W");
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::computeTBilinear( int iter )
{
    tic();
    // Feel::cout << "computeTBilinear" << std::endl;
    auto T = M_XhT->element();

    for( auto const& pairMat : M_matProps )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        if ( isDomain(marker) )
            {
                // Feel::cout << material << std::endl;
                auto alpha = material.getDouble("alpha");
                auto T0 = material.getDouble("T0");
                auto k0 = material.getDouble("k0");
                auto k = material.getScalar("k", {"T"}, {idv(M_T)}, {{"k0",k0},{"T0",T0},{"alpha",alpha}});

                if ( M_resolution == "picard" && iter > 0 )
                    {
                        form2( _test=M_XhT, _trial=M_XhT, _matrix=M_aT ) +=
                            integrate( markedelements(M_mesh, marker), k*inner(gradt(T),grad(T)) );
                    }
                else
                    {
                        form2( _test=M_XhT, _trial=M_XhT, _matrix=M_aT ) +=
                            integrate( markedelements(M_mesh, marker), k0*inner(gradt(T),grad(T)) );
                    }
            }
    }
    toc("computeTBilinear", M_verbosity >0);
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::computeTLinear( int iter )
{
    tic();
    // Feel::cout << "computeTLinear" << std::endl;
    auto T = M_XhT->element();

    for( auto const& pairMat : M_matProps )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;

        // check if marker is also defined in M_mesh_V
        if ( isDomain(marker) ) //nelements( markedelements(M_mesh_V, marker) ) > 0 )
            {
                // Feel::cout << "computeTLinear[" << material << std::endl;
                auto alpha = material.getDouble("alpha");
                auto T0 = material.getDouble("T0");
                auto sigma0 = material.getDouble("sigma0");
                auto sigma = material.getScalar("sigma", {"T"}, {idv(M_T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});

                if ( M_resolution == "picard" && iter > 0 )
                    {
                        form1( _test=M_XhT, _vector=M_rT ) +=
                            integrate( markedelements(M_mesh, marker), sigma*inner(gradv(M_V),gradv(M_V))*id(T) );
                    }
                else
                    {
                        form1( _test=M_XhT, _vector=M_rT ) +=
                            integrate( markedelements(M_mesh, marker), sigma0*inner(gradv(M_V),gradv(M_V))*id(T) );
                    }
            }
        else
            {
                form1( _test=M_XhT, _vector=M_rT ) +=
                    integrate( markedelements(M_mesh, marker), cst(0.)*id(T) );
            }
    }
    toc("computeTlinear", M_verbosity >0);
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::computeTBoundaryCond( int iter )
{
    tic();
    auto a = form2( _test=M_XhT, _trial=M_XhT, _matrix=M_aT);
    auto l = form1( _test=M_XhT, _vector=M_rT);
    auto T = M_XhT->element();

    auto itField = M_modelProps->boundaryConditions().find( "temperature");
    if ( itField != M_modelProps->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression1());
                std::string domain = exAtMarker.expression2();

                auto material = M_matProps[domain];
                // Feel::cout << marker << " " << material << "\n";
                auto alpha = material.getDouble("alpha");
                auto T0 = material.getDouble("T0");
                auto k0 = material.getDouble("k0");
                auto k = material.getScalar("k", {"T"}, {idv(M_T)}, {{"k0",k0},{"T0",T0},{"alpha",alpha}});

                if ( M_resolution == "picard" && iter > 0 )
                    {
                        a += integrate( markedfaces(M_mesh, marker),
                                        -k*(gradt(T)*N())*id(T)
                                        -k*(grad(T)*N())*idt(T)
                                        +k*M_penaldir*id(T)*idt(T)/hFace() );
                        l += integrate( markedfaces(M_mesh, marker),
                                        - k*grad(T)*N()*g
                                        + k*M_penaldir*id(T)*g/hFace() );
                    }
                else
                    {
                        a += integrate( markedfaces(M_mesh, marker),
                                        -k0*(gradt(T)*N())*id(T)
                                        -k0*(grad(T)*N())*idt(T)
                                        +k0*M_penaldir*id(T)*idt(T)/hFace() );
                        l += integrate( markedfaces(M_mesh, marker),
                                        - k0*grad(T)*N()*g
                                        + k0*M_penaldir*id(T)*g/hFace() );
                    }
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());
                l += integrate( markedfaces(M_mesh, marker), g*id(T) );
            }
        }
        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g1 = expr(exAtMarker.expression1());//h
                auto g2 = expr(exAtMarker.expression2());//Tw
                a += integrate( markedfaces(M_mesh, marker), g1*idt(T)*id(T) );
                l += integrate( markedfaces(M_mesh, marker), g1*g2*id(T) );
            }
        }
    }
    toc("computeTBoundaryCond", M_verbosity>0);
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::exportHOResults(double time)
{
    // CHECK( false ) << "exportHOResults not implement\n";
    if ( OrderT > 1)
        {
            auto m1 = lagrangeP1(_space=M_T.functionSpace())->mesh();
            auto XhVisu = Pch<1>(m1);

            auto opIVisu = opInterpolation(_domainSpace=M_T.functionSpace(),
                                           _imageSpace=XhVisu,
                                           _type=InterpolationNonConforme(false,true,false) );
            auto Visu = opIVisu->operator()(M_T);

            auto e = exporter( _mesh=m1, _name="HO_Thermics");
            e->step(0)->add( "T", Visu );
            e->save();
        }
    if ( OrderV > 1)
        {
            auto m1 = lagrangeP1(_space=M_V.functionSpace())->mesh();
            auto XhVisu = Pch<1>(m1);

            auto opIVisu = opInterpolation(_domainSpace=M_V.functionSpace(),
                                           _imageSpace=XhVisu,
                                           _type=InterpolationNonConforme(false,true,false) );
            auto Visu = opIVisu->operator()(M_V);

            auto e = exporter( _mesh=m1, _name="HO_Electrics");
            e->step(0)->add( "V", Visu );
            e->save();
        }
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::exportResults(double time)
{
#if defined(FEELPP_HAS_VTK)
    if (M_isHOVisu)
        {
            exportHOResults(time);
        }
#endif
    V_exportResults(time);
    T_exportResults(time);
}

template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::T_exportResults(double time)
{
#if defined(FEELPP_HAS_VTK)
    if (M_isHOVisu)
        {
            exportHOResults(time);
            //return;
        }
#endif

    auto e_T = exporter( _mesh=M_mesh, _name="Thermics");

    auto postProcess = M_modelProps->postProcess();
    auto itField = postProcess.find( "Fields");
    if ( itField != postProcess.end() )
    {
        for ( auto const& field : (*itField).second )
        {
            if ( field == "temperature" )
                e_T->step(time)->add( "temperature", M_T);
            if ( M_isNonLinear && field == "conductivity" )
                {
                    auto M_XhK = Pdh<OrderT>(M_mesh);
                    auto M_K = M_XhK->element();
                    for( auto const& pairMat : M_matProps )
                        {
                            auto marker = pairMat.first;
                            if ( isDomain(marker) )
                                {
                                    // Feel::cout << "K[" << marker << "] " << std::flush;
                                    auto material = pairMat.second;
                                    auto alpha = material.getDouble("alpha");
                                    auto T0 = material.getDouble("T0");
                                    auto k0 = material.getDouble("k0");
                                    auto k = material.getScalar("k", {"T"}, {idv(M_T)}, {{"k0",k0},{"T0",T0},{"alpha",alpha}});

                                    M_K += vf::project( _space=M_XhK, _range=markedelements(M_mesh, marker), _expr=k );
                                    // Feel::cout << "\n" << std::flush;
                                }
                        }
                    e_T->step(time)->add( "k", M_K);
                }
            // need first to adapt estimator class to new BCs
            std::string estimator_type  = soption(_name="thermoelectric.estimator_type");
            if ( field == "estimator" )
                {
                    // auto normL2 = normL2( elements(M_mesh), idv(M_T) );
                    // auto normH1 = normL2( elements(M_mesh), gradv(M_T) );

                    estimator estimator;
                    boost::tuple<double, double, p0_element_type> estimator_T;
                    if( estimator_type == "res")
                        {
                            estimator_T = estimator.residual_estimator_T(M_T, M_V, M_mesh, M_matProps, M_modelProps, M_weakdir);
                            e_T->step(time)->add( "estim_err", estimator_T.template get<2>());

                            Feel::cout << "\n ****** Residual estimators ****** \n";
                            Feel::cout << "            L2 norm\tH1 norm\n";
                            Feel::cout << "Temperature "
                                       << estimator_T.template get<0>() << "\t"
                                       << estimator_T.template get<1>() << "\n";
                                       // << math::abs( estimator_T.template get<0>() - normL2  )/normL2 << "\t"
                                       // << math::abs( estimator_T.template get<1>() - normH1  )/normH1 << "\n";
                            Feel::cout << " *********************************************" << std::endl;
                        }
                    else if( estimator_type == "zz" && (OrderV==1 && OrderT==1) )
                        {
                            estimator_T = estimator.zz_estimator(M_T, M_mesh);
                            e_T->step(time)->add( "estim_err", estimator_T.template get<2>());
                            Feel::cout<< "\n ****** ZZ estimators ****** \n";
                            Feel::cout << "            L2 norm\tH1 norm\n";
                            Feel::cout << "Temperature "
                                       << estimator_T.template get<0>() << "\t"
                                       << estimator_T.template get<1>() << "\n";
                            Feel::cout << " *********************************************" << std::endl;
                    }
                    else if( estimator_type == "anisotropic" && OrderT==1 )
                        {
#if defined( FEELPP_HAS_GMSH_ADAPT_H )
                            MeshAdapt mesh_adaptation;

                            std::vector<vectorN_type> measures = mesh_adaptation.measures();
                            std::vector<matrixN_type> directions = mesh_adaptation.directions();
                            boost::tuple<double, p0_element_type> estimator_T = estimator.anisotropic_estimator_T(M_T, M_V, M_mesh, M_matProps, M_modelProps, M_weakdir, M_penaldir, measures, directions);

                            Feel::cout << "\n ****** Anisotropic estimator ****** \n"
                                       << "              H1 semi-norm\n"
                                       << "Temperature " << estimator_T.template get<0>() << "\n"
                                       << " *********************************************** \n";
#else
                            Feel::cout << "ThermoElectricModel: Anisotropic estimator not available - requires gmsh to be patched \n";
#endif
                    }
                }
        }
    }

    e_T->save();
}
template<int Dim, int OrderV, int OrderT, int G_order>
void
Feel::thermoelectric<Dim, OrderV, OrderT, G_order>::V_exportResults(double time)
{
    auto e_V = exporter( _mesh=M_mesh_V, _name="Electrics");

    auto postProcess = M_modelProps->postProcess();
    auto itField = postProcess.find( "Fields");
    if ( itField != postProcess.end() )
    {
        for ( auto const& field : (*itField).second )
        {
            if ( field == "current" )
                {
                    e_V->step(time)->add( "current", M_C);
                }
            if ( field == "joules" )
                {
                    auto M_XhQ = Pdh<OrderV-1>(M_mesh_V);
                    auto M_Q = M_XhQ->element();
                    for( auto marker: M_mesh_V->markerNames() )
                        {
                            auto name = marker.first;
                            auto data = marker.second;
                            if ( data[1] == M_mesh_V->dimension() && isConductor(name) ) //&& nelements( markedelements(M_mesh_V, name) ) > 0 )
                                {
                                    // Feel::cout << "Qth[" << name << "] " << std::flush;M_isNonLinear
                                    auto material = M_matProps[name];
                                    auto sigma0 = material.getDouble("sigma0");

                                    if ( M_isNonLinear )
                                        {
                                            auto alpha = material.getDouble("alpha");
                                            auto T0 = material.getDouble("T0");
                                            auto sigma = material.getScalar("sigma", {"T"}, {idv(M_T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
                                            M_Q += vf::project(  _space=M_XhQ, _range=markedelements(M_mesh_V, name), _expr=sigma*inner(gradv(M_V)) );
                                        }
                                    else
                                        M_Q += vf::project(  _space=M_XhQ, _range=markedelements(M_mesh_V, name), _expr=cst(sigma0)*inner(gradv(M_V)) ); 
                                    // Feel::cout << "\n" << std::flush;
                                }
                        }

                    e_V->step(time)->add( "joules", M_Q);
                }
            if ( field == "potential" )
                e_V->step(time)->add( "potential", M_V);
            if ( M_isNonLinear && field == "conductivity" )
                {
                    auto M_XhS = Pdh<OrderV>(M_mesh_V);
                    auto M_S = M_XhS->element();
                    for( auto const& pairMat : M_matProps )
                        {
                            auto marker = pairMat.first;
                            if ( isConductor(marker) )
                                {
                                    auto material = pairMat.second;
                                    auto alpha = material.getDouble("alpha");
                                    auto T0 = material.getDouble("T0");
                                    auto sigma0 = material.getDouble("sigma0");

                                    auto sigma = material.getScalar("sigma", {"T"}, {idv(M_T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
                                    M_S += vf::project( _space=M_XhS, _range=markedelements(M_mesh_V, marker), _expr=sigma );
                                }
                        }
                    e_V->step(time)->add( "sigma", M_S);
                }
            // need first to adapt estimator class to new BCs
            std::string estimator_type  = soption(_name="thermoelectric.estimator_type");
            if ( field == "estimator" )
                {
                    // auto normL2 = normL2( elements(M_mesh_V), idv(M_V) );
                    // auto normH1 = normL2( elements(M_mesh_V), gradv(M_V) );

                    estimator estimator;
                    boost::tuple<double, double, p0_element_type> estimator_V;
                    if( estimator_type == "res")
                        {
                            estimator_V = estimator.residual_estimator_V(M_V, M_mesh_V, M_modelProps, M_weakdir, M_penaldir);
                            e_V->step(time)->add( "estim_err", estimator_V.template get<2>());

                            Feel::cout << "\n ****** Residual estimators ****** \n";
                            Feel::cout << "            L2 norm\tH1 norm\n";
                            Feel::cout << "Potential   "
                                       << estimator_V.template get<0>() << "\t"
                                       << estimator_V.template get<1>() << "\n";
                            Feel::cout << " *********************************************" << std::endl;
                        }
                    else if( estimator_type == "zz" && OrderV==1 )
                        {
                            estimator_V = estimator.zz_estimator(M_V, M_mesh_V);
                            e_V->step(time)->add( "estim_err", estimator_V.template get<2>());
                            Feel::cout<< "\n ****** ZZ estimators ****** \n";
                            Feel::cout << "            L2 norm\tH1 norm\n";
                            Feel::cout << "Potential   "
                                       << estimator_V.template get<0>() << "\t"
                                       << estimator_V.template get<1>() << "\n";
                                       // << math::abs( estimator_V.template get<0>() - normL2  )/normL2 << "\t"
                                       // << math::abs( estimator_V.template get<1>() - normH1  )/normH1 << "\n";
                            Feel::cout << " *********************************************" << std::endl;
                    }
                    else if( estimator_type == "anisotropic" && OrderV==1 )
                        {
#if defined( FEELPP_HAS_GMSH_ADAPT_H )
                            MeshAdapt mesh_adaptation;

                            std::vector<vectorN_type> measures = mesh_adaptation.measures();
                            std::vector<matrixN_type> directions = mesh_adaptation.directions();
                            boost::tuple<double, p0_element_type> estimator_V = estimator.anisotropic_estimator_V(M_V, M_mesh_V, M_matProps, M_modelProps, M_weakdir, M_penaldir, measures, directions);

                            Feel::cout << "\n ****** Anisotropic estimator ****** \n"
                                       << "              H1 semi-norm\n"
                                       << "Potential   " << estimator_V.template get<0>() << "\n"
                                       << " *********************************************** \n";
#else
                            Feel::cout << "ThermoElectricModel: Anisotropic estimator not available - requires gmsh to be patched \n";
#endif
                    }
                }
        }
    }

    e_V->save();
}

#endif
