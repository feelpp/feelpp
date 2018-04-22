/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Author(s): Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   Date: 2011-16-12

   Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)
   Copyright (C) CNRS

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file magnetostatic.cpp
   \author Cecile Daversin <daversin@math.unistra.fr>
   \date 2014-26-08
*/

#ifndef __MAGNETOSTATIC_HPP
#define __MAGNETOSTATIC_HPP 1

#ifndef FM_DIM
#define FM_DIM FEELPP_DIM
#endif

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>

#include <boost/spirit/include/qi.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/regex.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/range/join.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>

#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>

#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelpde/preconditionerblockms.hpp>
#include <feel/feelcore/removecomments.hpp>
#include <feel/feelcore/utility.hpp>

/** for high order visualisation **/
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>

/** include general hifimagnet applications options **/
#include <feel/feelmodels/hifimagnet/options.hpp>

/** include needed tools **/
#include <feel/feelmodels/hifimagnet/Tools/mesh_initializer.hpp>
#include <feel/feelmodels/hifimagnet/Tools/operators.hpp>

/** include Biot-Savart (TO BE MOVED) **/
#include <feel/feeldiscr/pchv.hpp>
#ifdef FEELPP_HAS_BIOTSAVART
#include <feel/feelmodels/hifimagnet/Magnetostatic/biot_savart.hpp>
#endif

#ifdef FEELPP_HAS_MAGNETTOOLS
/** include Bmap (TO BE MOVED) **/
#ifdef BMAP_IDF2
#include <feel/feelmodels/hifimagnet/Magnetostatic/bmap_new.hpp>
#else
#include <feel/feelmodels/hifimagnet/Magnetostatic/bmap_using_magnettools.hpp>
#endif
#endif


/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
Feel::po::options_description
MagnetostaticOptions()
{
    Feel::po::options_description magnetostatic_options("magnetostatic (linear) options");
    magnetostatic_options.add_options()
        ( "magnetostatic.A_guess_bmap",  Feel::po::value<bool>()->default_value( false ), "compute initial guess for  magnetic potential (A.m) using Bmap" )
        ( "magnetostatic.A_guess_filename", Feel::po::value<std::string>()->default_value( "" ), "initial guess for magnetic potential (h5 file)" )
        // for non-linear
        ("magnetostatic.eps-coeff", Feel::po::value<double>()->default_value( 0 ), "penalisation parameter in magnetostatic equation")
        ("magnetostatic.tolerance", Feel::po::value<double>()->default_value( 1e-14 ), "tolerance for nl computation (ferromagnetism)")
        ("magnetostatic.max_iterations", Feel::po::value<int>()->default_value( 1 ), "max iterations for nl computation (ferromagnetism)")
        // Dirichlet condition options
        ("magnetostatic.weakdir", Feel::po::value<bool>()->default_value( true ), "use weak Dirichlet condition" )
        ("magnetostatic.penaldir", Feel::po::value<double>()->default_value( 50 ), "penalisation parameter for the weak boundary Dirichlet formulation")
        // Print information on model
        ("magnetostatic.verbose", Feel::po::value<int>()->default_value( 1 ), "verbosity level")
        // Export results
        ("magnetostatic.export", Feel::po::value<bool>()->default_value( true ), "true = export results for magnetostatic")
        // json modelproperties
        ("magnetostatic.model_json", Feel::po::value<std::string>()->default_value( "model.json" ), "use json file for model properties")
        // background potential
        ( "magnetostatic.background_potential", Feel::po::value<bool>()->default_value( false ), "apply background potential")
        ;

    return magnetostatic_options
#ifdef FEELPP_HAS_BIOTSAVART
        .add( BiotSavartOptions() )
#endif
#ifdef FEELPP_HAS_MAGNETTOOLS
        .add( BMapOptions() )
#endif
        ;

}

inline
Feel::po::options_description
MagnetostaticOptionsLib()
{
    Feel::po::options_description magnetostatic_optionsLib("magnetostatic options");
    magnetostatic_optionsLib.add(Feel::backend_options("ms"))
        .add(Feel::blockms_options("blockms"))
        .add(Feel::ams_options("ams"));
    return magnetostatic_optionsLib;
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
template<int Dim, int Order, int G_order>
inline
Feel::AboutData
makeAboutMagnetostatic()
{
    Feel::AboutData about( (boost::format("magnetostatic_%1%D_P%2%_N%3%") % Dim % Order % G_order).str(),
                           "magnetostatic" ,
                           "0.1",
                           "nD(n=1,2,3) Potential and Temperature on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008-2009 Universite Joseph Fourier"
                           "Copyright (c) CNRS");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "");
    about.addAuthor("Cecile Daversin", "developer", "cecile.daversin@lncmi.cnrs.fr", "");
    about.addAuthor("Romain Hild", "developer", "romain.hild2@etu.unistra.fr", "");

    return about;

}

// Class magnetostatic
namespace Feel
{

    template<int Dim, int Order, int G_order>
    class magnetostatic
    {
    public:
        using self_type = magnetostatic<Dim, Order, G_order>;
        using self_ptrtype = boost::shared_ptr<self_type>;

        static const bool is_P1 = (Order==1);

        typedef Eigen::Matrix<double, Dim, 1> vectorN_type;
        typedef Eigen::Matrix<double, Dim, Dim> matrixN_type;

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
        typedef Simplex<FEELPP_DIM,G_order> convex_type;
        //! mesh type
        typedef Mesh<convex_type> mesh_type;
        //! mesh shared_ptr<> type
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

        //! Hcurl space
        typedef bases<Nedelec<0,NedelecKind::NED1> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type> space_type;
        typedef boost::shared_ptr<space_type> space_ptrtype;
        typedef typename space_type::element_type element_type;

        //! the Lagrange approximation function space type
        using lagrange_space_s_type = Pch_type<mesh_type,Order>; //P1 scalar space
        using lagrange_space_s_ptrtype = boost::shared_ptr<lagrange_space_s_type>;
        typedef typename lagrange_space_s_type::element_type lagrange_s_element_type;
        typedef std::map<std::string, lagrange_s_element_type> lagrange_s_map_type;

        using lagrange_space_v_type = Pchv_type<mesh_type,Order>; //P1 vectorial space
        using  lagrange_space_v_ptrtype = boost::shared_ptr<lagrange_space_v_type>;
        typedef typename lagrange_space_v_type::element_type lagrange_v_element_type;

        using P0_space_v_type = Pdhv_type<mesh_type,Order-1>; //P0 vectorial space
        using P0_space_v_ptrtype = boost::shared_ptr<P0_space_v_type>;
        typedef typename P0_space_v_type::element_type P0_v_element_type;

        //! Hdiv space
#if FEELPP_DIM==2
        typedef FunctionSpace<mesh_type, P0_v_type> B_space_type;
        typedef lagrange_space_s_ptrtype B_space_ptrtype;
        typedef lagrange_s_element_type B_element_type;
#else
        typedef bases<RaviartThomas<0> > RT_basis_type; //RT vectorial space
        typedef FunctionSpace<mesh_type, RT_basis_type> B_space_type;
        typedef boost::shared_ptr<B_space_type> B_space_ptrtype;
        typedef typename B_space_type::element_type B_element_type;
#endif
        //! Composite space Hcurl x Lagrange (for saddle point formulation)
        typedef bases< Nedelec<0,NedelecKind::NED1>, Lagrange<1,Scalar> > comp_basis_type;
        typedef FunctionSpace<mesh_type, comp_basis_type > comp_space_type;
        typedef boost::shared_ptr<comp_space_type> comp_space_ptrtype;
        typedef typename comp_space_type::element_type comp_element_type;

        typedef PreconditionerBlockMS<comp_space_type> blockms_preconditioner_type;
        typedef boost::shared_ptr<blockms_preconditioner_type> blockms_preconditioner_ptrtype;
        // standard preconditioner_type defined in feel/feelalg/

        typedef boost::function<void ( sparse_matrix_ptrtype& A,vector_ptrtype& F )> linearAssembly_function_type;

        //! the exporter factory type
        typedef Exporter<mesh_type> export_type;
        //! the exporter factory (shared_ptr<> type)
        typedef std::shared_ptr<export_type> export_ptrtype;

        //! Model properties type
        typedef ModelProperties model_prop_type;
        typedef std::shared_ptr<model_prop_type> model_prop_ptrtype;
        //! Material properties type
        typedef ModelMaterials material_map_type;
        //! Boundary conditions type
        typedef BoundaryConditions condition_map_type;

#ifdef FEELPP_HAS_BIOTSAVART
        // For biot_savart
        typedef BiotSavart<FEELPP_DIM, Order, G_order> biotsavart_type;
#endif

        typedef typename SolverLinear<double>::solve_return_type solve_return_type;

    private:

        int M_verbose;
        //! mesh characteristic size
        double M_meshSize; // Coil

        //! shape of the domain
        std::string M_shape;

        //! Geofiles input
        std::string M_geofile;
        std::string M_geo_depends;
        std::string M_geofile_path;

        //! Mesh adaptation method / type
        bool M_meshadapt;
        int M_meshadapt_itMax;
        double M_meshadapt_tol;
        std::string M_meshadapt_method;
        std::string M_meshadapt_type;

        /* mesh, pointers and spaces */
        mesh_ptrtype M_mesh;
        lagrange_space_s_ptrtype M_Xhl; //H1
        lagrange_space_v_ptrtype M_XhL; //H1 for vector
        space_ptrtype M_Xhn; //Hcurl
        B_space_ptrtype M_Xhr; // Hdiv
        comp_space_ptrtype M_Xh;

        bool M_weakdir;
        double M_penaldir;
        int M_max_iter;
        double M_tol;
        double M_eps;

        model_prop_ptrtype M_modelProps;
        material_map_type M_matProps;
        bool M_isNonLinear;

        bool M_isSaddlePoint;
        B_element_type M_B;
        element_type M_A;
        element_type M_A_bg;
        element_type M_AOld;
        // lagrange_s_element_type M_mu;
        // lagrange_s_map_type M_mu_map;

        backend_ptrtype M_backend;
        sparse_matrix_ptrtype M_a;
        vector_ptrtype M_r;

        void assembleVariationalForm();
        template<typename T> void updateBilinearForm( T );
        template<typename T> void updateLinearForm( T );
        template<typename T> void updateWeakBoundaryConditions( T );
        template<typename T> void updateStrongBoundaryConditions( T );

        void initOldA();
        void ferromagnetism();

    public:
        magnetostatic();
        static self_ptrtype New();
        void changeRepository();
        void init( mesh_ptrtype mesh = NULL );
        void solve();
        void exportResults(double time = 0 );

        mesh_ptrtype mesh() const { return M_mesh; }
        lagrange_space_s_ptrtype divergenceSpace() const { return M_Xhl; }
        space_ptrtype potentialSpace() const { return M_Xhn; }
        comp_space_ptrtype compositeSpace() const { return M_Xh; }

        B_space_ptrtype magneticFieldSpace() const { return M_Xhr; }
        B_element_type magneticField() const { return M_B; }
        // lagrange_s_map_type mu_Map() const { return M_mu_map; }

        element_type potentialField() const { return M_A; }
        model_prop_ptrtype modelProperties() const { return M_modelProps; }
        material_map_type materialProperties() const { return M_matProps; }
        bool isValidModel();
        bool isSaddlePoint() const { return M_isSaddlePoint; }
        bool isWeakDir() const { return M_weakdir; }
        double penalDir() const { return M_penaldir;}
        int nDof() const { return M_isSaddlePoint ? M_Xh->nDof() : M_Xhn->nDof(); }

        linearAssembly_function_type M_updateAssembly;
        linearAssembly_function_type M_updateRHS;
    };

    template<int Dim, int Order, int G_order> const bool magnetostatic<Dim, Order, G_order>::is_P1;
};

template<int Dim, int Order, int G_order>
Feel::magnetostatic<Dim, Order, G_order>::magnetostatic()
{
    M_verbose = ioption("magnetostatic.verbose");
    M_meshSize = doption( "hsize" );
    M_shape = soption( "shape" );
    M_geofile = soption( "geofile" );
    M_geo_depends = soption( "geo_depends" );
    M_geofile_path = soption( "geofile-path" );
    M_meshadapt = boption( "meshadapt" );
    M_meshadapt_itMax = ioption("meshadapt_maxiter");
    M_meshadapt_tol = doption("meshadapt_tol");
    M_meshadapt_method = soption( "meshadapt_method" );
    M_meshadapt_type = soption( "meshadapt_type" );
    M_weakdir = boption( "magnetostatic.weakdir" );
    M_penaldir = doption( "magnetostatic.penaldir" );
    M_tol = doption( "magnetostatic.tolerance" );
    M_max_iter = ioption( "magnetostatic.max_iterations" );
    M_eps = doption( "magnetostatic.eps-coeff" );
    M_modelProps = std::make_shared<model_prop_type>( Environment::expand( soption("magnetostatic.model_json") ) );
    if ( !isValidModel() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "Magnetostatic Model " << M_modelProps->model() << " not supported! Aborting" << std::endl;
        LOG(FATAL) << "Magnetostatic Model " << M_modelProps->model() << " not supported! Aborting";
    }
    M_isSaddlePoint = boost::icontains(M_modelProps->model(), "saddlepoint");
    M_isNonLinear = boost::icontains(M_modelProps->model(), "nonlinear");
    if ( M_isNonLinear )
        {
            M_tol = doption( "magnetostatic.tolerance" );
            M_max_iter = ioption( "magnetostatic.max_iterations" );
        }
    if ( M_isSaddlePoint && M_weakdir )
    {
        Feel::cout << "Weak Dirichlet conditions are not implemented in saddlepoint formulation! Aborting" << std::endl;
        LOG(FATAL) << "Weak Dirichlet conditions are not implemented in saddlepoint formulation! Aborting";
    }
    M_matProps = M_modelProps->materials();
}

template<int Dim, int Order, int G_order>
typename Feel::magnetostatic<Dim, Order, G_order>::self_ptrtype
Feel::magnetostatic<Dim, Order, G_order>::New()
{
    return boost::make_shared<self_type>();
}

template<int Dim, int Order, int G_order>
bool Feel::magnetostatic<Dim, Order, G_order>::isValidModel()
{
    Feel::cout << "MagnetostaticModel: isValidModel" << std::endl;

    bool b = false;
    b = b || boost::icontains(M_modelProps->model(),"magneto-saddlepoint")
        || boost::icontains(M_modelProps->model(),"magneto-regularized")
        || boost::icontains(M_modelProps->model(),"magneto-singular")
        || boost::icontains(M_modelProps->model(),"magnetostatic-saddlepoint")
        || boost::icontains(M_modelProps->model(),"magnetostatic-regularized")
        || boost::icontains(M_modelProps->model(),"magnetostatic-singular")
        || boost::icontains(M_modelProps->model(),"ferromagnetism-saddlepoint")
        || boost::icontains(M_modelProps->model(),"ferromagnetism-regularized")
        || boost::icontains(M_modelProps->model(),"coupled")
        || boost::icontains(M_modelProps->model(),"coupled-nonlinear");
    return b;
}


template<int Dim, int Order, int G_order>
void
Feel::magnetostatic<Dim, Order, G_order>::changeRepository( )
{
    /*** eventually change output dir***/
    std::string appname = "MagnetostaticModel";
    Environment::changeRepository( boost::format( "hifimagnet/%1%/%2%/%3%_P%4%_N%5%/" )
                                   % appname
                                   % M_geofile
                                   % M_modelProps->model()
                                   % Order
                                   % G_order );

}
/* Initialize mesh, spaces and elements and print info
 */
template<int Dim, int Order, int G_order>
void
Feel::magnetostatic<Dim, Order, G_order>::init( mesh_ptrtype mesh )
{
    tic();
    if ( !mesh )
    {
        mesh_initializer<Dim, G_order> mesh_init;
        M_mesh = mesh_init.initializeMesh(Environment::vm());
    } else
    {
        M_mesh = mesh;
    }
    toc("mesh", M_verbose > 1);

    // Loading conductors
    std::vector<std::string> conductor = vsoption(_name="conductor_volume");
    std::list<std::string> conductor_list;
    std::copy( conductor.begin(), conductor.end(), std::back_inserter( conductor_list ) );

    // check if conductor are properly declared as markers
    bool check_conductor_list = false;
    for( auto marker: M_mesh->markerNames() )
        {
            auto name = marker.first;
            auto data = marker.second;
            // Feel::cout << name << ", dim=" << data[1] << ", Dim=" << Dim << ",  M_mesh->dimension()=" <<  M_mesh->dimension() << std::endl;
            if ( data[1] == M_mesh->dimension() )
                {
                    check_conductor_list = check_conductor_list || std::find(std::begin(conductor_list), std::end(conductor_list), name) != std::end(conductor_list);
                    // Feel::cout << name << "[" ;
                    // Feel::cout << (std::find(std::begin(conductor_list), std::end(conductor_list), name) != std::end(conductor_list)) << ",";
                    // Feel::cout << (std::find(std::begin(conductor), std::end(conductor), name) != std::end(conductor)) << "]\n";
                }
        }
    // Feel::cout << std::endl;

    if (!check_conductor_list)
        {
            Feel::cout << "Magnetostatic Model : invalid conductor_list [";
            for (std::string item : conductor_list)
                cout << item << ",";
            Feel::cout << "]\n" << std::flush;
            LOG(FATAL) << "Magnetostatic Model : invalid conductor_list! Aborting";
        }

    tic();
    M_Xhl = lagrange_space_s_type::New( M_mesh );
    M_Xhn = space_type::New( M_mesh );
    M_Xh = comp_space_type::New( M_mesh ); //Hcurl x H1

#if FEELPP_DIM==3
    M_Xhr = B_space_type::New( M_mesh ); // Hdiv
#else
    M_Xhr = M_Xhl;
#endif

    M_B = M_Xhr->element();
    M_A = M_Xhn->element();
    toc("spaces", M_verbose > 1 );

    M_backend = backend(_name="ms", _rebuild=true);
    if ( M_isSaddlePoint )
    {
        M_a = M_backend->newMatrix( _test=M_Xh, _trial=M_Xh);
        M_r = M_backend->newVector( _test=M_Xh);
    }
    else
    {
        M_a = M_backend->newMatrix( _test=M_Xhn, _trial=M_Xhn);
        M_r = M_backend->newVector( _test=M_Xhn);
    }

    this->initOldA(); // should only be done once if ( M_meshadaptMethod != "no" )

    // Collect of options
    int proc_rank = Environment::worldComm().globalRank();

    if( M_verbose > 0)
    {
        Feel::cout << "model : " << M_modelProps->model() << std::endl;
        Feel::cout << "geofile name : " << M_geofile << std::endl;
        Feel::cout << "Number of dofs : " << M_Xh->nDof() << std::endl;
        Feel::cout << "Number of faces in mesh : " << M_mesh->numFaces() << std::endl;
        Feel::cout << "Number of edges in mesh : " << M_mesh->numEdges() << std::endl;
        Feel::cout << "Number of elts in mesh : " << M_mesh->numElements() << std::endl;

        if( M_weakdir )
            {
                Feel::cout << "Use weak Dirichlet conditions " << std::endl;
                Feel::cout << "penalisation coefficient : " << M_penaldir << "\n" << std::endl;
            }
        else
            Feel::cout << "Use strong Dirichlet conditions " << std::endl;

        if ( M_meshadapt )
            {
                Feel::cout << "Mesh adaptation method used : " << M_meshadapt_method  << std::endl;
                Feel::cout << "Mesh adaptation is : " << M_meshadapt_type  << std::endl;
                Feel::cout << "Mesh adaptation tolerance : " << M_meshadapt_tol << std::endl;
                Feel::cout << "iterations for mesh adaptation (max) : " << M_meshadapt_itMax << "\n" << std::endl;
            }
    }
}

/** Initialize old temperature and potential with initial guess or from a file */
template<int Dim, int Order, int G_order>
void
Feel::magnetostatic<Dim, Order, G_order>::initOldA()
{
    auto AInit = soption("magnetostatic.A_guess_filename");
    if ( !AInit.empty() )
    {
        if( M_verbose > 0 )
            cout << "Loading A_old from " << AInit << "\n";
        fs::path InitPath(AInit);
        if ( fs::exists( InitPath ) )
            M_AOld.loadHDF5( AInit );
        else
            throw std::logic_error( "No such file : " + AInit + " : magnetostatic.A_guess_filename" );
    }

    if ( boption("magnetostatic.A_guess_bmap") )
    {
        throw std::logic_error( "Not implemented yet : magnetostatic.A_guess_bmap" );
    }

    M_AOld = M_Xh->element(); // to init M_AOld to 0 vector
}

/* assemble the problem and solve it
 */
template<int Dim, int Order, int G_order>
void
Feel::magnetostatic<Dim, Order, G_order>::solve()
{
    Feel::cout << "Feel::magnetostatic solve:";
    Feel::cout << "M_isSaddlePoint=" << M_isSaddlePoint << ", ";
    Feel::cout << "Prec=" << soption("ms.pc-type") << std::endl;

    // saddlepoint formulation U=(u,p)
    auto U = M_Xh->element(); //magnetic potential (k)
    U = M_AOld;

    element_type u, u_old;
    // saddlepoint
    if ( M_isSaddlePoint )
        u = U.template element<0>();
    else
        u = M_Xhn->element();

    blockms_preconditioner_ptrtype blockms_prec;
    preconditioner_ptrtype ams_prec;
    if( M_isSaddlePoint && soption("ms.pc-type") == "blockms" )
        {
            Feel::cout << "BLOCKMS : preconditionner" << std::endl;
            blockms_prec = boost::make_shared<blockms_preconditioner_type>( M_Xh,*M_modelProps,"ms", M_a, M_eps );
        }
    if( !M_isSaddlePoint && soption("ms.pc-type") == "ams" )
        {
            Feel::cout << "AMS : preconditionner" << std::endl;
            Feel::cout << "M_eps=" << M_eps << "\n";
            Feel::cout << "ams.setAlphaBeta=" << boption(_name="ams.setAlphaBeta") << "\n";
            Feel::cout << "ams.useEdge=" << boption(_name="ams.useEdge") << std::endl;

            // Construct and hold the prec in a singleton. No need to define it here. It is identified with his prefix.
            auto prec = preconditioner(_pc=pcTypeConvertStrToEnum(soption("ms.pc-type")),
                                       _backend=backend(_name="ms"),
                                       _prefix="ms",
                                       _matrix=M_a);

            ams_prec = prec;
            auto Igrad = Grad( _domainSpace=M_Xhl, _imageSpace=M_Xhn);
            ams_prec->attachAuxiliarySparseMatrix("G",Igrad.matPtr());

            if(boption("ams.useEdge"))
                {
                    auto ozz = M_Xhn->element();
                    auto zoz = M_Xhn->element();
                    auto zzo = M_Xhn->element();
                    ozz.on(_range=elements(M_Xhn->mesh()),_expr=vec(cst(1),cst(0),cst(0)));
                    zoz.on(_range=elements(M_Xhn->mesh()),_expr=vec(cst(0),cst(1),cst(0)));
                    zzo.on(_range=elements(M_Xhn->mesh()),_expr=vec(cst(0),cst(0),cst(1)));
                    vector_ptrtype M_ozz = backend(_name="ms")->newVector(M_Xhn); *M_ozz = ozz; M_ozz->close();
                    vector_ptrtype M_zoz = backend(_name="ms")->newVector(M_Xhn); *M_zoz = zoz; M_zoz->close();
                    vector_ptrtype M_zzo = backend(_name="ms")->newVector(M_Xhn); *M_zzo = zzo; M_zzo->close();
                    ams_prec->attachAuxiliaryVector("Px",M_ozz);
                    ams_prec->attachAuxiliaryVector("Py",M_zoz);
                    ams_prec->attachAuxiliaryVector("Pz",M_zzo);
                }
            else
                {
                    auto x = M_Xhl->element();
                    auto y = M_Xhl->element();
                    auto z = M_Xhl->element();
                    x.on(_range=elements(M_Xhl->mesh()),_expr=Px());
                    y.on(_range=elements(M_Xhl->mesh()),_expr=Py());
                    z.on(_range=elements(M_Xhl->mesh()),_expr=Pz());
                    vector_ptrtype M_X = backend(_name="ms")->newVector(M_Xhl); *M_X = x; M_X->close();
                    vector_ptrtype M_Y = backend(_name="ms")->newVector(M_Xhl); *M_Y = y; M_Y->close();
                    vector_ptrtype M_Z = backend(_name="ms")->newVector(M_Xhl); *M_Z = z; M_Z->close();
                    ams_prec->attachAuxiliaryVector("X",M_X);
                    ams_prec->attachAuxiliaryVector("Y",M_Y);
                    ams_prec->attachAuxiliaryVector("Z",M_Z);
                }

            if ( M_eps == 0)
                {
                    Feel::cout << " set HYPRE_AMSSetBetaPoissonMatrix(solver, NULL)\n";
                    ams_prec->attachAuxiliarySparseMatrix("a_beta",NULL);
                }

            // to be checked
            if(boption(_name="ams.setAlphaBeta"))
                {
                    Feel::cout << " set HYPRE_AMSSetAlphaPoissonMatrix(solver, a_alpha)\n";
                    auto us = M_XhL->element();
                    auto vs = M_XhL->element();
                    auto a_alpha = form2(_test=M_XhL, _trial=M_XhL);
                    auto b_alpha = form1(_test=M_XhL);
                    for( auto const& pairMat : M_matProps )
                        {
                            auto marker = pairMat.first;
                            auto material = pairMat.second;
                            auto mu = material.getScalar("mu_mag");
                            a_alpha += integrate(_range=markedelements(M_XhL->mesh(),marker), _expr=1/mu*inner(gradt(us), grad(vs)) + M_eps*inner(idt(us),id(vs)));
                        }

                    auto itField = M_modelProps->boundaryConditions().find( "magnetic_potential");
                    if ( itField != M_modelProps->boundaryConditions().end() && !M_weakdir )
                        {
                            auto mapField = (*itField).second;
                            auto itType = mapField.find( "Dirichlet" );
                            if ( itType != mapField.end() )
                                {
                                    for ( auto const& exAtMarker : (*itType).second )
                                        {
                                            std::string marker = exAtMarker.marker();
                                            auto g = expr<FEELPP_DIM,1>(exAtMarker.expression());
                                            a_alpha += on(_range=markedfaces(M_XhL->mesh(),marker),_element=us, _expr=g, _rhs=b_alpha, _type="elimination_symmetric");
                                        }
                                }
#ifdef FEELPP_HAS_BIOTSAVART
                            itType = mapField.find( "Dirichlet_Biot" );
                            if ( itType != mapField.end() )
                                {
                                    for ( auto const& exAtMarker : (*itType).second )
                                        {
                                            std::string marker = exAtMarker.marker();
                                            std::string domain = exAtMarker.expression1();

                                            auto material = M_matProps[domain];
                                            auto mu = material.getScalar("mu_mag");
#if FEELPP_DIM == 2
                                            throw std::logic_error( "Computed 2D pbs Biot dirichlet BCs not implemented so far in AMSSetAlphaPoissonMatrix" );
#else
                                            auto subMesh = createSubmesh(M_mesh, markedfaces(M_mesh,marker));
                                            //auto Xh_bcs = Pchv<2>( subMesh );

                                            std::vector<std::string> conductor = vsoption(_name="conductor_volume");
                                            std::list<std::string> conductor_list;
                                            std::copy( conductor.begin(), conductor.end(), std::back_inserter( conductor_list ) );
                                            auto subMesh_cond = createSubmesh(M_mesh, markedelements(M_mesh, conductor_list));
                                            biotsavart_type biot_savart_model(Environment::vm(), M_mesh, marker);
                                            auto Xh_Js = P0_space_v_type::New( subMesh_cond );
                                            auto M_current = Xh_Js->element();

                                            std::vector<std::string> currentdensity = vsoption(_name="conductor_js");
                                            for (unsigned int i=0; i<conductor.size(); i++)
                                                {
                                                    if ( !currentdensity[i].empty() )
                                                        {
                                                            // get Js from expression
                                                            auto Js = expr<FEELPP_DIM,1>(currentdensity[i]);
                                                            M_current += vf::project(_space=Xh_Js, _range=markedelements(subMesh_cond, conductor[i]), _expr=Js);
                                                        }
                                                    else
                                                        {
                                                            // get Js from ???
                                                            Feel::cout << "Dirichlet_Biot[" << marker << "]: Js not know for " << conductor[i] << " : set to null" << std::endl;
                                                        }
                                                }

                                            auto A_any = biot_savart_model.biot_savart_A( _conductor_mesh=subMesh_cond,
                                                                                          _box_name=marker,
                                                                                          _current_density=M_current);
                                            auto A_BS = boost::any_cast<typename BiotSavartImpl<FEELPP_DIM,2,Order,G_order>::vec_box_element_type >( A_any );

                                            a_alpha += on( _range=markedfaces(M_mesh, marker), _rhs=b_alpha, _element=us, _expr=idv(A_BS), _type="elimination_symmetric");
                                            //throw std::logic_error( "Computed dirichlet BCs not implemented so far in AMSSetAlphaPoissonMatrix" );
#endif
                                        }
                                }
#endif
                        }
                    ams_prec->attachAuxiliarySparseMatrix("a_alpha",a_alpha.matrixPtr());
                    //a_alpha.matrixPtr()->printMatlab("Mag_Aalpha.m");

                    if (M_eps > 0.)
                        {
                            Feel::cout << " set HYPRE_AMSSetBetaPoissonMatrix(solver, a_beta)\n";
                            auto uu = M_Xhl->element();
                            auto a_beta = form2(_test=M_Xhl, _trial=M_Xhl);
                            auto b_beta = form1(_test=M_Xhl);
                            a_beta += integrate(_range=elements(M_Xhl->mesh()), _expr=M_eps*inner(grad(uu),gradt(uu)));

                            auto itField = M_modelProps->boundaryConditions().find( "magnetic_potential");
                            if ( itField != M_modelProps->boundaryConditions().end() && !M_weakdir )
                                {
                                    auto mapField = (*itField).second;
                                    auto itType = mapField.find( "Dirichlet" );
                                    if ( itType != mapField.end() )
                                        {
                                            for ( auto const& exAtMarker : (*itType).second )
                                                {
                                                    std::string marker = exAtMarker.marker();
                                                    a_beta += on(_range=markedfaces(M_Xhl->mesh(),marker),_element=uu, _expr=cst(0.), _rhs=b_beta, _type="elimination_symmetric");
                                                }
                                        }
#ifdef FEELPP_HAS_BIOTSAVART
                                    itType = mapField.find( "Dirichlet_Biot" );
                                    if ( itType != mapField.end() )
                                        {
                                            for ( auto const& exAtMarker : (*itType).second )
                                                {
                                                    std::string marker = exAtMarker.marker();
                                                    throw std::logic_error( "Computed dirichlet BCs not implemented so far in AMSSetAlphaPoissonMatrix" );
                                                }
                                        }
#endif
                                }
                            ams_prec->attachAuxiliarySparseMatrix("a_beta",a_beta.matrixPtr());
                            //a_beta.matrixPtr()->printMatlab("Mag_Abeta.m");
                        }
                }
        }

    // /******************************
    //  * Initialize mu_map with mu_0
    //  ******************************/
    // for( auto const& pairMat : M_matProps )
    // {
    //     auto marker = pairMat.first;
    //     auto material = pairMat.second;

    //     auto mu = material.getScalar("mu_mag");
    //     M_mu_map[marker] = vf::project( _space=M_Xhl, _range=markedelements(M_mesh,marker), _expr=mu );
    //     M_mu.add( M_mu_map[marker] );
    // }

    if( Environment::isMasterRank() && boost::icontains(M_modelProps->model(), "Ferromagnetism" ) )
        std::cout << "Iteration n 0" << std::endl;

    this->assembleVariationalForm();

    /***************************************************
     * Solve the problem using the precontionner or not
     ***************************************************/
    solve_return_type result;
    if ( M_isSaddlePoint )
    {
        tic();
        if(soption("ms.pc-type") == "blockms" )
            result = M_backend->solve(_matrix=M_a, _rhs=M_r, _solution=U, _prec=blockms_prec);
        else
            result = M_backend->solve(_matrix=M_a, _rhs=M_r, _solution=U);
        toc("saddle solve", M_verbose > 1);
    }
    else
    {
        tic();
        if(soption("ms.pc-type") == "ams" )
            result = M_backend->solve(_matrix=M_a, _rhs=M_r, _solution=u, _prec=ams_prec);
        else
            result = M_backend->solve(_matrix=M_a, _rhs=M_r, _solution=u);
        toc("regul solve", M_verbose > 1);
    }

    std::string msg = (boost::format("\n[MagnetostaticModel %1%] NbIter=%2% Residual=%3%\n") % soption("ms.pc-type") % result.nIterations() % result.residual()).str();
    if (result.isConverged())
        {
            Feel::cout << tc::green << msg << tc::reset << std::endl;
        }
    else
        {
            std::string errmsg = msg + " Failed to converge";
            throw std::logic_error( errmsg );
        }

    /***************************************
     * Get the potential and magnetic field
     ***************************************/
    if ( M_isSaddlePoint )
        M_A = U.template element<0>();
    else
        M_A = u;

#if defined(FEELPP_HAS_HDF5)
    tic();
    LOG(INFO) << " Saving results in HDF5 format (" << M_geofile << ")";
    std::string M_mshfile;
    if ( !M_geofile_path.empty() )
        M_mshfile = ( boost::format( "%1%/%2%" ) % Environment::expand( M_geofile_path ) % M_geofile ).str();
    fs::path M_mshfile_path(M_mshfile);
    if ( M_isSaddlePoint )
        M_A.saveHDF5(M_mshfile_path.stem().string()+"_A_phi.h5");
    else
        M_A.saveHDF5(M_mshfile_path.stem().string()+"_A.h5");
    LOG(INFO) << " HDF5 data saved";
    toc("magnetostatic [save A to hdf5]", M_verbose > 0);
#else
    cout << "Saving results in HDF5 not available - Feelpp is not been compiled with HDF5 support" << std::endl;
#endif

#if FEELPP_DIM == 2
    M_B = vf::project( _space=M_Xhl, _range=elements(M_mesh), _expr=curlxv(M_A) );
#else
    auto Icurl = Curl( _domainSpace=M_Xhn, _imageSpace=M_Xhr);
    M_B = Icurl(M_A);
#endif

    int proc_rank = Environment::worldComm().globalRank();

    auto A_m = minmax( _range=elements(M_mesh), _pset=_Q<4>(), _expr=inner(idv(M_A),idv(M_A)) );
    auto A_min = A_m.min();
    auto A_max = A_m.max();
    auto A_p_min = A_m.argmin();
    auto A_p_max = A_m.argmax();

    if( proc_rank == 0 )
        {
            std::cout << "\n MAGNETIC POTENTIAL VALUES :" << std::endl;
            Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "{", "}");
            std::cout << "A min : " << math::sqrt(A_min) << " T.m at " << A_p_min.format(CleanFmt) << "\n";
            std::cout << "A max : " << math::sqrt(A_max) << " T.m at " << A_p_max.format(CleanFmt) << "\n";
            std::cout << "******************************** \n" << std::endl;
        }

    auto B_m = minmax( _range=elements(M_mesh), _pset=_Q<4>(), _expr=inner(idv(M_B),idv(M_B)) );
    auto B_min = B_m.min();
    auto B_max = B_m.max();
    auto B_p_min = B_m.argmin();
    auto B_p_max = B_m.argmax();

    if( proc_rank == 0 )
        {
            std::cout << "\n MAGNETIC FIELD VALUES :" << std::endl;
            Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "{", "}");
            std::cout << "B min : " << math::sqrt(B_min) << " T at " << A_p_min.format(CleanFmt) << "\n";
            std::cout << "B max : " << math::sqrt(B_max) << " T at " << A_p_max.format(CleanFmt) << "\n";
            std::cout << "******************************** \n" << std::endl;
        }

    node_type pt(3);
    pt[0] = 0.; pt[1] = 0; pt[2] = 0;
    auto val = M_B(pt);
    auto bx = val(0,0,0); // evaluation de Bx
    auto by = val(1,0,0); // evaluation de By
    auto bz = val(2,0,0); // evaluation de Bz
    Feel::cout << "B(O) = {" << bx << "," << by << "," << bz << "}" << std::endl;

    if( boost::icontains(M_modelProps->model(), "ferromagnetism") )
        this->ferromagnetism();
}

template<int Dim, int Order, int G_order>
void
Feel::magnetostatic<Dim, Order, G_order>::assembleVariationalForm()
{
    tic();
    /**************************************************************
     * Compute linear form with the user or model defined function
     **************************************************************/
    if ( M_isSaddlePoint )
        this->updateLinearForm( M_Xh);
    else
        this->updateLinearForm( M_Xhn);

    /************************
     * Compute bilinear form
     ************************/
    if ( M_isSaddlePoint )
        this->updateBilinearForm( M_Xh);
    else
        this->updateBilinearForm( M_Xhn);
    /*************************
     * Add user defined terms
     *************************/
    if ( M_updateAssembly != NULL )
        this->M_updateAssembly( M_a, M_r );
    if ( M_updateRHS != NULL )
        this->M_updateRHS( M_a, M_r );

    // IF !USE_BIOT_SAVART
    /*******************************
     * add weak boundary conditions
     *******************************/
    if ( M_weakdir )
    {
        if ( M_isSaddlePoint )
            this->updateWeakBoundaryConditions( M_Xh);
        else
            this->updateWeakBoundaryConditions( M_Xhn);
    }

    /********************************************
     * add strong boundary conditions at the end
     ********************************************/
    if ( M_isSaddlePoint )
        this->updateStrongBoundaryConditions( M_Xh);
    else
        this->updateStrongBoundaryConditions( M_Xhn);

    toc("assemble var form", M_verbose > 1);
}

template<int Dim, int Order, int G_order>
template<typename T>
void
Feel::magnetostatic<Dim, Order, G_order>::updateBilinearForm( T Xh)
{
    tic();
    auto a = form2( _test=Xh, _trial=Xh, _matrix=M_a);
    auto u = M_Xhn->element();
    auto v = M_Xhn->element();

    // saddlepoint
    if ( typeid(T) == typeid(comp_space_ptrtype) )
    {
        auto phi = M_Xhl->element();
        auto psi = M_Xhl->element();
        a += integrate( elements(M_mesh), inner(trans(id(v)),gradt(phi)) );
        a += integrate( elements(M_mesh), inner(trans(idt(u)),grad(psi)) );
    }
    else
    {
        a += integrate( elements(M_mesh), M_eps*trans(idt(u))*id(v) );
    }

    for( auto const& pairMat : M_matProps )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;

        auto mu = material.getScalar("mu_mag");
        // if ( M_isNonLinear )
        //     {
        //         auto mu = material.getScalar("mu_mag"); // non linear expression here cf sigma in thermoelectricmodel
        //     }

#if FEELPP_DIM == 2
        a += integrate( markedelements(M_mesh, marker), 1/mu*curlxt(u)*curlx(v) );
#else
        a += integrate( markedelements(M_mesh, marker), 1/mu*trans(curlt(u))*curl(v) );
#endif
    }
    toc("update bilinear", M_verbose > 2);
}

template<int Dim, int Order, int G_order>
template<typename T>
void
Feel::magnetostatic<Dim, Order, G_order>::updateLinearForm( T Xh)
{
    tic();
    auto mu_0 = 4*boost::math::constants::pi<double>()*1e-7;
    if( soption(_name="units") == "mm" )
        {
            Feel::cout << "Units: mm -> change permeability\n" << std::flush;
            mu_0 *= 1000;
        }
    auto l = form1( _test=Xh, _vector=M_r);
    auto v = M_Xhn->element();

    std::vector<std::string> conductor = vsoption(_name="conductor_volume");
    std::vector<std::string> currentdensity = vsoption(_name="conductor_js");
    for (unsigned int i=0; i<conductor.size(); i++)
        {
            if ( currentdensity.size() != 0  && !currentdensity[i].empty() )
                {
                    auto Js = expr<Dim,1>( currentdensity[i] );
                    l += integrate( markedelements(M_mesh, conductor[i]), mu_0 * inner(Js,id(v)) );
                }
        }
    // auto itField = M_modelProps->boundaryConditions().find( "magnetic_potential");
    // if ( itField != M_modelProps->boundaryConditions().end() )
    // {
    //     auto mapField = (*itField).second;
    //     auto itType = mapField.find( "CurrentDensity" );
    //     if ( itType != mapField.end() )
    //     {
    //         for ( auto const& exAtMarker : (*itType).second )
    //         {
    //             auto marker = exAtMarker.marker();
    //             auto g = expr<FEELPP_DIM,1>(exAtMarker.expression());
    //             l += integrate( markedelements(M_mesh, marker), inner(g,id(v)) );
    //             // l += integrate( elements(M_mesh), mu_0*inner(g,id(v)) );
    //         }
    //     }
    // }
    toc("update linear", M_verbose > 2);
}

template<int Dim, int Order, int G_order>
template<typename T>
void
Feel::magnetostatic<Dim, Order, G_order>::updateWeakBoundaryConditions( T Xh)
{
    tic();
    auto a = form2( _test=Xh, _trial=Xh, _matrix=M_a);
    auto l = form1( _test=Xh, _vector=M_r);
    auto u = M_Xhn->element();
    auto v = M_Xhn->element();

    auto itField = M_modelProps->boundaryConditions().find( "magnetic_potential");
    if ( itField != M_modelProps->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                std::string domain = exAtMarker.expression2();

                auto material = M_matProps[domain];
                auto mu = material.getScalar("mu_mag");
#if FEELPP_DIM == 2
                auto g = expr(exAtMarker.expression1());
                l += integrate(markedfaces(M_mesh,marker), 1/mu*curlx(v)*g
                               + (cst(1.)/idv(mu))*M_penaldir*trans(g)*cross(id(v),N())/hFace() );
                a += integrate(markedfaces(M_mesh,marker), 1/mu*curlxt(u)*(cross(id(v),N()) )
                               + 1/mu*curlx(v)*(cross(idt(u),N()) )
                               + 1/mu*M_penaldir*trans(cross(idt(u),N()))*cross(id(v),N())/hFace() );

#else
                auto g = expr<3,1>(exAtMarker.expression());
                l += integrate(markedfaces(M_mesh,marker), 1/mu*trans(curl(v))*g
                               + 1/mu*M_penaldir*trans(g)*cross(id(v),N())/hFace() );
                a += integrate(markedfaces(M_mesh,marker), 1/mu*trans(curlt(u))*(cross(id(v),N()) )
                               + 1/mu*trans(curl(v))*(cross(idt(u),N()) )
                               + 1/mu*M_penaldir*trans(cross(idt(u),N()))*cross(id(v),N())/hFace() );
#endif
            }
        }
#ifdef FEELPP_HAS_BIOTSAVART
        itType = mapField.find( "Dirichlet_Biot" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                std::string domain = exAtMarker.expression1();
                auto material = M_matProps[domain];

#if FEELPP_DIM == 2
                Feel::cout << "Dirichlet_Biot[" << marker << "] not implemented for 2D problems\n";
#else
                auto subMesh = createSubmesh(M_mesh, markedfaces(M_mesh,marker));
                //auto Xh_bcs = Pchv<2>( subMesh );

                std::vector<std::string> conductor = vsoption(_name="conductor_volume");
                std::list<std::string> conductor_list;
                std::copy( conductor.begin(), conductor.end(), std::back_inserter( conductor_list ) );
                auto subMesh_cond = createSubmesh(M_mesh, markedelements(M_mesh, conductor_list));
                biotsavart_type biot_savart_model(Environment::vm(), M_mesh, marker);
                auto Xh_Js = P0_space_v_type::New( subMesh_cond );
                auto M_current = Xh_Js->element();
                std::vector<std::string> currentdensity = vsoption(_name="conductor_js");
                if ( currentdensity.size() != 0 )
                    {
                        for (unsigned int i=0; i<conductor.size(); i++)
                            {
                                if ( !currentdensity[i].empty() )
                                    {
                                        // get Js from expression
                                        auto Js = expr<FEELPP_DIM,1>(currentdensity[i]);
                                        M_current += vf::project(_space=Xh_Js, _range=markedelements(subMesh_cond, conductor[i]), _expr=Js);
                                    }
                                else
                                    {
                                        // get Js from ???
                                        throw std::logic_error( "Dirichlet_Biot[" + marker + "]: Js not know for " + conductor[i] );
                                    }
                            }

                        auto A_any = biot_savart_model.biot_savart_A( _conductor_mesh=subMesh_cond,
                                                              _box_name=marker,
                                                              _current_density=M_current);
                        auto A_BS = boost::any_cast<typename BiotSavartImpl<FEELPP_DIM,2,Order,G_order>::vec_box_element_type >( A_any );

                        auto mu = material.getScalar("mu_mag");
                        l += integrate(markedfaces(M_mesh,marker), 1/mu*trans(curl(v))*cross(idv(A_BS),N())
                                       + 1/mu*M_penaldir*trans(cross(idv(A_BS),N()))*cross(id(v),N())/hFace() );
                        a += integrate(markedfaces(M_mesh,marker), 1/mu*trans(curlt(u))*(cross(id(v),N()) )
                                       + 1/mu*trans(curl(v))*(cross(idt(u),N()) )
                                       + 1/mu*M_penaldir*trans(cross(idt(u),N()))*cross(id(v),N())/hFace() );

                        // what to do for saddlepoint??
                        if ( typeid(T) == typeid(comp_space_ptrtype) )
                            {
                                auto phi = M_Xhl->element();
                                auto psi = M_Xhl->element();
                                a += integrate(markedfaces(M_mesh,marker), M_penaldir/hFace() * idt(phi)*id(psi) );
                            }
                    }
#endif
            }
        }
#endif
    }
    toc("update weak bc", M_verbose > 2);
}

template<int Dim, int Order, int G_order>
template<typename T>
void
Feel::magnetostatic<Dim, Order, G_order>::updateStrongBoundaryConditions( T Xh)
{
    tic();
    auto a = form2( _test=Xh, _trial=Xh, _matrix=M_a);
    auto l = form1( _test=Xh, _vector=M_r);
    auto u = M_Xhn->element();
    // saddlepoint
    if ( typeid(T) == typeid(comp_space_ptrtype) )
    {
        auto phi = M_Xhl->element();
        a += on( _range=boundaryfaces(M_mesh), _rhs=l, _element=phi, _expr=cst(0.));
    }

    auto itField = M_modelProps->boundaryConditions().find( "magnetic_potential");
    if ( itField != M_modelProps->boundaryConditions().end() && !M_weakdir )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr<FEELPP_DIM,1>(exAtMarker.expression());
                a += on( _range=markedfaces(M_mesh, marker), _rhs=l, _element=u, _expr=g );
            }
        }
#ifdef FEELPP_HAS_BIOTSAVART
        itType = mapField.find( "Dirichlet_Biot" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                std::string domain = exAtMarker.expression1();

                auto material = M_matProps[domain];
                auto mu = material.getScalar("mu_mag");
#if FEELPP_DIM == 2
                cout << "Dirichlet_Biot not implemented for 2D problems\n";
#else
                auto subMesh = createSubmesh(M_mesh, markedfaces(M_mesh,marker));
                //auto Xh_bcs = Pchv<2>( subMesh );

                std::vector<std::string> conductor = vsoption(_name="conductor_volume");
                std::list<std::string> conductor_list;
                std::copy( conductor.begin(), conductor.end(), std::back_inserter( conductor_list ) );
                auto subMesh_cond = createSubmesh(M_mesh, markedelements(M_mesh, conductor_list));
                biotsavart_type biot_savart_model(Environment::vm(), M_mesh, marker);
                auto Xh_Js = P0_space_v_type::New( subMesh_cond );
                auto M_current = Xh_Js->element();

                std::vector<std::string> currentdensity = vsoption(_name="conductor_js");
                for (unsigned int i=0; i<conductor.size(); i++)
                    {
                        if ( !currentdensity[i].empty() )
                            {
                                // get Js from expression
                                auto Js = expr<FEELPP_DIM,1>(currentdensity[i]);
                                M_current += vf::project(_space=Xh_Js, _range=markedelements(subMesh_cond, conductor[i]), _expr=Js);
                            }
                        else
                            {
                                // get Js from ???
                                Feel::cout << "Dirichlet_Biot[" << marker << "]: Js not know for " << conductor[i] << " : set to null" << std::endl;
                            }
                    }

                auto A_any = biot_savart_model.biot_savart_A( _conductor_mesh=subMesh_cond,
                                                              _box_name=marker,
                                                              _current_density=M_current);
                auto A_BS = boost::any_cast<typename BiotSavartImpl<FEELPP_DIM,2,Order,G_order>::vec_box_element_type >( A_any );

                a += on( _range=markedfaces(M_mesh, marker), _rhs=l, _element=u, _expr=idv(A_BS) ); //_type="elimination_symmetric"??

                // what to do for saddlepoint??
                if ( typeid(T) == typeid(comp_space_ptrtype) )
                    {
                        auto phi = M_Xhl->element();
                        auto psi = M_Xhl->element();
                        a += on(markedfaces(M_mesh,marker), _rhs=l, _element=phi, _expr=cst(0.) ); //_type="elimination_symmetric"??
                    }
#endif
            }
        }
#endif
    }
    toc("update strong bc", M_verbose > 2);
}

template<int Dim, int Order, int G_order>
void
Feel::magnetostatic<Dim, Order, G_order>::ferromagnetism()
{
    tic();
    /**************************************************
     * Loop for ferromagnetism, update mu at each step
     **************************************************/
    int iter=1;

    comp_element_type U = M_Xh->element();
    element_type u = M_Xhn->element(), u_old;
    U = M_AOld;

    double UmUold_L2 = 0;
    double Uold_L2 = 1;
    double relative_U_norm = 1e+30;

    blockms_preconditioner_ptrtype blockms_prec;
    preconditioner_ptrtype ams_prec;
    if( M_isSaddlePoint && soption("ms.pc-type") == "blockms" )
        {
            Feel::cout << "BLOCKMS : preconditionner\n";
            blockms_prec = boost::make_shared<blockms_preconditioner_type>( M_Xh,*M_modelProps,"ms", M_a, M_eps );
        }
    if( !M_isSaddlePoint && soption("ms.pc-type") == "ams" )
        {
            Feel::cout << "AMS : preconditionner\n";
            if( Environment::isMasterRank() == 0)
                {
                    std::cout << "M_eps=" << M_eps << "\n";
                    std::cout << "ams.setAlphaBeta=" << boption(_name="ams.setAlphaBeta") << "\n";
                    std::cout << "ams.useEdge=" << boption(_name="ams.useEdge") << "\n";
                }
            throw std::logic_error( "AMS precond: not implemented yet for ferromagnetism" );
        }

    while ( relative_U_norm > M_tol )
    {
        if( Environment::isMasterRank() )
            std::cout << "Iteration n " << iter << std::endl;

        // /******************************************
        //  * update value of mu (for each materials)
        //  ******************************************/
        // auto B = vf::project( _space=M_Xhn, _range=elements(M_mesh), _expr=curlv(M_A) );
        // u_old = vf::project( M_Xhn, elements(M_mesh), idv(M_A) );
        // M_mu.zero();

        // for( auto const& pairMat : M_matProps )
        // {
        //     auto marker = pairMat.first;
        //     auto material = pairMat.second;
        //     auto mu = material.getScalar("mu_mag"); // non linear expression here cf sigma in thermoelectricmodel
        //     // auto mu_0 = material.getScalar("mu_mag");
        //     // auto Bs = material.getScalar("Bs");
        //     // auto uri = material.getScalar("kappa_ri");
        //     // auto mu = M_mu_map[marker]; //element_type
        //     // auto normH = vf::project(_space=M_Xhl, _range=elements(M_mesh), _expr=(cst(1.)/idv(mu))*inner(idv(B),idv(B)) );
        //     // if( material.getString("Bs") != "0.0" && material.getString("Bs") != "0." ) // Bs = 0 : this material is no ferromagnetic
        //     //     M_mu_map[marker] =
        //     //         vf::project( _space=M_Xhl, _range=markedelements(M_mesh,marker),
        //     //                      _expr=mu_0 + ( 2*Bs/(M_PI*idv(normH)) )*atan( (M_PI*(uri-1)*mu_0*idv(normH))/(2*Bs)  ));
        //     // M_mu += M_mu_map[marker];
        // }

        /***********************************
         * Assemble variational formulation
         ***********************************/
        this->assembleVariationalForm();

        /***************************************************
         * Solve the problem using the precontionner or not
         ***************************************************/
        solve_return_type result;
        if ( M_isSaddlePoint )
        {
            tic();
            if(soption("ms.pc-type") == "blockms" )
                M_backend->solve(_matrix=M_a, _rhs=M_r, _solution=U, _prec=blockms_prec);
            else
                M_backend->solve(_matrix=M_a, _rhs=M_r, _solution=U);
            toc("magneto_solve",FLAGS_v>0);
        }
        else
        {
            tic();
            if(soption("ms.pc-type") == "ams" )
                M_backend->solve(_matrix=M_a, _rhs=M_r, _solution=u, _prec=ams_prec);
            else
                M_backend->solve(_matrix=M_a, _rhs=M_r, _solution=u);
            toc("magneto_solve",FLAGS_v>0);
        }

        std::string msg = (boost::format("[MagnetostaticModel %1%] NbIter=%2% Residual=%3%\n") % soption("ms.pc-type") % result.nIterations() % result.residual()).str();
        if (result.isConverged())
            {
                Feel::cout << tc::green << msg << tc::reset << std::endl;
            }
        else
            {
                std::string errmsg = msg + " Failed to converge";
                throw std::logic_error( errmsg );
            }

        /******************************************
         * Update the potential and magnetic field
         ******************************************/
        if ( M_isSaddlePoint )
            M_A = U.template element<0>();
        else
            M_A = u;

#if defined(FEELPP_HAS_HDF5)
        tic();
        LOG(INFO) << " Saving results in HDF5 format (" << M_geofile << ")";
        std::string M_mshfile;
        if ( !M_geofile_path.empty() )
            M_mshfile = ( boost::format( "%1%/%2%" ) % Environment::expand( M_geofile_path ) % M_geofile ).str();
        fs::path M_mshfile_path(M_mshfile);
        if ( M_isSaddlePoint )
            M_A.saveHDF5(M_mshfile_path.stem().string()+"_A_phi.h5");
        else
            M_A.saveHDF5(M_mshfile_path.stem().string()+"_A.h5");
        LOG(INFO) << " HDF5 data saved";
        toc("magnetostatic [save A to hdf5]", M_verbose > 0);
#else
        cout << "Saving results in HDF5 not available - Feelpp is not been compiled with HDF5 support" << std::endl;
#endif

#if FEELPP_DIM == 2
        M_B = vf::project( _space=M_Xhl, _range=elements(M_mesh), _expr=curlxv(M_A) );
#else
        auto Icurl = Curl( _domainSpace=M_Xhn, _imageSpace=M_Xhr);
        M_B = Icurl(M_A);
#endif

        // Relative error calculation (u)
        //TODO : fix for saddle /stab

        UmUold_L2 = normL2( elements(M_mesh), idv(M_A) - idv(u_old) );
        Uold_L2 = normL2( elements(M_mesh), idv(u_old) );


        // Relative norm : ||U(k) - U(k-1)||_L2 / ||U(k-1)||_L2
        relative_U_norm = UmUold_L2/Uold_L2;

        if( Environment::worldComm().globalRank() == 0)
        {
            std::cout << "Relative error = " << relative_U_norm << std::endl;
        }

        iter++;
    } // while ferromagnetism
    toc("ferromagnetism", M_verbose > 1);
}

template<int Dim, int Order, int G_order>
void
Feel::magnetostatic<Dim, Order, G_order>::exportResults(double time)
{
    tic();
    auto e = exporter( _mesh=M_mesh,  _name="Magnetostics");

    auto postProcess = M_modelProps->postProcess();
    auto itField = postProcess.find( "Fields");
    if ( itField != postProcess.end() )
        {
            for ( auto const& field : (*itField).second )
            {
                if ( field == "magneticField" )
                    e->step(time)->add( "magneticField", M_B);
                if ( field == "magneticPotential" )
                    e->step(time)->add( "magneticPotential", M_A);
                if ( field == "currentDensity" )
                    {
                        auto M_XhJ = Pdhv<Order>(M_mesh);
                        auto M_J = M_XhJ->element();

                        std::vector<std::string> conductor = vsoption(_name="conductor_volume");
                        std::vector<std::string> currentdensity = vsoption(_name="conductor_js");
                        for (unsigned int i=0; i<conductor.size(); i++)
                            {
                                if ( currentdensity.size() && !currentdensity[i].empty() )
                                    {
                                        auto Js = expr<Dim,1>( currentdensity[i] );
                                        M_J += vf::project( _space=M_XhJ, _range=markedelements(M_mesh, conductor[i]), _expr=Js );
                                    }
                            }
                        e->step(time)->add( "currentDensity", M_J);
                    }
                if ( M_isNonLinear && field == "permeability" )
                    {
                        auto M_Xhmu = Pdh<Order>(M_mesh);
                        auto M_mu = M_Xhmu->element();
                        for( auto const& pairMat : M_matProps )
                            {
                                auto marker = pairMat.first;
                                Feel::cout << "mu[" << marker << "] " << std::flush;
                                auto material = pairMat.second;
                                auto mu = material.getScalar("mu_mag");

                                M_mu += vf::project( _space=M_Xhmu, _range=markedelements(M_mesh, marker), _expr=mu );
                                Feel::cout << "\n" << std::flush;
                            }
                        e->step(time)->add( "permeability", M_mu);
                    }
            }
        }

    e->save();
    LOG(INFO) << "exportResults done\n";
    toc("export", M_verbose > 1);
}

#endif /* __MAGNETOSTATIC_HPP 1 */
