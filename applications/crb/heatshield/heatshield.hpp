/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-13

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file heatshield.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-03-28
 */
#ifndef __HeatShield_H
#define __HeatShield_H 1

#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelalg/solvereigen.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/modelcrbbase.hpp>

#include <feel/feelts/bdf.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <feel/feelcore/pslogger.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>



namespace Feel
{

po::options_description
makeHeatShieldOptions()
{
    po::options_description heatshieldoptions( "HeatShield options" );
    heatshieldoptions.add_options()
    // mesh parameters
    ( "hsize", Feel::po::value<double>()->default_value( 1e-1 ), "first h value to start convergence" )
    ( "mshfile", Feel::po::value<std::string>()->default_value( "" ), "name of the gmsh file input")
    ( "do-not-use-operators-free", Feel::po::value<bool>()->default_value( true ), "never use operators free if true" )
    ( "load-mesh-already-partitioned", Feel::po::value<bool>()->default_value( "true" ), "load a mesh from mshfile that is already partitioned if true, else the mesh loaded need to be partitioned")
    ( "beta.A0", Feel::po::value<std::string>()->default_value( "" ), "expression of beta coefficients for A0" )
    ( "beta.A1", Feel::po::value<std::string>()->default_value( "" ), "expression of beta coefficients for A1" )
    ( "beta.A2", Feel::po::value<std::string>()->default_value( "" ), "expression of beta coefficients for A2" )
    ( "beta.F0.0", Feel::po::value<std::string>()->default_value( "" ), "expression of beta coefficients for F0" )
    ( "beta.F1.0", Feel::po::value<std::string>()->default_value( "" ), "expression of beta coefficients for F1" )
    ( "beta.M0", Feel::po::value<std::string>()->default_value( "" ), "expression of beta coefficients for M0" )
    ;
    return heatshieldoptions.add( Feel::feel_options() ).add( bdf_options( "heatshield" ) ).add( backend_options("backendl2") );
}
AboutData
makeHeatShieldAbout( std::string const& str = "heatShield" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "heat shield Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010,2011 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;
}


class ParameterDefinition
{
public :
    static const uint16_type ParameterSpaceDimension = 2;
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
};

template<int Order>
class FunctionSpaceDefinition
{
public :
    //static const uint16_type Order = 1;
    typedef double value_type;

    /*mesh*/
    typedef Simplex<2,1> entity_type; /*dim,order*/
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    static const bool is_time_dependent = true;
    static const bool is_linear = true;

};


/**
 * \class HeatShield
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
template<int Order>
class HeatShield : public ModelCrbBase< ParameterDefinition, FunctionSpaceDefinition<Order> >,
                   public boost::enable_shared_from_this< HeatShield<Order> >
{
public:


    typedef ModelCrbBase<ParameterDefinition,FunctionSpaceDefinition<Order> > super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;
    typedef typename super_type::beta_vector_light_type beta_vector_light_type;
    typedef typename super_type::affine_decomposition_light_type affine_decomposition_light_type;
    //operator free
    typedef typename super_type::operator_type operator_type;
    typedef typename super_type::operator_ptrtype operator_ptrtype;
    typedef typename super_type::operatorcomposite_type operatorcomposite_type;
    typedef typename super_type::operatorcomposite_ptrtype operatorcomposite_ptrtype;
    typedef typename super_type::functionalcomposite_type functionalcomposite_type;
    typedef typename super_type::functionalcomposite_ptrtype functionalcomposite_ptrtype;
    typedef typename super_type::functional_type functional_type;
    typedef typename super_type::functional_ptrtype functional_ptrtype;
    using super_type::computeBetaQm;
    /** @name Constants
     */
    //@{

    //static const uint16_type Order = 1;
    static const uint16_type ParameterSpaceDimension = 2;

    //@}

    /** @name Typedefs
     */
    //@{

    typedef double value_type;

    typedef typename FunctionSpaceDefinition<Order>::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename FunctionSpaceDefinition<Order>::basis_type basis_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    /*space*/
    typedef typename FunctionSpaceDefinition<Order>::space_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /*reduced basis space*/
    typedef ReducedBasisSpace<super_type, mesh_type, basis_type, value_type> rbfunctionspace_type;
    typedef boost::shared_ptr< rbfunctionspace_type > rbfunctionspace_ptrtype;

    /* parameter space */
    typedef ParameterDefinition::parameterspace_type parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;

    /* time discretization */
    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0, Scalar, Continuous > > > continuous_p0_space_type;
    typedef boost::shared_ptr<continuous_p0_space_type> continuous_p0_space_ptrtype;
    typedef typename continuous_p0_space_type::element_type continuous_p0_element_type;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0, Scalar> >, Discontinuous> p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    HeatShield();

    //! constructor from command line
    HeatShield( po::variables_map const& vm );


    //! copy constructor
    HeatShield( HeatShield const & );
    //! destructor
    ~HeatShield() {}

    //! initialization of the model
    void initModel();
    //@}

    /** @name Operator overloads
     */
    //@{

    //@}

    /** @name Accessors
     */
    //@{

    // \return the number of terms in affine decomposition of left hand
    // side bilinear form
    int Qa() const
    {
        return 3;
    }

    // \return the number of terms in affine decomposition of bilinear form
    // associated to mass matrix
    int Qm() const
    {
        return 1;
    }

    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     *
     * \return number of outputs associated to the model
     * in our case we have a compliant output and 2 others outputs : average temperature on boundaries
     */
    int Nl() const
    {
        return 2;
    }

    /**
     * \param l the index of output
     * \return number of terms  in affine decomposition of the \p q th output term
     * in our case no outputs depend on parameters
     */
    int Ql( int l ) const
    {
        return 1;
    }

    int mMaxA( int q )
    {
        if ( q < 3 )
            return 1;
        else
            throw std::logic_error( "[Model heatshield] ERROR : try to acces to mMaxA(q) with a bad value of q");
    }

    int mMaxM( int q )
    {
        if ( q < 1 )
            return 1;
        else
            throw std::logic_error( "[Model heatshield] ERROR : try to acces to mMaxM(q) with a bad value of q");
    }

    int mMaxF( int output_index, int q)
    {
        if ( q < 1 )
            return 1;
        else
            throw std::logic_error( "[Model heatshield] ERROR : try to acces to mMaxF(output_index,q) with a bad value of q");
    }


    /**
     * \brief Returns the function space
     */
    space_ptrtype functionSpace()
    {
        return Xh;
    }

    /**
     * \brief Returns the reduced basis function space
     */
    rbfunctionspace_ptrtype rBFunctionSpace()
    {
        return RbXh;
    }

    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
    {
        return M_Dmu;
    }

    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<beta_vector_light_type, beta_vector_light_type, std::vector<beta_vector_light_type>  >
    computeBetaQ( parameter_type const& mu , double time=1e30 )
    {

        if( M_use_ginac )
        {
            ginac_expressionA[0].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } } );
            auto projection=project(_space=continuous_p0, _expr=ginac_expressionA[0]);
            M_betaAq[0] = projection(0);
            ginac_expressionA[1].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } } );
            projection=project(_space=continuous_p0, _expr=ginac_expressionA[1]);
            M_betaAq[1] = projection(0);
            ginac_expressionA[2].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } } );
            projection=project(_space=continuous_p0, _expr=ginac_expressionA[2]);
            M_betaAq[2] = projection(0);

            projection=project(_space=continuous_p0, _expr=ginac_expressionM[0]);
            M_betaMq[0] = projection(0);

            ginac_expressionF[0].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } , { "surface", surface } } );
            projection=project(_space=continuous_p0, _expr=ginac_expressionF[0]);
            M_betaFq[0][0] = projection(0);
            ginac_expressionF[1].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } , { "surface", surface } } );
            projection=project(_space=continuous_p0, _expr=ginac_expressionF[1]);
            M_betaFq[1][0] = projection(0);

#if 0
            int idx=0;
            int nl = Nl();
            for(int i=0; i<nl; i++)
            {
                int ql=Ql(i);
                for(int j=0; j<ql; j++)
                {
                    ginac_expressionF[idx].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } , { "surface", surface } } );
                    auto projection=project(_space=continuous_p0, _expr=ginac_expressionF[idx]);
                    M_betaFq[i][j] = projection(0);
                    idx++;
                }
            }
#endif
#if 0
            LOG( INFO ) << "mu = "<<mu(0)<<" -- "<<mu(1);
            LOG( INFO ) <<"A0 : "<<M_betaAq[0]<<" -- should be 1";
            LOG( INFO ) <<"A1 : "<<M_betaAq[1]<<" -- should be "<<mu(0);
            LOG( INFO ) <<"A2 : "<<M_betaAq[2]<<" -- should be "<<mu(1);
            LOG( INFO ) <<"F0 : "<<M_betaFq[0][0]<<" -- should be "<<mu(0);
            LOG( INFO ) <<"F1 : "<<M_betaFq[1][0]<<" -- should be "<<1./surface;
            LOG( INFO ) <<"M0 : "<<M_betaMq[0]<<" -- should be 1";
#endif
        }//use ginac
        else
        {
            double biot_out   = mu( 0 );
            double biot_in    = mu( 1 );
            M_betaAq.resize( Qa() );
            M_betaAq[0] = 1 ;
            M_betaAq[1] = biot_out ;
            M_betaAq[2] = biot_in  ;

            M_betaMq.resize( Qm() );
            M_betaMq[0] = 1;

            M_betaFq.resize( Nl() );
            M_betaFq[0].resize( Ql(0) );
            M_betaFq[0][0] = biot_out;

            M_betaFq[1].resize( Ql(1) );
            M_betaFq[1][0] = 1./surface;
        }

        return boost::make_tuple( M_betaMq, M_betaAq, M_betaFq );
    }


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * run the convergence test
     */


    /**
     * \brief Returns the affine decomposition
     */
    affine_decomposition_light_type computeAffineDecompositionLight(){ return boost::make_tuple( M_Mq, M_Aq, M_Fq ); }

    void assemble();

    /**
     * inner product
     */
    sparse_matrix_ptrtype energyMatrix ( void )
    {
        return M;
    }

    /**
     * inner product for mass matrix
     */
    sparse_matrix_ptrtype massMatrix ( void )
    {
        return InnerMassMatrix;
    }

#if 0
    /**
     * inner product for the POD
     */
    sparse_matrix_ptrtype innerProductForPod ( void )
    {
        return Mpod;
    }
#endif
    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false, bool export_outputs=false );

    gmsh_ptrtype createGeo( double hsize );

    virtual operatorcomposite_ptrtype operatorCompositeLightA()
    {
        return M_compositeLightA;
    }
    virtual operatorcomposite_ptrtype operatorCompositeLightM()
    {
        return M_compositeLightM;
    }
    virtual std::vector< functionalcomposite_ptrtype > functionalCompositeLightF()
    {
        return M_compositeLightF;
    }

    bdf_ptrtype bdfModel(){ return M_bdf; }

    void initDataStructureForBetaCoeff();
    void buildGinacExpressions();

private:

    po::variables_map M_vm;

    backend_ptrtype backend;

    /* mesh parameters */
    double meshSize;

    int export_number;

    bool do_export;

    parameterspace_ptrtype M_Dmu;

    double surface;

    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;
    space_ptrtype Xh;
    continuous_p0_space_ptrtype continuous_p0;
    rbfunctionspace_ptrtype RbXh;

    sparse_matrix_ptrtype M,Mpod,InnerMassMatrix;

    std::vector<sparse_matrix_ptrtype> M_Aq;
    std::vector<sparse_matrix_ptrtype> M_Mq;
    std::vector<std::vector<vector_ptrtype> > M_Fq;

    std::vector<operator_ptrtype> M_Aq_free;
    std::vector<operator_ptrtype> M_Mq_free;
    std::vector<std::vector<functional_ptrtype> > M_Fq_free;

    operatorcomposite_ptrtype M_compositeLightA;
    operatorcomposite_ptrtype M_compositeLightM;
    std::vector< functionalcomposite_ptrtype > M_compositeLightF;

    beta_vector_light_type M_betaAq;
    beta_vector_light_type M_betaMq;
    std::vector<beta_vector_light_type> M_betaFq;

    bdf_ptrtype M_bdf;

    element_type u,v;

    bool M_use_ginac ;

    std::vector< Expr<GinacEx<2> > > ginac_expressionA;
    std::vector< Expr<GinacEx<2> > > ginac_expressionM;
    std::vector< Expr<GinacEx<2> > > ginac_expressionF;
};

template<int Order>
HeatShield<Order>::HeatShield()
    :
    backend( backend_type::build( BACKEND_PETSC ) ),
    meshSize( 2e-1 ),
    M_Dmu( new parameterspace_type )
{ }

template<int Order>
HeatShield<Order>::HeatShield( po::variables_map const& vm )
    :
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    meshSize( vm["hsize"].as<double>() ),
    M_Dmu( new parameterspace_type )
{
}

template<int Order>
void
HeatShield<Order>::initDataStructureForBetaCoeff()
{
    int qa = Qa();
    M_betaAq.resize( qa );
    int qm=Qm();
    M_betaMq.resize( qm );
    int nl = Nl();
    M_betaFq.resize( nl );
    for(int i=0; i<nl; i++)
    {
        int ql=Ql(i);
        M_betaFq[i].resize( ql );
    }

}

template<int Order>
void
HeatShield<Order>::buildGinacExpressions()
{

    int qa = Qa();
    std::vector< std::string > symbols_vec;
    symbols_vec.push_back("BiotOut");
    symbols_vec.push_back("BiotIn");
    symbols_vec.push_back("surface");
    for(int i=0; i<qa; i++)
    {
        std::string name = ( boost::format("beta.A%1%") %i ).str();
        std::string filename = ( boost::format("GinacA%1%") %i ).str();
        ginac_expressionA.push_back( expr( option(_name=name).template as<std::string>(), Symbols( symbols_vec ) , filename ) );
    }


    int qm=Qm();
    for(int i=0; i<qm; i++)
    {
        std::string name = ( boost::format("beta.M%1%") %i ).str();
        std::string filename = ( boost::format("GinacM%1%") %i ).str();
        ginac_expressionM.push_back( expr( option(_name=name).template as<std::string>(),  Symbols( symbols_vec ) , filename ) );
    }

    int nl = Nl();
    for(int i=0; i<nl; i++)
    {
        int ql=Ql(i);
        for(int j=0; j<ql; j++)
        {
            std::string name = ( boost::format("beta.F%1%.%2%") %i %j ).str();
            std::string filename = ( boost::format("GinacF%1%.%2%") %i %j ).str();
            ginac_expressionF.push_back( expr( option(_name=name).template as<std::string>(), Symbols( symbols_vec ) , filename ) );
        }
    }

}

template<int Order>
gmsh_ptrtype
HeatShield<Order>::createGeo( double hsize )
{
    gmsh_ptrtype gmshp( new Gmsh );
    std::ostringstream ostr;
    double H = hsize;
    double h = hsize*0.5;
    //double h = hsize*1;
    ostr <<"Point (1) = {0,  0, 0, "<<H<<"};\n"
         <<"Point (2) = {10, 0, 0, "<<H<<"};\n"
         <<"Point (3) = {10, 4, 0, "<<H<<"};\n"
         <<"Point (4) = {0,  4, 0, "<<H<<"};\n"
         <<"Point (10) = {1, 1, 0, "<<h<<"};\n"
         <<"Point (11) = {3, 1, 0, "<<h<<"};\n"
         <<"Point (12) = {3, 3, 0, "<<h<<"};\n"
         <<"Point (13) = {1, 3, 0, "<<h<<"};\n"
         <<"Point (20) = {4, 1, 0, "<<h<<"};\n"
         <<"Point (21) = {6, 1, 0, "<<h<<"};\n"
         <<"Point (22) = {6, 3, 0, "<<h<<"};\n"
         <<"Point (23) = {4, 3, 0, "<<h<<"};\n"
         <<"Point (30) = {7, 1, 0, "<<h<<"};\n"
         <<"Point (31) = {9, 1, 0, "<<h<<"};\n"
         <<"Point (32) = {9, 3, 0, "<<h<<"};\n"
         <<"Point (33) = {7, 3, 0, "<<h<<"};\n"
         <<"Line (101) = {1,2};\n"
         <<"Line (102) = {2,3};\n"
         <<"Line (103) = {3,4};\n"
         <<"Line (104) = {4,1};\n"
         <<"Line (110) = {10,11};\n"
         <<"Line (111) = {11,12};\n"
         <<"Line (112) = {12,13};\n"
         <<"Line (113) = {13,10};\n"
         <<"Line (120) = {20,21};\n"
         <<"Line (121) = {21,22};\n"
         <<"Line (122) = {22,23};\n"
         <<"Line (123) = {23,20};\n"
         <<"Line (130) = {30,31};\n"
         <<"Line (131) = {31,32};\n"
         <<"Line (132) = {32,33};\n"
         <<"Line (133) = {33,30};\n"
         <<"Line Loop (201) = {101, 102, 103, 104};\n"
         <<"Line Loop (210) = {110, 111, 112, 113};\n"
         <<"Line Loop (220) = {120, 121, 122, 123};\n"
         <<"Line Loop (230) = {130, 131, 132, 133};\n"
         <<"Plane Surface (300) = {201,-210,-220,-230};\n"
         <<"Physical Line (\"left\") = {104};\n"
         <<"Physical Line (\"right\") = {102};\n"
         <<"Physical Line (\"bottom\") = {101};\n"
         <<"Physical Line (\"top\") = {103};\n"
         <<"Physical Line (\"gamma_holes\") = {110,111,112,113, 120,121,122,123, 130,131,132,133};\n"
         <<"Physical Surface (\"Omega\") = {300};\n"
         ;
    std::ostringstream nameStr;
    nameStr.precision( 3 );
    nameStr << "heatshield_geo";
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}

template<int Order>
void HeatShield<Order>::initModel()
{

    using namespace Feel::vf;

    M_use_ginac = option(_name="crb.use-ginac-for-beta-expressions").template as<bool>();

    std::string mshfile_name = option("mshfile").as<std::string>();

    /*
     * First we create the mesh or load it if already exist
     */

    if( mshfile_name=="" )
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                               _desc = createGeo( meshSize ) );
    }
    else
    {

        bool load_mesh_already_partitioned=option(_name="load-mesh-already-partitioned").template as<bool>();
        if( ! load_mesh_already_partitioned )
        {
            int N = Environment::worldComm().globalSize();
            std::string mshfile = option("mshfile").as<std::string>();
            std::string mshfile_complete = option("mshfile").as<std::string>();
            auto pos = mshfile.find(".msh");
            mshfile.erase( pos , 4);
            std::string filename = (boost::format(mshfile+"-np%1%.msh") %N ).str();
            if( !fs::exists( filename ) )
            {
                super_type::partitionMesh( mshfile_complete, filename , 2 , 1 );
            }
            mesh = loadGMSHMesh( _mesh=new mesh_type,
                                 _filename=filename,
                                 _rebuild_partitions=false,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
        }
        else
        {
            mesh = loadGMSHMesh( _mesh=new mesh_type,
                                 _filename=option("mshfile").as<std::string>(),
                                 _rebuild_partitions=false,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
        }
    }


    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    continuous_p0 = continuous_p0_space_type::New( mesh );
    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << Xh->nLocalDof() << "\n";
        std::cout << "Number of dof " << Xh->nDof() << "\n";
    }
    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";

    RbXh = rbfunctionspace_type::New( _model=this->shared_from_this() , _mesh=mesh );

    surface = integrate( _range=elements( mesh ), _expr=cst( 1. ) ).evaluate()( 0,0 );
    //std::cout<<"surface : "<<surface<<std::endl;

    M_bdf = bdf( _space=Xh, _vm=M_vm, _name="heatshield" , _prefix="heatshield" );

    bool dont_use_operators_free = option(_name="do-not-use-operators-free").template as<bool>() ;
    if( dont_use_operators_free )
    {
        M_Aq.resize( this->Qa() );
        M_Mq.resize( this->Qm() );
        M_Fq.resize( this->Nl() );
        for(int l=0; l<Nl(); l++)
        {
            M_Fq[l].resize( Ql(l) );
        }
    }
    else
    {
        M_Aq_free.resize( this->Qa() );
        M_Mq_free.resize( this->Qm() );
        M_Fq_free.resize( this->Nl() );
        for(int l=0; l<Nl(); l++)
        {
            M_Fq_free[l].resize( Ql(l) );
        }
    }

    initDataStructureForBetaCoeff();
    if( M_use_ginac )
        buildGinacExpressions();

    typename Feel::ParameterSpace<ParameterSpaceDimension>::Element mu_min( M_Dmu );
    mu_min <<  /* Bi_out */ 1e-2 , /*Bi_in*/1e-3;
    M_Dmu->setMin( mu_min );
    typename Feel::ParameterSpace<ParameterSpaceDimension>::Element mu_max( M_Dmu );
    mu_max << /* Bi_out*/0.5   ,  /*Bi_in*/0.1;
    M_Dmu->setMax( mu_max );

    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";

    assemble();


} // HeatShield::init


template<int Order>
void HeatShield<Order>::assemble()
{

    using namespace Feel::vf;

    u = Xh->element();
    v = Xh->element();

    if( option(_name="do-not-use-operators-free").template as<bool>() )
    {
        M_Aq[0] = backend->newMatrix( Xh, Xh );
        M_Aq[1] = backend->newMatrix( Xh, Xh );
        M_Aq[2] = backend->newMatrix( Xh, Xh );
        M_Mq[0] = backend->newMatrix( Xh, Xh );
        M_Fq[0][0] = backend->newVector( Xh );
        M_Fq[1][0] = backend->newVector( Xh );
        form2(Xh, Xh, M_Aq[0]) = integrate( _range= elements( mesh ), _expr= gradt( u )*trans( grad( v ) ) );
        form2(Xh, Xh, M_Aq[1]) = integrate( _range= markedfaces( mesh, "left" ), _expr= idt( u )*id( v ) );
        form2(Xh, Xh, M_Aq[2]) = integrate( _range= markedfaces( mesh, "gamma_holes" ), _expr= idt( u )*id( v ) );
        form2(Xh, Xh, M_Mq[0]) = integrate ( _range=elements( mesh ), _expr=idt( u )*id( v ) );
        form1(Xh, M_Fq[0][0]) = integrate( _range=markedfaces( mesh,"left" ), _expr= id( v ) ) ;
        form1(Xh, M_Fq[1][0]) = integrate( _range=elements( mesh ), _expr= id( v ) ) ;
    }
    else
    {
        auto expr_a0 = integrate( _range= elements( mesh ), _expr= gradt( u )*trans( grad( v ) ) );
        auto operatorfree0=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a0 , _backend=backend );
        operatorfree0->setName("A0");
        M_Aq_free[0]=operatorfree0;

        auto expr_a1 = integrate( _range= markedfaces( mesh, "left" ), _expr= idt( u )*id( v ) );
        auto operatorfree1=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a1 , _backend=backend );
        operatorfree1->setName("A1");
        M_Aq_free[1]=operatorfree1;

        auto expr_a2 = integrate( _range= markedfaces( mesh, "gamma_holes" ), _expr= idt( u )*id( v ) );
        auto operatorfree2=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a2 , _backend=backend );
        operatorfree2->setName("A2");
        M_Aq_free[2]=operatorfree2;

        auto expr_f00 = integrate( _range=markedfaces( mesh,"left" ), _expr= id( v ) ) ;
        auto functionalfree00 = functionalLinearFree( _space=Xh , _expr=expr_f00 , _backend=backend );
        auto expr_f10 = integrate( _range=elements( mesh ), _expr= id( v ) ) ;
        auto functionalfree10 = functionalLinearFree( _space=Xh , _expr=expr_f10 , _backend=backend );
        M_Fq_free[0][0]=functionalfree00;
        M_Fq_free[1][0]=functionalfree10;

        //mass matrix
        auto expr_m0 = integrate ( _range=elements( mesh ), _expr=idt( u )*id( v ) );
        auto operatorfreeM0=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_m0 , _backend=backend );
        operatorfreeM0->setName("mass");
        M_Mq_free[0]=operatorfreeM0;
    }

    //for scalarProduct
    M = backend->newMatrix( _test=Xh, _trial=Xh );
    form2( Xh, Xh, M ) =
        integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) ) +
        integrate( _range= markedfaces( mesh, "left" ), _expr= 0.01 * idt( u )*id( v ) ) +
        integrate( _range= markedfaces( mesh, "gamma_holes" ), _expr= 0.001 * idt( u )*id( v ) )
        ;

    //scalar product used for mass matrix
    InnerMassMatrix = backend->newMatrix( _test=Xh, _trial=Xh );
    form2( Xh, Xh, InnerMassMatrix ) =
        integrate( _range=elements( mesh ), _expr=idt( u ) * id( v ) ) ;

    //scalar product used for the POD
    Mpod = backend->newMatrix( _test=Xh, _trial=Xh );
    form2( Xh, Xh, Mpod ) =
        integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) ) +
        integrate( _range= markedfaces( mesh, "left" ), _expr= 0.01 * idt( u )*id( v ) ) +
        integrate( _range= markedfaces( mesh, "gamma_holes" ), _expr= 0.001 * idt( u )*id( v ) )
        ;

    if( ! option(_name="do-not-use-operators-free").template as<bool>() )
    {
        M_compositeLightA = opLinearComposite( _domainSpace=Xh , _imageSpace=Xh );
        M_compositeLightA->addList( M_Aq_free );
        M_compositeLightM = opLinearComposite( _domainSpace=Xh , _imageSpace=Xh );
        M_compositeLightM->addList( M_Mq_free );
        M_compositeLightF.resize( this->Nl() );
        for(int output=0; output<this->Nl(); output++)
        {
            M_compositeLightF[output]=functionalLinearComposite( _space=Xh );
            M_compositeLightF[output]->addList( M_Fq_free[output] );
        }
    }

}


template<int Order>
double HeatShield<Order>::output( int output_index, parameter_type const& mu, element_type &u, bool need_to_solve , bool export_outputs )
{
    using namespace vf;

    CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    double s=0;

    bool dont_use_operators_free = option(_name="do-not-use-operators-free").template as<bool>() ;
    auto fqm = backend->newVector( Xh );
    if ( output_index<2 )
    {
        for ( int q=0; q<Ql( output_index ); q++ )
        {
            if( dont_use_operators_free )
            {
                s += M_betaFq[output_index][q]*dot( *M_Fq[output_index][q] , u );
            }
            else
            {
                M_Fq_free[output_index][q]->containerPtr( fqm );
                s += M_betaFq[output_index][q]*dot( *fqm , u );
            }
        }
    }
    else
    {
        throw std::logic_error( "[HeatShield::output] error with output_index : only 0 or 1 " );
    }

    return s ;
}

}

#endif /* __HeatShield_H */


