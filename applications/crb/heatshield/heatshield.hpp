/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-13

  Copyright (C) 2009-2014 Feel++ Consortium

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
#ifndef FEELPP_HEATSHIELD_HPP
#define FEELPP_HEATSHIELD_HPP 1

#include <boost/timer.hpp>
#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>

#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>

#include <feel/feelts/bdf.hpp>

//#include <feel/feelcore/pslogger.hpp>



namespace Feel
{

po::options_description
makeHeatShieldOptions()
{
    po::options_description heatshieldoptions( "HeatShield options" );
    heatshieldoptions.add_options()
    // mesh parameters
    ( "hsize", Feel::po::value<double>()->default_value( 1 ), "first h value to start convergence" )
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
    return heatshieldoptions.add( bdf_options( "heatshield" ) );
}
AboutData
makeHeatShieldAbout( std::string const& str = "heatShield" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "heat shield Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2009-2014 Feel++ Consortium");
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "", "" );
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

    typedef typename space_type::element_type element_type;
};


/**
 * \class HeatShield
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
template<int Order>
class HeatShield : public ModelCrbBase< ParameterDefinition, FunctionSpaceDefinition<Order>, TimeDependent >
{
public:

    typedef ModelCrbBase<ParameterDefinition, FunctionSpaceDefinition<Order> , TimeDependent > super_type;

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

    typedef typename super_type::element_type element_type;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::mesh_ptrtype mesh_ptrtype;

    typedef typename FunctionSpaceDefinition<Order>::space_type space_type;

    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    using super_type::computeBetaQm;

    //! initialization of the model
    void initModel();
    //@}

    gmsh_ptrtype createGeo( double hsize );

    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<beta_vector_light_type, beta_vector_light_type, std::vector<beta_vector_light_type>  >
    computeBetaQ( parameter_type const& mu , double time , bool only_terms_time_dependent=false )
    {
        if( M_use_ginac )
        {
            std::string symbol;

            ginac_expressionA[0].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } } );
            this->M_betaAq[0] = ginac_expressionA[0].evaluate();
            ginac_expressionA[1].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } } );
            this->M_betaAq[1] = ginac_expressionA[1].evaluate();
            ginac_expressionA[2].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } } );
            this->M_betaAq[2] = ginac_expressionA[2].evaluate();

            this->M_betaMq[0] = ginac_expressionM[0].evaluate();

            ginac_expressionF[0].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } , { "surface", surface } } );
            this->M_betaFq[0][0] = ginac_expressionF[0].evaluate();
            ginac_expressionF[1].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } , { "surface", surface } } );
            this->M_betaFq[1][0] = ginac_expressionF[1].evaluate();
#if 0
            int idx=0;
            int nl = M_Nl;
            for(int i=0; i<nl; i++)
            {
                int ql=Ql[i];
                for(int j=0; j<ql; j++)
                {
                    ginac_expressionF[idx].expression().setParameterValues( { { "BiotOut", mu(0) } , { "BiotIn", mu(1) } , { "surface", surface } } );
                    auto projection=project(_space=continuous_p0, _expr=ginac_expressionF[idx]);
                    this->M_betaFq[i][j] = projection(0);
                    idx++;
                }
            }
#endif
#if 0
            LOG( INFO ) << "mu = "<<mu(0)<<" -- "<<mu(1);
            LOG( INFO ) <<"A0 : "<<this->M_betaAq[0]<<" -- should be 1";
            LOG( INFO ) <<"A1 : "<<this->M_betaAq[1]<<" -- should be "<<mu(0);
            LOG( INFO ) <<"A2 : "<<this->M_betaAq[2]<<" -- should be "<<mu(1);
            LOG( INFO ) <<"F0 : "<<this->M_betaFq[0][0]<<" -- should be "<<mu(0);
            LOG( INFO ) <<"F1 : "<<this->M_betaFq[1][0]<<" -- should be "<<1./surface;
            LOG( INFO ) <<"M0 : "<<this->M_betaMq[0]<<" -- should be 1";
#endif
        }//use ginac
        else
        {
            double biot_out   = mu( 0 );
            double biot_in    = mu( 1 );
            if( ! only_terms_time_dependent )
            {
                this->M_betaAq[0] = 1 ;
                this->M_betaAq[1] = biot_out ;
                this->M_betaAq[2] = biot_in  ;
                this->M_betaMq[0] = 1;
            }
            this->M_betaFq[0][0] = biot_out;
            this->M_betaFq[1][0] = 1./surface;
        }

        return boost::make_tuple( this->M_betaMq, this->M_betaAq, this->M_betaFq );
    }


    /**
     * \brief Returns the affine decomposition
     */
    affine_decomposition_light_type computeAffineDecompositionLight(){ return boost::make_tuple( this->M_Mq, this->M_Aq, this->M_Fq ); }

    void assemble();

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    double output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false, bool export_outputs=false );

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

    void initDataStructureForBetaCoeff();
    void buildGinacExpressions();

    bdf_ptrtype bdfModel(){ return M_bdf; }

private:
    mesh_ptrtype mesh;

    int M_Qa;
    int M_Qm;
    std::vector< int > M_Ql;
    int M_Nl;

    double surface;
    std::vector<operator_ptrtype> M_Aq_free;
    std::vector<operator_ptrtype> M_Mq_free;
    std::vector<std::vector<functional_ptrtype> > M_Fq_free;

    operatorcomposite_ptrtype M_compositeLightA;
    operatorcomposite_ptrtype M_compositeLightM;
    std::vector< functionalcomposite_ptrtype > M_compositeLightF;

    element_type u,v;

    bool M_use_ginac ;

    std::vector< Expr<GinacEx<2> > > ginac_expressionA;
    std::vector< Expr<GinacEx<2> > > ginac_expressionM;
    std::vector< Expr<GinacEx<2> > > ginac_expressionF;

    bdf_ptrtype M_bdf;

};

template<int Order>
void
HeatShield<Order>::initDataStructureForBetaCoeff()
{
    this->M_betaAq.resize( M_Qa );
    this->M_betaMq.resize( M_Qm );
    this->M_betaFq.resize( M_Nl );
    for(int i=0; i<M_Nl; i++)
    {
        int ql=M_Ql[i];
        this->M_betaFq[i].resize( ql );
    }
}

template<int Order>
void
HeatShield<Order>::buildGinacExpressions()
{

    std::vector< std::string > symbols_vec;
    symbols_vec.push_back("BiotOut");
    symbols_vec.push_back("BiotIn");
    symbols_vec.push_back("surface");
    for(int i=0; i<M_Qa; i++)
    {
        std::string name = ( boost::format("beta.A%1%") %i ).str();
        std::string filename = ( boost::format("GinacA%1%") %i ).str();
        ginac_expressionA.push_back( expr( soption(_name=name), Symbols( symbols_vec ) , filename ) );
    }


    for(int i=0; i<M_Qm; i++)
    {
        std::string name = ( boost::format("beta.M%1%") %i ).str();
        std::string filename = ( boost::format("GinacM%1%") %i ).str();
        ginac_expressionM.push_back( expr( soption(_name=name),  Symbols( symbols_vec ) , filename ) );
    }

    for(int i=0; i<M_Nl; i++)
    {
        int ql=M_Ql[i];
        for(int j=0; j<ql; j++)
        {
            std::string name = ( boost::format("beta.F%1%.%2%") %i %j ).str();
            std::string filename = ( boost::format("GinacF%1%.%2%") %i %j ).str();
            ginac_expressionF.push_back( expr( soption(_name=name), Symbols( symbols_vec ) , filename ) );
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

    CHECK( this->is_linear && this->is_time_dependent ) << "Invalid model is_linear:" << this->is_linear << " is_time_dependent:" << this->is_time_dependent << "\n";
    LOG_IF( WARNING, ((this->Options&NonLinear) == NonLinear) ) << "Invalid model is_linear:" << this->is_linear << " is_time_dependent:" << this->is_time_dependent << "\n";

    std::string mshfile_name = option("mshfile").as<std::string>();
    double hsize = option("hsize").as<double>();

    M_Qa=3;
    M_Qm=1;
    M_Nl=2;
    M_Ql.resize( 2 );
    M_Ql[0]=1;
    M_Ql[1]=1;

    this->M_betaAq.resize( M_Qa );
    this->M_betaMq.resize( M_Qm );
    this->M_betaFq.resize( M_Nl );
    this->M_betaFq[0].resize( M_Ql[0] );
    this->M_betaFq[1].resize( M_Ql[1] );

    M_use_ginac = boption(_name="crb.use-ginac-for-beta-expressions");

    if( mshfile_name=="" )
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                               _desc = createGeo( hsize ) );
    }
    else
    {

        bool load_mesh_already_partitioned=boption(_name="load-mesh-already-partitioned");
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
    auto Xh = space_type::New( mesh );

    //this->setFunctionSpaces( Pch<Order>( loadMesh( _mesh=new Mesh<Simplex<2>> ) ) );
    this->setFunctionSpaces( Xh );

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << this->Xh->nLocalDof() << "\n";
        std::cout << "Number of dof " << this->Xh->nDof() << "\n";
    }

    surface = integrate( _range=elements( mesh ), _expr=cst( 1. ) ).evaluate()( 0,0 );

    M_bdf = bdf( _space=Xh, _vm=Environment::vm(), _name="heatshield" , _prefix="heatshield" );

    bool dont_use_operators_free = boption(_name="do-not-use-operators-free") ;
    if( dont_use_operators_free )
    {
        this->M_Aq.resize( M_Qa );
        this->M_Mq.resize( M_Qm );
        this->M_Fq.resize( M_Nl );
        for(int l=0; l<M_Nl; l++)
        {
            this->M_Fq[l].resize( M_Ql[l] );
        }
    }
    else
    {
        M_Aq_free.resize( M_Qa );
        M_Mq_free.resize( M_Qm );
        M_Fq_free.resize( M_Nl );
        for(int l=0; l<M_Nl; l++)
        {
            M_Fq_free[l].resize( M_Ql[l] );
        }
    }

    initDataStructureForBetaCoeff();
    if( M_use_ginac )
        buildGinacExpressions();

    //typename Feel::ParameterSpace<ParameterSpaceDimension>::Element mu_min( M_Dmu );
    auto mu_min = this->Dmu->element();
    mu_min <<  /* Bi_out */ 1e-2 , /*Bi_in*/1e-3;
    this->Dmu->setMin( mu_min );
    //typename Feel::ParameterSpace<ParameterSpaceDimension>::Element mu_max( M_Dmu );
    auto mu_max = this->Dmu->element();
    mu_max << /* Bi_out*/0.5   ,  /*Bi_in*/0.1;
    this->Dmu->setMax( mu_max );

    LOG(INFO) << "Number of dof " << this->Xh->nLocalDof() << "\n";

    assemble();


} // HeatShield::init


template<int Order>
void HeatShield<Order>::assemble()
{
    auto mesh = this->Xh->mesh();

    u = this->Xh->element();
    v = this->Xh->element();

    if( boption("do-not-use-operators-free") )
    {
        this->M_Aq[0] = backend()->newMatrix( this->Xh, this->Xh );
        this->M_Aq[1] = backend()->newMatrix( this->Xh, this->Xh );
        this->M_Aq[2] = backend()->newMatrix( this->Xh, this->Xh );
        this->M_Mq[0] = backend()->newMatrix( this->Xh, this->Xh );
        this->M_Fq[0][0] = backend()->newVector( this->Xh );
        this->M_Fq[1][0] = backend()->newVector( this->Xh );
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[0]) = integrate( _range= elements( mesh ), _expr= gradt( u )*trans( grad( v ) ) );
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[1]) = integrate( _range= markedfaces( mesh, "left" ), _expr= idt( u )*id( v ) );
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[2]) = integrate( _range= markedfaces( mesh, "gamma_holes" ), _expr= idt( u )*id( v ) );
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Mq[0]) = integrate ( _range=elements( mesh ), _expr=idt( u )*id( v ) );
        form1(_test=this->Xh, _vector=this->M_Fq[0][0]) = integrate( _range=markedfaces( mesh,"left" ), _expr= id( v ) ) ;
        form1(_test=this->Xh, _vector=this->M_Fq[1][0]) = integrate( _range=elements( mesh ), _expr= id( v ) ) ;
    }
    else
    {
        auto expr_a0 = integrate( _range= elements( mesh ), _expr= gradt( u )*trans( grad( v ) ) );
        auto operatorfree0=opLinearFree( _domainSpace=this->Xh , _imageSpace=this->Xh , _expr=expr_a0 );
        operatorfree0->setName("A0");
        M_Aq_free[0]=operatorfree0;

        auto expr_a1 = integrate( _range= markedfaces( mesh, "left" ), _expr= idt( u )*id( v ) );
        auto operatorfree1=opLinearFree( _domainSpace=this->Xh , _imageSpace=this->Xh , _expr=expr_a1 );
        operatorfree1->setName("A1");
        M_Aq_free[1]=operatorfree1;

        auto expr_a2 = integrate( _range= markedfaces( mesh, "gamma_holes" ), _expr= idt( u )*id( v ) );
        auto operatorfree2=opLinearFree( _domainSpace=this->Xh , _imageSpace=this->Xh , _expr=expr_a2 );
        operatorfree2->setName("A2");
        M_Aq_free[2]=operatorfree2;

        auto expr_f00 = integrate( _range=markedfaces( mesh,"left" ), _expr= id( v ) ) ;
        auto functionalfree00 = functionalLinearFree( _space=this->Xh , _expr=expr_f00  );
        auto expr_f10 = integrate( _range=elements( mesh ), _expr= id( v ) ) ;
        auto functionalfree10 = functionalLinearFree( _space=this->Xh , _expr=expr_f10 );
        M_Fq_free[0][0]=functionalfree00;
        M_Fq_free[1][0]=functionalfree10;

        //mass matrix
        auto expr_m0 = integrate ( _range=elements( mesh ), _expr=idt( u )*id( v ) );
        auto operatorfreeM0=opLinearFree( _domainSpace=this->Xh , _imageSpace=this->Xh , _expr=expr_m0 );
        operatorfreeM0->setName("mass");
        M_Mq_free[0]=operatorfreeM0;
    }

    //for scalarProduct
    auto M = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
    form2( _test=this->Xh, _trial=this->Xh, _matrix=M ) =
        integrate( elements( mesh ), gradt( u )*trans( grad( v ) ) ) +
        integrate( markedfaces( mesh, "left" ), 0.01 * idt( u )*id( v ) ) +
        integrate( markedfaces( mesh, "gamma_holes" ), 0.001 * idt( u )*id( v ) )
        ;
    this->addEnergyMatrix( M );

    //scalar product used for mass matrix
    auto InnerMassMatrix = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
    form2( _test=this->Xh, _trial=this->Xh, _matrix=InnerMassMatrix ) =
        integrate( _range=elements( mesh ), _expr=idt( u ) * id( v ) ) ;
    this->addMassMatrix(InnerMassMatrix);

#if 0
    //scalar product used for the POD
    Mpod = backend->newMatrix( _test=this->Xh, _trial=this->Xh );
    form2( _test=this->Xh, _trial=this->Xh, _matrix=Mpod ) =
        integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) ) +
        integrate( _range= markedfaces( mesh, "left" ), _expr= 0.01 * idt( u )*id( v ) ) +
        integrate( _range= markedfaces( mesh, "gamma_holes" ), _expr= 0.001 * idt( u )*id( v ) )
        ;
#endif

    if( ! boption(_name="do-not-use-operators-free") )
    {
        M_compositeLightA = opLinearComposite( _domainSpace=this->Xh , _imageSpace=this->Xh );
        M_compositeLightA->addList( M_Aq_free );
        M_compositeLightM = opLinearComposite( _domainSpace=this->Xh , _imageSpace=this->Xh );
        M_compositeLightM->addList( M_Mq_free );
        M_compositeLightF.resize( M_Nl );
        for(int output=0; output<M_Nl; output++)
        {
            M_compositeLightF[output]=functionalLinearComposite( _space=this->Xh );
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

    bool dont_use_operators_free = boption(_name="do-not-use-operators-free") ;
    auto fqm = backend()->newVector( this->Xh );
    if ( output_index<2 )
    {
        for ( int q=0; q<M_Ql[ output_index ]; q++ )
        {
            if( dont_use_operators_free )
            {
                s += this->M_betaFq[output_index][q]*dot( *this->M_Fq[output_index][q] , u );
            }
            else
            {
                M_Fq_free[output_index][q]->containerPtr( fqm );
                s += this->M_betaFq[output_index][q]*dot( *fqm , u );
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
