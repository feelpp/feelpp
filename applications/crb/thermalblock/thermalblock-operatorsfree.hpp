/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Abdoulaye Samake <abdoulaye.samake1@ujf-grenoble.fr>
   Date: 2011-06-02

   Copyright (C) 2009-2014 Feel++ Consortium

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
#ifndef FEELPP_ThermalBlockFree_HPP
#define FEELPP_ThermalBlockFree_HPP 1

#include <boost/timer.hpp>
#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>

#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>

namespace Feel
{

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */

po::options_description
makeThermalBlockFreeOptions()
{
    po::options_description thermalblockoptionsfree( "ThermalBlock options" );
    thermalblockoptionsfree.add_options()
    ( "hsize", po::value<double>()->default_value( 0.05 ), "mesh size" )
    ( "geofile", Feel::po::value<std::string>()->default_value( "" ), "name of the geofile input (used to store DB)")
    ( "mshfile", Feel::po::value<std::string>()->default_value( "" ), "name of the gmsh file input")
    ( "gamma_dir", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary Dirichlet formulation" )
    ( "load-mesh-already-partitioned", Feel::po::value<bool>()->default_value( "true" ), "load a mesh from mshfile that is already partitioned if true, else the mesh loaded need to be partitioned")
    ;
    return thermalblockoptionsfree;
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
AboutData
makeThermalBlockFreeAbout( std::string const& str = "thermalBlockfree" )
{
    AboutData about( "thermalblockfree" ,
                     "thermalblockfree" ,
                     "0.1",
                     "2D Heterogeneous Thermal Block Problem",
                     Feel::AboutData::License_GPL,
                     "Copyright (C) 2009-2014 Feel++ Consortium");
    about.addAuthor( "Abdoulaye Samake", "main developer", "abdoulaye.samake@ujf-grenoble.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "contributor", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "contributor", "", "" );
    return about;

}



class ParameterDefinition
{
public :
    static const uint16_type nx = 3;
    static const uint16_type ny = 3;
    static const uint16_type ParameterSpaceDimension = nx*ny-1;
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
};

class FunctionSpaceDefinition
{
public :
    typedef double value_type;

    static const uint16_type Order = 1;

    //! geometry entities type composing the mesh, here Simplex in Dimension 2 of Order 1
    typedef Simplex<2> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;

    //! the basis type of our approximation space
    typedef bases<Lagrange<Order,Scalar> > basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    typedef typename space_type::element_type element_type;
};


/**
 * \class ThermalBlock using operators free
 *
 * ThermalBlock Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
class ThermalBlockFree : public ModelCrbBase< ParameterDefinition , FunctionSpaceDefinition >
{

    static const uint16_type nx = 3;
    static const uint16_type ny = 3;

public:

    typedef ModelCrbBase<ParameterDefinition, FunctionSpaceDefinition> super_type;

    static const uint16_type ParameterSpaceDimension = 8;

    //! Polynomial order \f$P_1\f$
    static const uint16_type Order = 1;

    //! numerical type is double
    typedef double value_type;

    typedef typename super_type::element_type element_type;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::mesh_ptrtype mesh_ptrtype;

    //operator free
    typedef typename super_type::operator_type operator_type;
    typedef typename super_type::operator_ptrtype operator_ptrtype;
    typedef typename super_type::operatorcomposite_type operatorcomposite_type;
    typedef typename super_type::operatorcomposite_ptrtype operatorcomposite_ptrtype;
    typedef typename super_type::functionalcomposite_type functionalcomposite_type;
    typedef typename super_type::functionalcomposite_ptrtype functionalcomposite_ptrtype;
    typedef typename super_type::functional_type functional_type;
    typedef typename super_type::functional_ptrtype functional_ptrtype;

    //! initialisation of the model and definition of parameters values
    void initModel();

    std::string subdomainFromBoundary( std::string const& ) const;

    int subdomainId( std::string const& ) const;
    // \return the number of terms in affine decomposition of left hand
    // side bilinear form
    int Qa() const
    {
        return 10;
    }

    /**
     * \param l the index of output
     *
     * \return number of terms  in affine decomposition of the \p q th output term
     * in our case we have 1 term : 1 * \int_south v
     * but if u!=0 on the north we add (nx+1) terms (those on the north from dirichlet condition)
     */
    int Ql( int l ) const
    {
        return 1;
    }

    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     * in our case, the output is the right hand side of the primal problem
     *
     * \return number of outputs associated to the model
     */
    int Nl() const
    {
        return 1;
    }


    boost::tuple<beta_vector_light_type, std::vector<beta_vector_light_type> >
    computeBetaQ( parameter_type const& mu )
    {
        M_betaAq[0] = 1;

        for ( int i=1; i<nx*ny; i++ )
            M_betaAq[i] = mu( i-1 );
        int index=nx*ny;
        //penalization term
        M_betaAq[index]=1;

        //compliant output
        M_betaFq[0][0] = 1;

        return boost::make_tuple( M_betaAq, M_betaFq );
    }


    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);

    parameter_type refParameter()
    {
        return Dmu->min();
    }

    virtual operatorcomposite_ptrtype operatorCompositeLightA()
    {
        return M_compositeLightA;
    }
    virtual std::vector< functionalcomposite_ptrtype > functionalCompositeLightF()
    {
        return M_compositeLightF;
    }

private:
    value_type gamma_dir;

    mesh_ptrtype mmesh;


    std::map<std::string, std::pair<boost::timer, double> > timers;
    std::map<std::string,double> stats;

    std::vector<std::string> domainMarkers ;
    std::vector<std::string> northMarkers ;
    std::vector<std::string> westMarkers ;
    std::vector<std::string> eastMarkers ;
    std::vector<std::string> southMarkers ;

    //vectors contain index of subdomains along south and north boudaries
    //useful to implement computeBetaQm
    //they are filled in the initModel() function
    std::vector<int> south_subdomain_index;
    std::vector<int> north_subdomain_index;

    std::vector<operator_ptrtype> M_Aq_free;
    std::vector<std::vector<functional_ptrtype> > M_Fq_free;

    operatorcomposite_ptrtype M_compositeLightA;
    std::vector< functionalcomposite_ptrtype > M_compositeLightF;

    element_type u,v;

}; // ThermalBlockFree



gmsh_ptrtype
thermalBlockGeometryFree( int nx, int ny, double hsize )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;

    gmsh_ptrtype gmshp( new Gmsh );
    ostr << gmshp->preamble() <<" \n ";
    ostr << "nx=" << nx << ";\n"
         << "ny=" << ny << ";\n"
         << "hsize=" << hsize << ";\n"
         << "dx = 1/nx;\n"
         << "dy = 1/ny;\n"
         << "t1=0;\n"
         << "x=0;\n"
         << "y=0;\n"
         << "For j In {0:ny}\n"
         << "  For i In {0:nx}\n"
         << "  y = j*dy ;\n"
         << "  x = i*dx ;\n"
         << "  t1+=1;\n"
         << "  Point(t1) = {x, y, 0,  hsize} ;\n"
         << "  EndFor\n"
         << "EndFor\n"
         << "t2=0;\n"
         << "For j In {0:ny}\n"
         << "  For i In {0:nx-1}\n"
         << "  t2 +=1;\n"
         << "  Line(t2)={t2+j,t2+j+1};\n"
         << "//  Physical Line(t2+1)={t2};\n"
         << "  EndFor\n"
         << "EndFor\n";

    for ( int i = 1; i <= nx; ++i )
        ostr << "Physical Line(\"south_domain-"  << i << "\")={"<< i << "};\n";

    for ( int i = nx*ny+1, j=nx*( ny-1 )+1; i<= nx*( ny+1 ); ++i,++j )
        ostr << "  Physical Line(\"north_domain-" << j << "\") = {" << i << "};\n";

    ostr << "t3 = (ny+1)*nx;\n"
         << "t4 = 0;\n"
         << "For i In {0:nx}\n"
         << "  For j In {0:ny-1}\n"
         << "  t3 +=1;\n"
         << "  t4 +=1;\n"
         << " Line(t3)={t4,t4+nx+1};\n"
         << "// Physical Line(t3+1)={t3};\n"
         << " EndFor\n"
         << "EndFor\n";

    for ( int i = 1, j=0, k=1; i <= ny; ++i, j += nx+1, k += nx )
        ostr << "Physical Line(\"west_domain-"  << k << "\")={"<< nx*( ny+1 )+j+1 << "};\n";

    for ( int i = 1, j=0, k=nx; i <= ny; ++i, j += nx+1, k += nx )
        ostr << "Physical Line(\"east_domain-"  << k << "\")={"<< nx*( ny+2 )+j+1 << "};\n";

    ostr << "t5 = 0;\n"
         << "ne = (ny+1)*nx+1;\n"
         << "For j In {0:ny-1}\n"
         << "  For i In {0:nx-1}\n"
         << "  t5 +=1;\n"
         << "  Line Loop(t5)={t5,(ne+t5+j),-(nx+t5),-(ne+t5+j-1)};\n"
         << "  Plane Surface(t5)={t5};\n"
         << " EndFor\n"
         << "EndFor\n";

    for ( int i = 1, d= 1; i <= nx; ++i )
        for ( int j = 1; j <= ny; ++j, ++d )
        {
            ostr << "  Physical Surface(\"domain-"<< d << "\")={" << d << "};\n";
        }

    nameStr << "thermalblockfree";

    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}


void
ThermalBlockFree::initModel()
{

    gamma_dir=doption(_name="gamma_dir");

    std::string mshfile_name = option(_name="mshfile").as<std::string>();

    double hsize = doption(_name="hsize");

    if( mshfile_name=="" )
    {
        mmesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=thermalBlockGeometryFree( nx, ny, hsize ) );
    }
    else
    {

        bool load_mesh_already_partitioned=option(_name="load-mesh-already-partitioned").as<bool>();
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
            mmesh = loadGMSHMesh( _mesh=new mesh_type,
                                 _filename=filename,
                                 _rebuild_partitions=false,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
        }
        else
        {
            mmesh = loadGMSHMesh( _mesh=new mesh_type,
                                 _filename=option("mshfile").as<std::string>(),
                                 _rebuild_partitions=false,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
        }
    }

    auto names = mmesh->markerNames();
    southMarkers.clear();
    northMarkers.clear();
    eastMarkers.clear();
    westMarkers.clear();
    domainMarkers.clear();
    north_subdomain_index.clear();
    south_subdomain_index.clear();

    for ( auto it=names.begin(), en=names.end(); it!=en; ++it )
    {
        if ( it->first.find( "south_" ) != std::string::npos )
            southMarkers.push_back( it->first );

        if ( it->first.find( "north_" ) != std::string::npos )
            northMarkers.push_back( it->first );

        if ( it->first.find( "west_" ) != std::string::npos )
            westMarkers.push_back( it->first );

        if ( it->first.find( "east_" ) != std::string::npos )
            eastMarkers.push_back( it->first );

        if ( ( it->first.find( "domain-" ) != std::string::npos ) &&
                ( it->first.find( "_" ) == std::string::npos ) )
            domainMarkers.push_back( it->first );
    }

    /*
     * The function space and some associate elements are then defined
     */
    auto Xh = functionspace_type::New( mmesh );
    this->setFunctionSpaces( Xh );

    //  initialization
    M_Aq_free.resize( Qa() );
    M_Fq_free.resize( Nl() );
    for(int l=0; l<Nl(); l++)
    {
        M_Fq_free[l].resize( Ql(l) );
    }
    DVLOG(2) <<"[ThermalBlock::init] done allocating matrices/vectors \n";


    int qa = Qa();
    M_betaAq.resize( qa );

    int nl = Nl();
    M_betaFq.resize( nl );
    for(int i=0; i<nl; i++)
    {
        int ql=Ql(i);
        M_betaFq[i].resize( ql );
    }

    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();

    for ( int i=0; i<nx*ny-1; i++ )
    {
        mu_min[i]=0.1;
        mu_max[i]=10;
    }

    Dmu->setMin( mu_min );
    Dmu->setMax( mu_max );

    DVLOG(2) <<"[ThermalBlock::init] done with parameter space\n";
    u = Xh->element();
    v = Xh->element();

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << Xh->nLocalDof() << std::endl;
        std::cout << "Number of dof " << Xh->nDof() << std::endl;
    }

    int index=0;//index for M_Fq or M_Aq
    int subdomain_index; //index for subdomains along a boundary

    auto M = backend()->newMatrix( Xh , Xh );

    double mu_min_coeff=0.1;
    // on boundary north we have u=0 so term from weak dirichlet condition
    // vanish in the right hand side
    auto expr_f00 =
          integrate( _range=markedfaces( mmesh,"south_domain-1" ), _expr= id( v ) )
        + integrate( _range=markedfaces( mmesh,"south_domain-2" ), _expr= id( v ) )
        + integrate( _range=markedfaces( mmesh,"south_domain-3" ), _expr= id( v ) );
    auto functionalfree00 = functionalLinearFree( _space=Xh , _expr=expr_f00  );
    M_Fq_free[0][0]=functionalfree00;

    BOOST_FOREACH( auto marker, southMarkers )
    {
        subdomain_index = subdomainId( marker );
        south_subdomain_index.push_back( subdomain_index );
    }

    DVLOG(2) <<"[ThermalBlock::init] done with rhs\n";

    auto expr_a0 = integrate( markedelements( mmesh, "domain-1" ), gradt( u )*trans( grad( v ) ) );
    auto operatorfree0=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a0 );
    operatorfree0->setName("A0");
    M_Aq_free[0]=operatorfree0;

    auto expr_a1 = integrate( markedelements( mmesh, "domain-2" ), gradt( u )*trans( grad( v ) ) );
    auto operatorfree1=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a1 );
    operatorfree1->setName("A1");
    M_Aq_free[1]=operatorfree1;

    auto expr_a2 = integrate( markedelements( mmesh, "domain-3" ), gradt( u )*trans( grad( v ) ) );
    auto operatorfree2=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a2 );
    operatorfree2->setName("A2");
    M_Aq_free[2]=operatorfree2;

    auto expr_a3 = integrate( markedelements( mmesh, "domain-4" ), gradt( u )*trans( grad( v ) ) );
    auto operatorfree3=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a3 );
    operatorfree3->setName("A3");
    M_Aq_free[3]=operatorfree3;

    auto expr_a4 = integrate( markedelements( mmesh, "domain-5" ), gradt( u )*trans( grad( v ) ) );
    auto operatorfree4=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a4 );
    operatorfree4->setName("A4");
    M_Aq_free[4]=operatorfree4;

    auto expr_a5 = integrate( markedelements( mmesh, "domain-6" ), gradt( u )*trans( grad( v ) ) );
    auto operatorfree5=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a5 );
    operatorfree5->setName("A5");
    M_Aq_free[5]=operatorfree5;

    auto expr_a6 =
         integrate( markedelements( mmesh, "domain-7" ), gradt( u )*trans( grad( v ) ) )
        +integrate( markedfaces( mmesh, "north_domain-7" ),
                   -gradt( u )*vf::N()*id( v )
                   -grad( u )*vf::N()*idt( v )
                   );
    auto operatorfree6=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a6 );
    operatorfree6->setName("A6");
    M_Aq_free[6]=operatorfree6;

    auto expr_a7 =
         integrate( markedelements( mmesh, "domain-8" ), gradt( u )*trans( grad( v ) ) )
        +integrate( markedfaces( mmesh, "north_domain-8" ),
                    -gradt( u )*vf::N()*id( v )
                    -grad( u )*vf::N()*idt( v )
                    );
    auto operatorfree7=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a7 );
    operatorfree7->setName("A7");
    M_Aq_free[7]=operatorfree7;

    auto expr_a8 =
         integrate( markedelements( mmesh, "domain-9" ), gradt( u )*trans( grad( v ) ) )
        +integrate( markedfaces( mmesh, "north_domain-9" ),
                    -gradt( u )*vf::N()*id( v )
                    -grad( u )*vf::N()*idt( v )
                    );
    auto operatorfree8=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a8 );
    operatorfree8->setName("A8");
    M_Aq_free[8]=operatorfree8;

    auto expr_a9 =
         integrate( markedfaces( mmesh, "north_domain-7" ),gamma_dir*idt( u )*id( v )/h() )
        +integrate( markedfaces( mmesh, "north_domain-8" ),gamma_dir*idt( u )*id( v )/h() )
        +integrate( markedfaces( mmesh, "north_domain-9" ),gamma_dir*idt( u )*id( v )/h() );
    auto operatorfree9=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a9 );
    operatorfree9->setName("A9");
    M_Aq_free[9]=operatorfree9;


    form2( _test=Xh, _trial=Xh, _matrix=M ) = integrate( markedelements( mmesh, "domain-1" ), gradt( u )*trans( grad( v ) )  );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mmesh, "domain-2" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mmesh, "domain-3" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mmesh, "domain-4" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mmesh, "domain-5" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mmesh, "domain-6" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mmesh, "domain-7" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mmesh, "domain-8" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mmesh, "domain-9" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );

    form2( _test=Xh, _trial=Xh, _matrix=M ) +=  integrate( markedfaces( mmesh, "north_domain-7" ),
                                      -gradt( u )*vf::N()*id( v ) * mu_min_coeff
                                      -grad( u )*vf::N()*idt( v ) * mu_min_coeff
                                      );
    form2( _test=Xh, _trial=Xh, _matrix=M ) +=  integrate( markedfaces( mmesh, "north_domain-8" ),
                                      -gradt( u )*vf::N()*id( v ) * mu_min_coeff
                                      -grad( u )*vf::N()*idt( v ) * mu_min_coeff
                                      );
    form2( _test=Xh, _trial=Xh, _matrix=M ) +=  integrate( markedfaces( mmesh, "north_domain-9" ),
                                      -gradt( u )*vf::N()*id( v ) * mu_min_coeff
                                      -grad( u )*vf::N()*idt( v ) * mu_min_coeff
                                      );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedfaces( mmesh, "north_domain-7" ),gamma_dir*idt( u )*id( v )/h() );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedfaces( mmesh, "north_domain-8" ),gamma_dir*idt( u )*id( v )/h() );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedfaces( mmesh, "north_domain-9" ),gamma_dir*idt( u )*id( v )/h() );

    BOOST_FOREACH( auto marker, northMarkers )
    {
        std::string sid = subdomainFromBoundary( marker );
        subdomain_index = subdomainId( marker );
        north_subdomain_index.push_back( subdomain_index );
    }

    this->addEnergyMatrix( M );

    M_compositeLightA = opLinearComposite( _domainSpace=Xh , _imageSpace=Xh );
    M_compositeLightA->addList( M_Aq_free );
    M_compositeLightF.resize( 1 );
    int output=0;
    M_compositeLightF[output]=functionalLinearComposite( _space=Xh );
    M_compositeLightF[output]->addList( M_Fq_free[output] );

}//initModel()

double
ThermalBlockFree::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{

    double output=0;

    auto fq = backend()->newVector( this->Xh );
    if ( output_index==0 )
    {
        for ( int q=0; q<Ql( output_index ); q++ )
        {
            M_Fq_free[output_index][q]->containerPtr( fq );
            output += this->M_betaFq[output_index][q]*dot( *fq , u );
        }
    }
    else
    {
        throw std::logic_error( "[ThermalBlock::output] error with output_index : only 0 " );
    }

    return output;
}//output



std::string
ThermalBlockFree::subdomainFromBoundary( std::string const& boundary ) const
{
    typedef std::vector< std::string > split_vector_type;

    split_vector_type SplitVec; // #2: Search for tokens
    boost::split( SplitVec, boundary, boost::is_any_of( "_" ), boost::token_compress_on );
    return SplitVec[1];
}


int
ThermalBlockFree::subdomainId( std::string const& domain ) const
{
    typedef std::vector< std::string > split_vector_type;
    split_vector_type SplitVec; // #2: Search for tokens
    boost::split( SplitVec, domain, boost::is_any_of( "-" ), boost::token_compress_on );
    return boost::lexical_cast<int>( SplitVec[1] );
}


} // Feel

#endif /* __ThermalBlock_H */
