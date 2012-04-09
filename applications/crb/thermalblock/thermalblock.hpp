/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Abdoulaye Samake <abdoulaye.samake1@ujf-grenoble.fr>
   Date: 2011-06-02

   Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)

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
#ifndef __ThermalBlock_H
#define __ThermalBlock_H 1

#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/aitken.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelcrb/parameterspace.hpp>
#include <vector>

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
makeThermalBlockOptions()
{
    po::options_description thermalblockoptions( "ThermalBlock options" );
    thermalblockoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.05 ), "mesh size" )
    ( "nx", po::value<int>()->default_value( 3 ), "number of blocks in the x direction" )
    ( "ny", po::value<int>()->default_value( 3 ), "number of blocks in the y direction" )
    ( "gamma_dir", Feel::po::value<double>()->default_value( 10 ),
      "penalisation parameter for the weak boundary Dirichlet formulation" )
    ;
    return thermalblockoptions.add( Feel::feel_options() );
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
AboutData
makeThermalBlockAbout( std::string const& str = "thermalBlock" )
{
    AboutData about( "thermalblock" ,
                     "thermalblock" ,
                     "0.1",
                     "2D Heterogeneous Thermal Block Problem",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Abdoulaye Samake", "main developer", "abdoulaye.samake@ujf-grenoble.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "contributor", "christophe.prudhomme@ujf-grenoble.fr", "" );
    about.addAuthor( "Stephane Veys", "contributor", "stephane.veys@imag.fr", "" );
    return about;

}



using namespace vf;


/**
 * \class ThermalBlock
 *
 * ThermalBlock Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */


class ThermalBlock

{

    static const uint16_type nx = 3;
    static const uint16_type ny = 3;


public:
    static const uint16_type ParameterSpaceDimension = nx*ny;
    static const bool is_time_dependent = false;

    //! Polynomial order \f$P_2\f$
    static const uint16_type Order = 3;

    //! numerical type is double
    typedef double value_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    //! sparse matrix type associated with backend
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    //! vector type associated with backend
    typedef typename backend_type::vector_type vector_type;


    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef boost::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;


    //! geometry entities type composing the mesh, here Simplex in Dimension 2 of Order 1
    typedef Simplex<2> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! function space that holds piecewise constant (\f$P_0\f$) functions (e.g. to store material properties or partitioning
    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Scalar, Discontinuous> > > p0_space_type;
    //! an element type of the \f$P_0\f$ discontinuous function space
    typedef typename p0_space_type::element_type p0_element_type;

    //! the basis type of our approximation space
    typedef bases<Lagrange<Order,Scalar> > basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type, value_type> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;


    //! an element type of the approximation function space
    typedef typename functionspace_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    typedef  boost::numeric::ublas::vector<element_type> Vector_type ;


    /* parameter space */
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;
    typedef parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef parameterspace_type::sampling_type sampling_type;
    typedef parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    typedef boost::shared_ptr<element_type> element_ptrtype;

    typedef Eigen::VectorXd theta_vector_type;
    typedef boost::tuple<std::vector<sparse_matrix_ptrtype>, std::vector<std::vector<vector_ptrtype>  > > affine_decomposition_type;


    /**
     * Constructor
     *
     */

    ThermalBlock( po::variables_map const& vm )
        :
        M_vm ( vm ),
        M_backend( backend_type::build( vm ) ),
        meshSize( vm["hsize"].as<double>() ),
        gamma_dir( M_vm["gamma_dir"].as<double>() ),
        M_Dmu( new parameterspace_type ),
        timers()
    {
        this->init();
    }

    //! initialisation of the model and definition of parameters values
    void init();

    int numberOfBlocks() const
    {
        return nx*ny;
    }
    int numberOfBlocksInX() const
    {
        return nx;
    }
    int numberOfBlocksInY() const
    {
        return ny;
    }

    // std::map<std::string, boost::tuple<int, int> >
    // markerNames() const { return M_markername; }


    double l2Error( element_type& u );

    double h1Error( element_type& u );

    void exportResults( element_type& u );

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

    void checkout();

    void checkout( const double* X, unsigned long P, double* Y, unsigned long N );

    std::string subdomainFromBoundary( std::string const& ) const;

    int subdomainId( std::string const& ) const;
    // \return the number of terms in affine decomposition of left hand
    // side bilinear form
    int Qa() const
    {
        return nx*ny+nx+1;
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
     * \brief Returns the function space
     */
    functionspace_ptrtype functionSpace()
    {
        return Xh;
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

    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
    {
        return M_Dmu;
    }

    boost::tuple<theta_vector_type, std::vector<theta_vector_type> >
    computeThetaq( parameter_type const& mu , double time=0 )
    {
        M_thetaAq.resize( Qa() );

        //mu_i inside the domain (all subdomains index)
        for ( int i=0; i<nx*ny; i++ )
        {
            M_thetaAq( i ) = mu( i );
            //std::cout<<"[computeThetaAq] M_thetaAq("<<i<<") = mu ("<<i<<")"<<std::endl;
        }

        //IMPORTANT REMARK, subdomain indices begin at 1 and not 0
        int index_theta=nx*ny;

        for ( int i=0; i<north_subdomain_index.size(); i++ )
        {
            M_thetaAq( index_theta ) = mu ( north_subdomain_index[i]-1 );
            //std::cout<<"M_thetaAq("<<index_theta<<") = mu ("<<north_subdomain_index[i]-1<< ")"<<std::endl;
            index_theta++;
        }

        M_thetaAq( index_theta ) = 1;

        M_thetaFq.resize( Nl() );
        M_thetaFq[0].resize( Ql( 0 ) );

        for ( int i=0; i<Ql( 0 ); i++ ) M_thetaFq[0]( i ) = 1;

        return boost::make_tuple( M_thetaAq, M_thetaFq );
    }


    /**
     * \brief return the coefficient vector \p q component
     *
     */
    value_type thetaAq( int q ) const
    {
        return M_thetaAq( q );
    }



    /**
     * \brief return the coefficient vector
     */
    std::vector<theta_vector_type> const& thetaFq() const
    {
        return M_thetaFq;
    }


    /**
     * \return the \p q -th term of the \p l -th output
     */
    value_type thetaL( int l, int q ) const
    {
        return M_thetaFq[l]( q );
    }

    /**
     * set the mesh characteristic length to \p s
     */
    void setMeshSize( double s )
    {
        meshSize = s;
    }


    /**
     * \brief Returns the affine decomposition
     */
    affine_decomposition_type computeAffineDecomposition();

    /**
     * \brief solve the model for parameter \p mu
     * \param mu the model parameter
     * \param u the temperature field
     */
    void solve( parameter_type const& mu, element_ptrtype& u );

    /**
     * solve for a given parameter \p mu
     */
    void solve( parameter_type const& mu );

    /**
     * solve \f$ M u = f \f$
     */
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f );

    /**
     * returns the scalar product of the boost::shared_ptr vector x and
     * boost::shared_ptr vector y
     */
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y );

    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_type const& x, vector_type const& y );


    /**
     * update the PDE system with respect to \param mu
     */
    void update( parameter_type const& mu );

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu );

    /**
     * create a new matrix
     * \return the newly created matrix
     */
    sparse_matrix_ptrtype newMatrix() const;

private:

    po::variables_map M_vm;

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

    value_type gamma_dir;

    parameterspace_ptrtype M_Dmu;

    functionspace_ptrtype Xh;
    element_ptrtype pT;

    mesh_ptrtype mmesh;

    sparse_matrix_ptrtype D,M;
    vector_ptrtype F;



    std::map<std::string, std::pair<boost::timer, double> > timers;
    std::map<std::string,double> stats;

    std::vector<std::string> domainMarkers ;
    std::vector<std::string> northMarkers ;
    std::vector<std::string> westMarkers ;
    std::vector<std::string> eastMarkers ;
    std::vector<std::string> southMarkers ;

    //vectors contain index of subdomains along south and north boudaries
    //useful to implement computeThetaq
    //they are filled in the init() function
    std::vector<int> south_subdomain_index;
    std::vector<int> north_subdomain_index;

    std::vector<sparse_matrix_ptrtype> M_Aq;
    std::vector<std::vector<vector_ptrtype> > M_Fq;

    theta_vector_type M_thetaAq;
    std::vector<theta_vector_type> M_thetaFq;

}; // ThermalBlock


gmsh_ptrtype
thermalBlockGeometry( int nx, int ny, double hsize )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;

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

    nameStr << "thermalblock";

    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}


void
ThermalBlock::init()
{
    mmesh = createGMSHMesh( _mesh=new mesh_type,
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                            _desc=thermalBlockGeometry( nx, ny, meshSize ) );


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
    Xh = functionspace_type::New( mmesh );
    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    //  initialization
    M_Aq.resize( Qa() );

    for ( int i=0; i<Qa(); i++ )
    {
        M_Aq[i] = M_backend->newMatrix( Xh, Xh );
    }

    M_Fq.resize( 1 );
    M_Fq[0].resize( Ql( 0 ) );

    for ( int i=0; i<Ql( 0 ); i++ )
    {
        M_Fq[0][i] = M_backend->newVector( Xh );
    }

    Log() <<"[ThermalBlock::init] done allocating matrices/vectors \n";
    D = M_backend->newMatrix( Xh, Xh );
    F = M_backend->newVector( Xh );

    using namespace Feel::vf;
    Feel::ParameterSpace<nx*ny>::Element mu_min( M_Dmu );
    Feel::ParameterSpace<nx*ny>::Element mu_max( M_Dmu );

    for ( int i=0; i<nx*ny; i++ )
    {
        mu_min[i]=0.1;
        mu_max[i]=10;
    }

    M_Dmu->setMin( mu_min );
    M_Dmu->setMax( mu_max );

    Log() <<"[ThermalBlock::init] done with parameter space\n";
    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    form2( Xh, Xh, D,_init=true )=integrate( elements( mmesh ),0*idt( v )*id( v ), _Q<0>() );
    D->close();

    Log() << "Number of dof " << Xh->nLocalDof() << "\n";

    int index=0;//index for M_Fq or M_Aq
    int subdomain_index; //index for subdomains along a boundary

    // right hand side
    form1( Xh, M_Fq[0][0], _init=true ) ;
    BOOST_FOREACH( auto marker, southMarkers )
    {
        form1( Xh, M_Fq[0][0] ) += integrate( markedfaces( mmesh,marker ), id( v ) );
        subdomain_index = subdomainId( marker );
        south_subdomain_index.push_back( subdomain_index );
    }
    M_Fq[0][0]->close();
    Log() <<"[ThermalBlock::init] done with rhs\n";
    // on boundary north we have u=0 so term from weak dirichlet condition
    // vanish in the right hand side
    BOOST_FOREACH( auto domain, domainMarkers )
    {
        Log() <<"[ThermalBlock::init] domain " << domain << "\n";
        form2( Xh, Xh, M_Aq[index],_init=true ) =
            integrate( markedelements( mmesh, domain ), gradt( u )*trans( grad( v ) ) );
        M_Aq[index]->close();
        Log() <<"[ThermalBlock::init] done with Aq[" << index << "]\n";
        index++;

    }
    Log() <<"[ThermalBlock::init] done with domainMarkers\n";
    int last_index_Aq = Qa()-1;
    form2( Xh, Xh, M_Aq[last_index_Aq],_init=true );
    BOOST_FOREACH( auto marker, northMarkers )
    {
        std::string sid = subdomainFromBoundary( marker );
        form2( Xh, Xh, M_Aq[index], _init=true ) =  integrate( markedfaces( mmesh, marker ),
                -gradt( u )*vf::N()*id( v )
                -grad( u )*vf::N()*idt( v )
                                                             );
        M_Aq[index]->close();
        index++;

        form2( Xh, Xh, M_Aq[last_index_Aq] ) += integrate( markedfaces( mmesh, marker ),gamma_dir*idt( u )*id( v )/h() );

        subdomain_index = subdomainId( marker );
        north_subdomain_index.push_back( subdomain_index );
    }
    M_Aq[last_index_Aq]->close();
    Log() <<"[ThermalBlock::init] done with boundaryMarkers\n";

    M = M_backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, M, _init=true ) =
        integrate( elements( mmesh ), id( u )*idt( v ) + grad( u )*trans( gradt( u ) ) );
    M->close();
}//init()


const uint16_type ThermalBlock::Order;


typename ThermalBlock::affine_decomposition_type
ThermalBlock::computeAffineDecomposition()
{
    return boost::make_tuple( M_Aq, M_Fq );
}


void
ThermalBlock::update( parameter_type const& mu )
{
    D->zero();

    for ( size_type q = 0; q < M_Aq.size(); ++q )
    {
        //std::cout << "[affine decomp] scale q=" << q << " with " << M_thetaAq[q] << "\n";
        D->addMatrix( M_thetaAq[q], M_Aq[q] );
    }

    F->close();
    F->zero();

    for ( size_type q = 0; q < M_Fq[0].size(); ++q )
    {
        //std::cout << "[affine decomp] scale q=" << q << " with " << M_thetaFq[0][q] << "\n";
        F->add( M_thetaFq[0][q], M_Fq[0][q] );
    }
}


void
ThermalBlock::solve( parameter_type const& mu )
{
    //sd::cout<<"solve (mu) for mu = ["
    //int size = mu.size();
    //for(int i=0;i<size-1;i++) std::cout<<mu(i)<<" , ";
    //std::cout<< mu(size-1)<<"] "<<std::endl;
    element_ptrtype u( new element_type( Xh ) );
    this->solve( mu, u );
    this->exportResults( *u );
}


void
ThermalBlock::solve( parameter_type const& mu, element_ptrtype& u )
{
    this->computeThetaq( mu );
    this->update( mu );
    M_backend->solve( _matrix=D,  _solution=u, _rhs=F );
}


void
ThermalBlock::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    M_backend->solve( _matrix=M,  _solution=u, _rhs=f );
    //std::cout << "l2solve(u,f) done\n";
}




double
ThermalBlock::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}


double
ThermalBlock::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}


typename ThermalBlock::sparse_matrix_ptrtype
ThermalBlock::newMatrix() const
{
    return M_backend->newMatrix( Xh, Xh );
}



double
ThermalBlock::output( int output_index, parameter_type const& mu )
{
    using namespace vf;
    this->solve( mu, pT );
    vector_ptrtype U( M_backend->newVector( Xh ) );
    *U = *pT;

    // right hand side (compliant)
    if ( output_index == 0 )
    {
        return M_thetaFq[0]( 0 )*dot( M_Fq[0][0], U );
    }

    return 0;
}//output



std::string
ThermalBlock::subdomainFromBoundary( std::string const& boundary ) const
{
    typedef std::vector< std::string > split_vector_type;

    split_vector_type SplitVec; // #2: Search for tokens
    boost::split( SplitVec, boundary, boost::is_any_of( "_" ), boost::token_compress_on );
    return SplitVec[1];
}


int
ThermalBlock::subdomainId( std::string const& domain ) const
{
    typedef std::vector< std::string > split_vector_type;
    split_vector_type SplitVec; // #2: Search for tokens
    boost::split( SplitVec, domain, boost::is_any_of( "-" ), boost::token_compress_on );
    return boost::lexical_cast<int>( SplitVec[1] );
}



double
ThermalBlock::l2Error( element_type& u )
{
    auto Xh=u.functionSpace();
    auto mesh=Xh->mesh();
    value_type pi = M_PI;
    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );
    double L2error2 =integrate( elements( mesh ), ( idv( u )-g )*( idv( u )-g ) ).evaluate()( 0,0 );
    double error = math::sqrt( L2error2 );

    return error;
}


double
ThermalBlock::h1Error( element_type& u )
{
    auto Xh=u.functionSpace();
    auto mesh=Xh->mesh();
    value_type pi = M_PI;
    auto gradg = trans( +pi*cos( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() )*unitX()+
                        -pi*sin( pi*Px() )*sin( pi*Py() )*cos( pi*Pz() )*unitY()+
                        -pi*sin( pi*Px() )*cos( pi*Py() )*sin( pi*Pz() )*unitZ() );

    double semi_H1error =integrate( elements( mesh ),
                                    ( gradv( u )-gradg )*trans( ( gradv( u )-gradg ) ) ).evaluate()( 0,0 );
    double L2error2 = std::pow( l2Error( u ) , 2 );
    double error = math::sqrt( L2error2 + semi_H1error );

    return error;
}



void
ThermalBlock::exportResults( element_type& u )
{
    Log() << "exportResults starts\n";
    auto exporter = export_type::New( M_vm, "thermalblock" );
    exporter->step( 0 )->setMesh( u.mesh() );
    exporter->step( 0 )->add( "u", u );
    exporter->save();

    Log() << "exportResults done\n";
} // ThermalBlock::export



void
ThermalBlock::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute ThermalBlock\n";
    std::vector<double> X( 1 );
    X[0] = meshSize;
    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}


void
ThermalBlock::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    std::cout<<"[thermalblock.hpp] WARNING : function run is to implement, you can't use it"<<std::endl;
    exit( 0 );
} // ThermalBlock::run


void
ThermalBlock::checkout()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute ThermalBlock\n";
    std::vector<double> X( 1 );
    X[0] = meshSize;
    std::vector<double> Y( 3 );
    checkout( X.data(), X.size(), Y.data(), Y.size() );
}



void
ThermalBlock::checkout( const double* X, unsigned long P, double* Y, unsigned long N )
{
    meshSize=X[0];

    //std::srand(static_cast<unsigned>(std::time(0)));

    std::map<std::string,double> K;

    for ( auto it=domainMarkers.begin(), en=domainMarkers.end(); it!=en; ++it )
    {
        //        K.insert( std::make_pair( *it, std::exp( std::log(0.1)+(std::log(10)-std::log(0.1))*((double(std::rand())/double(RAND_MAX))) ) ) );
        K.insert( std::make_pair( *it, 1. ) );
    }

    value_type pi = M_PI;
    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );
    auto f = pi*pi*2*g;

    auto gradg = trans( +pi*cos( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() )*unitX()+
                        -pi*sin( pi*Px() )*sin( pi*Py() )*cos( pi*Pz() )*unitY()+
                        -pi*sin( pi*Px() )*cos( pi*Py() )*sin( pi*Pz() )*unitZ() );

    auto mygrad =  gradg*vf::N();

    element_type u( Xh, "u" );

    element_type v( Xh, "v" );

    auto B = M_backend->newVector( Xh );

    std::cout << "Building rhs restarts\n";
    form1( _test=Xh,_vector=B, _init=true ) ;

    BOOST_FOREACH( auto domain, domainMarkers )
    {
        std::cout << domain << " K[" << domain << "]=" << K[domain] << std::endl;
        form1( Xh, B ) += integrate( markedelements( mmesh, domain ),K[domain]*pi*pi*2*g*id( v ), _Q<10>() );
    }

    BOOST_FOREACH( auto marker, southMarkers )
    {
        // std::cout << "Building rhs for marker " << marker << std::endl;
        std::string sid = subdomainFromBoundary( marker );
        std::cout << marker << " K[" << sid << "]=" << K[sid] << std::endl;
        form1( Xh, B ) +=  integrate( markedfaces( mmesh, marker ),K[sid]*mygrad*id( v ), _Q<10>() )  ;
    }

    BOOST_FOREACH( auto marker, northMarkers )
    {
        // std::cout << "Building rhs for marker " << marker << std::endl;
        std::string sid = subdomainFromBoundary( marker );
        std::cout << marker << " K[" << sid << "]=" << K[sid] << std::endl;
        form1( Xh, B ) +=  integrate( markedfaces( mmesh, marker ), g*( -K[sid]*grad( v )*vf::N()+gamma_dir*id( v )/hFace() ) ), _Q<10>()  ;
    }

    BOOST_FOREACH( auto marker, westMarkers )
    {
        // std::cout << "Building rhs for marker " << marker << std::endl;
        std::string sid = subdomainFromBoundary( marker );
        std::cout << marker << " K[" << sid << "]=" << K[sid] << std::endl;
        form1( Xh, B ) +=  integrate( markedfaces( mmesh, marker ),K[sid]*mygrad*id( v ), _Q<10>() )  ;
    }

    BOOST_FOREACH( auto marker, eastMarkers )
    {
        // std::cout << "Building rhs for marker " << marker << std::endl;
        std::string sid = subdomainFromBoundary( marker );
        std::cout << marker << " K[" << sid << "]=" << K[sid] << std::endl;
        form1( Xh, B ) +=  integrate( markedfaces( mmesh, marker ), K[sid]*mygrad*id( v ),_Q<10>() )  ;
    }


    std::cout << "Building rhs done\n";

    auto A = M_backend->newMatrix( Xh, Xh );
    std::cout << "Building lhs restarts\n";
    form2( _test=Xh, _trial=Xh, _matrix=A, _init=true );

    BOOST_FOREACH( auto domain, domainMarkers )
    {
        // std::cout << "Building lhs for domain " << domain << std::endl;
        std::cout << domain << " K[" << domain << "]=" << K[domain] << std::endl;
        form2( Xh, Xh, A ) += integrate( markedelements( mmesh, domain ), K[domain]*gradt( u )*trans( grad( v ) ) );
    }
    BOOST_FOREACH( auto marker, northMarkers )
    {
        // std::cout << "Building lhs for marker " << marker << std::endl;
        std::string sid = subdomainFromBoundary( marker );
        std::cout << marker << " K[" << sid << "]=" << K[sid] << std::endl;
        // std::cout << "ident " << sid << std::endl;
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=  integrate( markedfaces( mmesh, marker ),
                -K[sid]*gradt( u )*vf::N()*id( v )
                -K[sid]*grad( u )*vf::N()*idt( v )
                +gamma_dir*idt( u )*id( v )/hFace() );
    }
    std::cout << "Building lhs done\n";

    std::cout<< "solve restarts " << "\n";

    M_backend->solve( _matrix=A, _solution=u, _rhs=B, _reuse_prec=false );

    std::cout<< "solve done " << "\n";

    double L2error2 = 0;

    double temp = 0;

    std::cout << "Compute error restarts\n";

    BOOST_FOREACH( auto domain, domainMarkers )
    {
        L2error2 += integrate( markedelements( mmesh, domain ), ( idv( u )-g )*( idv( u )-g ), _Q<15>() ).evaluate()( 0,0 );
        temp = integrate( markedelements( mmesh, domain ), ( idv( u )-g )*( idv( u )-g ), _Q<15>() ).evaluate()( 0,0 );
        std::cout << "L2error(" << domain <<")=" << math::sqrt( temp ) << "\n";
        temp = 0;
    }

    double L2error =   math::sqrt( L2error2 );

    std::cout << "||error||_L2 = " << L2error << "\n";

    std::cout << "Compute error done\n";
    //  double myerror = 0;

    //  BOOST_FOREACH( auto marker, eastMarkers )
    //  {
    //      myerror += integrate(markedfaces(mmesh, marker), cst(1.) ).evaluate()(0,0);
    //  }

    // std::cout << "myerror=" << myerror << "\n";

    export_ptrtype exporter( export_type::New( M_vm,
                             ( boost::format( "thermalblock-%1%-%2%" )
                               % Order
                               % 2 ).str() ) );

    if ( exporter->doExport() )
    {
        Log() << "exportResults starts\n";

        exporter->step( 0 )->setMesh( mmesh );

        exporter->step( 0 )->add( "u", u );
        // exporter->step(0)->add( "g", e );
        exporter->save();
        Log() << "exportResults done\n";
    }


} // ThermalBlock::checkout

} // Feel

#endif /* __ThermalBlock_H */



