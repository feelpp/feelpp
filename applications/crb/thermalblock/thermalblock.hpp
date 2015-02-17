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
#ifndef FEELPP_ThermalBlock_HPP
#define FEELPP_ThermalBlock_HPP 1

#include <boost/timer.hpp>
#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>

#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelcrb/crbmodelbase.hpp>

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
    ( "geofile", Feel::po::value<std::string>()->default_value( "" ), "name of the geofile input (used to store DB)")
    ( "mshfile", Feel::po::value<std::string>()->default_value( "" ), "name of the gmsh file input")
    ( "gamma_dir", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary Dirichlet formulation" )
    ( "load-mesh-already-partitioned", Feel::po::value<bool>()->default_value( "true" ), "load a mesh from mshfile that is already partitioned if true, else the mesh loaded need to be partitioned")
    ;
    return thermalblockoptions;
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
                     "Copyright (C) 2009-2014 Feel++ Consortium");
    about.addAuthor( "Abdoulaye Samake", "main developer", "abdoulaye.samake@ujf-grenoble.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "contributor", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "contributor", "", "" );
    return about;

}



using namespace vf;


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
 * \class ThermalBlock
 *
 * ThermalBlock Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
class ThermalBlock : public CRBModelBase< ParameterDefinition , FunctionSpaceDefinition >
{

    static const uint16_type nx = 3;
    static const uint16_type ny = 3;

public:

    typedef CRBModelBase<ParameterDefinition, FunctionSpaceDefinition> super_type;

    static const uint16_type ParameterSpaceDimension = nx*ny-1;

    using super_type::computeBetaQm;

    //! initialisation of the model and definition of parameters values
    void initModel();

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


    double l2Error( element_type& u );

    double h1Error( element_type& u );

    void checkout();

    void checkout( const double* X, unsigned long P, double* Y, unsigned long N );

    std::string subdomainFromBoundary( std::string const& ) const;

    int subdomainId( std::string const& ) const;

    steady_betaqm_type computeBetaQm( parameter_type const& mu )
        {
            setBetaA(0) = 1;

            for ( int i=1; i<nx*ny; i++ )
                setBetaA(i) = mu(i-1);

            int index=nx*ny;

            //penalization term
            setBetaA(index) = 1;

            //compliant output
            setBetaF(0,0) = 1;

            return boost::make_tuple( betaAqm(), betaFqm() );
     }


    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);

    bool referenceParametersGivenByUser() { return true; }
    parameter_type refParameter()
    {
        auto muref = Dmu->element();
        muref(0)=1; muref(1)=1; muref(2)=1;
        muref(3)=1; muref(4)=1; muref(5)=1;
        muref(6)=1; muref(7)=1;
        return muref;
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

}; // ThermalBlock

gmsh_ptrtype
thermalBlockGeometry( int nx, int ny, double hsize )
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

    nameStr << "thermalblock";

    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}


void
ThermalBlock::initModel()
{

    gamma_dir=doption(_name="gamma_dir");

    std::string mshfile_name = option(_name="mshfile").as<std::string>();

    double hsize = doption(_name="hsize");

    if( mshfile_name=="" )
    {
        mmesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=thermalBlockGeometry( nx, ny, hsize ) );
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

    DVLOG(2) <<"[ThermalBlock::init] done allocating matrices/vectors \n";

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
    auto u = Xh->element();
    auto v = Xh->element();

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << Xh->nLocalDof() << std::endl;
        std::cout << "Number of dof " << Xh->nDof() << std::endl;
    }

    int index=0;//index for M_Fq or M_Aq
    int subdomain_index; //index for subdomains along a boundary

    auto M = backend()->newMatrix( Xh , Xh );
    M->zero();
    auto muref=refParameter();
    double muref_coeff=muref(0);
    form1( _test=Xh, _vector=ptrF(0,0) );
    // on boundary north we have u=0 so term from weak dirichlet condition
    // vanish in the right hand side
    BOOST_FOREACH( auto marker, southMarkers )
    {
        form1( _test=Xh, _vector= ptrF(0,0) ) += integrate( markedfaces( mmesh,marker ), id( v ) );
        subdomain_index = subdomainId( marker );
        south_subdomain_index.push_back( subdomain_index );
    }


    DVLOG(2) <<"[ThermalBlock::init] done with rhs\n";
    BOOST_FOREACH( auto domain, domainMarkers )
    {
        DVLOG(2) <<"[ThermalBlock::init] domain " << domain << "\n";

        form2( _test=Xh, _trial=Xh, _matrix=ptrA(index) ) =
            integrate( markedelements( mmesh, domain ), gradt( u )*trans( grad( v ) ) );

        if( index == 0 )
        {
            form2( _test=Xh, _trial=Xh, _matrix=M ) +=
                integrate( markedelements( mmesh, domain ), gradt( u )*trans( grad( v ) ) );
        }
        else
        {
            form2( _test=Xh, _trial=Xh, _matrix=M ) +=
                integrate( markedelements( mmesh, domain ), gradt( u )*trans( grad( v ) ) * muref_coeff );
        }
        DVLOG(2) <<"[ThermalBlock::init] done with Aqm[" << index << "]\n";

        index++;
    }


    DVLOG(2) <<"[ThermalBlock::init] done with domainMarkers\n";
    int last_index = nx*ny;
    BOOST_FOREACH( auto marker, northMarkers )
    {
        std::string sid = subdomainFromBoundary( marker );
        subdomain_index = subdomainId( marker );

        form2( _test=Xh, _trial=Xh, _matrix=ptrA(subdomain_index-1) )
            += integrate( markedfaces( mmesh, marker ),
                          -gradt( u )*vf::N()*id( v )
                          -grad( u )*vf::N()*idt( v ) );

        form2( _test=Xh, _trial=Xh, _matrix=M )
            +=  integrate( markedfaces( mmesh, marker ),
                           -gradt( u )*vf::N()*id( v ) * muref_coeff
                           -grad( u )*vf::N()*idt( v ) * muref_coeff );

        form2( _test=Xh, _trial=Xh, _matrix=ptrA(last_index) )
            += integrate( markedfaces( mmesh, marker ),gamma_dir*idt( u )*id( v )/h() );
        form2( _test=Xh, _trial=Xh, _matrix=M )
            += integrate( markedfaces( mmesh, marker ),gamma_dir*idt( u )*id( v )/h() );

        north_subdomain_index.push_back( subdomain_index );
    }

    this->addEnergyMatrix( M );
    DVLOG(2) <<"[ThermalBlock::init] done with boundaryMarkers\n";
}//initModel()


double
ThermalBlock::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{
    double output=0;

    if ( output_index==0 )
    {
        for ( int q=0; q<Ql( output_index ); q++ )
        {
            element_ptrtype eltF( new element_type( Xh ) );
            *eltF = *Fqm(output_index,q);
            output += betaF(output_index,q)*dot( *eltF, u );
        }
    }
    else
    {
        throw std::logic_error( "[ThermalBlock::output] error with output_index : only 0 " );
    }

    return output;
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
ThermalBlock::checkout()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute ThermalBlock\n";
    std::vector<double> X( 1 );
    X[0] = option(_name="hsize").template as <double>();
    std::vector<double> Y( 3 );
    checkout( X.data(), X.size(), Y.data(), Y.size() );
}



void
ThermalBlock::checkout( const double* X, unsigned long P, double* Y, unsigned long N )
{
    double meshSize=X[0];

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

    auto B = backend()->newVector( Xh );

    std::cout << "Building rhs restarts\n";
    form1( _test=Xh,_vector=B, _init=true ) ;

    BOOST_FOREACH( auto domain, domainMarkers )
    {
        std::cout << domain << " K[" << domain << "]=" << K[domain] << std::endl;
        form1( _test=Xh, _vector=B ) += integrate( markedelements( mmesh, domain ),K[domain]*pi*pi*2*g*id( v ), _Q<10>() );
    }

    BOOST_FOREACH( auto marker, southMarkers )
    {
        // std::cout << "Building rhs for marker " << marker << std::endl;
        std::string sid = subdomainFromBoundary( marker );
        std::cout << marker << " K[" << sid << "]=" << K[sid] << std::endl;
        form1( _test=Xh, _vector=B ) +=  integrate( markedfaces( mmesh, marker ),K[sid]*mygrad*id( v ), _Q<10>() )  ;
    }

    BOOST_FOREACH( auto marker, northMarkers )
    {
        // std::cout << "Building rhs for marker " << marker << std::endl;
        std::string sid = subdomainFromBoundary( marker );
        std::cout << marker << " K[" << sid << "]=" << K[sid] << std::endl;
        form1( _test=Xh, _vector=B ) +=  integrate( markedfaces( mmesh, marker ), g*( -K[sid]*grad( v )*vf::N()+gamma_dir*id( v )/hFace() ) ), _Q<10>()  ;
    }

    BOOST_FOREACH( auto marker, westMarkers )
    {
        // std::cout << "Building rhs for marker " << marker << std::endl;
        std::string sid = subdomainFromBoundary( marker );
        std::cout << marker << " K[" << sid << "]=" << K[sid] << std::endl;
        form1( _test=Xh, _vector=B ) +=  integrate( markedfaces( mmesh, marker ),K[sid]*mygrad*id( v ), _Q<10>() )  ;
    }

    BOOST_FOREACH( auto marker, eastMarkers )
    {
        // std::cout << "Building rhs for marker " << marker << std::endl;
        std::string sid = subdomainFromBoundary( marker );
        std::cout << marker << " K[" << sid << "]=" << K[sid] << std::endl;
        form1( _test=Xh, _vector=B ) +=  integrate( markedfaces( mmesh, marker ), K[sid]*mygrad*id( v ),_Q<10>() )  ;
    }


    std::cout << "Building rhs done\n";

    auto A = backend()->newMatrix( Xh, Xh );
    std::cout << "Building lhs restarts\n";
    form2( _test=Xh, _trial=Xh, _matrix=A, _init=true );

    BOOST_FOREACH( auto domain, domainMarkers )
    {
        // std::cout << "Building lhs for domain " << domain << std::endl;
        std::cout << domain << " K[" << domain << "]=" << K[domain] << std::endl;
        form2( _test=Xh, _trial=Xh, _matrix=A ) += integrate( markedelements( mmesh, domain ), K[domain]*gradt( u )*trans( grad( v ) ) );
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

    backend()->solve( _matrix=A, _solution=u, _rhs=B, _reuse_prec=false );

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
} // ThermalBlock::checkout

} // Feel

#endif /* __ThermalBlock_H */
