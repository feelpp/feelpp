//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 14 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!

#include <feel/feelmor/crbplugin.hpp>
#include <ad.hpp>

namespace Feel {

po::options_description
makeAdvectionDiffusionOptions()
{
    po::options_description AdvectionDiffusionoptions( "AdvectionDiffusion options" );
    AdvectionDiffusionoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "mu1", po::value<double>()->default_value( 1 ), "lenght of the channel in [1;10]" )
    ( "mu2", po::value<double>()->default_value( 0.1 ), "Peclet number in [0.1;100]" )
    ( "no-export", "don't export results" )
    ;
    return AdvectionDiffusionoptions;
}
AboutData
makeAdvectionDiffusionAbout( std::string const& str )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "2D steady Advection-Diffusion",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2011 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

AdvectionDiffusion::AdvectionDiffusion()
    :
    super_type("ad"),
    meshSize( 0.01 ),
    M_do_export( true ),
    export_number( 0 )
{
    this->setPluginName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_NAME));
    this->setPluginLibName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_LIBNAME) );
}


AdvectionDiffusion::AdvectionDiffusion( po::variables_map const& vm )
    :
    super_type( "ad" ),
    M_vm( vm ),
    meshSize( vm["hsize"].as<double>() ),
    M_do_export( !vm.count( "no-export" ) ),
    export_number( 0 )
{
    this->setPluginName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_NAME));
    this->setPluginLibName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_LIBNAME) );
}
void
AdvectionDiffusion::initModel()
{
    using namespace Feel::vf;using Feel::vf::id;
    // geometry is a ]0,1[x]0,1[
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R( meshSize,"Omega",x1,x2 );
    R.setMarker( _type="line",_name="Inflow",_marker4=true );
    R.setMarker( _type="line",_name="Bottom",_marker1=true );
    R.setMarker( _type="line",_name="Top",_marker3=true );
    R.setMarker( _type="line",_name="Outflow",_marker2=true );
    R.setMarker( _type="surface",_name="Omega",_markerAll=true );
    mesh = R.createMesh( _mesh=new mesh_type, _name="Omega" );
    
    /*
     * The function space and some associate elements are then defined
     */
    auto Xh = space_type::New( mesh );
    this->setFunctionSpaces( Xh );
    if (Environment::isMasterRank() )
        std::cout << "Number of dof : "<< Xh->nDof() << std::endl;
    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    //  initialisation de A1 et A2
    auto mat_graph = stencil( _test=Xh,_trial=Xh)->graph();
    M_Aqm.resize( Qa() );
    for(int q=0; q<Qa(); q++)
    {
        M_Aqm[q].resize( 1 );
        M_Aqm[q][0] = backend()->newMatrix(0,0,0,0,mat_graph);
    }

    M_Fqm.resize( this->Nl() );
    for(int l=0; l<Nl(); l++)
    {
        M_Fqm[l].resize( Ql(l) );
        for(int q=0; q<Ql(l) ; q++)
        {
            M_Fqm[l][q].resize(1);
            M_Fqm[l][q][0] = backend()->newVector( Xh );
        }
    }

    D = backend()->newMatrix(0,0,0,0,mat_graph);
    F = backend()->newVector( Xh );
    Dmu->setDimension( 2 );
    auto mu_min = Dmu->element();
    mu_min << 1, 0.1;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 10,100;
    Dmu->setMax( mu_max );

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    std::cout << "Number of dof " << Xh->nLocalDof() << "\n";

    // right hand side
    form1( _test=Xh, _vector=M_Fqm[0][0][0] )
        = integrate( _range = markedfaces( mesh, "Bottom" ), _expr = id( v ) );
    M_Fqm[0][0][0]->close();

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        = integrate( _range = elements( mesh ), _expr = Py()*dxt( u )*id( v ) );
    M_Aqm[0][0]->close();

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[1][0] )
        = integrate( _range = elements( mesh ), _expr = dxt( u )*dx( v ) );
    M_Aqm[1][0]->close();

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )
        = integrate( _range = elements( mesh ), _expr = dyt( u )*dy( v ) );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )
        += integrate( _range = markedfaces( mesh,"Top" ),
                      _expr = - dyt( u )*Ny()*id( v ) - dy( u )*Ny()*idt( v ) + 20*idt( u )*id( v )/hFace() );
    M_Aqm[2][0]->close();

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[3][0] )
        = integrate( _range = markedfaces( mesh,"Inflow" ),
                     _expr = - dxt( u )*Nx()*id( v ) - dx( u )*Nx()*idt( v ) );
    M_Aqm[3][0]->close();

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[4][0] )
        = integrate( _range = markedfaces( mesh,"Inflow" ),
                     _expr = 20*idt( u )*id( v )/hFace() );
    M_Aqm[4][0]->close();

    M = backend()->newMatrix(0,0,0,0,mat_graph);

    form2( _test=Xh, _trial=Xh, _matrix=M )
        = integrate( _range = elements( mesh ), _expr = id( u )*idt( v ) + grad( u )*trans( gradt( u ) ) );
    M->close();

} // AdvectionDiffusion::run

// AdvectionDiffusion::sparse_matrix_ptrtype
// AdvectionDiffusion::newMatrix() const
// {
//     return backend()->newMatrix( Xh, Xh );
// }

// AdvectionDiffusion::vector_ptrtype
// AdvectionDiffusion::newVector() const
// {
//     return backend()->newVector( Xh );
// }

AdvectionDiffusion::affine_decomposition_type
AdvectionDiffusion::computeAffineDecomposition()
{
    return boost::make_tuple( M_Aqm, M_Fqm  );
}


void
AdvectionDiffusion::solve( sparse_matrix_ptrtype& D,
                           element_type& u,
                           vector_ptrtype& F )
{
    backend()->solve( _matrix=D, _solution=u, _rhs=F );
} // AdvectionDiffusion::solve


void
AdvectionDiffusion::exportResults( element_type& U , parameter_type const& mu )
{

    if ( M_do_export )
    {
        LOG(INFO) << "exportResults starts\n";

        std::string exp_name;
        export_ptrtype exporter;
        std::string mu_str;

        for ( int i=0; i<mu.size(); i++ )
        {
            mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
        }

        exp_name = "solution_with_parameters_" + mu_str;

        exporter = export_ptrtype( Exporter<mesh_type>::New( "ensight", exp_name  ) );
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->add( "u", U );
        exporter->save();
    }
} // AdvectionDiffusion::export

void
AdvectionDiffusion::update( parameter_type const& mu )
{
    *D = *M_Aqm[0][0];

    for ( size_type q = 1; q < M_Aqm.size(); ++q )
    {
        for ( size_type m = 0; m < mMaxA(q); ++m )
        {
            D->addMatrix( M_betaAqm[q][m] , M_Aqm[q][m] );
            //std::cout << "[affine decomp] scale q=" << q << " with " << M_betaAqm[q][0] << "\n";
        }
    }

    F->close();
    F->zero();

    for ( size_type q = 0; q < M_Fqm[0].size(); ++q )
    {
        for ( size_type m = 0; m < mMaxF(0,q); ++m )
        {
            //std::cout << "[affine decomp] scale q=" << q << " with " << M_betaFqm[0][q][0] << "\n";
            F->add( M_betaFqm[0][q][m], M_Fqm[0][q][m] );
        }
    }
}

AdvectionDiffusion::element_type
AdvectionDiffusion::solve( parameter_type const& mu )
{
    //std::cout << "solve(mu) for parameter " << mu << "\n";

    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
    this->exportResults( *T, mu );
    return *T;
    //this->exportResults( *T );

}

void
AdvectionDiffusion::solve( parameter_type const& mu, element_ptrtype& T )
{
    this->computeBetaQm( mu );
    this->update( mu );
    backend()->solve( _matrix=D,  _solution=T, _rhs=F );
    export_number++;
#if 0
    std::ofstream file;
    std::string mu_str;

    for ( int i=0; i<mu.size(); i++ )
    {
        mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
    }

    std::string number =  ( boost::format( "Exp_%1%" ) %export_number ).str();
    std::string name = "PFEMsolution" + mu_str + number;
    file.open( name,std::ios::out );

    for ( int i=0; i<T->size(); i++ ) file<<T->operator()( i )<<"\n";

    file.close();


    std::cout<<"pfem solution ok"<<std::endl;
    std::ofstream file_matrix;
    name = "PFEMmatrix" + mu_str + number;
    file_matrix.open( name,std::ios::out );
    file_matrix<<*D;
    file_matrix.close();

    std::cout<<"pfem matrix ok"<<std::endl;
    name = "PFEMrhs" + mu_str + number;
    std::ofstream file_rhs;
    file_rhs.open( name,std::ios::out );
    file_rhs<<*F;
    file_rhs.close();
    std::cout<<"pfem rhs ok"<<std::endl;
#endif

}

void
AdvectionDiffusion::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    backend()->solve( _matrix=M,  _solution=u, _rhs=f );
    //std::cout << "l2solve(u,f) done\n";
}

double
AdvectionDiffusion::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}
double
AdvectionDiffusion::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}

void
AdvectionDiffusion::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    //std::cout<<"RUN ::::::::::: "<<std::endl;
    using namespace vf;
    auto mu = Dmu->element();
    mu << X[0], X[1];
    static int do_init = true;

    if ( do_init )
    {
        meshSize = X[2];
        this->initModel();
        do_init = false;
    }

    this->solve( mu, pT );


    Y[0]=M_betaFqm[0][0][0]*integrate( _range = markedfaces( mesh, "Bottom" ), _expr = idv( *pT ) ).evaluate()( 0,0 );
}



double
AdvectionDiffusion::output( int output_index, parameter_type const& mu, element_type &u, bool need_to_solve )
{
    using namespace vf;
    if( need_to_solve )
        this->solve( mu, pT );
    else
        *pT = u;

    double output=0;

    // right hand side (compliant)
    if ( output_index == 0 )
    {
        //output = M_betaFqm[0][0][0]*dot( M_Fqm[0][0][0], U );
        for ( int q=0; q<Ql( output_index ); q++ )
        {
            for ( int m=0; m<mMaxF(output_index,q); m++ )
            {
                //element_ptrtype eltF( new element_type( Xh ) );
                //*eltF = *M_Fqm[output_index][q][m];
                //output += M_betaFqm[output_index][q][m]*dot( *eltF, *pT );
                output += M_betaFqm[output_index][q][m]*dot( *M_Fqm[output_index][q][m] , *pT );
            }
        }
    }

    else
    {
        throw std::logic_error( "[AdvectionDiffusion::output] error with output_index : only 0 " );
    }

    return output;
}

FEELPP_CRB_PLUGIN( AdvectionDiffusion, ad )

}

