/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2007-08-27

  Copyright (C) 2007 Unil

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file pbeq.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2007-08-27
 */

#include "pbeq.hpp"

//#include <boost/numeric/ublas/expression.hpp>

namespace Feel
{

void
Pbeq::loadFE()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    changeToMeshRepository();

    /*
     * logs will be in <feel repo>/<app name>/<entity>/P<p>/h_<h>
     */
    this->setLogs();

    /**
     * load mesh and spaces
     */

    // the reference compuational domain is composed by a dense sphere
    M_pbeqspace.setMeshSize ( this->vm()["hsize"].as<value_type>() );
    M_pbeqspace.setFarFactor( this->vm()["farfactor"].as<value_type>() );
    M_pbeqspace.setFarBnd   ( this->vm()["farBnd"].as<value_type>() );


    if ( Application::processId() == 0 )
    {
        if ( !( this->vm().count( "meshReuse" ) && M_pbeqspace.loadMesh( ) ) )
            M_pbeqspace.createMesh( this->vm().count( "geoOnly" ) );
    }

    application_type::barrier();

    if ( Application::processId() != 0 )
    {
        if ( !M_pbeqspace.loadMesh( ) )
            throw std::logic_error( "could not load mesh on processor != 0" );
    }
}

template<typename Mat, typename Vec1, typename Vec2>
void
Pbeq::solve( Mat const& D,
             Vec1& u,
             Vec2 const& F,
             bool is_sym  )
{
    timers["solver"].first.restart();

    //M_backend->set_symmetric( is_sym );

    vector_ptrtype U( M_backend->newVector(  M_pbeqspace.Xh() ) );
    M_backend->solve( D, D, U, F, false );
    u = *U;

    timers["solver"].second = timers["solver"].first.elapsed();
    LOG(INFO) << "[timer] solve: " << timers["solver"].second << "\n";
} // Pbeq::solve

template<typename Etype, typename Ktype>
void
Pbeq::getPBMatrix( sparse_matrix_ptrtype& PB,
                   element_type const& u,
                   element_type const& v,
                   Etype const& Epsilon,
                   Ktype const& kappa2 /*,
                                        vector_ptrtype& F */ )
{
    using namespace Feel::vf;

    im_type im;

    size_type pattern = Pattern::COUPLED;
    form2( M_pbeqspace.Xh(), M_pbeqspace.Xh(), PB, _init=true, _pattern=pattern ) =
        integrate( elements( *M_pbeqspace.mesh() ), im,
                   Epsilon*( M_jacInvStr2[0] * dxt( u ) * dx( v ) +
                             M_jacInvStr2[1] * dyt( u ) * dy( v ) +
                             M_jacInvStr2[2] * dzt( u ) * dz( v ) )
                   +  kappa2*idt( u )*id( v ) );

    PB->close();

    if ( this->vm().count( "export-matlab" ) )
        PB->printMatlab( "PBnobc.m" );

#ifdef PBEQ_USE_PETSC
    /*
    value_type penalisation_bc = this->vm()["penalbc"].as<value_type>();

    form2( M_pbeqspace.Xh(), M_pbeqspace.Xh(), PB ) +=
        integrate( boundaryfaces(*M_pbeqspace.mesh()), im,
                   ( - trans(id(v))*(gradt(u)*N())
                     - trans(idt(u))*(grad(v)*N())
                     + penalisation_bc*trans(idt(u))*id(v)/hFace()) );

    */

#else
    // input/output: vector_ptrtype& F
    /*
      form2( M_pbeqspace.Xh(), M_pbeqspace.Xh(), PB ) +=
        on( boundaryfaces(*M_pbeqspace.mesh()), u, *F, constant(0.0) );
    */
#endif

    if ( this->vm().count( "export-matlab" ) )
        PB->printMatlab( "PB.m" );

} // getPBMatrix





int
Pbeq::initMolecule( std::string const& moleculeName,
                    molecule_type& _molecule ) const
{
    std::ostringstream file;
    file << this->vm()["molDir"].as<std::string>()
         << "/"
         << moleculeName;
    int mol_read( 1 );

    switch ( M_filetype )
    {
    case molecule_type::PQR:
        mol_read = _molecule.readPQRfile( file.str() );
        break;

    case molecule_type::CRD:
        mol_read = _molecule.readCRDfile( file.str() );
        break;
    }

    if ( mol_read == 0 )
    {
        LOG(INFO) << "Molecule " << moleculeName << " has " << _molecule.size() << " atoms \n";
        _molecule.showMe();
        return _molecule.size();
    }

    LOG(INFO) << "problem in reading "
          << file.str()
          << "\n";

    //_molecule.showMe();

    return 0;


} // initMolecule

void
Pbeq::setDomain( molecule_type const& _molecule )
{

    uint16_type const dim( pbeqspace_type::Dim );
    node_type _min( dim ),  _max( dim );
    node_type  stretch( dim ), translation( dim );
    FEELPP_ASSERT( this->vm()["farBnd"].as<value_type>() > 0 )( dim ).error( "invalid farBnd" );

    _molecule.domainMinMax( _min,_max );

    stretch     = ( _max - _min )/2;
    translation = ( _max + _min )/2;

    LOG(INFO) << "i  \tmin,   \tmax, \tstretch,\ttranslation " << "\n";

    for ( int i( 0 ); i < dim; i++ )
        LOG(INFO) << i <<"\t"<<  _min[i] <<"\t"<<  _max[i] <<"\t"<< stretch[i] <<"\t\t"<< translation[i]
              << "\n";

    M_detJac = stretch[0];

    for ( int i = 1; i < dim; i++ )
        M_detJac *= stretch[i];

    FEELPP_ASSERT( M_detJac > 0 )( dim ).error( "invalid sign of the Jacobian" );

    M_jacInvStr2.resize( dim );

    for ( int i = 0; i < dim; i++ )
        M_jacInvStr2[i] = M_detJac/( stretch[i]*stretch[i] );

    M_pbeqspace.setStretch( stretch );
    M_pbeqspace.setTranslation( translation );

}

void
Pbeq::setCoeffs()
{
    //Parameter values

    Kb   = this->vm()["Kb"].  as<value_type>();
    pdie = this->vm()["pdie"].as<value_type>();
    sdie = this->vm()["sdie"].as<value_type>();
    temp = this->vm()["temp"].as<value_type>();
    Kappab   = 8*3.1415926*this->vm()["IonStr"].  as<value_type>() / ( Kb*temp*sdie );


    // rhsCoeff = 4 \pi e_c / (Kb*T) * detJac
    ec = 0.0854245;
    rhsCoeff = 4.* 3.1415926 * ec / ( Kb * temp );

    M_pbeqspace.setSmoothWindow( this->vm()["smooth"].as<value_type>() );

}


int
Pbeq::setReceptor( std::string const& receptorName )
{

    changeToSolutionRepository( receptorName );
    setCoeffs();

    initMolecule( receptorName,M_receptor );

    setDomain( M_receptor );

    deltaAndHeavisyde( M_receptor, M_F_receptor, M_H_receptor );

    return M_receptor.size();

}

int
Pbeq::setLigand( std::string const& ligandName )
{
    return initMolecule( ligandName,M_ligand );
}

int
Pbeq::setReclusterFile( std::string const& reclusterName )
{

    std::ostringstream file;
    file << this->vm()["molDir"].as<std::string>()
         << "/"
         << reclusterName
         << ".pdb.dock4";

    M_dock4file.open( file.str().c_str(), std::ifstream::in );

    if ( M_dock4file.fail() )
    {
        LOG(INFO) << "setReclusterFile::failed to open *" << reclusterName << "*"
              << "\n at " << file.str()
              << "\n";
        return -1;
    }

    if ( M_dock4file.eof() ) return -1;

    return 0;
}

int
Pbeq::recluster()
{
    return M_ligand.readDock4file( M_dock4file );

}

void
Pbeq::deltaAndHeavisyde( molecule_type const       &myMolecule,
                         vector_ptrtype            &F,
                         heavyside_element_ptrtype &H )
{
    using namespace Feel::vf;

    /*
     * Construction of the right hand side
     */

    timers["assembly"].first.restart();

    F = M_backend->newVector( M_pbeqspace.Xh() );

    M_pbeqspace.intvrho( myMolecule,F );
    F->close();

    // Checking consistency of the right hand side
    LOG(INFO) << myMolecule.name() << " total charge = " << myMolecule.totalCharge() << " == " << F->sum() << " = sum(F)\n";

    FEELPP_ASSERT(  std::abs( myMolecule.totalCharge() -  F->sum() )
                    < ( myMolecule.totalAbsCharge() +  F->l1Norm() )  * 1.e-9  ).error( "total charge and sum(F) are different" );


    F->scale( rhsCoeff* M_detJac );

    timers["assembly"].second = timers["assembly"].first.elapsed();
    LOG(INFO) << "[timer] rhs assembly: " << timers["assembly"].second << "\n";

    // right hand side done


    /*
     * Construction of the heavyside function
     */

    timers["heavyside"].first.restart();

    H.reset( new heavyside_element_type( M_pbeqspace.fastHeavyside( myMolecule ) ) );
    //H.reset( new heavyside_element_type( M_pbeqspace.heavyside( myMolecule )) );

    // Checking consistency (to be checked agains other hsize's)
    im_type im;
    double domain_size = integrate( elements( *M_pbeqspace.mesh() ), im, constant( 1.0 ) ).evaluate()( 0,0 );
    double intH        = integrate( elements( *M_pbeqspace.mesh() ), im, idv( *H ) ).evaluate()( 0,0 );

    //

    LOG(INFO) << "domain_size = " << M_detJac*domain_size
          << " intH = " << M_detJac*intH
          << "\n";



    timers["heavyside"].second += timers["heavyside"].first.elapsed();
    LOG(INFO) << "[timer] heavyside: " << timers["heavyside"].second << "\n";

}




template<typename Etype, typename Ktype>
value_type
Pbeq::buildSystemAndSolve( element_type        & u,
                           element_type   const& v,
                           Etype          const& Epsilon,
                           Ktype          const& kappa2,
                           vector_ptrtype const& F )
{

    /* vector_ptrtype  F( M_backend->newVector(M_pbeqspace.Xh()) );
    *F = *_F;
    */

    timers["assembly"].first.restart();

    sparse_matrix_ptrtype PB;
    PB = M_backend->newMatrix( M_pbeqspace.Xh(), M_pbeqspace.Xh() );

    getPBMatrix( PB, u, v, Epsilon, kappa2 );

    timers["assembly"].second = timers["assembly"].first.elapsed();
    LOG(INFO) << "[timer] assembly: " << timers["assembly"].second << "\n";


    /*
     * solving the system
     */

    this->solve( PB, u, F, false );

    if ( this->vm().count( "export-matlab" ) )
    {
        std::ostringstream ostr;
        ostr.precision( 3 );
        ostr << "F.m" ;
        F->printMatlab( ostr.str() );

        if ( Application::processId() == 0 )
            std::cout  << "Returning now, only checking matrices \n"
                       << std::flush;

        exit( 0 );

    }

    return postProcess( u, F, rhsCoeff );

}

value_type
Pbeq::solveReceptor( std::string const& receptorName )
{
    using namespace Feel::vf;

    setReceptor( receptorName );
    FEELPP_ASSERT(  M_receptor.size() > 0 ).error( "Receptor has zero length" );

    element_type u( M_pbeqspace.Xh(), "u" );
    element_type v( M_pbeqspace.Xh(), "v" );

    /*
     * Construction of the left hand side components
     */

    value_type jacKbsquare ( M_detJac*Kappab*Kappab );
    AUTO( Epsilon, pdie + ( sdie-pdie )* idv( *M_H_receptor ) );
    AUTO( kappa2, jacKbsquare* idv( *M_H_receptor ) );

    // unused. Why?
    value_type K ( Kb/std::sqrt( sdie ) );

    /*
     *  solvated solution
     */

    value_type solvatedEnergy = buildSystemAndSolve( u, v, Epsilon, kappa2, M_F_receptor );


    /*
     *  Vacuum solution
     */
    v=u;
    value_type  vacuumEnergy = buildSystemAndSolve( v, u, pdie, kappa2, M_F_receptor );

    /*
     *  export results
     */

    if ( this->vm().count( "export" ) )
        exportResults( *M_H_receptor,*M_H_receptor,u,v, M_F_receptor );

    return solvatedEnergy - vacuumEnergy;

} // end solveReceptor

value_type
Pbeq::solveLigand()
{

    FEELPP_ASSERT(  M_receptor.size() > 0 ).error( "Receptor has zero length" );

    using namespace Feel::vf;

    vector_ptrtype F;
    heavyside_element_ptrtype H;

    deltaAndHeavisyde( M_ligand, F, H );

    *F += *M_F_receptor;

    /*
     * Construction of the left hand side components
     */

    value_type jacKbsquare ( M_detJac*Kappab*Kappab );
    AUTO( Epsilon, pdie + ( sdie-pdie ) * idv( *H ) * idv( *M_H_receptor ) );
    AUTO( kappa2, jacKbsquare * idv( *H ) * idv( *M_H_receptor ) );

    // unused. Why?
    value_type K ( Kb/std::sqrt( sdie ) );


    // solutions

    element_type u( M_pbeqspace.Xh(), "u" );
    element_type v( M_pbeqspace.Xh(), "v" );

    /*
     *  solvated solution
     */

    value_type solvatedEnergy = buildSystemAndSolve( u, v, Epsilon, kappa2, F );

    /*
     *  Vacuum solution
     */

    v=u;
    value_type  vacuumEnergy = buildSystemAndSolve( v, u, pdie, kappa2, F );

    /*
     *  export results
     */

    if ( this->vm().count( "export" ) )
        exportResults( *M_H_receptor,*H,u,v,F );

    return solvatedEnergy - vacuumEnergy;

} // end solveLigand


value_type
Pbeq::postProcess( element_type const& u, vector_ptrtype const& F, value_type const& rhsCoeff ) const
{
    using namespace Feel::vf;

    // Computing Energy and norms

    timers["output"].first.restart();

    /*
    im_type im;
    value_type l2norm = integrate( elements(*M_pbeqspace.mesh()), im, val( idv(u)^2  ) ).evaluate()(0,0);
    // value_type global_l2norm = 0;
    // mpi::all_reduce( Application::comm(), l2norm, global_l2norm, std::plus<value_type>() );
    // LOG(INFO) << "L2norm = " << math::sqrt( global_l2norm ) << "\n";
    LOG(INFO) << "L2norm = " << math::sqrt( l2norm ) << "\n";
    */

#ifdef PBEQ_USE_SERIAL
    //value_type output(boost::numeric::ublas::inner_prod( *F, u ));
    value_type output = F->dot( u );
#else
    value_type output( inner_product( *F, u ) );
#endif

    //\Delta G_{ele} = 1/2 K_B*T /ec * \int_\Omega u \rho
    /*
      if ( Application::processId() == 0 )
        std::cout << "inner_prod ( F, u ) = " << output << "\n" << std::flush;
    */
    LOG(INFO) << "Energy = " << 1./( 8.* 3.1415926 *rhsCoeff*rhsCoeff ) << "*" << output << "\n";

    output = 1./( 8.* 3.1415926 *rhsCoeff*rhsCoeff ) * output;

    LOG(INFO) << "Energy = " << output << "\n";
    /*
    if ( Application::processId() == 0 )
        std::cout <<  " Energy = " << output << "\n"
                  << std::flush;
    */

    timers["output"].second += timers["output"].first.elapsed();
    LOG(INFO) << "[timer] Energy: " << timers["output"].second << "\n";

    return output;

} // end postProcess

void
Pbeq::exportResults( heavyside_element_type& H0, heavyside_element_type& H1,
                     element_type& u, element_type& v,
                     vector_ptrtype F )
{
    timers["export"].first.restart();

    timeset_type::step_ptrtype timeStep = timeSet->step( M_countExports++ );
    timeStep->setMesh( H0.functionSpace()->mesh() );
    timeStep->add( "H0", H0 );
    timeStep->add( "H1", H1 );
    timeStep->add( "u", u );
    timeStep->add( "v", v );

    element_type f( M_pbeqspace.Xh(), "f" );
    f = *F;
    timeStep->add( "rho", f );


    exporter->save();

    timers["export"].second = timers["export"].first.elapsed();
    LOG(INFO) << "[timer] export: " << timers["export"].second << "\n";

} // Pbeq::exportResults


void
Pbeq::changeToMeshRepository()
{

    this->changeRepository( boost::format( "%1%/h_%2%/farF_%3%/farB_%4%/" )
                            % this->about().appName()
                            % this->vm()["hsize"].as<value_type>()
                            % this->vm()["farfactor"].as<value_type>()
                            % this->vm()["farBnd"].as<value_type>()
                          );
}

void
Pbeq::changeToSolutionRepository( std::string const& receptorName )
{
    this->changeRepository( boost::format( "%1%/h_%2%/farF_%3%/farB_%4%/%5%/" )
                            % this->about().appName()
                            % this->vm()["hsize"].as<value_type>()
                            % this->vm()["farfactor"].as<value_type>()
                            % this->vm()["farBnd"].as<value_type>()
                            % receptorName
                          );

}

} // end namespace Feel

