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
   \file pbeq.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2007-08-27
 */

#ifndef __PBEQ_H
#define __PBEQ_H

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporterensight.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelvf/vf.hpp>

#include "pbeqapplication.hpp"
#include "pbeqspace.hpp"


namespace Feel
{

class Pbeq
    :
public application_type

{
    typedef application_type super;
public:

    // -- TYPEDEFS --

    /* derived from  PbeqSpace*/
    typedef PbeqSpace pbeqspace_type;
    typedef pbeqspace_type::heavyside_element_type    heavyside_element_type;
    typedef boost::shared_ptr<heavyside_element_type> heavyside_element_ptrtype;
    typedef pbeqspace_type::mesh_type                 mesh_type;
    typedef pbeqspace_type::element_type              element_type;
    typedef pbeqspace_type::im_type                   im_type;

    /* derived from molecule */
    typedef pbeqspace_type::molecule_type molecule_type;
    typedef molecule_type::node_type node_type;


    /*matrix*/
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef Exporter<mesh_type>::timeset_type timeset_type;

    Pbeq( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        exporter( new ExporterEnsight<mesh_type>( "pbeq" ) ),
        timeSet( new timeset_type( "pbeq" ) ),
        M_countExports( 0 ),
        timers(),
        stats()
    {
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( "pbeq" );
    }

    /**
     * alias for run()
     */
    /*
    void operator()()
    {
        run();
    }
    */

    /**
     * run the simulation
     */
    void loadFE();
    value_type solveReceptor( std::string const& receptorName );

    value_type solveLigand();

    int setReceptor( std::string const& receptorName );

    int setLigand  ( std::string const& ligandName );
    int setReclusterFile( std::string const& reclusterFile );
    int recluster();

    void setUsePQR()
    {
        M_filetype = molecule_type::PQR;
    }

    void setUseCRD()
    {
        M_filetype = molecule_type::CRD;
    }

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( heavyside_element_type& H0, heavyside_element_type& H1,
                        element_type& u, element_type& v,
                        vector_ptrtype F );

    template<typename Mat, typename Vec1, typename Vec2>
    void solve( Mat const& D,
                Vec1& u,
                Vec2 const& F,
                bool is_sym  );

    template<typename Etype, typename Ktype>
    void getPBMatrix( sparse_matrix_ptrtype& PB,
                      element_type const& u,
                      element_type const& v,
                      Etype const& Epsilon,
                      Ktype const& kappa2 /*,
                                           vector_ptrtype& F*/ );

    value_type postProcess( element_type const& u, vector_ptrtype const& F, value_type const& rhsCoeff ) const;

    void changeToMeshRepository();
    void changeToSolutionRepository( std::string const& receptorName );

    int initMolecule( std::string const& moleculeName,
                      molecule_type& _molecule ) const;

    void setDomain( molecule_type const& _molecule );
    void setCoeffs();


    void deltaAndHeavisyde( molecule_type const       &myMolecule,
                            vector_ptrtype            &F,
                            heavyside_element_ptrtype &H );

    template<typename Etype, typename Ktype>
    value_type buildSystemAndSolve( element_type         &u,
                                    element_type   const &v,
                                    Etype          const &Epsilon,
                                    Ktype          const &kappa2,
                                    vector_ptrtype const &F );

private:

    backend_ptrtype M_backend;
    boost::shared_ptr<export_type> exporter;
    export_type::timeset_ptrtype   timeSet;
    int                            M_countExports;


    mutable std::map<std::string,std::pair<boost::timer,value_type> > timers;
    std::map<std::string,value_type> stats;

    pbeqspace_type M_pbeqspace;

    molecule_type M_receptor;
    molecule_type M_ligand;

    value_type M_detJac;
    node_type  M_jacInvStr2;

    // set by setCoeffs()

    value_type Kb;
    value_type Kappab;
    value_type pdie;
    value_type sdie;
    value_type temp;

    value_type ec;
    value_type rhsCoeff;


    // left hand side

    vector_ptrtype            M_F_receptor;
    heavyside_element_ptrtype M_H_receptor;

    molecule_type::FILETYPE   M_filetype;

    std::ifstream M_dock4file;

}; // Pbeq

} // end namespace Feel


#endif // __PBEQ_H

