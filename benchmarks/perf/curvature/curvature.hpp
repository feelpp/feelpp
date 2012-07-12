/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2012-07-11

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file curvature.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-01-04
 */
#if !defined( __FEELPP_BENCH_CURVATURE_HPP)
#define __FEELPP_BENCH_CURVATURE_HPP 1

#include <boost/any.hpp>
#include <boost/utility.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/simget.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>



#include <feel/feelmesh/elements.hpp>

#include <feel/feelvf/vf.hpp>

namespace Feel
{
/**
 * \class Curvature class
 * \brief solves the curvature equations
 *
 */
template<int Dim,
         typename BasisU,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class Curvature
    :
public Simget
{
    typedef Simget super;
public:

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*mesh*/
    typedef Entity<Dim,1,Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    //# marker1 #
    typedef BasisU basis_u_type;
    typedef bases<basis_u_type> basis_type;
    //# endmarker1 #

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Curvature( std::string const& basis_name,
               po::variables_map const& vm, AboutData const& ad )
        :
        super( vm, ad ),
        M_backend(),
        M_basis_name( basis_name ),
        exporter()
    {
        mu = this->vm()["mu"].template as<value_type>();
        penalbc = this->vm()["bccoeff"].template as<value_type>();
    }


    std::string name() const
    {
        return M_basis_name;
    }

    /**
     * run the convergence test
     */
    void run();
    void run( const double* X, unsigned long P, double* Y, unsigned long N )
    {
        run();
    }

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type& v );

private:

    backend_ptrtype M_backend;
    std::string M_basis_name;
    double mu;
    double penalbc;

    boost::shared_ptr<export_type> exporter;
}; // Curvature


template<int Dim, typename BasisU, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Curvature<Dim, BasisU, Entity>::run()
{
    using namespace Feel::vf;

    int nparts = Environment::worldComm().size();
    bool prepare = this->vm()["benchmark.prepare"].template as<bool>();
    if ( prepare )
        nparts = this->vm()["benchmark.partitions"].template as<int>();

    if ( this->vm().count( "nochdir" ) == false )
    {
        this->changeRepository( boost::format( "perf/%1%/%2%/%3%/h_%4%/l_%5%/parts_%6%/" )
                                % this->about().appName()
                                % convex_type::name()
                                % M_basis_name
                                % meshSizeInit()
                                % level()
                                % nparts );
    }

    //! init backend
    M_backend = backend_type::build( this->vm() );
    exporter =  boost::shared_ptr<export_type>( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) );

    boost::mpi::timer t;

    double shear = this->vm()["shear"].template as<value_type>();
    bool recombine = this->vm()["recombine"].template as<bool>();

    /*
     * First we create the mesh, in the case of quads we wish to have
     * non-regular meshes to ensure that we don't have some super-convergence
     * associated to the mesh. Hence we use recombine=true to recombine
     * triangles generated by a Delaunay algorithm into quads
     */
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "circle" % Dim % 1 ).str() ,
                                              _shape="ellipsoid",
                                              _usenames=true,
                                              _dim=Dim,
                                              _h=meshSizeInit(),
                                              _shear=shear,
                                              _xmin=-1.,_xmax=1.,
                                              _ymin=-1.,_ymax=1. ),
                                _refine=level(),
                                _partitions=nparts );

    M_stats.put( "t.init.mesh",t.elapsed() );
    t.restart();
    if ( prepare ) return;

} // Curvature::run


template<int Dim, typename BasisU, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Curvature<Dim, BasisU, Entity>::exportResults( element_type& u, element_type& v )
{
    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( u.functionSpace()->mesh() );
        exporter->step( 0 )->add( "u", u );
        exporter->step( 0 )->add( "u_exact", v );
        exporter->save();
    }
} // Curvature::export
} // Feel

#endif // __FEELPP_BENCH_CURVATURE_HPP

