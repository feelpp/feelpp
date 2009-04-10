/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2006-02-01

  Copyright (C) 2006 EPFL

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file polyvis.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2006-02-01
 */
#ifndef POLYVIS_HPP
#define POLYVIS_HPP

#include <string>
#include <boost/plugin/virtual_constructors.hpp>

#include <life/lifepoly/polynomialset.hpp>
#include <life/lifepoly/operations.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifefilters/pointsettomesh.hpp>
#include <life/lifefilters/exporterensight.hpp>

namespace Life
{
class Polyvis {
public:

    virtual ~Polyvis() {}

    template<typename P, template<uint16_type> class PS>
    void
    save( std::string const& s, PolynomialSet<P,PS> const& pset )
    {
        using namespace Life;
        const uint16_type nDim = P::nDim;
        const int OrderPts = 10;
        PointSetToMesh<Simplex<nDim,OrderPts>, double> p2m;
        typedef typename PointSetToMesh<Simplex<nDim,1>, double>::mesh_type mesh_type;

        //GaussLobatto<Simplex<nDim,OrderPts>, OrderPts, double> GS;
        PointSetEquiSpaced<Simplex<nDim,OrderPts>, OrderPts, double> GS;
        //p2m.addBoundaryPoints( false );
        p2m.visit( &GS );

        typedef FunctionSpace<mesh_type, Lagrange<1, Scalar> > P1_type;
        typedef boost::shared_ptr<P1_type> P1_ptrtype;
        P1_ptrtype P1( P1_type::New( p2m.mesh() ) );

        std::vector<boost::shared_ptr<typename P1_type::element_type> > u( pset.coeff().size1()/pset.nComponents );
        ublas::matrix<double> evalpset( pset.evaluate( GS.points() ) );
        std::cout << "evalpset = " << evalpset << "\n";
        // evaluate the dubiner polynomials at the lobatto nodes
        for ( size_type i = 0;i < evalpset.size1()/pset.nComponents; ++i )
            {
                u[i] =  boost::shared_ptr<typename P1_type::element_type>( new typename P1_type::element_type( P1 ) );
                //u[i] = P1.newElement( "u" );
                //u[i].resize( evalpset.size2() );

                if ( pset.is_scalar )
                    {
                        std::cout << "u[" << i << "]=" << ublas::row( evalpset, i ) << "\n";
                        *u[i] = ublas::row( evalpset, i );
                    }
                else
                    {
                        int nC = pset.nComponents;
                        ublas::matrix<double> m( ublas::project( evalpset,
                                                                 ublas::range( i*nC, (i+1)*nC ),
                                                                 ublas::range( 0, evalpset.size2() ) ) );
                        std::cout << "m = " << m << "\n"
                                  << "v2m(m) = " << PolynomialSet<P,PS>::polyset_type::toMatrix( m ) << "\n";
                        *u[i] = ublas::row( PolynomialSet<P,PS>::polyset_type::toMatrix( m ), 0 );
                    }
                std::cout << s << "_u_" << i << *u[i] << "\n";
                std::cout << s << "_p_" << i << pset.polynomial( i ).evaluate( GS.points() ) << "\n";
            }


        typedef ExporterEnsight<mesh_type> export_type;
        typedef typename ExporterEnsight<mesh_type>::timeset_type timeset_type;
        export_type __ensight( s );
        typename export_type::timeset_ptrtype __ts( new timeset_type( s ) );
        __ts->setTimeIncrement( 1.0);
        __ensight.addTimeSet( __ts );

        typename timeset_type::step_ptrtype __step = __ts->step( 1.0 );

        __step->setMesh( p2m.mesh() );

        for ( size_type i = 0;i < u.size(); ++i )
            {
                std::ostringstream str;
                str << s << "_p_" << i;
                if ( pset.is_scalar )
                    __step->add( str.str(), *u[i]  );
                else
                    __step->add( str.str(), *u[i]  );
            }

        __ensight.save();
    }
};
}
namespace boost { namespace plugin {

    template<>
    struct virtual_constructors<Life::Polyvis>
    {
        typedef mpl::list<
            mpl::list<std::string>,
            mpl::list<std::string, int>,
            mpl::list<std::string, int, int>
        > type;
    };

}}

#endif

