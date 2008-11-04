/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-21

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file exportergmsh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-21
 */
#include <life/lifecore/life.hpp>

#include <life/lifefilters/exportergmsh.hpp>

namespace Life
{
namespace detail
{
template<typename TimesetType, typename FunctionType>
void
gmsh_save_value( std::ostream& out,
                 boost::shared_ptr<TimesetType> const& __ts,
                 FunctionType const& __u )
{
    typedef typename FunctionType::value_type value_type;
    typedef typename matrix_node<value_type>::type matrix_node_type;
    out << "$View\n";

    // view-name nb-time-steps
    out << __u.name() << " " << 1  << "\n";


    size_type nLines = 0;
    if ( __u.mesh()->element( 0 ).isALineShape() )
        nLines = __u.mesh()->elements().size();

    size_type nTriangles = 0;
    if ( __u.mesh()->element( 0 ).isATriangleShape() )
        nTriangles = __u.mesh()->elements().size();

    size_type nQuad = 0;
    if ( __u.mesh()->element( 0 ).isAQuadrangleShape() )
        nQuad = __u.mesh()->elements().size();

    size_type nTetra = 0;
    if ( __u.mesh()->element( 0 ).isATetrahedraShape() )
        nTetra = __u.mesh()->elements().size();

    size_type nHexa = 0;
    if ( __u.mesh()->element( 0 ).isAHexahedraShape() )
        nHexa = __u.mesh()->elements().size();


    Debug( 8007 ) << "    nLines = " << nLines << "\n";
    Debug( 8007 ) << "nTriangles = " << nTriangles << "\n";
    Debug( 8007 ) << "     nQuad = " << nQuad << "\n";
    Debug( 8007 ) << "    nTetra = " << nTetra << "\n";
    Debug( 8007 ) << "     nHexa = " << nHexa << "\n";

    // nb-scalar-points nb-vector-points nb-tensor-points
    out << "0 0 0 " << "\n";

    //nb-scalar-lines nb-vector-lines nb-tensor-lines
    if ( __u.is_scalar )
        out << nLines << " 0 0 " << "\n";
    else
        out << "0 " << nLines << " 0 " << "\n";

    //nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
    if ( __u.is_scalar )
        out << nTriangles  << " 0 0" << "\n";
    else
        out << "0 " << nTriangles << " 0" << "\n";

    //nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
    if ( __u.is_scalar )
        out << nQuad << " 0 0 " << "\n";
    else
        out << "0 " << nQuad << " 0 " << "\n";

    // nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
    if ( __u.is_scalar )
        out << nTetra << " 0 0"  << "\n";
    else
        out << "0 " << nTetra << " 0"  << "\n";

    // nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
    if ( __u.is_scalar )
        out << nHexa << " 0 0"  << "\n";
    else
        out << "0 " << nHexa << " 0"  << "\n";

    //nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
    out << "0 0 0 " << "\n";

    //nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
    out << "0 0 0 " << "\n";

    //nb-scalar-lines2 nb-vector-lines2 nb-tensor-lines2
    out << "0 0 0 " << "\n";

    //nb-scalar-triangles2 nb-vector-triangles2 nb-tensor-triangles2
    out << "0 0 0 " << "\n";

    //nb-scalar-quadrangles2 nb-vector-quadrangles2 nb-tensor-quadrangles2
    out << "0 0 0 " << "\n";

    //nb-scalar-tetrahedra2 nb-vector-tetrahedra2 nb-tensor-tetrahedra2
    out << "0 0 0 " << "\n";

    //nb-scalar-hexahedra2 nb-vector-hexahedra2 nb-tensor-hexahedra2
    out << "0 0 0 " << "\n";

    //nb-scalar-prisms2 nb-vector-prisms2 nb-tensor-prisms2
    out << "0 0 0 " << "\n";

    //nb-scalar-pyramids2 nb-vector-pyramids2 nb-tensor-pyramids2
    out << "0 0 0 " << "\n";

    //nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars
    out << "0 0 0 " << "\n";



    out.precision( 5 );
    out.setf( std::ios::scientific );
    std::for_each( __ts->beginStep(), __ts->endStep(),  out << lambda::bind( &TimesetType::step_type::time,
                                                                             *lambda::_1 ) << " " );
    out << "\n";


    //loop on all elements
    typedef typename FunctionType::mesh_type mesh_type;
    typedef typename mesh_type::element_const_iterator element_const_iterator;

    element_const_iterator __elit = __u.mesh()->beginElement();
    element_const_iterator __elen = __u.mesh()->endElement();
    for ( ; __elit != __elen; ++__elit )
        {
            uint16_type n_vert = __elit->nVertices();

            matrix_node_type const& __G = __elit->G();
            //std::cout << __G << "\n";
            typename matrix_node_type::const_iterator1 i1=__G.begin1();
            for ( uint16_type __i = 0;
                  i1!=__G.end1() && __i < n_vert; ++i1 )
                {
                    std::for_each( i1.begin(),i1.end(),
                                   out << lambda::_1 << " " );
                    out << "\n";
                }
            if ( __G.size1() == 1 )
                {
                    std::for_each( __G.begin1().begin(),
                                   __G.begin1().begin()+n_vert,
                                   ( out << 0 * lambda::_1 << " " ) );
                    out << "\n";
                    std::for_each( __G.begin1().begin(),
                                   __G.begin1().begin()+n_vert,
                                   out << 0 * lambda::_1 << " " );
                    out << "\n";
                }
            if ( __G.size1() == 2 )
                {
                    std::for_each( __G.begin1().begin(),
                                   __G.begin1().begin()+n_vert,
                                   out << 0 * lambda::_1 << " " );
                    out << "\n";
                }

            typedef typename mesh_type::element_type element_type;
            typedef typename element_type::point_const_iterator point_const_iterator;
            typedef typename FunctionType::id_type eval_fun_type;

            matrix_node_type M( __u.nComponents, n_vert );

            point_const_iterator __pit = __elit->beginPoint();
            point_const_iterator __pen = __elit->endPoint();
            for (uint16_type __i = 0 ;__pit != __pen && __i < n_vert; ++__pit, ++__i )
                {
                    eval_fun_type eval = __u( (*__pit)->node() );
                    for ( int __c = 0; __c < __u.nComponents; ++__c )
                        out << eval(__c, 0, 0) << " ";
                    for( int __c = __u.nComponents; __c < 3; ++__c )
                        out << " 0 ";
                    out << "\n";
                }
        }
    out << "$EndView\n";

}
template<typename MeshType>
void
gmsh_save_ascii( ExporterGmsh<MeshType> const& egmsh )
{
    Debug( 8007 ) << "[gmsh_save_ascii] saving in gmsh ascii file format\n";

    namespace lambda = boost::lambda;

    typedef ExporterGmsh<MeshType> exporter_type;
    typedef typename exporter_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::value_type value_type;
    typedef typename exporter_type::timeset_type timeset_type;
    typedef typename exporter_type::timeset_ptrtype timeset_ptrtype;
    typedef typename exporter_type::timeset_iterator timeset_iterator;
    typedef typename exporter_type::timeset_const_iterator timeset_const_iterator;

    typedef typename matrix_node<value_type>::type matrix_node_type;

    timeset_iterator __ts_it = egmsh.beginTimeSet();
    timeset_iterator __ts_en = egmsh.endTimeSet();
    while ( __ts_it != __ts_en )
        {
            timeset_ptrtype __ts = *__ts_it;
            std::string filename =  egmsh.prefix() + ".msh_data";

            std::ofstream out(filename.c_str());
            if (out.fail())
                {
                    Debug( 8007 ) << "cannot open " << filename.c_str() << "\n";
                    exit(0);
                }

            Debug( 8007 ) << "[ExporterGmsh] saving model " << __ts->name() << " at time step " << __ts->index() << " in " << filename << "\n";

            out << "$PostFormat\n"
                << "1.4 0 " << sizeof(double) << "\n"
                << "$EndPostFormat\n";

            //
            // write time step values
            //
            typename timeset_type::step_const_iterator __it = __ts->beginStep();
            typename timeset_type::step_const_iterator __end = __ts->endStep();

            typename timeset_type::step_type::mesh_ptrtype __m = ( *__it )->mesh();

            //
            // write data
            //
            __it = __ts->beginStep();
            while( __it != __end )
                {
                    typename timeset_type::step_ptrtype __step = *__it;;

                    __m = __step->mesh();
                    typedef typename timeset_type::step_type::mesh_type mesh_type;
                    typedef typename mesh_type::element_const_iterator element_const_iterator;

                    typename timeset_type::step_type::nodal_scalar_const_iterator __var = __step->beginNodalScalar();
                    typename timeset_type::step_type::nodal_scalar_const_iterator __varen = __step->endNodalScalar();
                    while( __var != __varen )
                        {
                            typename timeset_type::step_type::nodal_scalar_type const& __u = __var->second;
                            gmsh_save_value( out, __ts, __u );
                            ++__var;
                        }
                    typename timeset_type::step_type::element_scalar_const_iterator __varelit = __step->beginElementScalar();
                    typename timeset_type::step_type::element_scalar_const_iterator __varelen = __step->endElementScalar();
                    while( __varelit != __varelen )
                        {
                            typename timeset_type::step_type::element_scalar_type const& __u = __varelit->second;
                            gmsh_save_value( out, __ts, __u );
                            ++__varelit;
                        }

                    typename timeset_type::step_type::nodal_vector_const_iterator __vec = __step->beginNodalVector();
                    typename timeset_type::step_type::nodal_vector_const_iterator __vecen = __step->endNodalVector();
                    while( __vec != __vecen )
                        {
                            typename timeset_type::step_type::nodal_vector_type const& __u = __vec->second;
                            gmsh_save_value( out, __ts, __u );
                            ++__vec;
                        }
                    typename timeset_type::step_type::element_vector_const_iterator __vecelem = __step->beginElementVector();
                    typename timeset_type::step_type::element_vector_const_iterator __vecelemen = __step->endElementVector();
                    while( __vec != __vecen )
                        {
                            typename timeset_type::step_type::element_vector_type const& __u = __vecelem->second;
                            gmsh_save_value( out, __ts, __u );
                            ++__vec;
                        }
                    ++__it;
                }

            ++__ts_it;
        }
}
template<typename MeshType>
void
gmsh_save_binary( ExporterGmsh<MeshType> const& egmsh )
{
    Debug( 8007 ) << "[gmsh_save_ascii] saving in gmsh binary file format\n";

    namespace lambda = boost::lambda;

    typedef ExporterGmsh<MeshType> exporter_type;
    typedef typename exporter_type::mesh_type mesh_type;
    typedef typename mesh_type::value_type value_type;
    typedef typename exporter_type::timeset_type timeset_type;
    typedef typename exporter_type::timeset_ptrtype timeset_ptrtype;
    typedef typename exporter_type::timeset_iterator timeset_iterator;
    typedef typename exporter_type::timeset_const_iterator timeset_const_iterator;

    typedef typename matrix_node<value_type>::type matrix_node_type;

    timeset_iterator __ts_it = egmsh.beginTimeSet();
    timeset_iterator __ts_en = egmsh.endTimeSet();
    while ( __ts_it != __ts_en )
        {
            timeset_ptrtype __ts = *__ts_it;

            std::string filename =  egmsh.prefix() + ".msh_data";

            std::ofstream out(filename.c_str());
            if (out.fail())
                {
                    Debug( 8007 ) << "cannot open " << filename.c_str() << "\n";
                    exit(0);
                }

            Debug( 8007 ) << "[ExporterGmsh] saving model " << __ts->name() << " at time step " << __ts->index() << " in " << filename << "\n";

            out << "$PostFormat\n"
                << "1.4 0 " << sizeof(double) << "\n"
                << "$EndPostFormat\n"
                << "$View\n";

            // view-name nb-time-steps
            out << __ts->name() << " " << __ts->numberOfSteps()  << "\n";

            //
            // write time step values
            //
            typename timeset_type::step_const_iterator __it = __ts->beginStep();
            typename timeset_type::step_const_iterator __end = __ts->endStep();

            typename timeset_type::step_type::mesh_ptrtype __m = ( *__it )->mesh();

            size_type nLines = 0;
            if ( __m->element( 0 ).isALineShape() )
                nLines = __m->elements().size();

            size_type nTriangles = 0;
            if ( __m->element( 0 ).isATriangleShape() )
                nTriangles = __m->elements().size();

            size_type nQuad = 0;
            if ( __m->element( 0 ).isAQuadrangleShape() )
                nQuad = __m->elements().size();

            size_type nTetra = 0;
            if ( __m->element( 0 ).isATetrahedraShape() )
                nTetra = __m->elements().size();

            size_type nHexa = 0;
            if ( __m->element( 0 ).isAHexahedraShape() )
                nHexa = __m->elements().size();


            Debug( 8007 ) << "    nLines = " << nLines << "\n";
            Debug( 8007 ) << "nTriangles = " << nTriangles << "\n";
            Debug( 8007 ) << "     nQuad = " << nQuad << "\n";
            Debug( 8007 ) << "     nHexa = " << nHexa << "\n";

            // nb-scalar-points nb-vector-points nb-tensor-points
            out << "0 0 0 " << "\n";

            // nb-scalar-lines nb-vector-lines nb-tensor-lines
            out << nLines << " 0 0 " << "\n";

            // nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
            out << nTriangles  << " 0 0" << "\n";

            // nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
            out << nQuad << " 0 0 " << "\n";

            // nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
            out << nTetra << " 0 0"  << "\n";

            // nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
            out << nHexa << " 0 0 " << "\n";

            // nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
            out << "0 0 0 " << "\n";

            // nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
            out << "0 0 0 " << "\n";

            // nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars
            out << "0 0 0" << "\n";


            out.precision( 5 );
            out.setf( std::ios::scientific );
            std::for_each( __ts->beginStep(), __ts->endStep(),  out << lambda::bind( &timeset_type::step_type::time,
                                                                                     *lambda::_1 ) << " " );
            out << "\n";

            //
            // write data
            //
            __it = __ts->beginStep();
            while( __it != __end )
                {
                    typename timeset_type::step_ptrtype __step = *__it;;

                    //typename timeset_type::step_type::mesh_ptrtype __m = __step->mesh();
                    __m = __step->mesh();
                    typedef typename timeset_type::step_type::mesh_type mesh_type;
                    typedef typename mesh_type::element_const_iterator element_const_iterator;

                    // < scalar-point-value > ...

                    // < vector-point-value > ...
                    // < tensor-point-value > ...
                    // < scalar-line-value > ...
                    // < vector-line-value > ...
                    // < tensor-line-value > ...
                    // < scalar-triangle-value > ...
                    typename timeset_type::step_type::nodal_scalar_const_iterator __var = __step->beginNodalScalar();
                    typename timeset_type::step_type::nodal_scalar_const_iterator __varen = __step->endNodalScalar();
                    while( __var != __varen )
                        {
                            //loop on all elements
                            element_const_iterator __elit = __m->beginElement();
                            element_const_iterator __elen = __m->endElement();
                            for ( ; __elit != __elen; ++__elit )
                                {
                                    uint16_type n_vert = __elit->nVertices();

                                    matrix_node_type const& __G = __elit->G();
                                    //std::cout << __G << "\n";
                                    typename matrix_node_type::const_iterator1 i1=__G.begin1();
                                    for ( uint16_type __i = 0;
                                          i1!=__G.end1() && __i < n_vert; ++i1 )
                                        {
                                            std::for_each( i1.begin(),i1.end(),
                                                           out << lambda::_1 << " " );
                                            out << "\n";
                                        }
                                    if ( __G.size1() == 1 )
                                        {
                                            std::for_each( __G.begin1().begin(),
                                                           __G.begin1().begin()+n_vert,
                                                           ( out << 0 * lambda::_1 << " " ) );
                                            out << "\n";
                                            std::for_each( __G.begin1().begin(),
                                                           __G.begin1().begin()+n_vert,
                                                           out << 0 * lambda::_1 << " " );
                                            out << "\n";
                                        }
                                    if ( __G.size1() == 2 )
                                        {
                                            std::for_each( __G.begin1().begin(),
                                                           __G.begin1().begin()+n_vert,
                                                           out << 0 * lambda::_1 << " " );
                                            out << "\n";
                                        }

#warning TOBEFIXED
#if 0
                                    typename timeset_type::step_type::nodal_scalar_type const& __u = __var->second;


                                    typedef typename mesh_type::point_type point_type;
                                    typedef typename mesh_type::element_point_const_iterator point_const_iterator;

                                    point_const_iterator __pit = __elit->beginPoint();
                                    point_const_iterator __pen = __elit->endPoint();
                                    for (uint16_type __i = 0 ;__pit != __pen && __i < n_vert; ++__pit )
                                        {
                                            out << __u[ ( *__pit )->id() ] << "\n";
                                        }
#endif
                                }
                            ++__var;
                        }

                    typename timeset_type::step_type::nodal_vector_const_iterator __vec = __step->beginNodalVector();
                    typename timeset_type::step_type::nodal_vector_const_iterator __vecen = __step->endNodalVector();
                    while( __vec != __vecen )
                        {
                            ++__vec;
                        }
                    // < vector-triangle-value > ...
                    // < tensor-triangle-value > ...
                    // < scalar-quadrangle-value > ...
                    // < vector-quadrangle-value > ...
                    // < tensor-quadrangle-value > ...
                    // < scalar-tetrahedron-value > ..
                    // < vector-tetrahedron-value > ...
                    // < tensor-tetrahedron-value > ...
                    // < scalar-hexahedron-value > ...
                    // < vector-hexahedron-value > ...
                    // < tensor-hexahedron-value > ...
                    // < scalar-prism-value > ...
                    // < vector-prism-value > ...
                    // < tensor-prism-value > ...
                    // < scalar-pyramid-value > ...
                    // < vector-pyramid-value > ...
                    // < tensor-pyramid-value > ...
                    // < text2d > ... < text2d-chars > ...
                    // < text3d > ... < text3d-chars > ...

                    out << "$EndView\n";

                    ++__it;
                }
            out.close();

            ++__ts_it;
        }
}
} // detail



template<typename MeshType>
ExporterGmsh<MeshType>::ExporterGmsh( std::string const& __p, int freq )
    :
    super( "gmsh", __p, freq )
{

}
template<typename MeshType>
ExporterGmsh<MeshType>::ExporterGmsh( po::variables_map const& vm, std::string const& exp_prefix )
    :
    super( vm, exp_prefix )
{
}
template<typename MeshType>
ExporterGmsh<MeshType>::ExporterGmsh( ExporterGmsh const & __ex )
    :
    super( __ex )
{}
template<typename MeshType>
ExporterGmsh<MeshType>::~ExporterGmsh()
{}

template<typename MeshType>
void
ExporterGmsh<MeshType>::save() const
{
    static int freq = 0;

    Debug( 8007 ) << "[ExporterGmsh] checking if frequency is ok\n";

    if ( freq++ % this->freq()  )
        return;

    Debug( 8007 ) << "[ExporterGmsh] frequency is ok\n";

    Debug( 8007 ) << "[ExporterGmsh] save()...\n";

    if( this->fileType() == ASCII )
        detail::gmsh_save_ascii( *this );
    else if( this->fileType() == BINARY )
        detail::gmsh_save_binary( *this );

    Debug( 8007 ) << "[ExporterGmsh] saving done\n";
}

template<typename MeshType>
void
ExporterGmsh<MeshType>::visit( mesh_type* )
{
}

//
// explicit instances
//
template class ExporterGmsh<Mesh<GeoEntity<Simplex<1,1> > > >;
template class ExporterGmsh<Mesh<GeoEntity<Simplex<1,1,2> > > >;
template class ExporterGmsh<Mesh<GeoEntity<Simplex<2,1> > > >;
template class ExporterGmsh<Mesh<GeoEntity<Simplex<2,2> > > >;
template class ExporterGmsh<Mesh<GeoEntity<Simplex<2,1,3> > > >;
template class ExporterGmsh<Mesh<GeoEntity<Simplex<3,1> > > >;
template class ExporterGmsh<Mesh<GeoEntity<SimplexProduct<1,1> > > >;
template class ExporterGmsh<Mesh<GeoEntity<SimplexProduct<2,1> > > >;
template class ExporterGmsh<Mesh<GeoEntity<SimplexProduct<3,1> > > >;

}
