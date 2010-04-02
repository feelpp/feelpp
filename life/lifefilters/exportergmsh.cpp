/* -*- mode: c++ -*-

   This file is part of the Life library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2007-07-21

   Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file exportergmsh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-21
*/
#ifndef __EXPORTERGMSH_CPP
#define __EXPORTERGMSH_CPP 1

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
            typedef typename FunctionType::functionspace_type::mesh_type mesh_type;
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

            timeset_const_iterator __ts_it = egmsh.beginTimeSet();
            timeset_const_iterator __ts_en = egmsh.endTimeSet();
            while ( __ts_it != __ts_en )
                {
                    timeset_ptrtype __ts = *__ts_it;
                    std::string filename =  egmsh.prefix() + ".msh";

                    std::ofstream out(filename.c_str());
                    if (out.fail())
                        {
                            Debug( 8007 ) << "cannot open " << filename.c_str() << "\n";
                            exit(0);
                        }

                    Debug( 8007 ) << "[ExporterGmsh] saving model " << __ts->name() << " at time step " << __ts->index() << " in " << filename << "\n";

                    out << "$MeshFormat\n"
                        << "2 0 0 " << sizeof(double) << "\n"
                        << "$EndMeshFormat\n";

                    out << "$PhysicalNames\n";
                    // write Physical names here
                    out << "$EndPhysicalNames\n";

                    out << "$Nodes\n";
                    // Save mesh nodes here
                    out << "$EndNodes\n";

                    out << "$Elements\n";
                    // Save mesh elements here
                    out << "$EndElements\n";

#if 0
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
#endif
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

            timeset_const_iterator __ts_it = egmsh.beginTimeSet();
            timeset_const_iterator __ts_en = egmsh.endTimeSet();
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
            gmsh_save_ascii();
        //detail::gmsh_save_ascii( *this );
        else if( this->fileType() == BINARY )
            detail::gmsh_save_binary( *this );

        Debug( 8007 ) << "[ExporterGmsh] saving done\n";
    }

    template<typename MeshType>
    void
    ExporterGmsh<MeshType>::visit( mesh_type* )
    {
    }

    template<typename MeshType>
    void
    ExporterGmsh<MeshType>::gmsh_save_ascii() const
    {
        Debug( 8007 ) << "[gmsh_save_ascii] saving in gmsh ascii file format\n";

        timeset_const_iterator __ts_it = this->beginTimeSet();
        timeset_const_iterator __ts_en = this->endTimeSet();
        while ( __ts_it != __ts_en )
            {
                timeset_ptrtype __ts = *__ts_it;

                std::ostringstream __fname;

                __fname << __ts->name()  //<< this->prefix() //this->path()
                        << "-" << Application::nProcess() << "_" << Application::processId()
                        << ".msh";

                /*std::string filename =  this->prefix()
                  + __ts->name()
                  + "-" + Application::nProcess() + "_" + Application::processId()
                  + ".msh";*/
                std::ofstream out;

                typename timeset_type::step_const_iterator __stepIt = __ts->beginStep();
                typename timeset_type::step_const_iterator __stepIt_end = __ts->endStep();
                __stepIt = boost::prior( __stepIt_end );
                if( __stepIt != __stepIt_end )
                    {
                        step_ptrtype __step = *__stepIt;
                        if ( __step->isInMemory() )
                            {

                                if (__step->index()==1)
                                    {
                                        out.open(__fname.str().c_str(), ios::out);
                                    }
                                else
                                    {
                                        out.open(__fname.str().c_str(), ios::out | ios::app);
                                    }

                                if (out.fail())
                                    {
                                        Debug( 8007 ) << "cannot open " << __fname.str().c_str() << "\n";
                                        exit(0);
                                    }

                                Debug( 8007 ) << "[ExporterGmsh] saving model "
                                              << __ts->name() << " at time step "
                                              << __ts->index() << " in "
                                              << __fname.str() << "\n";


                                gmsh_save_Format( out );

                                // out << "$PhysicalNames\n";
                                // write Physical names here
                                //out << "$EndPhysicalNames\n";

                                gmsh_save_Nodes( out,__step);

                                gmsh_save_Elements( out, __step);

                                //gmsh_save_NodeData( out, __step);

                                gmsh_save_ElementNodeData( out, __step);
                            }
                    }
                ++__ts_it;
            }
    }

    template<typename MeshType>
    void
    ExporterGmsh<MeshType>::gmsh_save_file( std::ostream& out ) const
    {
    }

    template<typename MeshType>
    void
    ExporterGmsh<MeshType>::gmsh_save_Format( std::ostream& out ) const
    {
        out << "$MeshFormat\n"
            << "2.1 0 " << sizeof(double) << "\n"
            << "$EndMeshFormat\n";

    }


    template<typename MeshType>
    void
    ExporterGmsh<MeshType>::gmsh_save_Nodes( std::ostream& out,
                                             step_ptrtype __step ) const
    {
        out << "$Nodes\n";

        mesh_ptrtype mesh = __step->mesh();

        out << mesh->numPoints() << "\n";//number points

        point_const_iterator pt_it = mesh->beginPoint();
        point_const_iterator pt_en = mesh->endPoint();
        for ( ; pt_it!=pt_en ; ++pt_it )
            {
                out << pt_it->id()+1
                    <<" " << pt_it->node()[0];
                if ( mesh_type::nRealDim >= 2 )
                    out << " " << pt_it->node()[1];
                else
                    out << " 0";
                if ( mesh_type::nRealDim >= 3 )
                    out << " " << pt_it->node()[2];
                else
                    out << " 0";
                out << "\n";
            }

        out << "$EndNodes\n";

    }

    template<typename MeshType>
    void
    ExporterGmsh<MeshType>::gmsh_save_Elements( std::ostream& out,
                                                step_ptrtype __step ) const
    {
        out << "$Elements\n";

        mesh_ptrtype mesh = __step->mesh();

        out << mesh->numElements() << "\n";//number element

        typename mesh_type::element_const_iterator elt_it = mesh->beginElement();
        typename mesh_type::element_const_iterator elt_en = mesh->endElement();
        uint16_type nLocGeoPt;
        for ( ; elt_it != elt_en; ++elt_it )
            {
                out << elt_it->id()+1<<" ";
                nLocGeoPt = elt_it->nPoints();
                if (elt_it->isATriangleShape())
                    {
                        switch (mesh_type::nOrder) {
                        case 1 : //if (mesh_type::nOrder==1)
                            {
                                out << 2 ;//type triangle order 1
                                out<<" 2 99 2";
                                for (uint16_type p=0;p<nLocGeoPt;++p)
                                    out << " " << elt_it->point( p ).id()+1;
                                break;
                            }
                        case 2 : //else if (mesh_type::nOrder==2)
                            {
                                out << 9 ;//type triangle order 2
                                out<<" 2 99 2";
                                out << " " << elt_it->point( 0 ).id()+1;
                                out << " " << elt_it->point( 1 ).id()+1;
                                out << " " << elt_it->point( 2 ).id()+1;
                                out << " " << elt_it->point( 5 ).id()+1;
                                out << " " << elt_it->point( 3 ).id()+1;
                                out << " " << elt_it->point( 4 ).id()+1;
                                break;
                            }
                        case 3 : //else if (mesh_type::nOrder==3)
                            {
                                out << 21 ;//type triangle order 3
                                out<<" 2 99 2";
                                out << " " << elt_it->point( 0 ).id()+1;
                                out << " " << elt_it->point( 1 ).id()+1;
                                out << " " << elt_it->point( 2 ).id()+1;
                                out << " " << elt_it->point( 7 ).id()+1;
                                out << " " << elt_it->point( 8 ).id()+1;
                                out << " " << elt_it->point( 3 ).id()+1;
                                out << " " << elt_it->point( 4 ).id()+1;
                                out << " " << elt_it->point( 5 ).id()+1;
                                out << " " << elt_it->point( 6 ).id()+1;
                                out << " " << elt_it->point( 9 ).id()+1;
                                break;
                            }
                        case 4 : //else if (mesh_type::nOrder>3)
                            {
                                break;
#warning TOFILL
                            }
                        }
                    }
                else if ( elt_it->isAQuadrangleShape())
                    {
                        switch (mesh_type::nOrder) {
                        case 1 :
                            {
                                out << 4 ;//type quadrangle order 1
                                out<<" 2 99 2";
                                for (uint16_type p=0;p<nLocGeoPt;++p)
                                    out << " " << elt_it->point( p ).id()+1;
                                break;
                            }
                        case 2 :
                            {
                                out << 10 ;//type quadrangle order 2
                                out<<" 2 99 2";
                                out << " " << elt_it->point( 0 ).id()+1;
                                out << " " << elt_it->point( 1 ).id()+1;
                                out << " " << elt_it->point( 2 ).id()+1;
                                out << " " << elt_it->point( 3 ).id()+1;
                                out << " " << elt_it->point( 7 ).id()+1;
                                out << " " << elt_it->point( 4 ).id()+1;
                                out << " " << elt_it->point( 5 ).id()+1;
                                out << " " << elt_it->point( 6 ).id()+1;
                                break;
                            }
                        case 3 :
                            {
#if 0
                                out << 21 ;//???type quadrangle order 3
                                out<<" 2 99 2";
                                out << " " << elt_it->point( 0 ).id()+1;
                                out << " " << elt_it->point( 1 ).id()+1;
                                out << " " << elt_it->point( 2 ).id()+1;
                                out << " " << elt_it->point( 7 ).id()+1;
                                out << " " << elt_it->point( 8 ).id()+1;
                                out << " " << elt_it->point( 3 ).id()+1;
                                out << " " << elt_it->point( 4 ).id()+1;
                                out << " " << elt_it->point( 5 ).id()+1;
                                out << " " << elt_it->point( 6 ).id()+1;
                                out << " " << elt_it->point( 9 ).id()+1;
                                break;
#endif
                            }
                        case 4 : //else if (mesh_type::nOrder>3)
                            {
                                break;
#warning TOFILL
                            }
                        }

                    }
                out<<"\n";
            }
        out << "$EndElements\n";
    }



    template<typename MeshType>
    void
    ExporterGmsh<MeshType>::gmsh_save_NodeData( std::ostream& out, step_ptrtype __step ) const
    {
#if 0
        //!!!Not functionnal for curve element!!!

        typedef typename step_type::nodal_scalar_type nodal_scalar_type;
        typedef typename step_type::nodal_scalar_const_iterator nodal_scalar_const_iterator;

        //on parcourt le temp
        step_const_iterator __it = timeset->beginStep();
        step_const_iterator __end = timeset->endStep();
        for ( ; __it != __end ; ++__it)
            {
                nodal_scalar_const_iterator __var = (*__it)->beginNodalScalar();
                nodal_scalar_const_iterator __varen = (*__it)->endNodalScalar();

                out << "$NodeData\n";

                nodal_scalar_type const& __u = __var->second;
                //mesh_ptrtype mesh =__u.mesh();
                mesh_ptrtype mesh = (*__it)->mesh();

                out << "1\n";//number of string tag
                out << "a scalar node\n";
                out << "1\n";//number of real tag
                out << "0.0\n";
                out << "3\n";//number of integer tags:
                out << "0\n";//the time step (0; time steps always start at 0)
                out << "1\n";//n-component (1 is scalar) field
                out << mesh->numPoints() << "\n";//number associated nodal values

                point_const_iterator pt_it = mesh->beginPoint();
                point_const_iterator pt_en = mesh->endPoint();
                for ( ; pt_it!=pt_en ; ++pt_it )
                    {
                        out << pt_it->id()+1
                            <<" ";
                        out <<__u(pt_it->id());
                        //__u(pt_it->node());
                        out << "\n";
                    }
                out << "$EndNodeData\n";
            }
#endif
    }


    template<typename MeshType>
    void
    ExporterGmsh<MeshType>::gmsh_save_ElementNodeData( std::ostream& out,
                                                       step_ptrtype __step) const
    {

        typedef typename mesh_type::element_const_iterator element_mesh_const_iterator;

        typedef typename step_type::nodal_scalar_type nodal_scalar_type;
        typedef typename step_type::nodal_scalar_const_iterator nodal_scalar_const_iterator;
        typedef typename step_type::nodal_vector_type nodal_vectorial_type;
        typedef typename step_type::nodal_vector_const_iterator nodal_vectorial_const_iterator;


        mesh_ptrtype mesh = __step->mesh();

        uint16_type nLocalDof;
        size_type globaldof;

        nodal_scalar_const_iterator __varScal = __step->beginNodalScalar();
        nodal_scalar_const_iterator __varScal_end = __step->endNodalScalar();
        for ( ;__varScal!=__varScal_end ; ++__varScal)
            {
                out << "\n$ElementNodeData\n";

                nodal_scalar_type const& __u = __varScal->second;

                //uint __nbCompFieldGMSH;
                //    if (nodal_scalar_type::functionspace_type::is_scalar)         { __nbCompFieldGMSH=1; }
                //else if (nodal_scalar_type::functionspace_type::is_vectorial) { __nbCompFieldGMSH=3; }
                //else if (nodal_scalar_type::functionspace_type::is_tensor2)   { __nbCompFieldGMSH=9; }

                out << "1\n";//number of string tag
                out << "\"" << __varScal->first <<"\"\n";//a scalar node\n";
                out << "1\n";//number of real tag
                out << __step->time() << "\n";//"0.0\n";//the time value (0.0)
                out << "3\n";//number of integer tags:
                out << __step->index()-1 << "\n";//"0\n";//the time step (0; time steps always start at 0)
                out << "1\n";//n-component (1 is scalar) field
                out << mesh->numElements() << "\n";//number associated nodal values

                element_mesh_const_iterator elt_it = mesh->beginElement();
                element_mesh_const_iterator elt_en = mesh->endElement();
                if ( !__u.areGlobalValuesUpdated() )
                    __u.updateGlobalValues();

                for ( ; elt_it!=elt_en ; ++elt_it )
                    {
                        out << elt_it->id()+1;
                        //nLocGeoPt = elt_it->nPoints();
                        nLocalDof = mesh->numLocalVertices();
                        //nodal_scalar_type::functionspace_type::basis_type::nLocalDof;

                        out << " " << nLocalDof;
                        for ( uint16_type l = 0; l < nLocalDof; ++l )
                            {
                                out << " ";
                                globaldof = boost::get<0>(__u.functionSpace()->dof()->localToGlobal(elt_it->id(), l, 0 ));//l,c
                                out << __u( globaldof);
                            }
                        out << "\n";
                    }
                out << "$ElementEndNodeData\n";
            }

        nodal_vectorial_const_iterator __varVec = __step->beginNodalVector();
        nodal_vectorial_const_iterator __varVec_end = __step->endNodalVector();
        for ( ;__varVec!=__varVec_end ; ++__varVec)
            {
                out << "\n$ElementNodeData\n";

                nodal_vectorial_type const& __uVec = __varVec->second;

                uint16_type nComponents = __uVec.nComponents;

                out << "1\n";//number of string tag
                out << "\"" << __varVec->first <<"\"\n";//"a vectorial field\n";
                out << "1\n";//number of real tag
                out << __step->time() << "\n";//"0.0\n";//the time value (0.0)
                out << "3\n";//number of integer tags:
                out << __step->index() << "\n";//"0\n";//the time step (0; time steps always start at 0)
                out << "3\n";//n-component (3 is vectorial) field
                out << mesh->numElements() << "\n";//number associated nodal values

                element_mesh_const_iterator elt_it = mesh->beginElement();
                element_mesh_const_iterator elt_en = mesh->endElement();
                for ( ; elt_it!=elt_en ; ++elt_it )
                    {
                        out << elt_it->id()+1;
                        nLocalDof = mesh->numLocalVertices();
                        out << " " << nLocalDof;
                        for ( uint16_type l = 0; l < nLocalDof; ++l )
                            {
                                for( uint16_type c = 0; c < 3; ++c )
                                    {
                                        out << " ";
                                        if (c < nComponents)
                                            {
                                                globaldof = boost::get<0>(__uVec.functionSpace()->dof()->localToGlobal(elt_it->id(),
                                                                                                                       l,
                                                                                                                       c ));
                                                out << __uVec( globaldof);
                                            }
                                        else out << "0.0";
                                    }
                            }
                        out << "\n";
                    }
                out << "$ElementEndNodeData\n";
            }
    }


#if defined( LIFE_INSTANTIATION_MODE )
    //
    // explicit instances
    //
    template class ExporterGmsh<Mesh<Simplex<1,1> > >;
    template class ExporterGmsh<Mesh<Simplex<1,1,2> > >;
    template class ExporterGmsh<Mesh<Simplex<2,1> > >;
    template class ExporterGmsh<Mesh<Simplex<2,2> > >;
    template class ExporterGmsh<Mesh<Simplex<2,1,3> > >;
    template class ExporterGmsh<Mesh<Simplex<3,1> > >;
    template class ExporterGmsh<Mesh<SimplexProduct<1,1> > >;
    template class ExporterGmsh<Mesh<SimplexProduct<2,1> > >;
    template class ExporterGmsh<Mesh<SimplexProduct<3,1> > >;

    template class ExporterGmsh<Mesh<Simplex<2,3> > >;
    template class ExporterGmsh<Mesh<SimplexProduct<2,2> > >;
    template class ExporterGmsh<Mesh<SimplexProduct<2,3> > >;


#endif // LIFE_INSTANTIATION_MODE
}
#endif // __EXPORTERGMSH_CPP
