
#ifndef __geotool_H
#define __geotool_H 1


#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <map>
#include <life/lifealg/glas.hpp>
//#include <boost/parameter/keyword.hpp>
//#include <boost/parameter/preprocessor.hpp>
#include <life/lifecore/parameter.hpp>

#include <boost/preprocessor/tuple/elem.hpp>



/*_________________________________________________*/
/*_________________________________________________*/
/*_________________________________________________*/
//Bug si 0 argum√ent(a corriger)

# define GEOTOOL_SHAPE                                                  \
    ( 7, ( ( Rectangle     , 2, 1, 0, "rectangle"    , 2, RECTANGLE ),  \
           ( Quadrangle    , 2, 1, 0, "quadrangle"   , 4, QUADRANGLE ), \
           ( Circle        , 2, 1, 0, "circle"       , 2, CIRCLE    ),  \
           ( PartialDisque , 2, 1, 0, "partialdisque", 3, PARTIALDISQUE), \
           ( Special_1a    , 2, 2, 0, "special_1a"   , 1, SPECIAL_1A ), \
           ( Special_1b    , 2, 1, 0, "special_1b"   , 1, SPECIAL_1B ), \
           ( Hexaedre      , 3, 6, 1, "hexaedre"     , 8, HEXAEDRE  )   \
           )                                                            \
      )                                                                 \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_LINE_RECTANGLE          \
    ( 4, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_RECTANGLE       \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_LINE_QUADRANGLE         \
    ( 4, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_QUADRANGLE      \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/


/*_________________________________________________*/

# define GEOTOOL_MARKER_LINE_CIRCLE             \
    ( 1, ( ( 1, 2, ( 1,2 ) ) )                  \
      )                                         \
    /**/

# define GEOTOOL_MARKER_SURFACE_CIRCLE          \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_LINE_PARTIALDISQUE      \
    ( 2, ( ( 1, 4, ( 1,2,3,4 ) ),               \
           ( 2, 1, (    5    ) ) )              \
      )                                         \
    /**/

# define GEOTOOL_MARKER_SURFACE_PARTIALDISQUE   \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/


/*_________________________________________________*/

# define GEOTOOL_MARKER_LINE_SPECIAL_1A            \
    ( 4, ( ( 1, 2, ( 1,5 ) ),                      \
           ( 2, 2, ( 2,6 ) ),                      \
           ( 3, 2, ( 3,7 ) ),                      \
           ( 4, 2, ( 4,8 ) ) )                     \
      )                                            \
    /**/

# define GEOTOOL_MARKER_SURFACE_SPECIAL_1A      \
    ( 1, ( ( 1, 2, ( 1,2 ) ) )                  \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_LINE_SPECIAL_1B          \
    ( 3, ( ( 1, 2, ( 1,2 ) ),                    \
           ( 2, 1, ( 3   ) ),                      \
           ( 3, 1, ( 4   ) )                       \
           )                                     \
      )                                          \
    /**/

# define GEOTOOL_MARKER_SURFACE_SPECIAL_1B      \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_LINE_HEXAEDRE           \
    ( 12, ( (  1, 1, (  1 ) ),                  \
            (  2, 1, (  2 ) ),                  \
            (  3, 1, (  3 ) ),                  \
            (  4, 1, (  4 ) ),                  \
            (  5, 1, (  5 ) ),                  \
            (  6, 1, (  6 ) ),                  \
            (  7, 1, (  7 ) ),                  \
            (  8, 1, (  8 ) ),                  \
            (  9, 1, (  9 ) ),                  \
            ( 10, 1, ( 10 ) ),                  \
            ( 11, 1, ( 11 ) ),                  \
            ( 12, 1, ( 12 ) )                   \
            )                                   \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_HEXAEDRE        \
    ( 6, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ),                     \
           ( 5, 1, ( 5 ) ),                     \
           ( 6, 1, ( 6 ) )                      \
           )                                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_VOLUME_HEXAEDRE         \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_VOLUME_DEFAULT          \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

// Accessors

# define GEOTOOL_SHAPE_NAME_CLASS(i) BOOST_PP_TUPLE_ELEM(7, 0, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_DIM(i) BOOST_PP_TUPLE_ELEM(7, 1, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NBSURFACE(i) BOOST_PP_TUPLE_ELEM(7, 2, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NBVOLUME(i) BOOST_PP_TUPLE_ELEM(7, 3, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NAME_STR(i) BOOST_PP_TUPLE_ELEM(7, 4, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NBPARAM(i) BOOST_PP_TUPLE_ELEM(7, 5, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NAME_MACRO(i) BOOST_PP_TUPLE_ELEM(7, 6, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))

# define GEOTOOL_MARKER_LINE_INDICE(O) BOOST_PP_TUPLE_ELEM(3, 0, O)
# define GEOTOOL_MARKER_LINE_NBMARK(F,i) BOOST_PP_TUPLE_ELEM(3, 1, BOOST_PP_ARRAY_ELEM(i,F))
# define GEOTOOL_MARKER_LINE_ARRAYMARK(O) BOOST_PP_TUPLE_ELEM(3, 2, O)
# define GEOTOOL_MARKER_LINE_MARKVALUE(F,i,j)                           \
    BOOST_PP_TUPLE_ELEM( GEOTOOL_MARKER_LINE_NBMARK(F,i),j,GEOTOOL_MARKER_LINE_ARRAYMARK(BOOST_PP_ARRAY_ELEM(i, F)))

# define GEOTOOL_MARKER_SURFACE_INDICE(O) BOOST_PP_TUPLE_ELEM(3, 0, O)
# define GEOTOOL_MARKER_SURFACE_NBMARK(F,i) BOOST_PP_TUPLE_ELEM(3, 1, BOOST_PP_ARRAY_ELEM(i,F))
# define GEOTOOL_MARKER_SURFACE_ARRAYMARK(O) BOOST_PP_TUPLE_ELEM(3, 2, O)
# define GEOTOOL_MARKER_SURFACE_MARKVALUE(F,i,j)                        \
    BOOST_PP_TUPLE_ELEM( GEOTOOL_MARKER_SURFACE_NBMARK(F,i),j,GEOTOOL_MARKER_SURFACE_ARRAYMARK(BOOST_PP_ARRAY_ELEM(i, F)))

# define GEOTOOL_MARKER_VOLUME_INDICE(O) BOOST_PP_TUPLE_ELEM(3, 0, O)
# define GEOTOOL_MARKER_VOLUME_NBMARK(F,i) BOOST_PP_TUPLE_ELEM(3, 1, BOOST_PP_ARRAY_ELEM(i,F))
# define GEOTOOL_MARKER_VOLUME_ARRAYMARK(O) BOOST_PP_TUPLE_ELEM(3, 2, O)
# define GEOTOOL_MARKER_VOLUME_MARKVALUE(F,i,j)                         \
    BOOST_PP_TUPLE_ELEM( GEOTOOL_MARKER_VOLUME_NBMARK(F,i),j,GEOTOOL_MARKER_VOLUME_ARRAYMARK(BOOST_PP_ARRAY_ELEM(i, F)))

/*_________________________________________________*/
/*_________________________________________________*/
/*_________________________________________________*/


namespace Life {

    namespace GeoTool {

        typedef node<double>::type node_type;

        class GeoGMSHTool;
        typedef boost::shared_ptr< GeoGMSHTool> GeoGMSHTool_ptrtype;

        typedef std::map<uint,uint> map_data_type;
        typedef std::vector<map_data_type> vec_map_data_type;
        typedef boost::shared_ptr<vec_map_data_type> vec_map_data_ptrtype;

        typedef boost::tuple< GeoGMSHTool_ptrtype, vec_map_data_ptrtype ,std::string, std::string > data_geo_type;
        typedef boost::shared_ptr<data_geo_type> data_geo_ptrtype;


        void run(data_geo_ptrtype __dg);

        //A faire avec les BoostPP
        //void runRectangle(data_geo_ptrtype dg);
        //void runCircle(data_geo_ptrtype dg);
        //void runHexaedre(data_geo_ptrtype dg);



#define GEOTOOL_INSTANTIATES_FOR_COMP(r, state)                         \
        BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(2, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 1, state)) \
                            )                                           \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_INSTANTIATES_FOR_INCR(r, state)             \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 0, state)),	\
         BOOST_PP_TUPLE_ELEM(2, 1, state) )                 \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_INSTANTIATES_FOR(r,state)                               \
        void BOOST_PP_CAT(run,GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))) (data_geo_ptrtype dg); \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
        BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE(GEOTOOL_SHAPE),1) ),
                      GEOTOOL_INSTANTIATES_FOR_COMP,
                      GEOTOOL_INSTANTIATES_FOR_INCR,
                      GEOTOOL_INSTANTIATES_FOR );


    }

    namespace GeoTool {

        /*_________________________________________________*
         *_________________________________________________*
         *_________________________________________________*
         * GeoGMSHTool :                                   *
         *_________________________________________________*
         *_________________________________________________*
         *_________________________________________________*/


        class Node {
        public :

            Node()
                :
                _M_node(new node_type())
            {}

            Node(double __x) :
                _M_node(new node_type(1)) { (*_M_node)(0)=__x;}

            Node(double __x, double __y) :
                _M_node(new node_type(2))
            {
                (*_M_node)(0)=__x;
                (*_M_node)(1)=__y;
            }

            Node(double __x, double __y, double __z) :
                _M_node(new node_type(3))
            {
                (*_M_node)(0)=__x;
                (*_M_node)(1)=__y;
                (*_M_node)(2)=__z;
            }

            Node(Node const & m)
                :
                _M_node(m._M_node)
            {}

            Node operator=(Node const & m)
            {
                _M_node.reset(new node_type(*(m._M_node)));
                return *this;
            }

            double operator()(uint n)
            {
                return this->getNode()(n);
            }

            /*boost::shared_ptr<*/node_type
            getNode() const { return *_M_node;}

            boost::shared_ptr<node_type> _M_node;
        };

        /*_________________________________________________*/

        class Loop {
        public :

            Loop(Loop const & L) : _M_loop(L._M_loop) {}

            Loop() {_M_loop.clear();}

            void  operator=(Loop m) { this->_M_loop=m._M_loop; }
            Loop  operator>>(int __n) { _M_loop.push_back(__n); return *this; }

            std::list<int>::const_iterator begin() const { return _M_loop.begin();}
            std::list<int>::const_iterator end() const { return _M_loop.end(); }

            std::list<int> _M_loop;
        };


        /*_________________________________________________*
         *_________________________________________________*
         *_________________________________________________*
         * GeoGMSHTool :                                   *
         *_________________________________________________*
         *_________________________________________________*
         *_________________________________________________*/


        class GeoGMSHTool
        {
        public:

            typedef node<double>::type node_type;
            // list de < nameMesh, meshSize >
            typedef boost::tuple<std::string,double> names_base_type;
            typedef std::list< names_base_type > names_type;
            typedef std::map< std::string, names_type > map_shape_names_type;
            typedef names_type::const_iterator names_const_iterator_type;
            typedef map_shape_names_type::const_iterator map_shape_names_const_iterator_type;

            typedef boost::tuple<std::string,std::string,uint> marker_base_type;
            typedef std::map<std::string,std::list<marker_base_type > > marker_markerName_type;
            typedef std::map< std::string, marker_markerName_type > marker_type_type;
            typedef std::map< std::string, marker_type_type > marker_name_type;
            typedef std::map< std::string, marker_type_type > marker_shape_type;

            typedef marker_markerName_type::const_iterator marker_markerName_const_iterator_type;
            typedef marker_type_type::const_iterator marker_type_const_iterator_type;
            typedef marker_name_type::const_iterator marker_name_const_iterator_type;
            typedef marker_shape_type::const_iterator marker_shape_const_iterator_type;

            typedef std::vector<node_type> parameter_rectangle_type;
            typedef std::map<std::string, parameter_rectangle_type > parameter_name_type;
            typedef std::map<std::string, parameter_name_type > parameter_shape_type;
            typedef parameter_name_type::const_iterator parameter_name_const_iterator_type;
            typedef parameter_shape_type::const_iterator parameter_shape_const_iterator_type;

            // gestion des surfaces : shape,name,value
            typedef boost::tuple<std::string,std::string, uint > surface_type;
            typedef std::list< surface_type > surface_type_type;
            typedef std::list< surface_type_type > surface_name_type;
            typedef surface_type_type::const_iterator surface_type_const_iterator_type;
            typedef surface_name_type::const_iterator surface_name_const_iterator_type;

            // gestion des volumes : shape,name,value
            typedef boost::tuple<std::string,std::string, uint > volume_type;
            typedef std::list< volume_type > volume_type_type;
            typedef std::list< volume_type_type > volume_name_type;
            typedef volume_type_type::const_iterator volume_type_const_iterator_type;
            typedef volume_name_type::const_iterator volume_name_const_iterator_type;

            GeoGMSHTool(std::string __shape="NO_SHAPE", std::string __name="NO_NAME", double __meshSize=0.1)
                :
                _M_cptPt(1),
                _M_cptLine(1),
                _M_cptLineLoop(1),
                _M_cptSurface(1),
                _M_cptSurfaceLoop(1),
                _M_cptVolume(1),
                _M_surfaceList(new surface_name_type()),
                _M_volumeList(new volume_name_type()),
                _M_ostrSurfaceLoop( new std::ostringstream()),
                _M_map_Shape( new map_shape_names_type()),
                _M_paramShape( new parameter_shape_type()),
                //_M_markShape( new marker_shape_type()),
                _M_markShape( new marker_type_type()),
                _M_ostr( new std::ostringstream())
            {
            }

            GeoGMSHTool( GeoGMSHTool const & m )
                :
                _M_cptPt(m._M_cptPt),
                _M_cptLine(m._M_cptLine),
                _M_cptLineLoop(m._M_cptLineLoop),
                _M_cptSurface(m._M_cptSurface),
                _M_cptSurfaceLoop(m._M_cptSurfaceLoop),
                _M_cptVolume(m._M_cptVolume),

                _M_surfaceList(new surface_name_type(*(m._M_surfaceList))),
                _M_volumeList(new volume_name_type(*(m._M_volumeList))),

                _M_ostrSurfaceLoop(new std::ostringstream()),

                _M_map_Shape(new map_shape_names_type(*(m._M_map_Shape))),
                _M_paramShape(new parameter_shape_type(*(m._M_paramShape))),
                //_M_markShape(new marker_shape_type(*(m._M_markShape))),
                _M_markShape(new marker_type_type(*(m._M_markShape))),
                _M_ostr(new std::ostringstream())
            {
                updateOstr((m._M_ostr)->str());
                *_M_ostrSurfaceLoop << (m._M_ostrSurfaceLoop)->str();
            }

            void zeroCpt()
            {
                _M_cptPt=1;
                _M_cptLine=1;
                _M_cptLineLoop=1;
                _M_cptSurface=1;
                _M_cptSurfaceLoop=1;
                _M_cptVolume=1;

                surface_name_type::iterator itSurf = this->_M_surfaceList->begin();
                surface_name_type::iterator itSurf_end = this->_M_surfaceList->end();
                for ( ; itSurf != itSurf_end; ++itSurf)
                    {
                        surface_type_type::iterator itSurf2 = itSurf->begin();
                        surface_type_type::iterator itSurf2_end = itSurf->end();
                        for ( ; itSurf2 != itSurf2_end; ++itSurf2)
                            {
                                boost::get<2>(*itSurf2)=0;
                            }
                    }


                _M_ostrSurfaceLoop.reset(new std::ostringstream());

                volume_name_type::iterator itVol = this->_M_volumeList->begin();
                volume_name_type::iterator itVol_end = this->_M_volumeList->end();
                for ( ; itVol != itVol_end; ++itVol)
                    {
                        volume_type_type::iterator itVol2 = itVol->begin();
                        volume_type_type::iterator itVol2_end = itVol->end();
                        for ( ; itVol2 != itVol2_end; ++itVol2)
                            {
                                boost::get<2>(*itVol2)=0;
                            }
                    }


            }

            void
            operator=( GeoGMSHTool const & m )
            {
                _M_cptPt = m._M_cptPt;
                _M_cptLine = m._M_cptLine;
                _M_cptLineLoop = m._M_cptLineLoop;
                _M_cptSurface = m._M_cptSurface;
                _M_cptSurfaceLoop = m._M_cptSurfaceLoop;
                _M_cptVolume = m._M_cptVolume;

                _M_surfaceList.reset(new surface_name_type(*(m._M_surfaceList)));
                _M_volumeList.reset(new volume_name_type(*(m._M_volumeList)));

                _M_ostrSurfaceLoop.reset(new std::ostringstream());
                *_M_ostrSurfaceLoop << (m._M_ostrSurfaceLoop)->str();

                _M_map_Shape.reset(new map_shape_names_type(*(m._M_map_Shape)));
                _M_paramShape.reset(new parameter_shape_type(*(m._M_paramShape)));
                //_M_markShape.reset(new marker_shape_type(*(m._M_markShape)));
                _M_markShape.reset(new marker_type_type(*(m._M_markShape)));
                _M_ostr.reset(new std::ostringstream());
                updateOstr((m._M_ostr)->str());
            }

            GeoGMSHTool operator+(const GeoGMSHTool & m);
            GeoGMSHTool operator-(const GeoGMSHTool & m);

            GeoGMSHTool opFusion(const GeoGMSHTool & m,int __typeop);

            void init();

            /*
             *
             */
            void initData(std::string __shape,
                          std::string __name,
                          double __meshSize,
                          std::vector<GeoTool::Node> & __param,
                          uint dim,
                          uint __nbsurface,
                          uint __nbvolume)
            {
                boost::tuple<std::string,double> __id = boost::make_tuple(__name, __meshSize);

                (*(_M_map_Shape))[__shape].push_back(__id);

                (*(_M_paramShape))[__shape][__name].resize(__param.size());
                for (uint n=0;n<__param.size();++n)
                    {
                        (*(_M_paramShape))[__shape][__name][n] = __param[n].getNode();
                    }

                if (dim>=2)
                    {
                        //Attention 0 par defaut pour dire que ce n'est pas initialiser
                        for(uint n=0;n<__nbsurface;++n)
                            {
                                std::list< boost::tuple<std::string,std::string, uint >	>__listTemp;
                                __listTemp.push_back( boost::make_tuple(__shape,__name,0));
                                _M_surfaceList->push_back( __listTemp);
                            }
                    }
                if (dim==3)
                    {
                        //Attention 0 par defaut pour dire que ce n'est pas initialiser
                        for(uint n=0;n<__nbvolume;++n)
                            {
                                std::list< boost::tuple<std::string,std::string, uint >	>__listTemp;__listTemp.clear();
                                __listTemp.push_back( boost::make_tuple(__shape,__name,0));
                                _M_volumeList->push_back( __listTemp);
                            }

                    }

            }

            /*
             *Utile pour la fct geoStr()
             *Pas de maj pour cptSurface et cptVolume car traitement different
             */
            void updateData(GeoGMSHTool const & m)
            {
                _M_cptPt = m._M_cptPt;
                _M_cptLine = m._M_cptLine;
                _M_cptLineLoop = m._M_cptLineLoop;
                _M_cptSurfaceLoop = m._M_cptSurfaceLoop;

                _M_cptSurface = m._M_cptSurface;
                _M_cptVolume = m._M_cptVolume;

                _M_map_Shape = m._M_map_Shape;
                _M_paramShape = m._M_paramShape;
                _M_markShape = m._M_markShape;

                _M_surfaceList.reset(new surface_name_type(*(m._M_surfaceList)));
                _M_volumeList.reset(new volume_name_type(*(m._M_volumeList)));

                _M_ostrSurfaceLoop.reset(new std::ostringstream());
                *_M_ostrSurfaceLoop << (m._M_ostrSurfaceLoop)->str();
            }

            /*
             * Update the output stringstream wich generate the gmsh code
             */
            void updateOstr( std::string __str)
            {
                *_M_ostr << __str;
            }

            /*
             * Generate the gmsh code
             */
            std::string geoStr();

            /*
             * Clean
             */
            void cleanOstr() { _M_ostr.reset(new std::ostringstream()); }

            /*_________________________________________________*
             *_________________________________________________*
             * Accessor                                        *
             *_________________________________________________*
             *_________________________________________________*/

            uint cptPt() const { return _M_cptPt;}
            uint cptLine() const { return _M_cptLine;}
            uint cptLineLoop() const { return _M_cptLineLoop;}
            uint cptSurface() const { return _M_cptSurface;}
            uint cptSurfaceLoop() const { return _M_cptSurfaceLoop;}
            uint cptVolume() const { return _M_cptVolume;}

            /*_________________________________________________*
             * Parameter
             *_________________________________________________*/

            parameter_shape_const_iterator_type
            paramShapeBegin() const
            {
                return _M_paramShape->begin();
            }

            parameter_shape_const_iterator_type paramShapeEnd() const { return _M_paramShape->end();}

            parameter_name_const_iterator_type
            paramNameBegin(std::string __shape) const
            {
                return _M_paramShape->find(__shape)->second.begin();
            }

            parameter_name_const_iterator_type
            paramNameEnd(std::string __shape) const
            {
                return _M_paramShape->find(__shape)->second.end();
            }

            parameter_rectangle_type
            getParameter(std::string __shape, std::string __name) const
            {
                return _M_paramShape->find(__shape)->second.find(__name)->second;
            }

            /*_________________________________________________*
             * Marker
             *_________________________________________________*/
            /*
            marker_shape_const_iterator_type
            markShapeBegin() const
            {
                return _M_markShape->begin();
            }

            marker_shape_const_iterator_type
            markShapeEnd() const
            {
                return _M_markShape->end();
            }*/


            marker_type_const_iterator_type
            markerTypeBegin(/*std::string __shape*/) const
            {
                //return _M_markShape->find(__shape)->second.begin();
                return _M_markShape->begin();
            }

            marker_type_const_iterator_type
            markerTypeEnd(/*std::string __shape*/) const
            {
                //return _M_markShape->find(__shape)->second.end();
                return _M_markShape->end();
            }
            /*
            marker_type_type
            markerType(std::string __shape) const
            {
                //return _M_markShape->find(__shape)->second;
                }*/

            marker_markerName_const_iterator_type
            markerMarkerNameBegin(/*std::string __shape,*/ std::string __type) const
            {
                //return _M_markShape->find(__shape)->second.find(__type)->second.begin();
                return _M_markShape->find(__type)->second.begin();
            }

            marker_markerName_const_iterator_type
            markerMarkerNameEnd(/*std::string __shape,*/ std::string __type) const
            {
                //return _M_markShape->find(__shape)->second.find(__type)->second.end();
                return _M_markShape->find(__type)->second.end();
            }

            marker_markerName_type
            markerMarkerName(/*std::string __shape,*/ std::string __type) const
            {
                //return _M_markShape->find(__shape)->second.find(__type)->second;
                return _M_markShape->find(__type)->second;
            }

            std::list<marker_base_type>::const_iterator
            markerListIndiceBegin(/*std::string __shape,*/ std::string __type ,std::string __markerName) const
            {
                return _M_markShape->find(__type)->second.find(__markerName)->second.begin();
            }

            std::list<marker_base_type>::const_iterator
            markerListIndiceEnd(/*std::string __shape,*/ std::string __type ,std::string __markerName) const
            {
                //return _M_markShape->find(__shape)->second.find(__type)->second.find(__markerName)->second.end();
                return _M_markShape->find(__type)->second.find(__markerName)->second.end();
            }


            std::list<marker_base_type>
            getMarkerName(/*std::string __shape,*/ std::string __type ,std::string __markerName) const
            {
                //return _M_markShape->find(__shape)->second.find(__type)->second.find(__markerName)->second;
                return _M_markShape->find(__type)->second.find(__markerName)->second;
            }

            /*_________________________________________________*
             *_________________________________________________*
             * Members                                         *
             *_________________________________________________*
             *_________________________________________________*/


            // memory
            uint _M_cptPt;
            uint _M_cptLine;
            uint _M_cptLineLoop;
            uint _M_cptSurface;
            uint _M_cptSurfaceLoop;
            uint _M_cptVolume;

            // gestion des surface : shape,name,value
            // value is the marker associated to the planeSurface (init to 0 and to use when call geoStr())
            //std::list< std::list< boost::tuple<std::string,std::string, uint > > > _M_surfaceList;
            boost::shared_ptr<surface_name_type> _M_surfaceList;
            boost::shared_ptr<volume_name_type> _M_volumeList;

            boost::shared_ptr<std::ostringstream> _M_ostrSurfaceLoop;


            // data containers
            boost::shared_ptr<map_shape_names_type> _M_map_Shape;
            boost::shared_ptr<parameter_shape_type> _M_paramShape;
            //boost::shared_ptr<marker_shape_type> _M_markShape;
            boost::shared_ptr<marker_type_type> _M_markShape;

            // output string
            boost::shared_ptr<std::ostringstream> _M_ostr;

        };

        GeoGMSHTool
        GeoGMSHTool::operator-(const GeoGMSHTool & m)
        {
            return this->opFusion(m,2);
        }

        GeoGMSHTool
        GeoGMSHTool::operator+(const GeoGMSHTool & m)
        {
            return this->opFusion(m,1);
        }



        GeoGMSHTool
        GeoGMSHTool::opFusion(const GeoGMSHTool & m,int __typeOp)
        {

            GeoGMSHTool __geoTool;

            //Add Plane Surface for operator + : (((rect,u1,_)))+(((circ,u2,_))) -> (((rect,u1,_)),((circ,u2,_)))
            if (__typeOp==1)
                {
                    __geoTool._M_surfaceList.reset(new surface_name_type(*(this->_M_surfaceList)));
                    surface_name_const_iterator_type itSurf = m._M_surfaceList->begin();
                    surface_name_const_iterator_type itSurf_end = m._M_surfaceList->end();
                    for ( ; itSurf != itSurf_end; ++itSurf)
                        {
                            surface_type_const_iterator_type itSurf2 = itSurf->begin();
                            surface_type_const_iterator_type itSurf2_end = itSurf->end();
                            surface_type_type __listTemp;

                            for ( ; itSurf2 != itSurf2_end; ++itSurf2)
                                {
                                    __listTemp.push_back(*itSurf2);
                                }
                            __geoTool._M_surfaceList->push_back(__listTemp);
                            }
                }
            // Add Plane Surface for operator - : (((rect,u1,_)))-(((circ,u2,_))) -> (((rect,u1,_),(circ,u2,_)))
            else if (__typeOp==2)
                {
                    __geoTool._M_surfaceList.reset(new surface_name_type(*(this->_M_surfaceList)));
                    surface_name_const_iterator_type itSurf = m._M_surfaceList->begin();
                    surface_name_const_iterator_type itSurf_end = m._M_surfaceList->end();
                    for ( ; itSurf != itSurf_end; ++itSurf)
                        {
                            surface_type_const_iterator_type itSurf2 = itSurf->begin();
                            surface_type_const_iterator_type itSurf2_end = itSurf->end();
                            surface_type_type __listTemp;
                            for ( ; itSurf2 != itSurf2_end; ++itSurf2)
                                {
                                    __geoTool._M_surfaceList->begin()->push_back(*itSurf2);
                                }
                        }
                }

            //get data from this (easy)
            __geoTool._M_map_Shape.reset( new map_shape_names_type( *this->_M_map_Shape) );
            __geoTool._M_paramShape.reset( new parameter_shape_type( *this->_M_paramShape) );
            //__geoTool._M_markShape.reset( new marker_shape_type (*this->_M_markShape) );
            __geoTool._M_markShape.reset( new marker_type_type (*this->_M_markShape) );

            //get data from (more hard because no duplication)
            map_shape_names_const_iterator_type itShape = m._M_map_Shape->begin();
            map_shape_names_const_iterator_type itShape_end = m._M_map_Shape->end();
            while (itShape!=itShape_end)
                {
                    names_const_iterator_type itName = itShape->second.begin();
                    names_const_iterator_type itName_end = itShape->second.end();
                    while (itName!=itName_end)
                        {
                            //verify that should not duplicate the name of shape
                            bool __find=false;
                            names_base_type __temp;
                            names_const_iterator_type itLs = (*(__geoTool._M_map_Shape))[itShape->first].begin();
                            names_const_iterator_type itLs_end = (*(__geoTool._M_map_Shape))[itShape->first].end();
                            while (itLs!=itLs_end) {
                                if ( boost::get<0>(*itLs) == boost::get<0>(*itName)) __find=true;
                                ++itLs;
                            }
                            if (!__find) (*(__geoTool._M_map_Shape))[itShape->first].push_back(*itName);

                            (*(__geoTool. _M_paramShape))[itShape->first][ boost::get<0>(*itName)] = m.getParameter(itShape->first,boost::get<0>(*itName));

                            ++itName;
                        }
                    ++itShape;
                }

            //update marker
            //marker_shape_const_iterator_type itShape2 = m.markShapeBegin();
            //marker_shape_const_iterator_type itShape2_end = m.markShapeEnd();
            //while (itShape2!=itShape2_end)
            //    {
                    marker_type_const_iterator_type itMarkType = m.markerTypeBegin(/*itShape2->first*/);
                    marker_type_const_iterator_type itMarkType_end = m.markerTypeEnd(/*itShape2->first*/);
                        while (itMarkType!=itMarkType_end)
                        {
                            if (__typeOp==1 || itMarkType->first=="line")
                                {
                                    marker_markerName_const_iterator_type itMarkName = m.markerMarkerNameBegin(/*itShape2->first,*/
                                                                                                               itMarkType->first);
                                    marker_markerName_const_iterator_type itMarkName_end = m.markerMarkerNameEnd(/*itShape2->first,*/
                                                                                                                 itMarkType->first);
                                    while (itMarkName!=itMarkName_end)
                                        {
                                            //Est-ce UTILE?????
                                            if ( !m.getMarkerName(/*itShape2->first,*/itMarkType->first, itMarkName->first).empty() )
                                                {
                                                    std::list<marker_base_type>::const_iterator itLRef=m.markerListIndiceBegin(/*itShape2->first,*/
                                                                                                                               itMarkType->first,
                                                                                                                               itMarkName->first);
                                                    std::list<marker_base_type>::const_iterator itLRef_end=m.markerListIndiceEnd(/*itShape2->first,*/
                                                                                                                                 itMarkType->first,
                                                                                                                                 itMarkName->first);
                                                    while(itLRef != itLRef_end )
                                                        {
                                                            (*(__geoTool._M_markShape))/*["rectangle"itShape2->first]*/[itMarkType->first][itMarkName->first].push_back(*itLRef);
                                                            ++itLRef;
                                                        }
                                                }
                                            ++itMarkName;
                                        }
                                }
                            ++itMarkType;
                        }
                        //++itShape2;
                        // }

            return __geoTool;
        }


        void
        GeoGMSHTool::init()
        {
            *_M_ostr << "Mesh.MshFileVersion = " << 2 << ";\n";
        }



        std::string
        GeoGMSHTool::geoStr()
        {
            // clean
            _M_ostr.reset( new std::ostringstream());
            this->zeroCpt();

            // version of gmsh
            init();

            //data memory ( type->shape->name )
            std::vector<std::map<std::string,std::map<std::string, std::map<uint,uint> > > > __dataMemGlob(6);

            //iterate on the shape
            map_shape_names_const_iterator_type itShape =  _M_map_Shape->begin();
            map_shape_names_const_iterator_type itShape_end = _M_map_Shape->end();
            while (itShape!=itShape_end)
                {
                    //iterate on the name of shape
                    names_const_iterator_type itName = itShape->second.begin();
                    names_const_iterator_type itName_end = itShape->second.end();
                    while (itName!=itName_end)
                        {

                            *_M_ostr << "h=" << boost::get<1>(*itName) << ";\n";
                            //data memory
                            vec_map_data_ptrtype __dataMem(new vec_map_data_type(6));

                            GeoGMSHTool_ptrtype __geoTool(new GeoGMSHTool());
                            __geoTool->updateData(*this);
                            __geoTool->cleanOstr();
                            GeoTool::data_geo_ptrtype __data_geoTool(new GeoTool::data_geo_type(boost::make_tuple( __geoTool,
                                                                                                                   __dataMem,
                                                                                                                   itShape->first,
                                                                                                                   boost::get<0>(*itName)
                                                                                                                   )));
                            // generate the code for the geometry
                            run(__data_geoTool);

                            __dataMemGlob[0][itShape->first][boost::get<0>(*itName)] = (*(boost::get<1>(*__data_geoTool)))[0];//pts
                            __dataMemGlob[1][itShape->first][boost::get<0>(*itName)] = (*(boost::get<1>(*__data_geoTool)))[1];//lines
                            __dataMemGlob[2][itShape->first][boost::get<0>(*itName)] = (*(boost::get<1>(*__data_geoTool)))[2];//lineLoop
                            __dataMemGlob[3][itShape->first][boost::get<0>(*itName)] = (*(boost::get<1>(*__data_geoTool)))[3];//Surface
                            __dataMemGlob[4][itShape->first][boost::get<0>(*itName)] = (*(boost::get<1>(*__data_geoTool)))[4];//SurfaceLoop
                            __dataMemGlob[5][itShape->first][boost::get<0>(*itName)] = (*(boost::get<1>(*__data_geoTool)))[5];//Volume

                            // get infos
                            this->updateData( *boost::get<0>(*__data_geoTool));
                            this->updateOstr( boost::get<0>(*__data_geoTool)->_M_ostr->str());

                            ++itName;
                        }
                    ++itShape;
                }


            //Write the planes surfaces
            //Fait ici a cause des opertateurs (+,-)
            surface_name_const_iterator_type itSurf = this->_M_surfaceList->begin();
            surface_name_const_iterator_type itSurf_end = this->_M_surfaceList->end();
            std::ostringstream __surface_str;
            //counter of surface
            uint __surfnumber=1;
            for ( ; itSurf != itSurf_end; ++itSurf)
                {
                    surface_type_const_iterator_type itSurf2 = itSurf->begin();
                    //Attention : On fait a cause des op - : sinon les markers surfaces sont incorrectes(l'idee est de marquer la 1ere sous-surface)
                    __dataMemGlob[3][boost::get<0>(*itSurf2)][boost::get<1>(*itSurf2)][__surfnumber]=__surfnumber;

                    surface_type_const_iterator_type itSurf2_end = --itSurf->end();

                    __surface_str << "Plane Surface(" << __surfnumber << ") = {" ;

                    for ( ; itSurf2 != itSurf2_end; ++itSurf2)
                        {
                            __surface_str << boost::get<2>(*itSurf2) << ",";
                        }
                    __surface_str << boost::get<2>(*itSurf2);

                    __surface_str << "};\n";

                    ++__surfnumber;

                }
            this->updateOstr(__surface_str.str());

            //Write the surfaces loops
            this->updateOstr(_M_ostrSurfaceLoop->str());


            //Write the volumes
            //Fait ici a cause des opertateurs (+,-)
            volume_name_const_iterator_type itVol = this->_M_volumeList->begin();
            volume_name_const_iterator_type itVol_end = this->_M_volumeList->end();
            std::ostringstream __volume_str;
            //counter of surface
            uint __volnumber=1;
            for ( ; itVol != itVol_end; ++itVol)
                {
                    volume_type_const_iterator_type itVol2 = itVol->begin();
                    volume_type_const_iterator_type itVol2_end = --itVol->end();

                    __volume_str << "Volume(" << __volnumber << ") = {" ;

                    for ( ; itVol2 != itVol2_end; ++itVol2)
                        {
                            __volume_str << boost::get<2>(*itVol2) << ",";
                        }
                    __volume_str << boost::get<2>(*itVol2);

                    __volume_str << "};\n";

                    ++__volnumber;

                }
            this->updateOstr(__volume_str.str());


            //write marker
            //iterate on the shape
            //itShape =  _M_map_Shape->begin();
            //itShape_end = _M_map_Shape->end();
            //while (itShape!=itShape_end)
            //{
                    // generate the code for the marker
                    marker_type_const_iterator_type itMarkType= (*(_M_markShape))/*[itShape->first]*/.begin();
                    marker_type_const_iterator_type itMarkType_end=(*(_M_markShape))/*[itShape->first]*/.end();
                    while (itMarkType!=itMarkType_end)
                        {
                            marker_markerName_const_iterator_type itMarkName = (*(_M_markShape))/*[itShape->first]*/[itMarkType->first].begin();
                            marker_markerName_const_iterator_type itMarkName_end=(*(_M_markShape))/*[itShape->first]*/[itMarkType->first].end();
                            while ( itMarkName!=itMarkName_end)
                                {
                                    if (itMarkType->first=="line")
                                        {
                                            *_M_ostr << "Physical Line(\"" << itMarkName->first << "\") = {";

                                            std::list<marker_base_type>::const_iterator itMark = itMarkName->second.begin();
                                            std::list<marker_base_type>::const_iterator itMark_end = --itMarkName->second.end();
                                            while (itMark!=itMark_end)
                                                {
                                                    *_M_ostr << __dataMemGlob[1][boost::get<0>(*itMark)][boost::get<1>(*itMark)][boost::get<2>(*itMark)]<<",";
                                                    ++itMark;
                                                }
                                            *_M_ostr << __dataMemGlob[1][boost::get<0>(*itMark)][boost::get<1>(*itMark)][boost::get<2>(*itMark)] << "};\n";
                                        }
                                    else if (itMarkType->first=="surface")
                                        {
                                            *_M_ostr << "Physical Surface(\"" << itMarkName->first << "\") = {";

                                            std::list<marker_base_type>::const_iterator itMark = itMarkName->second.begin();
                                            std::list<marker_base_type>::const_iterator itMark_end = --itMarkName->second.end();
                                            while (itMark!=itMark_end)
                                                {
                                                    *_M_ostr << __dataMemGlob[3][boost::get<0>(*itMark)][boost::get<1>(*itMark)][boost::get<2>(*itMark)]<<",";
                                                    ++itMark;
                                                }
                                            *_M_ostr << __dataMemGlob[3][boost::get<0>(*itMark)][boost::get<1>(*itMark)][boost::get<2>(*itMark)] << "};\n";

                                        }
                                    ++itMarkName;
                                }
                            ++itMarkType;
                        }
                    //++itShape;
                    //}


            return _M_ostr->str();

        }



        /*_________________________________________________*
         *_________________________________________________*
         *_________________________________________________*
         * Function on the namespace                       *
         *_________________________________________________*
         *_________________________________________________*
         *_________________________________________________*/


#define GEOTOOL_GENERATE_RUN(r,state)                                   \
        if( boost::get<2>(*__dg) ==  GEOTOOL_SHAPE_NAME_STR(BOOST_PP_TUPLE_ELEM(2,0,state)) ) \
            {                                                           \
                BOOST_PP_CAT(run,GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))) (__dg); \
            }                                                           \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/


        // A refaire avec les boost pp
        void run( data_geo_ptrtype __dg)
        {
            BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE(GEOTOOL_SHAPE),1) ),
                          GEOTOOL_INSTANTIATES_FOR_COMP,
                          GEOTOOL_INSTANTIATES_FOR_INCR,
                          GEOTOOL_GENERATE_RUN );
        }

        /*_________________________________________________*
         *_________________________________________________*
         * Function on the namespace                       *
         *_________________________________________________*
         *_________________________________________________*/

        template <uint Numero>
        node_type
        param(data_geo_ptrtype __dg)
        {
            std::string __shape = boost::get<2>(*__dg);
            std::string __name = boost::get<3>(*__dg);
            //node_type __node = boost::get<Numero>(boost::get<0>(*__dg)->_M_paramShape->find(__shape)->second[__name]);
            node_type __node =(boost::get<0>(*__dg)->_M_paramShape->find(__shape)->second[__name])[Numero];

            __node.resize(3);
            if (__node.size()<2)
                __node(1)=0.0;
            if (__node.size()<3)
                __node(2)=0.0;

            return __node;
        }

        /*_________________________________________________*
         *_________________________________________________*
         * Function on the namespace                       *
         *_________________________________________________*
         *_________________________________________________*/


        void writePoint(uint __numLoc, data_geo_ptrtype __dg ,double __x1,double __x2, double __x3=0)
        {
            (*(boost::get<1>(*__dg)))[0][__numLoc] = boost::get<0>(*__dg)->cptPt(); //            __mapPt[0][__numLoc] = boost::get<0>(*__dg)->cptPt();
            std::ostringstream __ostr;
            __ostr << "Point(" << boost::get<0>(*__dg)->cptPt()
                   << ") = {"
                   << __x1 << ","
                   << __x2 << ","
                   << __x3 <<", h};\n";
            boost::get<0>(*__dg)->updateOstr(__ostr.str());
            ++(boost::get<0>(*__dg)->_M_cptPt);
        }

        /*_________________________________________________*/

        void
        writeLine(uint __numLoc, data_geo_ptrtype __dg ,uint __n1, uint __n2)
        {
            (*(boost::get<1>(*__dg)))[1][__numLoc] = boost::get<0>(*__dg)->cptLine();
            std::ostringstream __ostr;
            __ostr << "Line(" << boost::get<0>(*__dg)->cptLine()
                   << ") = {"
                   << (*(boost::get<1>(*__dg)))[0][__n1] << ","
                   << (*(boost::get<1>(*__dg)))[0][__n2] << "};\n";

            boost::get<0>(*__dg)->updateOstr(__ostr.str());

            ++boost::get<0>(*__dg)->_M_cptLine;
        }

        /*_________________________________________________*/

        void
        writeCircle(uint __numLoc, data_geo_ptrtype __dg ,uint __n1, uint __n2, uint __n3)
        {
            (*(boost::get<1>(*__dg)))[1][__numLoc] = boost::get<0>(*__dg)->cptLine();

            std::ostringstream __ostr;
            __ostr << "Circle(" << boost::get<0>(*__dg)->cptLine() //cptCircle
                   << ") = {"
                   << (*(boost::get<1>(*__dg)))[0][__n1] << ","
                   << (*(boost::get<1>(*__dg)))[0][__n2] << ","
                   << (*(boost::get<1>(*__dg)))[0][__n3] << "};\n";

            boost::get<0>(*__dg)->updateOstr(__ostr.str());

            //++boost::get<0>(*__dg)->_M_cptCircle;
            ++boost::get<0>(*__dg)->_M_cptLine;
        }
        /*_________________________________________________*/

        void
        writeSpline(uint __numLoc, data_geo_ptrtype __dg ,Loop __loop)
        {
            (*(boost::get<1>(*__dg)))[1][__numLoc] = boost::get<0>(*__dg)->cptLine();

            std::ostringstream __ostr;
            __ostr << "Spline(" << boost::get<0>(*__dg)->cptLine() //cptCircle
                   << ") = {";

            std::list<int>::const_iterator it= __loop.begin();
            std::list<int>::const_iterator it_end= --__loop.end();
            while (it!=it_end)
                {
                    __ostr << (*(boost::get<1>(*__dg)))[0][*it] <<"," ;
                    ++it;
                }
            __ostr << (*(boost::get<1>(*__dg)))[0][*it] << "};\n";

            boost::get<0>(*__dg)->updateOstr(__ostr.str());

            ++boost::get<0>(*__dg)->_M_cptLine;
        }

        /*_________________________________________________*/

        void
        writeLineLoop(uint __numLoc, data_geo_ptrtype __dg , Loop /*const*/ __loop )
        {
             (*(boost::get<1>(*__dg)))[2][__numLoc] = boost::get<0>(*__dg)->cptLineLoop();

            std::ostringstream __ostr;
            __ostr << "Line Loop(" << boost::get<0>(*__dg)->cptLineLoop()
                   << ") = {" ;
            std::list<int>::const_iterator it= __loop.begin();
            std::list<int>::const_iterator it_end= --__loop.end();
            while (it!=it_end)
                {
                    if (*it>0)
                        __ostr << (*(boost::get<1>(*__dg)))[1][*it] <<"," ;
                    else
                        __ostr << "-" << (*(boost::get<1>(*__dg)))[1][-*it] <<"," ;
                    ++it;
                }

            if (*it>0)
                __ostr << (*(boost::get<1>(*__dg)))[1][*it] << "};\n";
            else
                __ostr << "-" << (*(boost::get<1>(*__dg)))[1][-*it] << "};\n";

            boost::get<0>(*__dg)->updateOstr(__ostr.str());

            ++boost::get<0>(*__dg)->_M_cptLineLoop;
        }

        /*_________________________________________________*/


        //ici on n'ecrit pas, on memorise cause des operations de difference
        //l'ecriture est realise dans geoStr()
        void
        writePlaneSurface(uint __numLoc, data_geo_ptrtype __dg , uint __ind)
        {
            (*(boost::get<1>(*__dg)))[3][__numLoc] = boost::get<0>(*__dg)->cptSurface();

            bool __find=false;
            //Memorize in surfaceList
            GeoGMSHTool::surface_name_type::iterator itSurf = boost::get<0>(*__dg)->_M_surfaceList->begin();
            GeoGMSHTool::surface_name_type::iterator itSurf_end = boost::get<0>(*__dg)->_M_surfaceList->end();
            for( ; itSurf !=itSurf_end;++itSurf)
                {
                    GeoGMSHTool::surface_type_type::iterator itSurf2 = itSurf->begin();
                    GeoGMSHTool::surface_type_type::iterator itSurf2_end = itSurf->end();
                    while (itSurf2 !=itSurf2_end)
                        {
                            if (boost::get<0>(*itSurf2) == boost::get<2>(*__dg))
                                {
                                    if (boost::get<1>(*itSurf2) == boost::get<3>(*__dg))
                                        {
                                            //on cherche la 1ere surface non init
                                            if (boost::get<2>(*itSurf2) == 0 && !__find)
                                                {
                                                    boost::get<2>(*itSurf2) = (*(boost::get<1>(*__dg)))[2][__ind];
                                                    __find=true;
                                                }
                                        }
                                }
                            ++itSurf2;
                        }
                }
            ++boost::get<0>(*__dg)->_M_cptSurface;
        }

        void
        writeSurfaceLoop(uint __numLoc, data_geo_ptrtype __dg , Loop /*const*/ __loop )
        {
            (*(boost::get<1>(*__dg)))[4][__numLoc] = boost::get<0>(*__dg)->cptSurfaceLoop();
            std::ostringstream __ostr;
            __ostr << "Surface Loop(" << boost::get<0>(*__dg)->cptSurfaceLoop()
                   << ") = {" ;
            std::list<int>::const_iterator it= __loop.begin();
            std::list<int>::const_iterator it_end= --__loop.end();
            while (it!=it_end)
                {
                    if (*it>0)
                        __ostr << (*(boost::get<1>(*__dg)))[3][*it] <<"," ;
                    else
                        __ostr << "-" << (*(boost::get<1>(*__dg)))[3][-*it] <<"," ;
                    ++it;
                }

            if (*it>0)
                __ostr << (*(boost::get<1>(*__dg)))[3][*it] << "};\n";
            else
                __ostr << "-" << (*(boost::get<1>(*__dg)))[3][-*it] << "};\n";

            //boost::get<0>(*__dg)->updateOstr(__ostr.str());
            *(boost::get<0>(*__dg)->_M_ostrSurfaceLoop) << __ostr.str();

            ++boost::get<0>(*__dg)->_M_cptSurfaceLoop;
        }

        //ici on n'ecrit pas, on memorise cause des operations de difference
        //l'ecriture est realise dans geoStr()
        void
        writeVolume(uint __numLoc, data_geo_ptrtype __dg , uint __ind)
        {
            (*(boost::get<1>(*__dg)))[5][__numLoc] = boost::get<0>(*__dg)->cptVolume();

            bool __find=false;
            //Memorize in volumeList
            GeoGMSHTool::volume_name_type::iterator itSurf = boost::get<0>(*__dg)->_M_volumeList->begin();
            GeoGMSHTool::volume_name_type::iterator itSurf_end = boost::get<0>(*__dg)->_M_volumeList->end();
            for( ; itSurf !=itSurf_end;++itSurf)
                {
                    GeoGMSHTool::volume_type_type::iterator itSurf2 = itSurf->begin();
                    GeoGMSHTool::volume_type_type::iterator itSurf2_end = itSurf->end();
                    while (itSurf2 !=itSurf2_end)
                        {
                            if (boost::get<0>(*itSurf2) == boost::get<2>(*__dg))
                                {
                                    if (boost::get<1>(*itSurf2) == boost::get<3>(*__dg))
                                        {
                                            //on cherche la 1ere surface non init
                                            if (boost::get<2>(*itSurf2) == 0 && !__find)
                                                {
                                                    boost::get<2>(*itSurf2) = (*(boost::get<1>(*__dg)))[4][__ind];
                                                    __find=true;
                                                }
                                        }
                                }
                            ++itSurf2;
                        }
                }
            ++boost::get<0>(*__dg)->_M_cptVolume;
        }


        /*_________________________________________________*/


        /*_________________________________________________*
         *_________________________________________________*
         *_________________________________________________*
         * PREPROCESSOR METHODS                            *
         *_________________________________________________*
         *_________________________________________________*
         *_________________________________________________*/




#define GEOTOOL_FOR_COMP2(r, state)                                     \
        BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(4, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(4, 1, state)) \
                            )                                           \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_FOR_INCR2(r, state)                         \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(4, 0, state)),    \
         BOOST_PP_TUPLE_ELEM(4, 1, state),                  \
         BOOST_PP_TUPLE_ELEM(4, 2, state),                  \
         BOOST_PP_TUPLE_ELEM(4, 3, state) )                 \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_FOR_MARKER_LINE_MACRO2(r, state)                        \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_LINE_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_LINE_, \
                                                                                             GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) )), \
                                                                                BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_FOR_MARKER_SURFACE_MACRO2(r, state)                     \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_SURFACE_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                                GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) )), \
                                                                                   BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                   BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_FOR_MARKER_VOLUME_MACRO2(r, state)                      \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_VOLUME_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                               BOOST_PP_IF(BOOST_PP_GREATER(GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(4,3,state)) ,0), \
                                                                                                           GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state)),\
                                                                                                           DEFAULT)), \
                                                                                  BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                  BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_FOR_MARKER_LINE_MACRO(r, state)                         \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_LINE_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_LINE_, \
                                                                                    GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state))), \
                                                                       BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_LINE_MACRO2) \
                    }                                                   \

        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_FOR_MARKER_SURFACE_MACRO(r, state)                      \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_SURFACE_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                       GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state)) ), \
                                                                          BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_SURFACE_MACRO2) \
                    }                                                   \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_FOR_MARKER_VOLUME_MACRO(r, state)                       \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_VOLUME_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                      BOOST_PP_IF(BOOST_PP_GREATER(GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(3,2,state)) ,0), \
                                                                                                  GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state)),\
                                                                                                  DEFAULT )), \
                                                                          BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_VOLUME_MACRO2) \
                    }                                                   \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/

#define GEOTOOL_FOR_COMP1(r, state)                                     \
        BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(3, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(3, 1, state)) \
                            )                                           \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_FOR_INCR1(r, state)                         \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(3, 0, state)),	\
         BOOST_PP_TUPLE_ELEM(3, 1, state),                  \
         BOOST_PP_TUPLE_ELEM(3, 2, state) )                 \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_FOR_COMP(r, state)                                      \
        BOOST_PP_LESS( BOOST_PP_TUPLE_ELEM(2, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 1, state)) \
                            )                                           \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_FOR_INCR(r, state)                          \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 0, state)),	\
         BOOST_PP_TUPLE_ELEM(2, 1, state) )                 \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_SHAPE_PARAM(r, state)                                   \
        _M_param[BOOST_PP_TUPLE_ELEM(2,0,state)] = BOOST_PP_CAT( __param, \
                                                                 BOOST_PP_TUPLE_ELEM(2,0,state) ); \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_SHAPE_FOR_PARAM_SIGNATURE(r, state)                     \
        Node BOOST_PP_CAT( __param, BOOST_PP_TUPLE_ELEM(2,0,state) ) BOOST_PP_COMMA() \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_SHAPE_PARAM_SIGNATURE(state)                           \
        BOOST_PP_FOR( (0, BOOST_PP_SUB(GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state)),1) ), \
                      GEOTOOL_FOR_COMP,                                 \
                      GEOTOOL_FOR_INCR,                                 \
                      GEOTOOL_SHAPE_FOR_PARAM_SIGNATURE)                \
        /**/
        /*_________________________________________________*/
        /*                                                 */
        /**/
#define GEOTOOL_SHAPE_CLASS(r,state)                                    \
        class GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state)) : public GeoGMSHTool \
        {                                                               \
        public :                                                        \
                                                                        \
            typedef GeoGMSHTool::node_type node_type;                   \
            typedef GeoTool::Node Node;                                 \
                                                                        \
                                                                        \
            GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))(double __meshSize, \
                                                                     std::string __name, \
                                                                     GEOTOOL_SHAPE_PARAM_SIGNATURE(state) \
                                                                     uint type = 0 ) /*Ne sert a rien, juste a cause de la virgule au dessus)*/ \
            /*Node __param0,                                            \
              Node __param1 )*/                                         \
                :                                                       \
                GeoGMSHTool( shape(), __name, __meshSize),              \
                _M_name(__name)                                         \
                {                                                       \
                    _M_param.resize( GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state))); \
                    BOOST_PP_FOR( (0, BOOST_PP_SUB(GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state)),1) ), \
                                  GEOTOOL_FOR_COMP,                     \
                                  GEOTOOL_FOR_INCR,                     \
                                  GEOTOOL_SHAPE_PARAM);                 \
                                                                        \
                    initData(shape(),                                   \
                             __name,                                    \
                             __meshSize,                                \
                             _M_param,                                  \
                             GEOTOOL_SHAPE_DIM(BOOST_PP_TUPLE_ELEM(2,0,state)),\
                             GEOTOOL_SHAPE_NBSURFACE(BOOST_PP_TUPLE_ELEM(2,0,state)), \
                             GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(2,0,state))); \
                }                                                       \
                                                                        \
                                                                        \
                                                                        \
                                                                        \
            GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))(const GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state)) & m) \
                :                                                       \
                GeoGMSHTool(m),                                         \
                _M_param(m._M_param)                                    \
                    {}                                                  \
                                                                        \
            BOOST_PARAMETER_MEMBER_FUNCTION(                            \
                                            (void),                     \
                                            setMarker,                  \
                                            tag,                        \
                                            (required                   \
                                             ( type, (std::string))		\
                                             ( name, (std::string)) )   \
                                            (optional                   \
                                             (markerAll, (bool), false) \
                                             (marker1, (bool), false)   \
                                             (marker2, (bool), false)   \
                                             (marker3, (bool), false)   \
                                             (marker4, (bool), false)   \
                                             (marker5, (bool), false)   \
                                             (marker6, (bool), false)   \
                                             (marker7, (bool), false)   \
                                             (marker8, (bool), false)   \
                                             (marker9, (bool), false)   \
                                             (marker10, (bool), false)   \
                                             (marker11, (bool), false)   \
                                             (marker12, (bool), false)   \
                                             ))                         \
                {                                                       \
                                                                        \
                    if (markerAll) {                                    \
                        marker1=true;                                   \
                        marker2=true;                                   \
                        marker3=true;                                   \
                        marker4=true;                                   \
                        marker5=true;                                   \
                        marker6=true;                                   \
                        marker7=true;                                   \
                        marker8=true;                                   \
                        marker9=true;                                   \
                        marker10=true;                                   \
                        marker11=true;                                   \
                        marker12=true;                                   \
                    }                                                   \
                                                                        \
                    std::list<marker_base_type > __listMarker = (*(_M_markShape))/*[this->shape()]*/[type][name]; \
                                                                        \
                                                                        \
                    if (type=="line")                                   \
                        {                                               \
                            BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE( BOOST_PP_CAT(GEOTOOL_MARKER_LINE_, \
                                                                                             GEOTOOL_SHAPE_NAME_MACRO( BOOST_PP_TUPLE_ELEM(2,0,state)))), \
                                                           1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
                                          GEOTOOL_FOR_COMP1,            \
                                          GEOTOOL_FOR_INCR1,            \
                                          GEOTOOL_FOR_MARKER_LINE_MACRO) \
                                }                                       \
                    else if (type=="surface")                           \
                        {                                               \
                            BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE( BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                             GEOTOOL_SHAPE_NAME_MACRO( BOOST_PP_TUPLE_ELEM(2,0,state)))), \
                                                           1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
                                          GEOTOOL_FOR_COMP1,            \
                                          GEOTOOL_FOR_INCR1,            \
                                          GEOTOOL_FOR_MARKER_SURFACE_MACRO) \
                                }                                       \
                    else if (type=="volume")                            \
                        {                                               \
                            BOOST_PP_FOR( (0, BOOST_PP_SUB(  GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(2,0,state)), \
                                                           1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
                                          GEOTOOL_FOR_COMP1,            \
                                          GEOTOOL_FOR_INCR1,            \
                                          GEOTOOL_FOR_MARKER_VOLUME_MACRO) \
                                }                                       \
                                                                        \
                    (*(_M_markShape))/*[this->shape()]*/[type][name] = __listMarker; \
                }                                                       \
                                                                        \
                                                                        \
            std::string _M_name;                                        \
                                                                        \
            static const std::string shape() { return GEOTOOL_SHAPE_NAME_STR(BOOST_PP_TUPLE_ELEM(2,0, state));} \
            const std::string name() const {return _M_name;}			\
                                                                        \
                                                                        \
            std::vector<GeoTool::Node> _M_param;                        \
                                                                        \
        };                                                              \
        /**/
        /*_________________________________________________*/
        /*_________________________________________________*/
        /*                                                 */
        /**/



        //creation des classes representants les objets geotool
        BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE(GEOTOOL_SHAPE),1) ),
                      GEOTOOL_FOR_COMP,
                      GEOTOOL_FOR_INCR,
                      GEOTOOL_SHAPE_CLASS );





        /*_________________________________________________*
         *_________________________________________________*
         * Function user                                   *
         *_________________________________________________*
         *_________________________________________________*/

        void
        runRectangle(data_geo_ptrtype dg)
        {

            node_type PtA = param<0>(dg);
            node_type PtB = param<1>(dg);

            writePoint( 1, dg , PtA(0), PtA(1) );
            writePoint( 2, dg , PtB(0), PtA(1) );
            writePoint( 3, dg , PtB(0), PtB(1) );
            writePoint( 4, dg , PtA(0), PtB(1) );

            writeLine( 1, dg , 1 , 2);
            writeLine( 2, dg , 2 , 3);
            writeLine( 3, dg , 3 , 4);
            writeLine( 4, dg , 4 , 1);

            writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4);

            writePlaneSurface( 1, dg, 1);

        }

        void
        runQuadrangle(data_geo_ptrtype dg)
        {

            node_type PtA = param<0>(dg);
            node_type PtB = param<1>(dg);
            node_type PtC = param<2>(dg);
            node_type PtD = param<3>(dg);

            writePoint( 1, dg , PtA(0), PtA(1) );
            writePoint( 2, dg , PtB(0), PtB(1) );
            writePoint( 3, dg , PtC(0), PtC(1) );
            writePoint( 4, dg , PtD(0), PtD(1) );

            writeLine( 1, dg , 1 , 2);
            writeLine( 2, dg , 2 , 3);
            writeLine( 3, dg , 3 , 4);
            writeLine( 4, dg , 4 , 1);

            writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4);

            writePlaneSurface( 1, dg, 1);

        }


        void
        runCircle(data_geo_ptrtype dg)
        {

            node_type PtA = param<0>(dg);
            node_type PtB = param<1>(dg);

            writePoint( 1, dg , PtB(0), PtB(1) );
            writePoint( 2, dg , PtA(0), PtA(1) );
            writePoint( 3, dg , 2*PtA(0)-PtB(0), 2*PtA(1)-PtB(1) );

            writeCircle( 1, dg, 1, 2, 3);
            writeCircle( 2, dg, 3, 2, 1);

            writeLineLoop( 1, dg, Loop()>>1>>2);

            writePlaneSurface( 1, dg, 1);
        }



        void
        runPartialDisque(data_geo_ptrtype dg)
        {

            node_type center = param<0>(dg);
            node_type rayon = param<1>(dg);
            node_type angle1 = param<2>(dg);
            node_type angle2 = param<3>(dg);

#if 0
            node base1;

            writePoint( 1, dg , center(0), center(1) );
            writePoint( 2, dg , ptdeb(0), ptdeb(1) );
            writePoint( 3, dg , ptfin(0), ptfin(1) );

            if (angle1<90)
                {
                    if (angle2>90)
                        {
                            writePoint( 4, dg , ref1(0), ref1(1) );
                            if (angle2>180)
                                {
                                    writePoint( 5, dg , ref2(0), ref2(1) );
                                    if (angle2>270)
                                        {
                                            writePoint( 6, dg , ref3(0), ref3(1) );
                                        }
                                }
                        }
                }
            else if (angle1<180)
                {
                    if (angle2>180)
                        {
                            writePoint( 6, dg , ref2(0), ref2(1) );
                            if (angle2>270)
                                {
                                    writePoint( 6, dg , ref3(0), ref3(1) );
                                }
                        }
                }
            else if (angle1<270)
                {
                    if (angle2>270)
                        {
                            writePoint( 6, dg , ref3(0), ref3(1) );
                        }
                }

#endif
            writeCircle( 1, dg, 1, 2, 3);
            writeCircle( 2, dg, 3, 2, 1);

            writeLineLoop( 1, dg, Loop()>>1>>2);

            writePlaneSurface( 1, dg, 1);


        }

        void
        runSpecial_1a(data_geo_ptrtype dg)
        {
            double yh=1.0;
            Node a1( 0.0, yh);
            Node a2( 3.0, yh);
            Node a3( 6.0, yh+0.7);
            Node a4( 7.3, yh-0.5);
            Node a5( 8.5, yh);
            Node a6(11.0, yh);
            double ep=0.3;
            //_______________________________________________//
            writePoint( 1, dg , a1(0), a1(1) );
            writePoint( 2, dg , a2(0), a2(1) );
            writePoint( 3, dg , a3(0), a3(1) );
            writePoint( 4, dg , a4(0), a4(1) );
            writePoint( 5, dg , a5(0), a5(1) );
            writePoint( 6, dg , a6(0), a6(1) );

            writeSpline( 1, dg, Loop()>>1>>2>>3>>4>>5>>6);

            writePoint( 7, dg , a1(0), a1(1)+ep );
            writePoint( 8, dg , a2(0), a2(1)+ep );
            writePoint( 9, dg , a3(0), a3(1)+ep );
            writePoint( 10, dg , a4(0), a4(1)+ep );
            writePoint( 11, dg , a5(0), a5(1)+ep );
            writePoint( 12, dg , a6(0), a6(1)+ep );

            writeSpline( 2, dg, Loop()>>7>>8>>9>>10>>11>>12);

            writeLine( 3, dg, 1,7);
            writeLine( 4, dg, 6,12);

            writeLineLoop( 1, dg, Loop()>>1>>4>>-2>>-3);
            writePlaneSurface( 1, dg, 1);

            writePoint( 13, dg , a1(0), -a1(1) );
            writePoint( 14, dg , a2(0), -a2(1) );
            writePoint( 15, dg , a3(0), -a3(1) );
            writePoint( 16, dg , a4(0), -a4(1) );
            writePoint( 17, dg , a5(0), -a5(1) );
            writePoint( 18, dg , a6(0), -a6(1) );

            writeSpline( 5, dg, Loop()>>13>>14>>15>>16>>17>>18);

            writePoint( 19, dg , a1(0), -a1(1)-ep );
            writePoint( 20, dg , a2(0), -a2(1)-ep );
            writePoint( 21, dg , a3(0), -a3(1)-ep );
            writePoint( 22, dg , a4(0), -a4(1)-ep );
            writePoint( 23, dg , a5(0), -a5(1)-ep );
            writePoint( 24, dg , a6(0), -a6(1)-ep );

            writeSpline( 6, dg, Loop()>>19>>20>>21>>22>>23>>24);

            writeLine( 7, dg, 13,19);
            writeLine( 8, dg, 18,24);

            writeLineLoop( 2, dg, Loop()>>5>>8>>-6>>-7);
            writePlaneSurface( 2, dg, 2);

        }

        void
        runSpecial_1b(data_geo_ptrtype dg)
        {
            double yh=1.0;
            Node a1( 0.0, yh);
            Node a2( 3.0, yh);
            Node a3( 6.0, yh+0.7);
            Node a4( 7.3, yh-0.5);
            Node a5( 8.5, yh);
            Node a6(11.0, yh);
            double ep=0.3;
            //_______________________________________________//
            writePoint( 1, dg , a1(0), a1(1) );
            writePoint( 2, dg , a2(0), a2(1) );
            writePoint( 3, dg , a3(0), a3(1) );
            writePoint( 4, dg , a4(0), a4(1) );
            writePoint( 5, dg , a5(0), a5(1) );
            writePoint( 6, dg , a6(0), a6(1) );

            writeSpline( 1, dg, Loop()>>1>>2>>3>>4>>5>>6);

            writePoint( 7, dg , a1(0), -a1(1) );
            writePoint( 8, dg , a2(0), -a2(1) );
            writePoint( 9, dg , a3(0), -a3(1) );
            writePoint( 10, dg , a4(0), -a4(1) );
            writePoint( 11, dg , a5(0), -a5(1) );
            writePoint( 12, dg , a6(0), -a6(1) );

            writeSpline( 2, dg, Loop()>>7>>8>>9>>10>>11>>12);

            writeLine( 3, dg, 1,7);
            writeLine( 4, dg, 6,12);

            writeLineLoop( 1, dg, Loop()>>1>>4>>-2>>-3);

            writePlaneSurface( 1, dg, 1);

        }

        void
        runHexaedre(data_geo_ptrtype dg)
        {

            node_type Pt1 = param<0>(dg);
            node_type Pt2 = param<1>(dg);
            node_type Pt3 = param<2>(dg);
            node_type Pt4 = param<3>(dg);
            node_type Pt5 = param<4>(dg);
            node_type Pt6 = param<5>(dg);
            node_type Pt7 = param<6>(dg);
            node_type Pt8 = param<7>(dg);

            writePoint( 1, dg , Pt1(0), Pt1(1), Pt1(2) );
            writePoint( 2, dg , Pt2(0), Pt2(1), Pt2(2) );
            writePoint( 3, dg , Pt3(0), Pt3(1), Pt3(2) );
            writePoint( 4, dg , Pt4(0), Pt4(1), Pt4(2) );
            writePoint( 5, dg , Pt5(0), Pt5(1), Pt5(2) );
            writePoint( 6, dg , Pt6(0), Pt6(1), Pt6(2) );
            writePoint( 7, dg , Pt7(0), Pt7(1), Pt7(2) );
            writePoint( 8, dg , Pt8(0), Pt8(1), Pt8(2) );

            writeLine( 1, dg , 1 , 2);
            writeLine( 2, dg , 2 , 3);
            writeLine( 3, dg , 3 , 4);
            writeLine( 4, dg , 4 , 1);
            writeLine( 5, dg , 5 , 6);
            writeLine( 6, dg , 6 , 7);
            writeLine( 7, dg , 7 , 8);
            writeLine( 8, dg , 8 , 5);
            writeLine( 9, dg , 1 , 5);
            writeLine( 10, dg , 2 , 6);
            writeLine( 11, dg , 3 , 7);
            writeLine( 12, dg , 4 , 8);

            writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4);
            writePlaneSurface( 1, dg, 1);
            writeLineLoop( 2, dg, Loop()>>5>>6>>7>>8);
            writePlaneSurface( 2, dg, 2);
            writeLineLoop( 3, dg, Loop()>>1>>10>>-5>>-9);
            writePlaneSurface( 3, dg, 3);
            writeLineLoop( 4, dg, Loop()>>10>>6>>-11>>-2);
            writePlaneSurface( 4, dg, 4);
            writeLineLoop( 5, dg, Loop()>>11>>7>>-12>>-3);
            writePlaneSurface( 5, dg, 5);
            writeLineLoop( 6, dg, Loop()>>9>>-8>>-12>>4);
            writePlaneSurface( 6, dg, 6);

            writeSurfaceLoop( 1, dg, Loop()>>1>>2>>3>>4>>5>>6);

            writeVolume(1, dg, 1);
        }





    }//GeoTool

} //Life
#endif /* __geotool_H */
