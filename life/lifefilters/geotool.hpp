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


# define GEOTOOL_SHAPE					\
  ( 2, ( ( Rectangle, "rectangle", 2, RECTANGLE ),	\
       	 ( Circle   , "circle"   , 2, CIRCLE    ) )	\
    )							\
  /**/

# define GEOTOOL_SHAPE_NAME_CLASS(i) BOOST_PP_TUPLE_ELEM(4, 0, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NAME_STR(i) BOOST_PP_TUPLE_ELEM(4, 1, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NBPARAM(i) BOOST_PP_TUPLE_ELEM(4, 2, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NAME_MACRO(i) BOOST_PP_TUPLE_ELEM(4, 3, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))



# define GEOTOOL_MARKER_RECTANGLE		\
  ( 4, ( ( 1, 1, ( 1 ) ),			\
	 ( 2, 1, ( 2 ) ),			\
	 ( 3, 1, ( 3 ) ),			\
	 ( 4, 1, ( 4 ) ) )			\
    )						\
  /**/



# define GEOTOOL_MARKER_CIRCLE			\
  ( 1, ( ( 1, 2, ( 1,2 ) ) )			\
    )						\
  /**/


// Accessors for the operator datatype.


# define GEOTOOL_MARKER_INDICE(O) BOOST_PP_TUPLE_ELEM(3, 0, O)

# define GEOTOOL_MARKER_NBMARK(F,i) BOOST_PP_TUPLE_ELEM(3, 1, BOOST_PP_ARRAY_ELEM(i,F))

# define GEOTOOL_MARKER_ARRAYMARK(O) BOOST_PP_TUPLE_ELEM(3, 2, O)
# define GEOTOOL_MARKER_MARKVALUE(F,i,j) BOOST_PP_TUPLE_ELEM( GEOTOOL_MARKER_NBMARK(F,i),j,GEOTOOL_MARKER_ARRAYMARK(BOOST_PP_ARRAY_ELEM(i, F)))


namespace Life {

  BOOST_PARAMETER_NAME(markerName)
  BOOST_PARAMETER_NAME(marker1)
  BOOST_PARAMETER_NAME(marker2)
  BOOST_PARAMETER_NAME(marker3)
  BOOST_PARAMETER_NAME(marker4)


  
  namespace GeoTool {
    
    typedef node<double>::type node_type;

    class GeoGMSHTool;
    typedef boost::shared_ptr< GeoGMSHTool> GeoGMSHTool_ptrtype;

    typedef boost::tuple< GeoGMSHTool_ptrtype  ,std::string, std::string > data_geo_type;
    typedef boost::shared_ptr<data_geo_type> data_geo_ptrtype;

    typedef std::map<uint,uint> map_data_type;
    typedef std::vector<map_data_type> vec_map_data_type;

    void run(data_geo_ptrtype __dg, vec_map_data_type & __dataMem);

    void runRectangle(data_geo_ptrtype dg, vec_map_data_type & __dataMem);
    void runCircle(data_geo_ptrtype dg, vec_map_data_type & __dataMem);

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
	(*_M_node)(2)=__y;
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


      /*boost::shared_ptr<*/node_type
      getNode() const { return *_M_node;}

      boost::shared_ptr<node_type> _M_node;
    };

    /*_________________________________________________*/

    /*    template <uint Numero>
	  node_type
	  param(data_geo_ptrtype __dg)
	  {
	  std::string __shape = boost::get<1>(*__dg);
	  std::string __name = boost::get<2>(*__dg);
	  node_type __node = boost::get<Numero>(boost::get<0>(*__dg)->_M_paramShape->find(__shape)->second[__name]);
	  __node.resize(3);
	  if (__node.size()<2)
	  __node(1)=0.0;
	  if (__node.size()<3)
	  __node(2)=0.0;

	  return __node;
	  }*/

    /*_________________________________________________*/

    class Loop {
    public :
      
      Loop(Loop const & L) : _M_loop(L._M_loop) {}

      Loop() {_M_loop.clear();}

      void  operator=(Loop m) { this->_M_loop=m._M_loop; }
      Loop  operator>>(uint __n) { _M_loop.push_back(__n); return *this; }
      
      std::list<uint>::const_iterator begin() const { return _M_loop.begin();}
      std::list<uint>::const_iterator end() const { return _M_loop.end(); }

      std::list<uint> _M_loop;
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
      typedef std::list< boost::tuple<std::string,double> > names_type;
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



      GeoGMSHTool(std::string __shape="NO_SHAPE", std::string __name="NO_NAME", double __meshSize=0.1)
	:
	_M_cptPt(1),
	_M_cptLine(1),
	//_M_cptCircle(1),
	_M_cptLineLoop(1),
	_M_cptSurface(1),
	_M_map_Shape( new map_shape_names_type()),
	_M_paramShape( new parameter_shape_type()),
	_M_markShape( new marker_shape_type()),
	_M_ostr( new std::ostringstream())
      {
      }

      GeoGMSHTool( GeoGMSHTool const & m ) 
	:
	_M_cptPt(m._M_cptPt),
	_M_cptLine(m._M_cptLine),
	//_M_cptCircle(m._M_cptCircle),
	_M_cptLineLoop(m._M_cptLineLoop),
	_M_cptSurface(m._M_cptSurface),
	_M_surfaceList(m._M_surfaceList),
	_M_map_Shape(new map_shape_names_type(*(m._M_map_Shape))),
	_M_paramShape(new parameter_shape_type(*(m._M_paramShape))),
	_M_markShape(new marker_shape_type(*(m._M_markShape))),
	_M_ostr(new std::ostringstream())
      {
	updateOstr((m._M_ostr)->str());
      }

      void zeroCpt() 
      { 
	_M_cptPt=1;
	_M_cptLine=1;
	//_M_cptCircle=1;
	_M_cptLineLoop=1;
	_M_cptSurface=1;
      }

      void
      operator=( GeoGMSHTool const & m ) 
      {
	_M_cptPt = m._M_cptPt;
	_M_cptLine = m._M_cptLine;
	//_M_cptCircle = m._M_cptCircle;
	_M_cptLineLoop = m._M_cptLineLoop;
	_M_cptSurface = m._M_cptSurface;
	_M_surfaceList = m._M_surfaceList;
	_M_map_Shape.reset(new map_shape_names_type(*(m._M_map_Shape)));
	_M_paramShape.reset(new parameter_shape_type(*(m._M_paramShape)));
	_M_markShape.reset(new marker_shape_type(*(m._M_markShape)));
	_M_ostr.reset(new std::ostringstream());
	updateOstr((m._M_ostr)->str());
      }

      GeoGMSHTool operator+(const GeoGMSHTool & m);
      GeoGMSHTool operator-(const GeoGMSHTool & m);

      GeoGMSHTool opFusion(const GeoGMSHTool & m,int __typeop);

      void init();

      void initData(std::string __shape, 
		    std::string __name, 
		    double __meshSize, 
		    std::vector<GeoTool::Node> & __param) 
      {
	boost::tuple<std::string,double> __id = boost::make_tuple(__name, __meshSize);

	(*(_M_map_Shape))[__shape].push_back(__id);
	(*(_M_markShape))[__shape][__name];
	//(*(_M_paramShape))[__shape][__name] = boost::make_tuple(*__param1.getNode(),*__param2.getNode());
	(*(_M_paramShape))[__shape][__name].resize(__param.size());
	(*(_M_paramShape))[__shape][__name][0] = __param[0].getNode();
	(*(_M_paramShape))[__shape][__name][1] = __param[1].getNode();

	//Attention 0 par defaut
	std::list< boost::tuple<std::string,std::string, uint >	>__listTemp;
	__listTemp.push_back( boost::make_tuple(__shape,__name,0));
	_M_surfaceList.push_back( __listTemp);
      }
      
      void updateData(GeoGMSHTool const & m) {
	_M_cptPt = m._M_cptPt;
	_M_cptLine = m._M_cptLine;
	//_M_cptCircle = m._M_cptCircle;
	_M_cptLineLoop = m._M_cptLineLoop;
	_M_cptSurface = m._M_cptSurface;
	
	_M_map_Shape = m._M_map_Shape;
	_M_paramShape = m._M_paramShape;
	_M_markShape = m._M_markShape;

	//Est-ce à faire?
	_M_surfaceList = m._M_surfaceList;

	//cleanOstr();
	//updateOstr((m._M_ostr)->str());
	//_M_ostr.reset(new std::ostringstream());
	// *_M_ostr<<(m._M_ostr)->str();
      }


      void updateOstr( std::string __str)
      {
	*_M_ostr << __str;
      }

      std::string geoStr();

      void cleanOstr() { _M_ostr.reset(new std::ostringstream()); }

      //void runRectangle(node_type xbg, node_type xhd);


      /*_________________________________________________*
       *_________________________________________________*
       * Accessor                                        *
       *_________________________________________________*
       *_________________________________________________*/


      double cptPt() const { return _M_cptPt;}
      double cptLine() const { return _M_cptLine;}
      //double cptCircle() const { return _M_cptCircle;}
      double cptLineLoop() const { return _M_cptLineLoop;}
      double cptSurface() const { return _M_cptSurface;}

  
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

      marker_shape_const_iterator_type
      markRectangleBegin() const
      { 
	return _M_markShape->begin();
      }
 
      marker_shape_const_iterator_type 
      markRectangleEnd() const 
      { 
	return _M_markShape->end();
      }

      /*      marker_name_const_iterator_type
	      markNameBegin(std::string __shape) const
	      { 
	      return _M_markShape->find(__shape)->second.begin();
	      }*/

      marker_type_const_iterator_type
      markerTypeBegin(std::string __shape) const 
      { 
	return _M_markShape->find(__shape)->second.begin();
      }

      marker_type_const_iterator_type 
      markerTypeEnd(std::string __shape) const 
      { 
	return _M_markShape->find(__shape)->second.end();
      }

      marker_type_type 
      markerType(std::string __shape) const 
      { 
	return _M_markShape->find(__shape)->second;
      }

      marker_markerName_const_iterator_type
      markerMarkerNameBegin(std::string __shape, std::string __type) const 
      { 
	return _M_markShape->find(__shape)->second.find(__type)->second.begin();
      }

      marker_markerName_const_iterator_type 
      markerMarkerNameEnd(std::string __shape, std::string __type) const 
      { 
	return _M_markShape->find(__shape)->second.find(__type)->second.end();
      }

      marker_markerName_type 
      markerMarkerName(std::string __shape, std::string __type) const 
      { 
	return _M_markShape->find(__shape)->second.find(__type)->second;
      }


      std::list<marker_base_type>::const_iterator
      markerListIndiceBegin(std::string __shape, std::string __type ,std::string __markerName) const
      { 
	return _M_markShape->find(__shape)->second.find(__type)->second.find(__markerName)->second.begin();
      }

      std::list<marker_base_type>::const_iterator
      markerListIndiceEnd(std::string __shape, std::string __type ,std::string __markerName) const
      { 
	return _M_markShape->find(__shape)->second.find(__type)->second.find(__markerName)->second.end();
      }


      std::list<marker_base_type>
      getMarkerName(std::string __shape, std::string __type ,std::string __markerName) const
      { 
	return _M_markShape->find(__shape)->second.find(__type)->second.find(__markerName)->second;
      }

  
      /*_________________________________________________*
       *_________________________________________________*
       * Members                                         *
       *_________________________________________________*
       *_________________________________________________*/


      // memory 
      uint _M_cptPt;
      uint _M_cptLine;
      //uint _M_cptCircle;
      uint _M_cptLineLoop;
      uint _M_cptSurface;

      // gestion des surface : shape,name,value
      std::list< std::list< boost::tuple<std::string,std::string, uint > > > _M_surfaceList;

      // data containers
      boost::shared_ptr<map_shape_names_type> _M_map_Shape;
      boost::shared_ptr<parameter_shape_type> _M_paramShape;
      boost::shared_ptr<marker_shape_type> _M_markShape;

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
    
      __geoTool._M_cptPt = this->_M_cptPt + m.cptPt()-1;
      __geoTool._M_cptLine = this->_M_cptLine + m.cptLine()-1;
      //__geoTool._M_cptCircle = this->_M_cptCircle + m.cptCircle()-1;
      __geoTool._M_cptLineLoop = this->_M_cptLineLoop + m.cptLineLoop()-1;
      __geoTool._M_cptSurface = this->_M_cptSurface + m.cptSurface()-1;


      //NEW!!!!!!!!!
      if (__typeOp==1)
	{
	  __geoTool._M_surfaceList = this->_M_surfaceList;
	  std::list<std::list< boost::tuple<std::string,std::string,uint> > >::const_iterator itSurf = m._M_surfaceList.begin();
	  std::list<std::list< boost::tuple<std::string,std::string,uint> > >::const_iterator itSurf_end = m._M_surfaceList.end();
	  for ( ; itSurf != itSurf_end; ++itSurf) 
	    { 
	      std::list< boost::tuple<std::string,std::string,uint> >::const_iterator itSurf2 = itSurf->begin();
	      std::list< boost::tuple<std::string,std::string,uint> >::const_iterator itSurf2_end = itSurf->end();
	      std::list< boost::tuple<std::string,std::string,uint> > __listTemp;

	      for ( ; itSurf2 != itSurf2_end; ++itSurf2) 
		{ 
		  __listTemp.push_back(*itSurf2);
		}	
	      //__geoTool._M_surfaceList.push_back(*itSurf);  
	      __geoTool._M_surfaceList.push_back(__listTemp);  
	    }
	}
      else
	{
	  //NEW!!!!!!!!!
	  __geoTool._M_surfaceList = this->_M_surfaceList;
	  std::list<std::list< boost::tuple<std::string,std::string,uint> > >::const_iterator itSurf = m._M_surfaceList.begin();
	  std::list<std::list< boost::tuple<std::string,std::string,uint> > >::const_iterator itSurf_end = m._M_surfaceList.end();
	  for ( ; itSurf != itSurf_end; ++itSurf) 
	    { 
	      std::list< boost::tuple<std::string,std::string,uint> >::const_iterator itSurf2 = itSurf->begin();
	      std::list< boost::tuple<std::string,std::string,uint> >::const_iterator itSurf2_end = itSurf->end();
	      std::list< boost::tuple<std::string,std::string,uint> > __listTemp;

	      for ( ; itSurf2 != itSurf2_end; ++itSurf2) 
		{ 
		  __geoTool._M_surfaceList.begin()->push_back(*itSurf2);  

		  //__listTemp.push_back(*itSurf2);
		  //itSu.push_back(*itSurf);
		}	
	    }
	}

      __geoTool._M_ostr.reset(new std::ostringstream());
      __geoTool.updateOstr(this->_M_ostr->str());
      __geoTool.updateOstr(m._M_ostr->str());
     
      //get data from this (easy)
      __geoTool._M_map_Shape.reset( new map_shape_names_type( *this->_M_map_Shape) );
      __geoTool._M_paramShape.reset( new parameter_shape_type( *this->_M_paramShape) );
      __geoTool._M_markShape.reset( new marker_shape_type (*this->_M_markShape) );

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
	      std::list<boost::tuple<std::string,double> >::const_iterator itLs = (*(__geoTool._M_map_Shape))[itShape->first].begin();
	      std::list<boost::tuple<std::string,double> >::const_iterator itLs_end = (*(__geoTool._M_map_Shape))[itShape->first].end();
	      bool __find=false;
	      boost::tuple<std::string,double> __temp;
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

      map_shape_names_const_iterator_type itShape2 = m._M_map_Shape->begin();
      map_shape_names_const_iterator_type itShape2_end = m._M_map_Shape->end();
      while (itShape2!=itShape2_end)
	{
	  marker_type_const_iterator_type itMarkType = m.markerTypeBegin(itShape2->first);
	  marker_type_const_iterator_type itMarkType_end = m.markerTypeEnd(itShape2->first);
	  while (itMarkType!=itMarkType_end)
	    {
	      marker_markerName_const_iterator_type itMarkName = m.markerMarkerNameBegin(itShape2->first,itMarkType->first);
	      marker_markerName_const_iterator_type itMarkName_end = m.markerMarkerNameEnd(itShape2->first,itMarkType->first);
	      while (itMarkName!=itMarkName_end)
		{
		  if ( !m.getMarkerName(itShape2->first,itMarkType->first,itMarkName->first).empty() ) 
		    {
		      std::list<marker_base_type>::const_iterator itLRef = m.markerListIndiceBegin(itShape2->first,itMarkType->first,itMarkName->first);
		      std::list<marker_base_type>::const_iterator itLRef_end = m.markerListIndiceEnd(itShape2->first,itMarkType->first,itMarkName->first);
		      while(itLRef != itLRef_end ) {
			(*(__geoTool._M_markShape))[itShape2->first][itMarkType->first][itMarkName->first].push_back(*itLRef); 
			++itLRef;
		      }
		    }
		  ++itMarkName;    
		}
	      ++itMarkType;
	    }
	  ++itShape2;
	}
      
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

      // gmsh : part init
      init();

      //data memory ( type->shape->name )
      std::vector<std::map<std::string,std::map<std::string, std::map<uint,uint> > > > __dataMemGlob(5);

      //iterate on the shape
      map_shape_names_const_iterator_type itShape =  _M_map_Shape->begin();
      map_shape_names_const_iterator_type itShape_end = _M_map_Shape->end();
      while (itShape!=itShape_end)
	{
	  //iterate on the name of shape
	  std::list<boost::tuple<std::string,double> >::const_iterator itName = itShape->second.begin();
	  std::list<boost::tuple<std::string,double> >::const_iterator itName_end = itShape->second.end();
	  while (itName!=itName_end)
	    {    

	      *_M_ostr << "h=" << boost::get<1>(*itName) << ";\n";
	      //data memory
	      vec_map_data_type __dataMem(8);
      
	      GeoGMSHTool_ptrtype __geoTool(new GeoGMSHTool());     
	      __geoTool->updateData(*this);
	      __geoTool->cleanOstr();
	      GeoTool::data_geo_ptrtype __data_geoTool(new GeoTool::data_geo_type(boost::make_tuple( __geoTool,
												     itShape->first,
												     boost::get<0>(*itName)
												     )));
	      // generate the code for the geometry
	      run(__data_geoTool,__dataMem);

	      __dataMemGlob[0][itShape->first][boost::get<0>(*itName)] = __dataMem[0];//pts
	      __dataMemGlob[1][itShape->first][boost::get<0>(*itName)] = __dataMem[1];//lines
	      __dataMemGlob[2][itShape->first][boost::get<0>(*itName)] = __dataMem[2];//lineLoop
	      __dataMemGlob[3][itShape->first][boost::get<0>(*itName)] = __dataMem[3];//Surface

	      // get infos
	      this->updateData( *boost::get<0>(*__data_geoTool));
	      this->updateOstr( boost::get<0>(*__data_geoTool)->_M_ostr->str());

	      ++itName;
	    }
	  ++itShape;
	}


      //NEW!!!!!!!!!
      //Ecriture des planes surfaces
      //Fait ici a cause des opertateurs (+,-)
      std::list<std::list< boost::tuple<std::string,std::string,uint> > >::const_iterator itSurf = this->_M_surfaceList.begin();
      std::list<std::list< boost::tuple<std::string,std::string,uint> > >::const_iterator itSurf_end = this->_M_surfaceList.end();
      std::ostringstream __surface_str;

      for ( ; itSurf != itSurf_end; ++itSurf)
	{
	  std::list< boost::tuple<std::string,std::string,uint> >::const_iterator itSurf2 = itSurf->begin();
	  std::list< boost::tuple<std::string,std::string,uint> >::const_iterator itSurf2_end = --itSurf->end();

	  //on suppose une seule surface par objet <=> numero local=1
	  __dataMemGlob[3][boost::get<0>(*itSurf2)][boost::get<1>(*itSurf2)][1]=this->cptSurface();
	  __surface_str << "Plane Surface(" << this->cptSurface() << ") = {" ;

	  for ( ; itSurf2 != itSurf2_end; ++itSurf2)
	    {
	      __surface_str << boost::get<2>(*itSurf2) << ",";
	    }	  
	  __surface_str << boost::get<2>(*itSurf2);

	  __surface_str << "};\n";
	  ++this->_M_cptSurface;

	}


      this->updateOstr(__surface_str.str());

      //iterate on the shape
      itShape =  _M_map_Shape->begin();
      itShape_end = _M_map_Shape->end();
      while (itShape!=itShape_end)
	{
	  // generate the code for the marker
	  marker_type_const_iterator_type itMarkType= (*(_M_markShape))[itShape->first].begin();
	  marker_type_const_iterator_type itMarkType_end=(*(_M_markShape))[itShape->first].end();
	  while (itMarkType!=itMarkType_end)
	    {
	      marker_markerName_const_iterator_type itMarkName = (*(_M_markShape))[itShape->first][itMarkType->first].begin();
	      marker_markerName_const_iterator_type itMarkName_end=(*(_M_markShape))[itShape->first][itMarkType->first].end();
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
	  ++itShape;
	}

     
      return _M_ostr->str();

    }



    /*_________________________________________________*
     *_________________________________________________*
     *_________________________________________________*
     * Function on the namespace                       *
     *_________________________________________________*
     *_________________________________________________*
     *_________________________________________________*/



    void run( data_geo_ptrtype __dg, vec_map_data_type & __dataMem)
    {

      if( boost::get<1>(*__dg) =="rectangle")
	{
	  runRectangle(__dg,__dataMem);
	}
      if( boost::get<1>(*__dg) =="circle")
	{
	  runCircle(__dg,__dataMem);
	}

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
      std::string __shape = boost::get<1>(*__dg);
      std::string __name = boost::get<2>(*__dg);
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


    void writePoint(uint __numLoc,       vec_map_data_type & __mapPt, data_geo_ptrtype __dg ,double __x1,double __x2, double __x3=0)
    {
      __mapPt[0][__numLoc] = boost::get<0>(*__dg)->cptPt();
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
    writeLine(uint __numLoc, vec_map_data_type & __mapLine, data_geo_ptrtype __dg ,uint __n1, uint __n2)
    {
      __mapLine[1][__numLoc] = boost::get<0>(*__dg)->cptLine();
      std::ostringstream __ostr;
      __ostr << "Line(" << boost::get<0>(*__dg)->cptLine()
	     << ") = {" 
	     << __mapLine[0][__n1] << "," 
	     << __mapLine[0][__n2] << "};\n";

      boost::get<0>(*__dg)->updateOstr(__ostr.str());

      ++boost::get<0>(*__dg)->_M_cptLine;
    }

    /*_________________________________________________*/

    void
    writeCircle(uint __numLoc, vec_map_data_type & __mapLine, data_geo_ptrtype __dg ,uint __n1, uint __n2, uint __n3)
    {

      //      __mapLine[4][__numLoc] = boost::get<0>(*__dg)->cptCircle();
      __mapLine[1][__numLoc] = boost::get<0>(*__dg)->cptLine();

      std::ostringstream __ostr;
      __ostr << "Circle(" << boost::get<0>(*__dg)->cptLine() //cptCircle
	     << ") = {" 
	     << __mapLine[0][__n1] << "," 
	     << __mapLine[0][__n2] << "," 
	     << __mapLine[0][__n3] << "};\n";

      boost::get<0>(*__dg)->updateOstr(__ostr.str());

      //++boost::get<0>(*__dg)->_M_cptCircle;
      ++boost::get<0>(*__dg)->_M_cptLine;
    }

    /*_________________________________________________*/

    void
    writeLineLoop(uint __numLoc, vec_map_data_type & __mapLine, data_geo_ptrtype __dg , Loop /*const*/ __loop, uint __typeBaseShape=1 )
    {
      __mapLine[2][__numLoc] = boost::get<0>(*__dg)->cptLineLoop();
      std::ostringstream __ostr;
      __ostr << "Line Loop(" << boost::get<0>(*__dg)->cptLineLoop()
	     << ") = {" ;
      std::list<uint>::const_iterator it= __loop.begin();
      std::list<uint>::const_iterator it_end= --__loop.end();
      while (it!=it_end)
	{
	  __ostr << __mapLine[__typeBaseShape][*it] <<"," ;
	  ++it;
	}

      __ostr << __mapLine[__typeBaseShape][*it] << "};\n";

      boost::get<0>(*__dg)->updateOstr(__ostr.str());

      ++boost::get<0>(*__dg)->_M_cptLineLoop;
    }

    /*_________________________________________________*/


    //ici on n'ecrit pas, on memorise cause des operations de difference
    //l'ecriture est realise dans geoStr()
    void
    writePlaneSurface(uint __numLoc, vec_map_data_type & __mapLine, data_geo_ptrtype __dg , uint __ind)
    {
      __mapLine[3][__numLoc] = boost::get<0>(*__dg)->cptSurface();
      //std::ostringstream __ostr;
      //__ostr << "Plane Surface(" << boost::get<0>(*__dg)->cptSurface()
      //<< ") = {" ;

      /*      std::list<uint>::const_iterator it= __loop.begin();
	      std::list<uint>::const_iterator it_end= --__loop.end();
	      while (it!=it_end)
	      {
	      __ostr << __mapLine[1][*it] <<"," ;
	      ++it;
	      }*/

      //__ostr << __mapLine[2][__ind] << "};\n";
      //boost::get<0>(*__dg)->updateOstr(__ostr.str());

      std::list<std::list<boost::tuple<std::string,std::string,uint> > >::iterator itSurf = boost::get<0>(*__dg)->_M_surfaceList.begin();
      std::list<std::list<boost::tuple<std::string,std::string,uint> > >::iterator itSurf_end = boost::get<0>(*__dg)->_M_surfaceList.end();
      for( ; itSurf !=itSurf_end;++itSurf)
	{
	  std::list<boost::tuple<std::string,std::string,uint> >::iterator itSurf2 = itSurf->begin();
	  std::list<boost::tuple<std::string,std::string,uint> >::iterator itSurf2_end = itSurf->end();
	  while (itSurf2 !=itSurf2_end) 
	    {
	      if (boost::get<0>(*itSurf2) == boost::get<1>(*__dg))
		{
		  if (boost::get<1>(*itSurf2) == boost::get<2>(*__dg))
		    {
		      boost::get<2>(*itSurf2) = __mapLine[2][__ind];
		    }
		}
	      ++itSurf2;
	    }
	}

      //      ++boost::get<0>(*__dg)->_M_cptSurface;
    }



    /*_________________________________________________*
     *_________________________________________________*
     *_________________________________________________*
     * PREPROCESSOR METHODS                            *
     *_________________________________________________*
     *_________________________________________________*
     *_________________________________________________*/




#define GEOTOOL_FOR_COMP2(r, state)					\
    BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(4, 0, state),		\
			BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(4, 1, state))	\
			)						\
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/
#define GEOTOOL_FOR_INCR2(r, state)				\
    (								\
     BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(4, 0, state)),		\
     BOOST_PP_TUPLE_ELEM(4, 1, state),				\
     BOOST_PP_TUPLE_ELEM(4, 2, state),				\
     BOOST_PP_TUPLE_ELEM(4, 3, state) )				\
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/
#define GEOTOOL_FOR_MARKER_MACRO2(r, state)				\
    __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
					     GEOTOOL_MARKER_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_, \
										    GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) )), \
								       BOOST_PP_TUPLE_ELEM(4, 2, state), \
								       BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
					     )				\
			   );						\
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/
#define GEOTOOL_FOR_MARKER_MACRO(r, state)				\
    if (BOOST_PP_CAT(marker,						\
		     BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
				  1 ) ) )				\
      {									\
	BOOST_PP_FOR( (0,						\
		       BOOST_PP_SUB(GEOTOOL_MARKER_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_,\
								       GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state))), \
							  BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
		       BOOST_PP_TUPLE_ELEM(3, 0, state),		\
		       BOOST_PP_TUPLE_ELEM(3, 2, state)			\
		       ),						\
		      GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_MACRO2) \
	  }								\
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/
#define GEOTOOL_FOR_COMP1(r, state)					\
    BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(3, 0, state),		\
			BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(3, 1, state))	\
			)						\
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/
#define GEOTOOL_FOR_INCR1(r, state)			\
    (							\
     BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(3, 0, state)),	\
     BOOST_PP_TUPLE_ELEM(3, 1, state),			\
     BOOST_PP_TUPLE_ELEM(3, 2, state) )			\
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/    
#define GEOTOOL_FOR_COMP(r, state)					\
    BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(2, 0, state),		\
			BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 1, state))	\
			)						\
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/
#define GEOTOOL_FOR_INCR(r, state)			\
    (							\
     BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 0, state)),	\
     BOOST_PP_TUPLE_ELEM(2, 1, state) )			\
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/
#define GEOTOOL_SHAPE_PARAM(r, state)					\
    _M_param[BOOST_PP_TUPLE_ELEM(2,0,state)] = BOOST_PP_CAT( __param,	\
							     BOOST_PP_TUPLE_ELEM(2,0,state) ); \
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/
#define GEOTOOL_SHAPE_CLASS(r,state)					\
    class GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state)) : public GeoGMSHTool \
    {									\
    public :								\
      									\
      typedef GeoGMSHTool::node_type node_type;				\
      typedef GeoTool::Node Node;					\
      									\
      									\
      GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))(double __meshSize, \
							       std::string __name, \
							       Node __param0, \
							       Node __param1 ) \
	:								\
	GeoGMSHTool( shape(), __name, __meshSize),			\
	_M_name(__name)							\
	{								\
	  _M_param.resize( GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state))); \
	  BOOST_PP_FOR( (0, BOOST_PP_SUB(GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state)),1) ), \
			GEOTOOL_FOR_COMP,				\
			GEOTOOL_FOR_INCR,				\
			GEOTOOL_SHAPE_PARAM);				\
	  								\
	  initData(shape(),__name, __meshSize, _M_param);		\
	}								\
									\
									\
									\
									\
      GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))(const GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state)) & m) \
	:								\
	GeoGMSHTool(m),							\
	_M_param(m._M_param)						\
	  {}								\
									\
      BOOST_PARAMETER_MEMBER_FUNCTION(					\
				      (void),				\
				      setMarker,			\
				      tag,				\
				      (required				\
				       ( type, (std::string))		\
				       ( name, (std::string)) )		\
				      (optional				\
				       (marker1, (bool), false)		\
				       (marker2, (bool), false)		\
				       (marker3, (bool), false)		\
				       (marker4, (bool), false)))	\
	{								\
									\
	  std::list<marker_base_type > __listMarker = (*(_M_markShape))[this->shape()][type][name]; \
  									\
									\
	  BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE( BOOST_PP_CAT(GEOTOOL_MARKER_, \
									   GEOTOOL_SHAPE_NAME_MACRO( BOOST_PP_TUPLE_ELEM(2,0,state)))), \
					 1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
			GEOTOOL_FOR_COMP1,				\
			GEOTOOL_FOR_INCR1,				\
			GEOTOOL_FOR_MARKER_MACRO)			\
									\
	    (*(_M_markShape))[this->shape()][type][name] = __listMarker; \
	}								\
									\
									\
      std::string _M_name;						\
									\
      static const std::string shape() { return GEOTOOL_SHAPE_NAME_STR(BOOST_PP_TUPLE_ELEM(2,0, state));} \
      const std::string name() const {return _M_name;}			\
      									\
									\
      std::vector<GeoTool::Node> _M_param;				\
									\
    };									\
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
    runRectangle(data_geo_ptrtype dg, vec_map_data_type & __dataMem)
    {

      node_type PtA = param<0>(dg);
      node_type PtB = param<1>(dg);

      writePoint( 1, __dataMem , dg , PtA(0), PtA(1) );
      writePoint( 2, __dataMem , dg , PtB(0), PtA(1) );
      writePoint( 3, __dataMem , dg , PtB(0), PtB(1) );
      writePoint( 4, __dataMem , dg , PtA(0), PtB(1) );

      writeLine( 1, __dataMem , dg , 1 , 2);
      writeLine( 2, __dataMem , dg , 2 , 3);
      writeLine( 3, __dataMem , dg , 3 , 4);
      writeLine( 4, __dataMem , dg , 4 , 1);

      writeLineLoop( 1, __dataMem , dg, Loop()>>1>>2>>3>>4);

      writePlaneSurface( 1, __dataMem , dg, 1);
      
    }


    void
    runCircle(data_geo_ptrtype dg, vec_map_data_type & __dataMem)
    {

      node_type PtA = param<0>(dg);
      node_type PtB = param<1>(dg);

      writePoint( 1, __dataMem , dg , PtB(0), PtB(1) );
      writePoint( 2, __dataMem , dg , PtA(0), PtA(1) );
      writePoint( 3, __dataMem , dg , 2*PtA(0)-PtB(0), 2*PtA(1)-PtB(1) );

      writeCircle( 1,__dataMem , dg, 1, 2, 3);
      writeCircle( 2,__dataMem , dg, 3, 2, 1);

      //Artung le 4 c'est pour identifier la loop faite d arc2cercle
      //sinon 1 par defaut pour Line -> A simplifier 
      writeLineLoop( 1, __dataMem , dg, Loop()>>1>>2);

      writePlaneSurface( 1, __dataMem , dg, 1);
      
    }





  }//GeoTool

} //Life
