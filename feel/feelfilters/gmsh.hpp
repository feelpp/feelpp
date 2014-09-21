/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4 tw=0

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-10

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2014 Feel++ Consortium
  

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
   \file gmsh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-02-10
 */
#ifndef __gmsh_H
#define __gmsh_H 1

#include <boost/type_traits.hpp>

#include <boost/filesystem.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/icl/type_traits/is_map.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelgmsh.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>
#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelfilters/gmshenums.hpp>
#include <feel/feelfilters/periodicentities.hpp>


#if defined( FEELPP_HAS_GMSH_H )
class GModel;
#endif

// there is a macro called sign in Gmsh that conflicts with
// at least one member function sign() from DofTable.
// hence we undefine the macro sign after including Gmsh headers
#undef sign

namespace Feel
{
extern const char* FEELPP_GMSH_FORMAT_VERSION;
}

namespace Feel
{
/**
 * \class Gmsh
 * \brief Gmsh Mesh Generator
 *
 * The \c Gmsh class helps with the generation of meshes using the
 * program \c gmsh. Typically one would generate a \c .geo std::string
 * for example and pass it to \c generate() along with the prefix \c
 * name of the \c .geo and \c .msh files to create/generate.
 *
 * \code
 * std::ostringstream ostr;
 * std::string fname;
 * ostr << "h=" << 0.1 << ";\n"
 * << "Point(1) = {0,0,0,h};\n"
 * << "Point(2) = {1,0,0,h};\n"
 * << "Point(3) = {1,1,0,h};\n"
 * << "Point(4) = {0,1,0,h};\n"
 * << "Line(1) = {1,2};\n"
 * << "Line(2) = {2,3};\n"
 * << "Line(3) = {3,4};\n"
 * << "Line(4) = {4,1};\n"
 * << "Line Loop(4) = {1,2,3,4};\n"
 * << "Plane Surface(5) = {4};\n"
 * << "Physical Line(10) = {1,3};\n"
 * << "Physical Line(20) = {2,4};\n"
 * << "Physical Surface(6) = {5};\n";
 * Gmsh gmsh;
 * // generate the \c .geo file using \c ostr and
 * // call gmsh to generate the \c msh file from the \c .geo file
 * fname = gmsh.generate( "triangle", ostr.str() );
 * // fname can now be passed to \c ImporterGmsh for example
 * mesh_type mesh;
 * ImporterGmsh<mesh_type> import( fname );
 * mesh->accept( import );
 * \endcode
 * \ingroup Importer
 * \author Christophe Prud'homme
 * \see http://www.geuz.org/gmsh/
 */
class Gmsh
{
public:

    enum DomainType { GMSH_REFERENCE_DOMAIN = 0, GMSH_REAL_DOMAIN };

    struct Factory
    {
        typedef Feel::Singleton< Feel::Factory< Gmsh, std::string > > type;
    };

    /** \name Constructors, destructor
     */
    //@{

    Gmsh( int nDim = 1, int nOrder = GMSH_ORDER_ONE, WorldComm const& worldComm=Environment::worldComm() );
    Gmsh( Gmsh const & __g );
    virtual ~Gmsh();

    //@}

    /** \name Operator overloads
     */
    //@{

    /**
     * assignment operator
     *
     * \param __g another Gmsh object instance
     *
     * \return the newly generated Gmsh object
     */
    Gmsh& operator=( Gmsh const& __g );

    static boost::shared_ptr<Gmsh> New( po::variables_map const& vm );
    static boost::shared_ptr<Gmsh> New( std::string const& shape, uint16_type d = 2,
                                        uint16_type o = 1, std::string const& ct = "simplex",
                                        WorldComm const& worldComm = Environment::worldComm() );

    //@}

    /** \name Accessors
     */
    //@{

    /**
     * \return mesh dimension
     */
    int dimension() const
        {
            return M_dimension;
        }

    /**
     * get the order of the elements of the mesh
     * \return the order of the elements of the mesh
     */
    GMSH_ORDER order() const
        {
            return ( GMSH_ORDER ) M_order;
        }

    /**
     * @return the file format version
     */
    std::string version() const
        {
            return M_version;
        }

    /**
     * @return file format
     */
    GMSH_FORMAT format() const
        {
            return M_format;
        }

    /**
     * @brief get the geometry
     * @return the gmsh geo description
     */
    std::pair<std::string,std::string> geo() const { return M_geo; }

    /**
     * @return true if gmsh format is ascii
     */
    bool isASCIIFormat() const { return M_format == GMSH_FORMAT_ASCII; }

    /**
     * @return true if gmsh format is binary
     */
    bool isBinaryFormat() const { return M_format == GMSH_FORMAT_BINARY; }

    /**
     * \return the name of the file
     */
    std::string prefix() const
        {
            return M_name;
        }

    PeriodicEntities const& periodic() const { return M_periodic; }

    /**
     * \return bounding box
     */
    std::vector<std::pair<double,double> > const& boundingBox() const
        {
            return M_I;
        }
    double xmin() const
        {
            return M_I[0].first;
        }
    double xmax() const
        {
            return M_I[0].second;
        }
    double ymin() const
        {
            return M_I[1].first;
        }
    double ymax() const
        {
            return M_I[1].second;
        }
    double zmin() const
        {
            return M_I[2].first;
        }
    double zmax() const
        {
            return M_I[2].second;
        }

    /**
     * \return characteristic length
     */
    double const& h() const
        {
            return M_h;
        }

    /**
     * \return the geometry description
     */
    std::string description() const
        {
            std::string d = this->getDescription();

            if ( !d.empty() ) M_desc = d;

            return M_desc;
        }

    //! \brief Get the value of a GMSH geometry parameter.
    //! If the parameter does not match any parameter, the function throws
    //! an out_of_range exception.
    //!     \param _name Geo parameter name.
    //! \return Return the geo parameter value.
    double geoParameter( std::string const& _name )
        {
            return boost::lexical_cast<double>( M_geoParamMap.at( _name ) );
        }

    //! \brief Get all GMSH geometry parameters.
    //! \return Return a map containing the geo gmsh geometry parameters as {par,value}.
    std::map<std::string, std::string> geoParameters()
        {
            return M_geoParamMap;
        }

    /**
     * add the mid point of the domain
     */
    bool addMidPoint() const
        {
            return M_addmidpoint;
        }

    /**
     * \return true if use the physical name, 'false' otherwise
     */
    bool usePhysicalNames() const
        {
            return M_usePhysicalNames;
        }

    //! \return the world comm
    WorldComm const& worldComm() const
        {
            return M_worldComm;
        }

    //! \return the nnumber of partitions
    int numberOfPartitions() const
        {
            return M_partitions;
        }

    //! \return true if save msh file by partitions, false otherwise
    bool mshFileByPartition() const
        {
            return M_partition_file;
        }

    //! \return the partitioner
    GMSH_PARTITIONER partitioner() const
        {
            return M_partitioner;
        }

    //! get the shear
    double shear() const
        {
            return M_shear;
        }

    //! return true if recombine, false otherwise
    bool recombine() const
        {
            return M_recombine;
        }
    int structuredMesh() const
        {
            return M_structured;
        }
    int refinementLevels() const
        {
            return M_refine_levels;
        }

    /**
     * @brief get the Gmsh GModel data structure
     * @return the Gmsh GModel data structure
     */
    GModel* gModel() const { return M_gmodel; }

    //@}

    /** \name  Mutators
     */
    //@{

    /**
     * set the dimension
     */
    Gmsh& setDimension( int dim )
        {
            M_dimension = dim;
            return *this;
        }

    /**
     * the gmsh generator to generate a reference domain
     * \return the mesh generator
     *
     * \code
     * Gmsh gmsh;
     * gmsh = gmsh.ref();
     * \endcode
     */
    Gmsh& ref()
        {
            this->setReferenceDomain();
            return *this;
        }

    /**
     * set the characteristic length
     * \param h the characteristic length
     * \return the mesh generator
     */
    Gmsh& h( double _h )
        {
            this->setCharacteristicLength( _h );
            return *this;
        }

    /**
     * set the order of the elements of the mesh it can be either
     * GMSH_ORDER_ONE (order 1/linear) or GMSH_ORDER_TWO(order
     * 2/quadratic)
     *
     * \param o order of the elements
     */
    void setOrder( int o )
        {
            M_order = ( GMSH_ORDER ) o;
            M_geoParamMap["ElementOrder"]=boost::lexical_cast<std::string>(o);
        }

    /**
     * set the file format \p version in ascii or binary \p format
     */
    void setVersion( std::string version, GMSH_FORMAT format = GMSH_FORMAT_ASCII )
        {
            if ( version != "1" && version != "2" && version != FEELPP_GMSH_FORMAT_VERSION )
                throw std::invalid_argument( "invalid gmsh file format version" );

            M_version = version;
            M_format = format;
        }

    /**
     * if \p in is true, read in-memory geo files
     */
    void setInMemory( bool in )
        {
            M_in_memory = in;
        }

    /**
     * set file \p format: ascii or binary
     */
    void setFileFormat( GMSH_FORMAT format )
        {
            M_format = format;
        }

    /**
     * @brief set the gmsh geo
     * @param g gmsh geo
     */
    void setGeo( std::pair<std::string,std::string> const& g ) { M_geo = g; }

    /**
     * set the description of the geometry
     */
    void setDescription( std::string const& desc )
        {
            M_desc = desc;
        }
    void setSubStructuring( bool substruct )
        {
            M_substructuring = substruct;
        }
    bool subStructuring() const
        {
            return M_substructuring;
        }
    /**
     * set the prefix of the Gmsh files
     */
    void setPrefix( std::string const& name )
        {
            M_name = name;
        }
    virtual void setX( std::pair<double,double> const& x )
        {
            FEELPP_ASSERT( dimension() >= 1 )( dimension() ).error( "invalid dimension" );
            M_I[0] = x;
            M_geoParamMap["xmin"]=boost::lexical_cast<std::string>(x.first);
            M_geoParamMap["xmax"]=boost::lexical_cast<std::string>(x.second);
        }
    virtual void setY( std::pair<double,double> const& y )
        {
            FEELPP_ASSERT( dimension() >= 2 )( dimension() ).warn( "invalid dimension" );

            if ( dimension() >= 2 )
            {
                M_I[1] = y;
                M_geoParamMap["ymin"]=boost::lexical_cast<std::string>(y.first);
                M_geoParamMap["ymax"]=boost::lexical_cast<std::string>(y.second);
            }
        }
    virtual void setZ( std::pair<double,double> const& z )
        {
            FEELPP_ASSERT( dimension() >= 3 )( dimension() ).warn( "invalid dimension" );

            if ( dimension() >= 3 )
            {
                M_I[2] = z;
                M_geoParamMap["zmin"]=boost::lexical_cast<std::string>(z.first);
                M_geoParamMap["zmax"]=boost::lexical_cast<std::string>(z.second);
            }
        }

    //! the gmsh generator to generate a reference domain
    virtual void setReferenceDomain()
        {
            if ( dimension() >= 1 )
                M_I[0] = std::make_pair( -1, 1 );

            if ( dimension() >= 2 )
                M_I[1] = std::make_pair( -1, 1 );

            if ( dimension() >= 3 )
                M_I[2] = std::make_pair( -1, 1 );
        }

    //! set the characteristic length to \p h
    virtual void setCharacteristicLength( double _h )
        {
            M_h = _h;
        }

    //! set number of subdivison in x-direction
    virtual void setNx( double _nx )
        {
            M_nx = _nx;
        }

    //! set number of subdivison in y-direction
    virtual void setNy( double _ny )
        {
            M_ny = _ny;
        }

    //! set number of subdivison in z-direction
    virtual void setNz( double _nz )
        {
            M_nz = _nz;
        }


    /**
     * if add is true, set M_addmidpoint to true, false otherwise
     */
    void setAddMidPoint( bool add )
        {
            M_addmidpoint = add;
        }

    //! \brief Modify an existing geo parameter.
    //! If the parameter does not match any parameter, the function throws
    //! an out_of_range exception.
    //!     \param _name Geo parameter name.
    //!     \param _value Geo parameter value.
    //! \return Return the current Gmsh object.
    void setGeoParameter( std::string const& _name, double _value )
        {
            M_geoParamMap[ _name ] = boost::lexical_cast<std::string>( _value );
        }

    //! \brief Modify geo gmsh geometry parameters from a map of parameters.
    //! If the parameter does not match any parameter, the function throws
    //! an out_of_range exception.
    //!     \param geomap A map containing the geo parameters (param,value).
    void setGeoParameters( std::map<std::string, std::string> const& geomap, bool _update=1 )
        {
            if( _update )
            {
                for( const auto& iter : geomap)
                    M_geoParamMap[iter.first] = iter.second;
            }
            else
                M_geoParamMap = geomap;
        }

    /**
     * Set the use of physical names to describe the boundaries of the domain: if \p option
     * is set to true then the generator will generate a PhysicalNames Section and replace
     * numerical id by strings for the Physical boundaries
     */
    void usePhysicalNames( bool option )
        {
            M_usePhysicalNames = option;
        }

    //! set the communicator
    void setWorldComm( WorldComm const& _worldcomm )
        {
            M_worldComm = _worldcomm;
            M_partitions = M_worldComm.size();
        }

    //! set the number of partitions
    void setNumberOfPartitions( int n )
        {
            M_partitions = n;
        }

    //! set save msh file by partitions
    void setMshFileByPartition( bool p )
        {
            M_partition_file = p;
        }

    void setStructuredMesh( int s )
        {
            M_structured = s;
        }
    void setRefinementLevels( int levels )
        {
            M_refine_levels = levels;
        }
    //! set the partitioner
    void setPartitioner( GMSH_PARTITIONER const& p )
        {
            M_partitioner = p;
        }

    //! shear the domain
    void setShear( double _shear )
        {
            M_shear = _shear;
        }


    //! recombine simplices into quads
    void setRecombine( bool _recombine )
        {
            M_recombine = _recombine;
        }

    void setPeriodic( PeriodicEntities const& p )
        {
            M_periodic = p;
        }
    //@}

    /** \name  Methods
     */
    //@{

    /**
     * generate a Gmsh msh file from \p name
     */
    std::string generate( std::string const& name ) const;

    /**
     * \param name  filename prefix to create the \c geo and \c msh file from \p geo
     * \param geo gmsh geometry description
     * \param forceRebuild if true, rebuild the mesh even if geofile is unchanged
     *        if false, rebuild only if geo file has changed.
     *        Useful if generateGeo has been called outside or if gmsh lybrary has changed.
     * \return the name of the mesh file generate by \c gmsh (with the \c .msh extension)
     */
    boost::tuple<std::string,bool>
    generate( std::string const& name,
              std::string const& geo,
              bool const forceRebuild = false,
              bool const parametric = false,
              bool const modifGeo = true ) const;

    /**
     * refine the mesh uniformly by splitting
     * - \param name  name of the gmsh mesh file
     * - \param level the number of refinements
     */
    std::string refine( std::string const& name, int level = 1, bool const parametric = false ) const;

    //! \return the preamble for gmsh geometries
    std::string preamble() const;

    /**
     * \return the content of the geo file \p file in a \p std::string
     */
    std::string getDescriptionFromFile( std::string const& file ) const;

    /**
     * load mesh and generate a new partion of this mesh
     */
    void rebuildPartitionMsh( std::string const& nameMshInput,std::string const& nameMshOutput ) const;

    //! Extract all parameters from a geo gmsh geometry description and store them into a map.
    //! \param geo Gmsh geometry description.
    //! \return Geo parameter map containing each parameter and its value.
    std::map<std::string, std::string>  retrieveGeoParameters( std::string const& geo ) const;

    //! \brief Create a map from a list of geometry parameters string and separated
    //! by a character `:`.
    //!     \param geopars List of parameters as `key=value`. Each new parameter
    //! is separated by a char `:`.
    //! \return Return a map of GMSH geometry parameters and their values. If the string
    //! is empty, it returns an empty map.
    static std::map<std::string, std::string> gpstr2map( std::string const& geopars );

    //@}

protected:

    /**
     * sublass must provide the geo description
     */
    virtual std::string getDescription() const
        {
            return std::string();
        }

    /**
     * \param name  filename prefix to create the \c geo
     * \param geo gmsh geometry description
     * \return returns whether geo file has changed or not. Usually called inside generate, but may be used to
     *         just generate the geo file.
     *         Note: if you use it alone, generate will call this routine again, hence generate needs to know
     *         whether to regenerate the mesh or not
     */
    bool generateGeo( std::string const& name, std::string const& geo,bool const modifGeo=true ) const;

private:

    void generate( std::string const& __name, uint16_type dim, bool parametric = false ) const;

    std::string  prefix( std::string const& __name, uint16_type dim ) const;

protected:
    //! communicator
    WorldComm M_worldComm;

    //! mesh dimension
    int M_dimension;

    //! mesh order
    int M_order;
    // gmsh
    std::string M_version;

    // gmsh file format (ascii or binary)
    GMSH_FORMAT M_format;

    // name of the file
    std::string M_name;

    // description of the geometry
    mutable std::string M_desc;

    // geometry parameters map
    std::map< std::string, std::string > M_geoParamMap;

    bool M_in_memory;

    //! bounding box
    std::vector<std::pair<double,double> > M_I;
    //! characteristic length
    double M_h;
    //! number of discretization in X direction
    double M_nx;
    //! number of discretization in Y direction
    double M_ny;
    //! number of discretization in Z direction
    double M_nz;
    //! mid point
    bool M_addmidpoint;
    //! add physical names to msh files
    bool M_usePhysicalNames;

    //! partitioner type
    GMSH_PARTITIONER M_partitioner;
    //! number of partitions
    int M_partitions;
    //! save msh file by partition
    bool M_partition_file;
    //! shear
    double M_shear;
    //! recombine simplices into hypercubes
    bool M_recombine;
    // build structured mesh
    int M_structured;
    //! number of refinement levels
    int M_refine_levels;

    bool M_substructuring;

    PeriodicEntities M_periodic;

    mutable std::pair<std::string,std::string> M_geo;

#if defined( FEELPP_HAS_GMSH_H )
    mutable GModel*  M_gmodel;
#endif
};

///! \typedef gmsh_type Gmsh
typedef Gmsh gmsh_type;
///! \typedef gmsh_ptrtype boost:shared_ptr<gmsh_type>
typedef boost::shared_ptr<gmsh_type> gmsh_ptrtype;

} // Feel

#endif /* __Gmsh_H */
