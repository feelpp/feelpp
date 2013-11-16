/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4 tw=0

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-10

  Copyright (C) 2005,2006 EPFL
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
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>
#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelfilters/gmshenums.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelmesh/meshmover.hpp>

namespace Feel
{
extern const char* FEELPP_GMSH_FORMAT_VERSION;
}

#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/exportergmsh.hpp>

namespace Feel
{
class PeriodicEntities: public std::map<int,std::pair<int,int> > {};

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
    static boost::shared_ptr<Gmsh> New( std::string const& shape, uint16_type d = 2, uint16_type o = 1, std::string const& ct = "simplex" );

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
     * set file \p format: ascii or binary
     */
    void setFileFormat( GMSH_FORMAT format )
        {
            M_format = format;
        }

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
            M_geoParamMap.at( _name ) = boost::lexical_cast<std::string>( _value );
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
                    M_geoParamMap.at(iter.first) = iter.second;
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

    //! bounding box
    std::vector<std::pair<double,double> > M_I;
    //! characteristic length
    double M_h;
    //! number of discretization in X direction
    double M_nx;
    //! number of discretization in Y direction
    double M_ny;
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
};

///! \typedef gmsh_type Gmsh
typedef Gmsh gmsh_type;
///! \typedef gmsh_ptrtype boost:shared_ptr<gmsh_type>
typedef boost::shared_ptr<gmsh_type> gmsh_ptrtype;

/// \cond DETAIL
namespace detail
{
template<typename Args, typename Tag=tag::mesh>
struct mesh
{
    typedef typename boost::remove_pointer<
        typename boost::remove_const<
            typename boost::remove_reference<
                typename parameter::binding<Args, Tag>::type
                >::type
            >::type
    >::type _type;
typedef typename mpl::if_<is_shared_ptr<_type>,
                          mpl::identity<typename _type::element_type>,
                          mpl::identity<_type> >::type::type type;
typedef boost::shared_ptr<type> ptrtype;
};

template<typename Args, typename Tag=tag::geoentity>
struct meshFromGeoEntity
{
    typedef typename boost::remove_pointer<
        typename boost::remove_const<
            typename boost::remove_reference<
                typename parameter::binding<Args, Tag>::type
                >::type
            >::type
    >::type _type;

    typedef typename _type::GeoShape GeoShape;
    typedef typename mpl::if_< mpl::bool_<GeoShape::is_simplex>,
                               mpl::identity< Mesh< Simplex< GeoShape::nDim,GeoShape::nOrder,GeoShape::nRealDim> > >,
                               mpl::identity< Mesh< Hypercube< GeoShape::nDim,GeoShape::nOrder,GeoShape::nRealDim> > >
                               >::type::type type;

    typedef PointSet<GeoShape, typename type::value_type> pointset_type;
};


template <typename ElementSpaceType>
void
straightenMeshUpdateEdgesOnBoundaryIsolated( ElementSpaceType & straightener, mpl::int_<0> /**/ )
{}
template <typename ElementSpaceType>
void
straightenMeshUpdateEdgesOnBoundaryIsolated( ElementSpaceType & straightener, mpl::int_<1> /**/ )
{}
template <typename ElementSpaceType>
void
straightenMeshUpdateEdgesOnBoundaryIsolated( ElementSpaceType & straightener, mpl::int_<2> /**/ )
{}
template <typename ElementSpaceType>
void
straightenMeshUpdateEdgesOnBoundaryIsolated( ElementSpaceType & straightener, mpl::int_<3> /**/ )
{
    typedef typename ElementSpaceType::functionspace_type space_type;
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::dof_type::fe_type fe_type;

    auto const ncdof = space_type::dof_type::nComponents;
    auto const dofshift = fe_type::nDofPerVertex*mesh_type::element_type::numVertices;

    auto mesh = straightener.mesh();
    auto const myrank = mesh->worldComm().localRank();

    std::set<size_type> edgeIdFoundToUpdate;

    auto itedge = mesh->beginEdgeOnBoundary();
    auto const enedge = mesh->endEdgeOnBoundary();
    for ( ; itedge!=enedge ; ++itedge )
    {
        if (itedge->processId()!=myrank || itedge->numberOfProcGhost()==0 ) continue;

        auto const theedgeid = itedge->id();

        std::set<size_type> ghostFaceIdFoundOnBoundary;

        auto itprocghost=itedge->elementsGhost().begin();
        auto const enprocghost=itedge->elementsGhost().end();
        for ( ; itprocghost!=enprocghost ; ++itprocghost)
        {
            auto iteltghost = itprocghost->second.begin();
            auto const eneltghost = itprocghost->second.end();
            for ( ; iteltghost!=eneltghost ; ++iteltghost )
            {
                auto const& eltGhost = mesh->element(*iteltghost,itprocghost->first);
                for ( uint16_type f = 0 ; f < mesh_type::element_type::numTopologicalFaces ; ++f )
                {
                    auto const& theface = eltGhost.face(f);
                    if ( theface.isOnBoundary() )
                    {
                        bool findEdge=false;
                        for ( uint16_type e = 0; e < mesh_type::face_type::numEdges && !findEdge ; ++e )
                        {
                            if ( theface.edge(e).id() == theedgeid) { findEdge=true; ghostFaceIdFoundOnBoundary.insert(theface.id());}
                        }
                    }
                }
            }
        } // for ( ; itprocghost!=enprocghost ; ++itprocghost)

        // if 2 faces are find then the edge must be straigten
        if (ghostFaceIdFoundOnBoundary.size()==2) edgeIdFoundToUpdate.insert(theedgeid);

    } // for ( ; itedge!=enedge ; ++itedge )


    if (edgeIdFoundToUpdate.size() > 0)
        {
            auto iteltactif = mesh->beginElementOnBoundary();
            auto const eneltactif = mesh->endElementOnBoundary();
            for ( ; iteltactif!=eneltactif ; ++iteltactif )
            {
                //if (iteltactif->processId()!=myrank) continue;

                for ( uint16_type e = 0; e < mesh_type::element_type::numEdges ; ++e )
                {
                    if ( edgeIdFoundToUpdate.find(iteltactif->edge(e).id()) != edgeIdFoundToUpdate.end())
                    {
                        //std::cout << "find edge " << std::endl;
                        auto const idEltFind = iteltactif->id();
                        for ( uint16_type locdof = 0 ; locdof<fe_type::nDofPerEdge ; ++locdof )
                            {
                                auto const local_id = dofshift + e*fe_type::nDofPerEdge + locdof;

                                for ( uint16_type comp = 0; comp < ncdof; ++comp )
                                    {
                                        auto const globdof = straightener.functionSpace()->dof()->localToGlobal( idEltFind, local_id, comp ).template get<0>();
                                        //std::cout << straightener.functionSpace()->dof()->dofPoint( globdof ).template get<0>() << std::endl;
                                        straightener(globdof) = 0;
                                    }
                            }
                    }
                }

            } // for ( ; iteltactif!=eneltactif ; ++iteltactif )
        } // if (edgeIdFoundToUpdate.size() > 0)


} // straightenMeshUpdateEdgesOnBoundaryIsolated

} // namespace detail
/// \endcond

/**
 *
 * \brief straighten the internal faces of a high order mesh
 *
 * \arg mesh mesh data structure
 */
BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::mesh<Args>::ptrtype ), // return type
    straightenMesh,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( mesh, * )
        )
    ( optional
      ( refine,          *( boost::is_integral<mpl::_> ), 0 )
      ( save,          *( boost::is_integral<mpl::_> ), 0 )
      ( worldcomm, ( WorldComm ), Environment::worldComm() )
        ) )
{
    typedef typename Feel::detail::mesh<Args>::type _mesh_type;
    typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

    VLOG(1) << "straighten mesh of order " <<  _mesh_type::nOrder << " start";

    _mesh_ptrtype _mesh( mesh );

    using namespace vf;
    typedef FunctionSpace<_mesh_type,bases<Lagrange<_mesh_type::nOrder,Vectorial> > > space_t;
#if defined(FEELPP_ENABLE_MPI_MODE)
    auto Xh = space_t::New( _mesh=_mesh, _worldscomm=std::vector<WorldComm>(1,worldcomm) );
#else
    auto Xh = space_t::New( _mesh=_mesh );
#endif

    auto xHo = vf::project( _space=Xh, _range=elements( mesh ), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_HO );
    auto xLo = vf::project( _space=Xh, _range=elements( mesh ), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_O1 );
    auto xHoBdy = vf::project( _space=Xh, _range=boundaryfaces( mesh ), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_HO );
    auto xLoBdy = vf::project( _space=Xh, _range=boundaryfaces( mesh ), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_O1 );

    auto straightener = Xh->element();
    straightener=( xLo-xHo )-( xLoBdy-xHoBdy );

    if (worldcomm.localSize()>1)
        Feel::detail::straightenMeshUpdateEdgesOnBoundaryIsolated( straightener,mpl::int_<_mesh_type::nDim>() );

    double norm_mean_value = integrate( _range=boundaryfaces( _mesh ), _expr=idv( straightener ) ).evaluate(true,worldcomm).norm();

    if ( norm_mean_value > 1e-12 )
        std::cout << "the straightening process may have failed\n"
                  << "norm of component-wise mean value of displacement on the boundary should be 0"
                  << "norm_mean_value: "  << norm_mean_value << "\n"
                  << "you should consider not using straightenMesh()\n"
                  << "\n";

    boost::shared_ptr<Exporter<_mesh_type,_mesh_type::nOrder> > exporter;

    if ( save )
    {
        exporter = Exporter<_mesh_type,_mesh_type::nOrder>::New( "gmsh"/*test_app->vm()*/, "straightener" );
        exporter->step( 0 )->setMesh( _mesh );
        exporter->step( 0 )->add( "xHo", xHo );
        exporter->step( 0 )->add( "xLo", xLo );
        exporter->step( 0 )->add( "xHoBdy", xHoBdy );
        exporter->step( 0 )->add( "xLoBdy", xLoBdy );
        exporter->step( 0 )->add( "straightener", straightener );
        exporter->save();
    }

    MeshMover<_mesh_type> meshmove;
    meshmove.apply( _mesh, straightener );

    VLOG(1) << "straighten mesh of order " <<  _mesh_type::nOrder << " finish";

    return _mesh;
}

/**
 *
 * \brief load a mesh data structure (hold in a shared_ptr<>) using GMSH
 *
 * \arg mesh mesh data structure
 * \arg filename filename string (with extension)
 * \arg refine optionally refine with \p refine levels the mesh (default: 0)
 * \arg update update the mesh data structure (build internal faces and edges) (default : true)
 * \arg physical_are_elementary_regions boolean to load specific meshes formats (default : false)
 */
BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::mesh<Args>::ptrtype ), // return type
    loadGMSHMesh,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( mesh, * )
      ( filename, * )
        ) // 4. one required parameter, and

    ( optional
      ( straighten,          *( boost::is_integral<mpl::_> ), option(_name="gmsh.straighten").template as<bool>() )
      ( refine,          *( boost::is_integral<mpl::_> ), option(_name="gmsh.refine").template as<int>() )
      ( update,          *( boost::is_integral<mpl::_> ), 0 )
      ( physical_are_elementary_regions,		   *, option(_name="gmsh.physical_are_elementary_regions").template as<bool>() )
      ( worldcomm,       *, Environment::worldComm() )
      ( rebuild_partitions,	(bool), option(_name="gmsh.partition").template as<bool>() )
      ( rebuild_partitions_filename,	*, filename )
      ( partitions,      *( boost::is_integral<mpl::_> ), worldcomm.globalSize() )
      ( partitioner,     *( boost::is_integral<mpl::_> ), option(_name="gmsh.partitioner").template as<int>() )
      ( partition_file,   *( boost::is_integral<mpl::_> ), 0 )
        )
    )
{
    typedef typename Feel::detail::mesh<Args>::type _mesh_type;
    typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

    _mesh_ptrtype _mesh( mesh );
    _mesh->setWorldComm( worldcomm );

    std::string filename_with_path = Environment::findFile( filename );
    if ( filename_with_path.empty() )
    {
        std::vector<std::string> plist = Environment::geoPathList();
        std::ostringstream ostr;
        std::for_each( plist.begin(), plist.end(), [&ostr]( std::string s ) { ostr << " - " << s << "\n"; } );
        CHECK( !filename_with_path.empty() ) << "File " << filename << " cannot be found in the following paths list:\n " << ostr.str();
    }

    Gmsh gmsh( _mesh_type::nDim,_mesh_type::nOrder, worldcomm );
    gmsh.setRefinementLevels( refine );
    gmsh.setNumberOfPartitions( partitions );
    gmsh.setPartitioner( (GMSH_PARTITIONER)partitioner );
    gmsh.setMshFileByPartition( partition_file );


    // refinement if option is enabled to a value greater or equal to 1
    if ( refine )
    {
        filename_with_path = gmsh.refine( filename_with_path, refine );
    }
    else if ( rebuild_partitions )
    {
        gmsh.rebuildPartitionMsh(filename_with_path,rebuild_partitions_filename);
        filename_with_path=rebuild_partitions_filename;
    }

    ImporterGmsh<_mesh_type> import( filename_with_path, FEELPP_GMSH_FORMAT_VERSION, worldcomm );

    // need to replace physical_region by elementary_region while reading
    if ( physical_are_elementary_regions )
    {
        import.setElementRegionAsPhysicalRegion( physical_are_elementary_regions );
    }

    _mesh->accept( import );

    if ( update )
    {
        _mesh->components().reset();
        _mesh->components().set( update );
        _mesh->updateForUse();
    }

    else
    {
        _mesh->components().reset();
    }

    if ( straighten && _mesh_type::nOrder > 1 )
        return straightenMesh( _mesh=_mesh,
                               _worldcomm=worldcomm.subWorldComm() );

    return _mesh;
}

/**
 *
 * \brief save a mesh data structure (hold in a shared_ptr<>) in the GMSH format
 *
 * \arg mesh mesh data structure
 * \arg filename filename string (with extension)
 */
BOOST_PARAMETER_FUNCTION(
    ( void ),          // return type
    saveGMSHMesh,    // 2. function name
    tag,             // 3. namespace of tag types
    ( required
      ( mesh, * )
      ( filename, * ) ) // 4. one required parameter, and
    ( optional
      ( parametricnodes,          *( boost::is_integral<mpl::_> ), 0 ) )
    )
{
    typedef typename Feel::detail::mesh<Args>::type _mesh_type;
    typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

#if BOOST_FILESYSTEM_VERSION == 3
    ExporterGmsh<_mesh_type,1> exporter( fs::path( filename ).stem().string(), 1,  mesh->worldComm() );
#elif BOOST_FILESYSTEM_VERSION == 2
    ExporterGmsh<_mesh_type,1> exporter( fs::path( filename ).stem(), 1, mesh->worldComm() );
#endif
    exporter.saveMesh( filename, mesh, parametricnodes );

}

BOOST_PARAMETER_FUNCTION(
    ( void ),  // return type
    saveGeoEntityAsGMSHMesh,    // 2. function name
    tag,             // 3. namespace of tag types
    ( required
      ( geoentity, * )
      ( filename, * ) ) // 4. one required parameter, and
    ( optional
      ( pointset, *, typename Feel::detail::meshFromGeoEntity<Args>::pointset_type() ) )
    )
{
    typedef typename Feel::detail::meshFromGeoEntity<Args>::type _mesh_type;

#if BOOST_FILESYSTEM_VERSION == 3
    ExporterGmsh<_mesh_type,1> exporter( fs::path( filename ).stem().string(), 1,  Environment::worldComm().subWorldCommSeq() );
#elif BOOST_FILESYSTEM_VERSION == 2
    ExporterGmsh<_mesh_type,1> exporter( fs::path( filename ).stem(), 1, Environment::worldComm().subWorldCommSeq() );
#endif
    exporter.gmshSaveOneElementAsMesh( filename, geoentity, pointset );
}


/**
 *
 * \brief create a mesh data structure (hold in a shared_ptr<>) using GMSH
 *
 * \arg mesh mesh data structure
 * \arg descprition
 * \arg h (float, optional, default = 0.1)
 * \arg order (integer, optional, default = 1)
 * \arg parametricnodes (boolean, optional, default = 0)
 * \arg refine (boolean, optional, default = 0)
 * \arg update (boolean, optional, default = 0)
 * \arg force_rebuild boolean (boolean, optional, default = 0)
 * \arg physical_are_elementary_regions change file format (optional, default = false)
 */
BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::mesh<Args>::ptrtype ), // return type
    createGMSHMesh,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( mesh, * )
      ( desc, * )
        ) // 4. one required parameter, and

    ( optional
      ( format,         *, option(_name="gmsh.format").template as<int>() )
      ( h,              *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.hsize").template as<double>() )
      ( geo_parameters,  *( boost::icl::is_map<mpl::_> ), Gmsh::gpstr2map("") )
      ( parametricnodes, *( boost::is_integral<mpl::_> ), 0 )
      ( straighten,      *( boost::is_integral<mpl::_> ), option(_name="gmsh.straighten").template as<bool>() )
      ( refine,          *( boost::is_integral<mpl::_> ), option(_name="gmsh.refine").template as<int>() )
      ( structured,          *( boost::is_integral<mpl::_> ), option(_name="gmsh.structured").template as<int>() )
      ( update,          *( boost::is_integral<mpl::_> ), MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK )
      ( force_rebuild,   *( boost::is_integral<mpl::_> ), 0 )
      ( physical_are_elementary_regions,           *,false )
      ( periodic,        *, PeriodicEntities() )
      ( rebuild_partitions,	(bool), option(_name="gmsh.partition").template as<bool>() )
      ( rebuild_partitions_filename, *( boost::is_convertible<mpl::_,std::string> )	, desc->prefix()+".msh" )
      ( worldcomm,      *, Environment::worldComm() )
      ( partitions,   *( boost::is_integral<mpl::_> ), worldcomm.globalSize() )
      ( partition_file,   *( boost::is_integral<mpl::_> ), 0 )
      ( partitioner,   *( boost::is_integral<mpl::_> ), GMSH_PARTITIONER_CHACO )
        )
    )
{
    typedef typename Feel::detail::mesh<Args>::type _mesh_type;
    typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

    _mesh_ptrtype _mesh( mesh );
    _mesh->setWorldComm( worldcomm );

    if ( worldcomm.isActive() )
    {
        desc->setDimension( mesh->nDim );
        desc->setOrder( mesh->nOrder );
        desc->setWorldComm( worldcomm );
        desc->setNumberOfPartitions( partitions );
        desc->setPartitioner( (GMSH_PARTITIONER) partitioner );
        desc->setMshFileByPartition( partition_file );
        desc->setRefinementLevels( refine );
        desc->setFileFormat( (GMSH_FORMAT)format );
        desc->setStructuredMesh( structured );
        desc->setPeriodic( periodic );

        std::string fname;
        bool generated_or_modified;
        boost::tie( fname, generated_or_modified ) = desc->generate( desc->prefix(), desc->description(), force_rebuild, parametricnodes );

        // refinement if option is enabled to a value greater or equal to 1
        // do not refine if the mesh/geo file was previously generated or modified
        if ( refine && !generated_or_modified )
        {
            VLOG(1) << "Refine mesh ( level: " << refine << ")\n";
            fname = desc->refine( fname, refine, parametricnodes );
        }

        if ( rebuild_partitions )
        {
            desc->rebuildPartitionMsh(fname,rebuild_partitions_filename);
            fname=rebuild_partitions_filename;
        }

        ImporterGmsh<_mesh_type> import( fname, FEELPP_GMSH_FORMAT_VERSION, worldcomm );

        // need to replace physical_regions by elementary_regions for specific meshes
        if ( physical_are_elementary_regions )
        {
            import.setElementRegionAsPhysicalRegion( physical_are_elementary_regions );
        }

        _mesh->accept( import );

        if ( update )
        {
            _mesh->components().reset();
            if ( desc->subStructuring() )
            {
                _mesh->components().set( update|MESH_PROPAGATE_MARKERS );
                _mesh->setSubStructuring( true );
            }
            else
                _mesh->components().set( update );
            _mesh->updateForUse();
        }

        else
        {
            _mesh->components().reset();
        }

        if ( straighten && _mesh_type::nOrder > 1 )
            return straightenMesh( _mesh=_mesh,
                                   _worldcomm=worldcomm.subWorldComm() );
    }
    return _mesh;
}


/**
 * \brief generate a simple geometrical domain from required and optional parameters
 *
 * List of required parameters:
 *  - \param _name name of the file that will ge generated without extension
 *  - \param _shape shape of the domain to be generated (simplex or hypercube)
 * List of optional parameters:
 *  - \param _dim dimension of the domain (default: 2)
 *  - \param _order order of the geometry (default: 1)
 *  - \param _h characteristic size of the mesh (default: 0.1)
 *  - \param _convex type of convex used to mesh the domain (default: simplex) (simplex or hypercube)
 *  - \param _addmidpoint add middle point (default: true )
 *  - \param _xmin minimum x coordinate (default: 0)
 *  - \param _xmax maximum x coordinate (default: 1)
 *  - \param _ymin minimum y coordinate (default: 0)
 *  - \param _ymax maximum y coordinate (default: 1)
 *  - \param _zmin minimum z coordinate (default: 0)
 *  - \param _zmax maximum z coordinate (default: 1)
 *
 * \attention this function uses the Boost.Parameter library that allows to
 * enter the parameter in any order.
 *
 */
BOOST_PARAMETER_FUNCTION(
    ( gmsh_ptrtype ), // return type
    domain,    // 2. function name
    tag,           // 3. namespace of tag types
    ( required
      ( name,           *( boost::is_convertible<mpl::_,std::string> ) )
      )
    ( optional
      ( shape,          *( boost::is_convertible<mpl::_,std::string> ),  option(_name="gmsh.domain.shape").template as<std::string>() )
      ( shear,          *( boost::is_arithmetic<mpl::_> )    ,  option(_name="gmsh.domain.shear").template as<double>() )
      ( recombine,      *( boost::is_integral<mpl::_> )    , option(_name="gmsh.domain.recombine").template as<bool>() )
      ( dim,              *( boost::is_integral<mpl::_> ), 3 )
      ( order,              *( boost::is_integral<mpl::_> ), 1 )
      ( geo_parameters,  *( boost::icl::is_map<mpl::_> ), Gmsh::gpstr2map("") )
      ( h,              *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.hsize").template as<double>() )
      ( convex,         *( boost::is_convertible<mpl::_,std::string> ), option(_name="gmsh.domain.convex").template as<std::string>() )
      ( addmidpoint,    *( boost::is_integral<mpl::_> ), option(_name="gmsh.domain.addmidpoint").template as<bool>() )
      ( usenames,       *( boost::is_integral<mpl::_> ), option(_name="gmsh.domain.usenames").template as<bool>() )
      ( xmin,           *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.domain.xmin").template as<double>() )
      ( xmax,           *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.domain.xmax").template as<double>())
      ( ymin,           *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.domain.ymin").template as<double>() )
      ( ymax,           *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.domain.ymax").template as<double>() )
      ( zmin,           *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.domain.zmin").template as<double>() )
      ( zmax,           *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.domain.zmax").template as<double>() )
      ( nx,             *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.domain.nx").template as<double>() )
      ( ny,             *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.domain.ny").template as<double>() )
      ( substructuring, *( boost::is_integral<mpl::_> ), option(_name="gmsh.domain.substructuring").template as<bool>() ) ) )
{
    gmsh_ptrtype gmsh_ptr = Gmsh::New( shape, 3, 1, convex );

    gmsh_ptr->setPrefix( name );
    gmsh_ptr->setGeoParameters( gmsh_ptr->retrieveGeoParameters( gmsh_ptr->description() ), 0 );
    gmsh_ptr->setGeoParameters( geo_parameters );
    gmsh_ptr->setCharacteristicLength( h );
    gmsh_ptr->setAddMidPoint( addmidpoint );
    gmsh_ptr->usePhysicalNames( usenames );
    gmsh_ptr->setShear( shear );
    gmsh_ptr->setRecombine( recombine );
    gmsh_ptr->setX( std::make_pair( xmin, xmax ) );
    gmsh_ptr->setY( std::make_pair( ymin, ymax ) );
    gmsh_ptr->setZ( std::make_pair( zmin, zmax ) );
    gmsh_ptr->setNx( nx );
    gmsh_ptr->setNy( ny );
    gmsh_ptr->setSubStructuring( substructuring );
    return gmsh_ptr;
}

/**
 * \brief geo return a gmsh_ptrtype of a .geo mesh
 *
 * \arg filename
 * \arg dimension
 * \arg order (optional, default = 1)
 * \arg h (optional, default = 0.1 )
 */
BOOST_PARAMETER_FUNCTION(
    ( gmsh_ptrtype ), // return type
    geo,    // 2. function name
    tag,           // 3. namespace of tag types
    ( required
      ( filename,       *( boost::is_convertible<mpl::_,std::string> ) ) )
    ( optional
      ( h,              *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.hsize").template as<double>() )
      ( geo_parameters,    *( boost::icl::is_map<mpl::_> ), Gmsh::gpstr2map( option(_name="gmsh.geo-variables-list").template as<std::string>() ) )
      ( dim,              *( boost::is_integral<mpl::_> ), 3 )
      ( order,              *( boost::is_integral<mpl::_> ), 1 )
      ( files_path, *( boost::is_convertible<mpl::_,std::string> ), Environment::localGeoRepository() )
      ( depends, *( boost::is_convertible<mpl::_,std::string> ), option(_name="gmsh.depends").template as<std::string>() )
      ( worldcomm,       (WorldComm), Environment::worldComm() ) )
    )

{
    gmsh_ptrtype gmsh_ptr( new Gmsh( 3, 1, worldcomm ) );

    gmsh_ptr->setCharacteristicLength( h );
#if BOOST_FILESYSTEM_VERSION == 3
    gmsh_ptr->setPrefix( fs::path( filename ).stem().string() );
#elif BOOST_FILESYSTEM_VERSION == 2
    gmsh_ptr->setPrefix( fs::path( filename ).stem() );
#endif

    std::string filename_with_path = Environment::findFile( filename );
    if ( filename_with_path.empty() )
    {
        std::vector<std::string> plist = Environment::geoPathList();
        std::ostringstream ostr;
        std::for_each( plist.begin(), plist.end(), [&ostr]( std::string s ) { ostr << " - " << s << "\n"; } );
        CHECK( !filename_with_path.empty() ) << "File " << filename << " cannot be found in the following paths list:\n " << ostr.str();
    }

    gmsh_ptr->setDescription( gmsh_ptr->getDescriptionFromFile( filename_with_path ) );
    gmsh_ptr->setGeoParameters( gmsh_ptr->retrieveGeoParameters( gmsh_ptr->description() ), 0 );
    gmsh_ptr->setGeoParameters( geo_parameters );

    if( worldcomm.globalRank() == worldcomm.masterRank() )
    {
        fs::path cp = fs::current_path();
        std::vector<std::string> depends_on_files;
        if ( !depends.empty() )
            algorithm::split( depends_on_files, depends, algorithm::is_any_of( ":,; " ), algorithm::token_compress_on );
        // copy include/merged files needed by geometry file
        boost::for_each( depends_on_files,
                         [&cp, &files_path]( std::string const& _filename )
                         {
                             fs::path file_path( files_path );
                             file_path /= _filename;

                             try
                             {
                                 boost::system::error_code ec;

                                 if ( !( fs::exists( file_path ) && fs::is_regular_file( file_path ) ) )
                                     std::cout << "File : " << file_path << " doesn't exist or is not a regular file" << std::endl;

                                 else if ( !fs::exists( cp / _filename )  )
                                     fs::copy_file( file_path, fs::path( _filename ), fs::copy_option::none );

                             }

                             catch ( const fs::filesystem_error& e )
                             {
                                 std::cerr << "Error: " << e.what() << std::endl;
                             }
                         } );
    }
    worldcomm.barrier();


    return gmsh_ptr;

}



/**
 * \brief convert to gmsh format
 *
 * \arg filename
 * \arg dim (optional, default = 3)
 * \arg order (optional, default = 1)
 */
BOOST_PARAMETER_FUNCTION(
    ( gmsh_ptrtype ), // return type
    mshconvert,    // 2. function name
    tag,           // 3. namespace of tag types
    ( required
      ( filename,       *( boost::is_convertible<mpl::_,std::string> ) ) )
    ( optional
      ( dim,              *( boost::is_integral<mpl::_> ), 3 )
      ( order,              *( boost::is_integral<mpl::_> ), 1 ) )
    )
{
    gmsh_ptrtype gmsh_ptr( new Gmsh( 3, 1 ) );
#if BOOST_FILESYSTEM_VERSION == 3
    gmsh_ptr->setPrefix( fs::path( filename ).stem().string() );
#elif BOOST_FILESYSTEM_VERSION == 2
    gmsh_ptr->setPrefix( fs::path( filename ).stem() );
#endif

    // first try in the current path
    if ( fs::exists( filename ) )
        gmsh_ptr->setDescription( ( boost::format( "Merge \"%1%\";\n" ) % filename ).str() );

    else if ( fs::exists( fs::path( Environment::localGeoRepository() ) / filename ) )
        gmsh_ptr->setDescription( ( boost::format( "Merge \"%1%\";\n" ) % ( fs::path( Environment::localGeoRepository() ) / filename ).string() ).str() );

    else if ( Environment::systemGeoRepository().template get<1>()  &&
              fs::exists( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ) )
        gmsh_ptr->setDescription( ( boost::format( "Merge \"%1%\";\n" ) % ( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ).string() ).str() );

    else
    {
        std::ostringstream ostr;
        ostr << "File " << filename << " was not found neither in current directory or in " << Environment::localGeoRepository() << " or in " << Environment::systemGeoRepository();
        throw std::invalid_argument( ostr.str() );
    }

    return gmsh_ptr;
}

/**
 * \brief convert to msh format
 *
 * \arg filename
 */
BOOST_PARAMETER_FUNCTION(
    ( std::string ), // return type
    img2msh,    // 2. function name
    tag,           // 3. namespace of tag types
    ( required
      ( filename,       *( boost::is_convertible<mpl::_,std::string> ) ) )
    ( optional
      ( prefix,       *( boost::is_convertible<mpl::_,std::string> ), fs::path( filename ).stem() ) )
    )
{
    gmsh_ptrtype gmsh_ptr( new Gmsh( 2, 1 ) );
    gmsh_ptr->setPrefix( prefix );
    std::string meshname = ( boost::format( "%1%-0.msh" ) % prefix ).str();

    // first try in the current path
    if ( fs::exists( filename ) )
        gmsh_ptr->setDescription( ( boost::format( "Merge \"%1%\";\nSave View [0] \"%2%\";\n" ) % filename % meshname ).str() );

    else if ( fs::exists( fs::path( Environment::localGeoRepository() ) / filename ) )
        gmsh_ptr->setDescription( ( boost::format( "Merge \"%1%\";\nSave View [0] \"%2%\";\n" ) % ( fs::path( Environment::localGeoRepository() ) / filename ).string() % meshname ).str() );

    else if ( Environment::systemGeoRepository().template get<1>()  &&
              fs::exists( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ) )
        gmsh_ptr->setDescription( ( boost::format( "Merge \"%1%\";\nSave View [0] \"%2%\";\n" ) % ( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ).string() % meshname ).str() );

    else
    {
        std::ostringstream ostr;
        ostr << "File " << filename << " was not found neither in current directory or in " << Environment::localGeoRepository() << " or in " << Environment::systemGeoRepository();
        throw std::invalid_argument( ostr.str() );
    }

    gmsh_ptr->generate( gmsh_ptr->prefix(), gmsh_ptr->description() );
    return meshname;
}



/**
 * build a mesh of the unit segment [0,1]
 */
boost::shared_ptr<Mesh<Simplex<1> > > unitSegment( double h = option(_name="gmsh.hsize").as<double>() );

/**
 * build a mesh of the unit square [0,1]^2 using triangles
 */
boost::shared_ptr<Mesh<Simplex<2> > > unitSquare( double h = option(_name="gmsh.hsize").as<double>(),
                                                  PeriodicEntities pe = PeriodicEntities() );

/**
 * build a mesh of the unit circle using triangles
 */
template<int Ngeo=1>
inline
boost::shared_ptr<Mesh<Simplex<2,Ngeo> > >
unitCircle( double h = option(_name="gmsh.hsize").template as<double>() )
{
    return createGMSHMesh(_mesh=new Mesh<Simplex<2,Ngeo> >,
                          _desc=domain( _name="square",
                                        _shape="ellipsoid",
                                        _dim=2,
                                        _xmin=-1,
                                        _ymin=-1,
                                        _h=h ) );
}

/**
 * build a mesh of the unit circle using triangles
 */
template<int Ngeo=1>
inline
boost::shared_ptr<Mesh<Simplex<3,Ngeo> > >
unitSphere( double h = option(_name="gmsh.hsize").template as<double>() )
{
    return createGMSHMesh(_mesh=new Mesh<Simplex<3,Ngeo> >,
                          _desc=domain( _name="sphere",
                                        _shape="ellipsoid",
                                        _dim=3,
                                        _xmin=-1,
                                        _ymin=-1,
                                        _zmin=-1,
                                        _h= h ) );
}


/**
 * build a mesh of the unit square [0,1]^3 using tetrahedrons
 */
boost::shared_ptr<Mesh<Simplex<3> > > unitCube( double h = option(_name="gmsh.hsize").as<double>() );

template<int Dim, typename Convex=Simplex<Dim>>
inline
boost::shared_ptr<Mesh<Convex> >
unitHypercube( double h = option(_name="gmsh.hsize").template as<double>() )
{
    return createGMSHMesh(_mesh=new Mesh<Convex>,
                          _desc=domain( _name="hypercube",
                                        _shape="hypercube",
                                        _convex=Convex::type(),
                                        _dim=Dim,
                                        _h=h ) );
}



/**
 *
 * \brief load a mesh data structure (hold in a shared_ptr<>) using GMSH
 *
 * \arg mesh mesh data structure
 * \arg filename filename string (with extension)
 * \arg refine optionally refine with \p refine levels the mesh (default: 0)
 * \arg update update the mesh data structure (build internal faces and edges) (default : true)
 * \arg physical_are_elementary_regions boolean to load specific meshes formats (default : false)
 */
BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::mesh<Args>::ptrtype ), // return type
    loadMesh,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( mesh, *)

        ) // 4. one required parameter, and

    ( optional
      ( filename, *( boost::is_convertible<mpl::_,std::string> ), option(_name="gmsh.filename").template as<std::string>() )
      ( h,              *( boost::is_arithmetic<mpl::_> ), option(_name="gmsh.hsize").template as<double>() )
      ( straighten,          (bool), option(_name="gmsh.straighten").template as<bool>() )
      ( refine,          *( boost::is_integral<mpl::_> ), option(_name="gmsh.refine").template as<int>() )
      ( update,          *( boost::is_integral<mpl::_> ), MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES )
      ( physical_are_elementary_regions,		   (bool), option(_name="gmsh.physical_are_elementary_regions").template as<bool>() )
      ( worldcomm,       (WorldComm), Environment::worldComm() )
      ( force_rebuild,   *( boost::is_integral<mpl::_> ), option(_name="gmsh.rebuild").template as<bool>() )
      ( rebuild_partitions,	(bool), option(_name="gmsh.partition").template as<bool>() )
      ( rebuild_partitions_filename, *( boost::is_convertible<mpl::_,std::string> )	, filename )
      ( partitions,      *( boost::is_integral<mpl::_> ), worldcomm.globalSize() )
      ( partitioner,     *( boost::is_integral<mpl::_> ), option(_name="gmsh.partitioner").template as<int>() )
      ( partition_file,   *( boost::is_integral<mpl::_> ), 0 )
      ( depends, *( boost::is_convertible<mpl::_,std::string> ), option(_name="gmsh.depends").template as<std::string>() )
        )
    )
{
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsequenced"
#endif
    typedef typename Feel::detail::mesh<Args>::type _mesh_type;
    typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

    // look for mesh_name in various directories (executable directory, current directory. ...)
    // return an empty string if the file is not found

    fs::path mesh_name=fs::path(Environment::findFile(filename));
    LOG_IF( WARNING, mesh_name.extension() != ".geo" && mesh_name.extension() != ".msh" )
        << "Invalid filename " << filename << " it should have either the .geo or .msh extension\n";


    if ( mesh_name.extension() == ".geo" )
    {
        return createGMSHMesh( _mesh=mesh,
                               _desc=geo( _filename=mesh_name.string(),
                                          _h=h,_depends=depends,
                                          _worldcomm=worldcomm ),
                               _h=h,
                               _straighten=straighten,
                               _refine=refine,
                               _update=update,
                               _physical_are_elementary_regions=physical_are_elementary_regions,
                               _force_rebuild=force_rebuild,
                               _worldcomm=worldcomm,
                               _rebuild_partitions=rebuild_partitions,
                               _rebuild_partitions_filename=rebuild_partitions_filename,
                               _partitions=partitions,
                               _partitioner=partitioner,
                               _partition_file=partition_file
            );
    }

    if ( mesh_name.extension() == ".msh"  )
    {
        return loadGMSHMesh( _mesh=mesh,
                             _filename=mesh_name.string(),
                             _straighten=straighten,
                             _refine=refine,
                             _update=update,
                             _physical_are_elementary_regions=physical_are_elementary_regions,
                             _worldcomm=worldcomm,
                             _rebuild_partitions=rebuild_partitions,
                             _rebuild_partitions_filename=rebuild_partitions_filename,
                             _partitions=partitions,
                             _partitioner=partitioner,
                             _partition_file=partition_file
            );

    }

    LOG(WARNING) << "File " << mesh_name << " not found, generating instead an hypercube in " << _mesh_type::nDim << "D geometry and mesh...";
    return createGMSHMesh(_mesh=mesh,
                          _desc=domain( _name=option(_name="gmsh.domain.shape").template as<std::string>(), _h=h ),
                          _h=h,
                          _refine=refine,
                          _update=update,
                          _physical_are_elementary_regions=physical_are_elementary_regions,
                          _force_rebuild=force_rebuild,
                          _worldcomm=worldcomm,
                          _rebuild_partitions=rebuild_partitions,
                          _rebuild_partitions_filename=rebuild_partitions_filename,
                          _partitions=partitions,
                          _partitioner=partitioner,
                          _partition_file=partition_file );
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
}

} // Feel

#endif /* __Gmsh_H */
