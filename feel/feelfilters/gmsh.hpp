
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-02-10
 */
#ifndef __gmsh_H
#define __gmsh_H 1

#include <boost/type_traits.hpp>

#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/range/algorithm/for_each.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>
#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelfilters/gmshenums.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feelfilters/exporterquick.hpp>

namespace Feel
{
extern const char* FEELPP_GMSH_FORMAT_VERSION;
}

#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/exportergmsh.hpp>

namespace Feel
{
enum GMSH_PARTITIONER
{
    GMSH_PARTITIONER_CHACO = 1,
    GMSH_PARTITIONER_METIS = 2

};
enum GMSH_ORDER
{
    GMSH_ORDER_ONE = 1,
    GMSH_ORDER_TWO = 2,
    GMSH_ORDER_THREE = 3,
    GMSH_ORDER_FOUR = 4,
    GMSH_ORDER_FIVE = 5
};
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

    Gmsh( int nDim = 1, int nOrder = GMSH_ORDER_ONE );
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
    Gmsh& operator=( Gmsh const& __g )
        {
            if (  this != &__g )
            {
                M_dimension = __g.M_dimension;
                M_order = __g.M_order;
                M_version = __g.M_version;
                M_addmidpoint = __g.M_addmidpoint;
                M_usePhysicalNames = __g.M_usePhysicalNames;
                M_shear = __g.M_shear;
            }
            return *this;
        }

    static boost::shared_ptr<Gmsh> New( po::variables_map const& vm );
    static boost::shared_ptr<Gmsh> New( std::string const& shape, uint16_type d = 2, uint16_type o = 1, std::string const& ct = "simplex" );

    //@}

    /** \name Accessors
     */
    //@{

    /**
     * \return mesh dimension
     */
    int dimension() const { return M_dimension; }

    /**
     * get the order of the elements of the mesh
     * \return the order of the elements of the mesh
     */
    GMSH_ORDER order() const { return (GMSH_ORDER) M_order; }

    /**
     * @return the file format version
     */
    std::string version() const { return M_version; }

    /**
     * \return the name of the file
     */
    std::string prefix() const { return M_name; }

    /**
     * \return bounding box
     */
    std::vector<std::pair<double,double> > const& boundingBox() const {return M_I;}
    double xmin() const { return M_I[0].first; }
    double xmax() const { return M_I[0].second; }
    double ymin() const { return M_I[1].first; }
    double ymax() const { return M_I[1].second; }
    double zmin() const { return M_I[2].first; }
    double zmax() const { return M_I[2].second; }

    /**
     * \return characteristic length
     */
    double const& h() const {return M_h; }

    /**
     * \return the geometry description
     */
    std::string description() const { std::string d = this->getDescription(); if ( !d.empty() ) M_desc = d; return M_desc; }

    /**
     * add the mid point of the domain
     */
    bool addMidPoint() const { return M_addmidpoint; }

    /**
     * \return true if use the physical name, 'false' otherwise
     */
    bool usePhysicalNames() const { return M_usePhysicalNames; }

    //! \return the world comm
    WorldComm const& worldComm() const { return M_worldComm; }

    //! \return the nnumber of partitions
    int numberOfPartitions() const { return M_partitions; }

    //! \return true if save msh file by partitions, false otherwise
    bool mshFileByPartition() const { return M_partition_file; }

    //! \return the partitioner
    GMSH_PARTITIONER partitioner() const { return M_partitioner; }

    //! get the shear
    double shear() const { return M_shear; }

    //! return true if recombine, false otherwise
    bool recombine() const { return M_recombine; }

    //@}

    /** \name  Mutators
     */
    //@{

    /**
     * set the dimension
     */
    Gmsh& setDimension( int dim ) { M_dimension = dim; return *this; }

    /**
     * the gmsh generator to generate a reference domain
     * \return the mesh generator
     *
     * \code
     * Gmsh gmsh;
     * gmsh = gmsh.ref();
     * \endcode
     */
    Gmsh& ref() { this->setReferenceDomain(); return *this; }

    /**
     * set the characteristic length
     * \param h the characteristic length
     * \return the mesh generator
     */
    Gmsh& h( double _h ) { this->setCharacteristicLength( _h ); return *this; }

    /**
     * set the order of the elements of the mesh it can be either
     * GMSH_ORDER_ONE (order 1/linear) or GMSH_ORDER_TWO(order
     * 2/quadratic)
     *
     * \param o order of the elements
     */
    void setOrder( int o )
        {
            M_order = (GMSH_ORDER) o;
        }

    /**
     * set the file format version
     */
    void setVersion( std::string version )
        {
            if ( version != "1" && version != "2" && version != FEELPP_GMSH_FORMAT_VERSION )
                throw std::invalid_argument( "invalid gmsh file format version" );
            M_version = version;
        }

    /**
     * set the description of the geometry
     */
    void setDescription( std::string const& desc )
        {
            M_desc = desc;
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
        }
    virtual void setY( std::pair<double,double> const& y )
        {
            FEELPP_ASSERT( dimension() >= 2 )( dimension() ).warn( "invalid dimension" );
            if ( dimension() >= 2 )
                M_I[1] = y;
        }
    virtual void setZ( std::pair<double,double> const& z )
        {
            FEELPP_ASSERT( dimension() >= 3 )( dimension() ).warn( "invalid dimension" );
            if ( dimension() >= 3 )
                M_I[2] = z;
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
    virtual void setCharacteristicLength( double _h ) { M_h = _h; }

    /**
     * if add is true, set M_addmidpoint to true, false otherwise
     */
    void setAddMidPoint( bool add ) { M_addmidpoint = add; }

    /**
     * Set the use of physical names to describe the boundaries of the domain: if \p option
     * is set to true then the generator will generate a PhysicalNames Section and replace
     * numerical id by strings for the Physical boundaries
     */
    void usePhysicalNames( bool option ) { M_usePhysicalNames = option; }

    //! set the communicator
    void setWorldComm(WorldComm const& _worldcomm) { M_worldComm = _worldcomm; }

    //! set the number of partitions
    void setNumberOfPartitions( int n ) { M_partitions = n; }

    //! set save msh file by partitions
    void setMshFileByPartition( bool p ) { M_partition_file = p; }

    //! set the partitioner
    void setPartitioner( GMSH_PARTITIONER const& p ) {  M_partitioner = p; }

    //! shear the domain
    void setShear( double _shear ) { M_shear = _shear; }

    //! recombine simplices into quads
    void setRecombine( bool _recombine ) { M_recombine = _recombine; }

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
    std::string generate( std::string const& name,
                          std::string const& geo,
                          bool const forceRebuild = false,
                          bool const parametric = false,
                          bool const modifGeo = true) const;

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
    //@}

protected:

    /**
     * sublass must provide the geo description
     */
    virtual std::string getDescription() const { return std::string(); }

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

    // name of the file
    std::string M_name;

    // description of the geometry
    mutable std::string M_desc;

    //! bounding box
    std::vector<std::pair<double,double> > M_I;
    //! characteristic length
    double M_h;
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
                              mpl::identity<typename _type::value_type>,
                              mpl::identity<_type> >::type::type type;
    typedef boost::shared_ptr<type> ptrtype;
};
}
/// \endcond

/**
 *
 * \brief straighten the internal faces of a high order mesh
 *
 * \arg mesh mesh data structure
 */
BOOST_PARAMETER_FUNCTION(
    (typename detail::mesh<Args>::ptrtype), // return type
    straightenMesh,    // 2. function name

    tag,           // 3. namespace of tag types

    (required
     (mesh, *)
        )
    (optional
     (refine,          *(boost::is_integral<mpl::_>), 0 )
     (save,          *(boost::is_integral<mpl::_>), 0 )
        ))
{
    typedef typename detail::mesh<Args>::type _mesh_type;
    typedef typename detail::mesh<Args>::ptrtype _mesh_ptrtype;

    _mesh_ptrtype _mesh( mesh );

    using namespace vf;
    typedef FunctionSpace<_mesh_type,bases<Lagrange<_mesh_type::nOrder,Vectorial> > > space_t;
    auto Xh = space_t::New( _mesh );
    auto xHo = vf::project( _space=Xh, _range=elements(mesh), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_HO );
    auto xLo = vf::project( _space=Xh, _range=elements(mesh), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_O1 );
    auto xHoBdy = vf::project( _space=Xh, _range=boundaryfaces(mesh), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_HO );
    auto xLoBdy = vf::project( _space=Xh, _range=boundaryfaces(mesh), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_O1 );
    auto straightener = Xh->element();
    straightener=(xLo-xHo)-(xLoBdy-xHoBdy);
    double norm_mean_value = integrate( _range=boundaryfaces(_mesh), _expr=idv(straightener) ).evaluate().norm();
    if ( norm_mean_value > 1e-12 )
        std::cout << "the straightening process may have failed\n"
                  << "norm of component-wise mean value of displacement on the boundary should be 0"
                  << "norm_mean_value: "  << norm_mean_value << "\n"
                  << "you should consider not using straightenMesh()\n"
                  << "\n";

    boost::shared_ptr<Exporter<_mesh_type,_mesh_type::nOrder> > exporter;
    if ( save )
    {
        exporter.reset( Exporter<_mesh_type,_mesh_type::nOrder>::New( "gmsh"/*test_app->vm()*/, "straightener") );
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
    (typename detail::mesh<Args>::ptrtype), // return type
    loadGMSHMesh,    // 2. function name

    tag,           // 3. namespace of tag types

    (required
     (mesh, *)
     (filename, *)
        ) // 4. one required parameter, and

    (optional
     (straighten,          *(boost::is_integral<mpl::_>), 1 )
     (refine,          *(boost::is_integral<mpl::_>), 0 )
     (update,          *(boost::is_integral<mpl::_>), 0 )
	 (physical_are_elementary_regions,		   *,false)
        )
    )
{
    typedef typename detail::mesh<Args>::type _mesh_type;
    typedef typename detail::mesh<Args>::ptrtype _mesh_ptrtype;

    _mesh_ptrtype _mesh( mesh );

    // refinement if option is enabled to a value greater or equal to 1
    if ( refine )
    {
        Gmsh gmsh;
        gmsh.refine( filename, refine );
    }

    ImporterGmsh<_mesh_type> import( filename );

    // need to replace physical_region by elementary_region while reading
    if (physical_are_elementary_regions)
    {
        import.setElementRegionAsPhysicalRegion(physical_are_elementary_regions);
    }

    _mesh->accept(import);

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
        return straightenMesh( _mesh );
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
    (void),            // return type
    saveGMSHMesh,    // 2. function name
    tag,             // 3. namespace of tag types
    (required
     (mesh, *)
     (filename, *))  // 4. one required parameter, and
    (optional
     (parametricnodes,          *(boost::is_integral<mpl::_>), 0 ) )
    )
{
    typedef typename detail::mesh<Args>::type _mesh_type;
    typedef typename detail::mesh<Args>::ptrtype _mesh_ptrtype;

#if BOOST_FILESYSTEM_VERSION == 3
    ExporterGmsh<_mesh_type,1> exporter( fs::path(filename).stem().string() );
#elif BOOST_FILESYSTEM_VERSION == 2
    ExporterGmsh<_mesh_type,1> exporter( fs::path(filename).stem() );
#endif
    exporter.saveMesh( filename, mesh, parametricnodes );
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
    (typename detail::mesh<Args>::ptrtype), // return type
    createGMSHMesh,    // 2. function name

    tag,           // 3. namespace of tag types

    (required
     (mesh, *)
     (desc, *)
        ) // 4. one required parameter, and

    (optional
     (h,              *(boost::is_arithmetic<mpl::_>), 0.1 )
     (parametricnodes,*(boost::is_integral<mpl::_>), 0 )
     (straighten,     *(boost::is_integral<mpl::_>), 1 )
     (refine,          *(boost::is_integral<mpl::_>), 0 )
     (update,          *(boost::is_integral<mpl::_>), 0 )
     (force_rebuild,   *(boost::is_integral<mpl::_>), 0 )
     (physical_are_elementary_regions,           *,false)
     (partitions,   *(boost::is_integral<mpl::_>), 1 )
     (partition_file,   *(boost::is_integral<mpl::_>), 0 )
     (partitioner,   *(boost::is_integral<mpl::_>), GMSH_PARTITIONER_CHACO )
     (worldcomm,      *, WorldComm() )
        )
    )
{
    typedef typename detail::mesh<Args>::type _mesh_type;
    typedef typename detail::mesh<Args>::ptrtype _mesh_ptrtype;

    _mesh_ptrtype _mesh( mesh );
    _mesh->setWorldComm(worldcomm);
    if (worldcomm.isActive())
        {
            desc->setDimension(mesh->nDim);
            desc->setOrder(mesh->nOrder);
            desc->setWorldComm(worldcomm);
            desc->setNumberOfPartitions( partitions );
            desc->setPartitioner( partitioner );
            desc->setMshFileByPartition( partition_file );

            std::string fname = desc->generate( desc->prefix(), desc->description(), force_rebuild, parametricnodes );
            // refinement if option is enabled to a value greater or equal to 1
            if ( refine )
                {
                    Debug() << "Refine mesh ( level: " << refine << ")\n";
                    Gmsh gmsh;
                    fname = gmsh.refine( fname, refine, parametricnodes );
                }

            ImporterGmsh<_mesh_type> import( fname, FEELPP_GMSH_FORMAT_VERSION, worldcomm );
            // need to replace physical_regions by elementary_regions for specific meshes
            if (physical_are_elementary_regions)
                {
                    import.setElementRegionAsPhysicalRegion(physical_are_elementary_regions);
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
                return straightenMesh( _mesh );
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
    (gmsh_ptrtype), // return type
    domain,    // 2. function name
    tag,           // 3. namespace of tag types
    (required
     (name,           *(boost::is_convertible<mpl::_,std::string>))
     (shape,          *(boost::is_convertible<mpl::_,std::string>)))
    (optional
     (shear,          *(boost::is_arithmetic<mpl::_>)    , 0)
     (recombine,      *(boost::is_integral<mpl::_>)    , 0)
     (dim,              *(boost::is_integral<mpl::_>), 3 )
     (order,              *(boost::is_integral<mpl::_>), 1 )
     (h,              *(boost::is_arithmetic<mpl::_>), double(0.1) )
     (convex,         *(boost::is_convertible<mpl::_,std::string>), "Simplex")
     (addmidpoint,    *(boost::is_integral<mpl::_>), true )
     (usenames,       *(boost::is_integral<mpl::_>), false )
     (xmin,           *(boost::is_arithmetic<mpl::_>), 0. )
     (xmax,           *(boost::is_arithmetic<mpl::_>), 1 )
     (ymin,           *(boost::is_arithmetic<mpl::_>), 0. )
     (ymax,           *(boost::is_arithmetic<mpl::_>), 1 )
     (zmin,           *(boost::is_arithmetic<mpl::_>), 0. )
     (zmax,           *(boost::is_arithmetic<mpl::_>), 1 )))
{
    gmsh_ptrtype gmsh_ptr = Gmsh::New( shape, 3, 1, convex );

    gmsh_ptr->setPrefix( name );
    gmsh_ptr->setCharacteristicLength( h );
    gmsh_ptr->setAddMidPoint( addmidpoint );
    gmsh_ptr->usePhysicalNames( usenames );
    gmsh_ptr->setShear( shear );
    gmsh_ptr->setRecombine( recombine );
    gmsh_ptr->setX( std::make_pair( xmin, xmax ) );
    gmsh_ptr->setY( std::make_pair( ymin, ymax ) );
    gmsh_ptr->setZ( std::make_pair( zmin, zmax ) );
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
    (gmsh_ptrtype), // return type
    geo,    // 2. function name
    tag,           // 3. namespace of tag types
    (required
     (filename,       *(boost::is_convertible<mpl::_,std::string>)))
    (optional
     (h,              *(boost::is_arithmetic<mpl::_>), double(0.1) )
     (dim,              *(boost::is_integral<mpl::_>), 3 )
     (order,              *(boost::is_integral<mpl::_>), 1 )
     (files_path, *(boost::is_convertible<mpl::_,std::string>), std::string("."))
     (depends, *(boost::is_convertible<mpl::_,std::list<std::string> >), std::list<std::string>() ))
    )

{
    gmsh_ptrtype gmsh_ptr( new Gmsh( 3, 1 ) );

    gmsh_ptr->setCharacteristicLength( h );
#if BOOST_FILESYSTEM_VERSION == 3
    gmsh_ptr->setPrefix( fs::path(filename).stem().string() );
#elif BOOST_FILESYSTEM_VERSION == 2
    gmsh_ptr->setPrefix( fs::path(filename).stem() );
#endif

    fs::path cp;
    try
    {
        fs::current_path(cp);
    }
    catch(...)
    {

    }
    // first try in the current path
    if ( fs::exists( cp / filename ) )
      {
          gmsh_ptr->setDescription(gmsh_ptr->getDescriptionFromFile((cp/filename).string()));
      }
    else if ( fs::exists( fs::path(Environment::localGeoRepository()) / filename ) )
      {
          gmsh_ptr->setDescription( gmsh_ptr->getDescriptionFromFile((fs::path(Environment::localGeoRepository()) / filename).string()) );
      }
    else if ( Environment::systemGeoRepository().template get<1>()  &&
              fs::exists( fs::path(Environment::systemGeoRepository().get<0>()) / filename ) )
      {
          gmsh_ptr->setDescription( gmsh_ptr->getDescriptionFromFile((fs::path(Environment::systemGeoRepository().get<0>()) / filename).string()) );
      }
    else
    {
        std::ostringstream ostr;
        ostr << "File " << filename << " was not found neither in current directory or in " << Environment::localGeoRepository() << " or in " << Environment::systemGeoRepository();
        throw std::invalid_argument( ostr.str() );
    }

    // copy include/merged files needed by geometry file
    boost::for_each( depends,
                     [&cp, &files_path]( std::string const& _filename)
                     {
                         fs::path file_path( files_path );
                         file_path /= _filename;
                         try
                         {
                             boost::system::error_code ec;
                             if( !( fs::exists(file_path) && fs::is_regular_file( file_path ) ) )
                                 std::cout << "File : " << _filename << " doesn't exist or is not a regular file" << std::endl;
                             else if ( !fs::exists(cp / _filename)  )
                                 fs::copy_file(file_path, fs::path(_filename), fs::copy_option::none );

                         }
                         catch (const fs::filesystem_error& e)
                         {
                             std::cerr << "Error: " << e.what() << std::endl;
                         }
                     });

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
    (gmsh_ptrtype), // return type
    mshconvert,    // 2. function name
    tag,           // 3. namespace of tag types
    (required
     (filename,       *(boost::is_convertible<mpl::_,std::string>)))
    (optional
     (dim,              *(boost::is_integral<mpl::_>), 3 )
     (order,              *(boost::is_integral<mpl::_>), 1 ))
    )
{
    gmsh_ptrtype gmsh_ptr( new Gmsh( 3, 1 ) );
#if BOOST_FILESYSTEM_VERSION == 3
    gmsh_ptr->setPrefix( fs::path(filename).stem().string() );
#elif BOOST_FILESYSTEM_VERSION == 2
    gmsh_ptr->setPrefix( fs::path(filename).stem() );
#endif

    // first try in the current path
    if ( fs::exists( filename ) )
        gmsh_ptr->setDescription((boost::format( "Merge \"%1%\";\n" ) % filename ).str());
    else if ( fs::exists( fs::path(Environment::localGeoRepository()) / filename ) )
        gmsh_ptr->setDescription((boost::format( "Merge \"%1%\";\n" ) % (fs::path(Environment::localGeoRepository()) / filename).string()).str() );
    else if ( Environment::systemGeoRepository().template get<1>()  &&
              fs::exists( fs::path(Environment::systemGeoRepository().get<0>()) / filename ) )
        gmsh_ptr->setDescription( (boost::format( "Merge \"%1%\";\n" ) % (fs::path(Environment::systemGeoRepository().get<0>()) / filename).string()).str() );
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
    (std::string), // return type
    img2msh,    // 2. function name
    tag,           // 3. namespace of tag types
    (required
     (filename,       *(boost::is_convertible<mpl::_,std::string>)))
    (optional
     (prefix,       *(boost::is_convertible<mpl::_,std::string>), fs::path(filename).stem()))
    )
{
    gmsh_ptrtype gmsh_ptr( new Gmsh( 2, 1 ) );
    gmsh_ptr->setPrefix( prefix );
    std::string meshname = (boost::format( "%1%-0.msh" ) % prefix ).str();

    // first try in the current path
    if ( fs::exists( filename ) )
        gmsh_ptr->setDescription((boost::format( "Merge \"%1%\";\nSave View [0] \"%2%\";\n" ) % filename % meshname ).str());
    else if ( fs::exists( fs::path(Environment::localGeoRepository()) / filename ) )
        gmsh_ptr->setDescription((boost::format( "Merge \"%1%\";\nSave View [0] \"%2%\";\n" ) % (fs::path(Environment::localGeoRepository()) / filename).string() % meshname ).str() );
    else if ( Environment::systemGeoRepository().template get<1>()  &&
              fs::exists( fs::path(Environment::systemGeoRepository().get<0>()) / filename ) )
        gmsh_ptr->setDescription( (boost::format( "Merge \"%1%\";\nSave View [0] \"%2%\";\n" ) % (fs::path(Environment::systemGeoRepository().get<0>()) / filename).string() % meshname ).str() );
    else
    {
        std::ostringstream ostr;
        ostr << "File " << filename << " was not found neither in current directory or in " << Environment::localGeoRepository() << " or in " << Environment::systemGeoRepository();
        throw std::invalid_argument( ostr.str() );
    }
    gmsh_ptr->generate(gmsh_ptr->prefix(), gmsh_ptr->description());
    return meshname;
}


} // Feel

#endif /* __Gmsh_H */
