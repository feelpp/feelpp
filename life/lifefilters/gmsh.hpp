/* -*- mode: c++; coding: utf-8 -*-

  This file is part of the Life library

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

#include <life/lifecore/life.hpp>
#include <life/lifecore/parameter.hpp>
#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>
#include <life/lifefilters/gmshenums.hpp>
namespace Life
{
const std::string LIFE_GMSH_FORMAT_VERSION="2.1";
}

#include <life/lifefilters/importergmsh.hpp>

namespace Life
{

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
        typedef Life::Singleton< Life::Factory< Gmsh, std::string > > type;
    };

    /** \name Constructors, destructor
     */
    //@{

    Gmsh( int nDim = 1, int nOrder = GMSH_ORDER_ONE )
        :
        M_dimension( nDim ),
        M_order( nOrder ),
        M_version( LIFE_GMSH_FORMAT_VERSION ),
        M_I( nDim ),
        M_h( 0.1 ),
        M_addmidpoint( true ),
        M_usePhysicalNames( false )
        {}
    Gmsh( Gmsh const & __g )
        :
        M_dimension( __g.M_dimension ),
        M_order( __g.M_order ),
        M_version( __g.M_version ),
        M_addmidpoint( __g.M_addmidpoint ),
        M_usePhysicalNames( __g.M_usePhysicalNames )
        {}
    ~Gmsh()
        {}

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
    std::string description() const { return this->getDescription(); }

    /**
     * add the mid point of the domain
     */
    bool addMidPoint() const { return M_addmidpoint; }

    /**
     * \return true if use the physical name, 'false' otherwise
     */
    bool usePhysicalNames() const { return M_usePhysicalNames; }


    //@}

    /** \name  Mutators
     */
    //@{



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
            if ( version != "1" && version != "2" && version != LIFE_GMSH_FORMAT_VERSION )
                throw std::invalid_argument( "invalid gmsh file format version" );
            M_version = version;
        }

    virtual void setX( std::pair<double,double> const& x )
        {
            LIFE_ASSERT( dimension() >= 1 )( dimension() ).error( "invalid dimension" );
            M_I[0] = x;
        }
    virtual void setY( std::pair<double,double> const& y )
        {
            LIFE_ASSERT( dimension() >= 2 )( dimension() ).warn( "invalid dimension" );
            if ( dimension() >= 2 )
                M_I[1] = y;
        }
    virtual void setZ( std::pair<double,double> const& z )
        {
            LIFE_ASSERT( dimension() >= 3 )( dimension() ).warn( "invalid dimension" );
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

    //@}

    /** \name  Methods
     */
    //@{

    std::string generateLine( std::string const& __name, double __h );

    std::string generateCube( std::string const& __name, double __h );

    std::string generateCircle( std::string const& __name, double __h );

    std::string generateSquare( std::string const& __name, double __h );

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
    std::string generate( std::string const& name, std::string const& geo, bool const forceRebuild = false) const;

    /**
     * refine the mesh uniformly by splitting
     * - \param name  name of the gmsh mesh file
     * - \param level the number of refinements
     */
    void refine( std::string const& name, int level = 1 ) const;

    //! \return the preamble for gmsh geometries
    std::string preamble() const;

    //@}

protected:

    /**
     * sublass must provide the geo description
     */
    virtual std::string getDescription() const {}

    /**
     * \param name  filename prefix to create the \c geo
     * \param geo gmsh geometry description
     * \return returns whether geo file has changed or not. Usually called inside generate, but may be used to
     *         just generate the geo file.
     *         Note: if you use it alone, generate will call this routine again, hence generate needs to know
     *         whether to regenerate the mesh or not
     */
    bool generateGeo( std::string const& name, std::string const& geo ) const;

private:

    void generate( std::string const& __name, uint16_type dim ) const;

    std::string  prefix( std::string const& __name, uint16_type dim ) const;

protected:
    mpi::communicator M_comm;

    //! mesh dimension
    int M_dimension;

    //! mesh order
    int M_order;
    // gmsh
    std::string M_version;
    //! bounding box
    std::vector<std::pair<double,double> > M_I;
    //! characteristic length
    double M_h;
    //! mid point
    bool M_addmidpoint;
    //! add physical names to msh files
    bool M_usePhysicalNames;
};

///! \typedef gmsh_type Gmsh
typedef Gmsh gmsh_type;
///! \typedef gmsh_ptrtype boost:shared_ptr<gmsh_type>
typedef boost::shared_ptr<gmsh_type> gmsh_ptrtype;

/// \cond DETAIL
namespace detail
{
template<typename Args>
struct mesh
{
    typedef typename boost::remove_pointer<
        typename boost::remove_const<
            typename boost::remove_reference<
                typename parameter::binding<Args, tag::mesh>::type
                >::type
            >::type
        >::type type;
    typedef boost::shared_ptr<type> ptrtype;
};
}
/// \endcond

/**
 *
 * \brief load a mesh data structure (hold in a shared_ptr<>) using GMSH
 *
 * \arg mesh mesh data structure
 * \arg filename filename string (with extension)
 * \arg refine optionally refine with \p refine levels the mesh (default: 0)
 * \arg update update the mesh data structure (build internal faces and edges) (default : true)
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
     (refine,          *(boost::is_integral<mpl::_>), 0 )
     (update,          *(boost::is_integral<mpl::_>), 0 )
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
    _mesh->accept( import );

    if ( update )
    {
        _mesh->components().set( update );
        _mesh->updateForUse();
    }
    return _mesh;
}

/**
 *
 * \brief create a mesh data structure (hold in a shared_ptr<>) using GMSH
 *
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
     (h,              *(boost::is_floating_point<mpl::_>), 0.1 )
     (order,          *(boost::is_integral<mpl::_>), 1 )
     (refine,          *(boost::is_integral<mpl::_>), 0 )
        )
    )
{
    typedef typename detail::mesh<Args>::type _mesh_type;
    typedef typename detail::mesh<Args>::ptrtype _mesh_ptrtype;

    _mesh_ptrtype _mesh( mesh );
    std::string fname = boost::get<2>( desc )->generate( boost::get<0>( desc ), boost::get<1>( desc ) );

    // refinement if option is enabled to a value greater or equal to 1
    if ( refine )
    {
        Gmsh gmsh;
        gmsh.refine( fname, refine );
    }

    ImporterGmsh<_mesh_type> import( fname );
    _mesh->accept( import );

    return _mesh;
}



BOOST_PARAMETER_FUNCTION(
    (boost::tuple<std::string,std::string,gmsh_ptrtype>), // return type
    domain,    // 2. function name
    tag,           // 3. namespace of tag types
    (required
     (name,           *(boost::is_convertible<mpl::_,std::string>))
     (shape,          *(boost::is_convertible<mpl::_,std::string>)))
    (optional
     (dim,            *(boost::is_integral<mpl::_>)      , 2)
     (order,          *(boost::is_integral<mpl::_>)      , 1)
     (h,              *(boost::is_floating_point<mpl::_>), double(0.1) )
     (convex,         *(boost::is_convertible<mpl::_,std::string>), "simplex")
     (addmidpoint,    *(boost::is_integral<mpl::_>), true )
     (usenames,       *(boost::is_integral<mpl::_>), false )
     (xmin,           *(boost::is_floating_point<mpl::_>), 0. )
     (xmax,           *(boost::is_floating_point<mpl::_>), 1 )
     (ymin,           *(boost::is_floating_point<mpl::_>), 0. )
     (ymax,           *(boost::is_floating_point<mpl::_>), 1 )
     (zmin,           *(boost::is_floating_point<mpl::_>), 0. )
     (zmax,           *(boost::is_floating_point<mpl::_>), 1 )))
{
    gmsh_ptrtype gmsh_ptr = Gmsh::New( shape, dim, order, convex );

    gmsh_ptr->setCharacteristicLength( h );
    gmsh_ptr->setAddMidPoint( addmidpoint );
    gmsh_ptr->usePhysicalNames( usenames );
    gmsh_ptr->setX( std::make_pair( xmin, xmax ) );
    if ( dim >= 2 )
        gmsh_ptr->setY( std::make_pair( ymin, ymax ) );
    if ( dim >= 3 )
        gmsh_ptr->setZ( std::make_pair( zmin, zmax ) );
    return boost::make_tuple( name, gmsh_ptr->description(), gmsh_ptr );
}

/**
 * \fn boost::tuple<std::string,std::string,gmsh_ptrtype> domain( std::string _name, std::string _shape, int _dim, int _order, double _h, std::string _convex)
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
} // Life

#endif /* __Gmsh_H */
