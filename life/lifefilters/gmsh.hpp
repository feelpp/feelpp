/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-02-10

  Copyright (C) 2005,2006 EPFL
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
   \file gmsh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-02-10
 */
#ifndef __gmsh_H
#define __gmsh_H 1

#include <life/lifecore/life.hpp>
#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>

namespace Life
{

enum GMSH_ORDER
{
    GMSH_ORDER_ONE = 1,
    GMSH_ORDER_TWO = 2,
    GMSH_ORDER_THREE = 3,
    GMSH_ORDER_FOUR = 4
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

    struct Factory
    {
        typedef Life::Singleton< Life::Factory< Gmsh, std::string > > type;
    };

    /** \name Constructors, destructor
     */
    //@{

    Gmsh()
        :
        _M_order( GMSH_ORDER_ONE ),
        M_version( 2 )
        {}
    Gmsh( Gmsh const & __g )
        :
        _M_order( __g._M_order ),
        M_version( __g.M_version )
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
                _M_order = __g._M_order;
                M_version = __g.M_version;
            }
            return *this;
        }

    //@}

    /** \name Accessors
     */
    //@{

    /**
     * get the order of the elements of the mesh
     * \return the order of the elements of the mesh
     */
    GMSH_ORDER order() const { return _M_order; }

    /**
     * @return the file format version
     */
    int version() const { return M_version; }
    //@}

    /** \name  Mutators
     */
    //@{

    /**
     * set the order of the elements of the mesh it can be either
     * GMSH_ORDER_ONE (order 1/linear) or GMSH_ORDER_TWO(order
     * 2/quadratic)
     *
     * \param o order of the elements
     */
    void setOrder( int o )
        {
            _M_order = (GMSH_ORDER) o;
        }

    /**
     * set the file format version
     */
    void setVersion( int version )
    {
        if ( version != 1 && version != 2 )
            throw std::invalid_argument( "invalid gmsh file format version" );
        M_version = version;
    }

    //@}

    /** \name  Methods
     */
    //@{

    std::string generateLine( std::string const& __name, double __h );

    std::string generateCube( std::string const& __name, double __h );

    std::string generateCircle( std::string const& __name, double __h );

    std::string generateSquare( std::string const& __name, double __h );

    /**
     * \param name  filename prefix to create the \c geo and \c msh file from \p geo
     * \param geo gmsh geometry description
     * \param forceRebuild if true, rebuild the mesh even if geofile is unchanged
     *        if false, rebuild only if geo file has changed.
     *        Useful if generateGeo has been called outside or if gmsh lybrary has changed.
     * \return the name of the mesh file generate by \c gmsh (with the \c .msh extension)
     */
    std::string generate( std::string const& name, std::string const& geo, bool const forceRebuild = false) const;

    //@}

protected:

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

private:
    GMSH_ORDER _M_order;
    int M_version;
};

}
#endif /* __Gmsh_H */
