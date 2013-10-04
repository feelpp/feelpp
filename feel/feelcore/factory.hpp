/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-10-04

  Copyright (C) 2009 Université de Grenoble 1
  Copyright (C) 2004 EPFL

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
   \file factory.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-10-04
 */
#ifndef __factory_H
#define __factory_H 1

#include <map>
#include <stdexcept>
#include <string>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/typeinfo.hpp>

namespace Feel
{
/*! \struct FactoryDefaultError
  Manages the "Unknown Type" error in an object Factory.
*/
template <
typename IdentifierType,
         class AbstractProduct
         >
struct FactoryDefaultError
{
    class Exception
        :
    public std::exception
    {
    public:
        Exception( IdentifierType id )
            :
            std::exception(),
            M_ex()
        {
            M_ex = this->getEx( id );

        }
        ~Exception() throw()
        {}
        const char* what() const throw ()
        {
            return M_ex.c_str();
        }
        std::string getEx( std::string const& id )
        {
            std::ostringstream __ex_str;
            __ex_str << "[Factory] Unknown Type : " << id;
            return __ex_str.str();
        }
        template<typename T>
        std::string getEx( T const& id )
        {
            std::ostringstream __ex_str;
            __ex_str << "[Factory] Unknown Type : ";
            return __ex_str.str();
        }
    private:
        std::string M_ex;
    };

    static AbstractProduct* onUnknownType( IdentifierType id )
    {
        throw Exception( id );
    }
};

/*!
  \class Factory
 *\ingroup Core
 *\brief Implements a generic object Factory

  \sa FactoryDefaultError, FactoryClone, TypeInfo

  @author Christophe Prud'homme
*/
template
<
class AbstractProduct,
      typename IdentifierType,
      typename ProductCreator = boost::function<AbstractProduct*()>,
      template<typename, class> class FactoryErrorPolicy = FactoryDefaultError
      >
class Factory
    :
public FactoryErrorPolicy<IdentifierType,AbstractProduct>
{
public:


    /** @name Typedefs
     */
    //@{
    typedef IdentifierType identifier_type;
    typedef AbstractProduct product_type;
    typedef ProductCreator creator_type;

    typedef FactoryErrorPolicy<identifier_type,product_type> super;

    //@}


    /** @name  Methods
     */
    //@{

    /**
     * Register a product.
     *
     * A product is composed of an identifier (typically a
     * std::string) and a functor that will create the associated
     * object.
     *
     * @param id identifier for the object to be registered
     * @param creator the functor that will create the registered
     * object
     *
     * @return true if registration went fine, false otherwise
     */
    bool registerProduct( const identifier_type& id, creator_type creator )
    {
        DVLOG(2) << "Registered type with id : " << id << "\n";
        return M_associations.insert( typename id_to_product_type::value_type( id, creator ) ).second;
    }

    /**
     * Unregister a product
     *
     * @param id
     * @sa registerProduct
     * @return true if unregistration went fine, false otherwise
     */
    bool unregisterProduct( const identifier_type& id )
    {
        DVLOG(2) << "Unregistered type with id : " << id << "\n";
        return M_associations.erase( id ) == 1;
    }

    /**
     * Create an object from a product registered in the Factory using
     * identifier \c id
     *
     * @param id identifier of the product to instantiate
     *
     * @return the object associate with \c id
     */
    product_type* createObject( const identifier_type& id )
    {
        typename id_to_product_type::const_iterator i = M_associations.find( id );

        if ( i != M_associations.end() )
        {
            DVLOG(2) << "Creating type with id : " << id << "\n";
            return ( i->second )();
        }

        DVLOG(2) << "Unknown type with id : " << id << "\n";
        return super::onUnknownType( id );
    }



    //@}
private:
    typedef std::map<identifier_type, creator_type> id_to_product_type;
    id_to_product_type M_associations;

};

/*!
  \class FactoryClone
 *\ingroup Core
 *\brief Implements a generic cloning object Factory

  \sa Factory, FactoryDefaultError

  \author Christophe Prud'homme
*/
template <
class AbstractProduct,
      class ProductCreator = boost::function<AbstractProduct* ( const AbstractProduct* )>,
      template<typename, class> class FactoryErrorPolicy = FactoryDefaultError
      >
class FactoryClone
    :
public FactoryErrorPolicy<TypeInfo, AbstractProduct>
{
public:


    /** @name Typedefs
     */
    //@{

    typedef FactoryErrorPolicy<TypeInfo,AbstractProduct> super;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    bool registerProduct( const TypeInfo& id, ProductCreator creator )
    {
        return M_associations.insert( typename id_to_product_type::value_type( id, creator ) ).second;
    }

    bool unregisterProduct( const TypeInfo& id )
    {
        return M_associations.erase( id ) == 1;
    }

    AbstractProduct* createObject( const AbstractProduct* model )
    {
        if ( model == 0 ) return 0;

        typename id_to_product_type::const_iterator i = M_associations.find( typeid( *model ) );

        if ( i != M_associations.end() )
        {
            return ( i->second )( model );
        }

        return super::onUnknownType( typeid( *model ) );
    }

    //@}

private:
    typedef std::map<TypeInfo, ProductCreator> id_to_product_type;
    id_to_product_type M_associations;
};


}
#endif /* __Factory_H */
