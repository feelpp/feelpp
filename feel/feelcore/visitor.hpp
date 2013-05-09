/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-03-21

  Copyright (C) 2007,2009 Université de Grenoble 1
  Copyright (C) 2005,2006,2009 EPFL

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
   \file visitor.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-03-21
 */
#ifndef __Visitor_H
#define __Visitor_H 1

#include <boost/mpl/assert.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/pop_front.hpp>

#include <feel/feelcore/feel.hpp>

namespace Feel
{
namespace mpl=boost::mpl;
/**
 * \class VisitorBase
 *\ingroup Core
 *\brief The base class of any Acyclic Visitor
 *
 *  @author Christophe Prud'homme <Christophe.Prudhomme@ann.jussieu.fr>
 *  @see Gamma, Helm, Johnson, Vlissides,   Design Patterns  Pub: Addison Wesley
 *  @see Moderner C++ Design by Andrei Alexandrescu Pub: Addison Wesley
 */
class VisitorBase
{
public:
    virtual ~VisitorBase() {}
};

/**
 * \class Visitor
 *
 *  This class is the base class to implement the Visitor Pattern.
 *  Let's A be a visitor class to a class B:
 *  \code
 *  class A: public Visitor<B>
 *  {
 *  public:
 *  ...
 *  void visit(B * b)
 *  {
 *  ..implement the functionnality to be added to a B class..
 *  }
 *  };
 *  class B:
 *  {
 *  public:
 *  void accept(Visitor<B>* b){  b->visit(this); }
 *  };
 *  \endcode
 *
 * @author Christophe Prud'homme <Christophe.Prudhomme@ann.jussieu.fr>
 * @see Gamma, Helm, Johnson, Vlissides,
 *  Design Patterns
 *  Pub: Addison Wesley
 */
template<
class T,
      typename R = void
      >
class Visitor
{
public:

    /** @name Typedefs
     */
    //@{
    typedef R return_type;

    //@}

    /** @name Constructors, Destructors and methods
     */
    //@{

    //! virtual base destructor
    virtual ~Visitor()
    {
        // do nothing here
    }

    //! visit a data structure
    virtual return_type visit( T* ) = 0;

    //! visit a data structure
    return_type visit( T& __t )
    {
        return visit( &__t );
    }

    //@}
};

template<
class TList,
      typename R = void
      >
class VisitorList
    :
public Visitor<typename mpl::front<TList>::type, R>,
public mpl::if_<mpl::greater<mpl::size<TList>,mpl::long_<1l> >,
    mpl::identity<VisitorList<mpl::pop_front<TList>,R> >,
    mpl::identity<mpl::identity<VisitorBase> > >::type::type
{
#if 0
    typedef typename mpl::if_<mpl::equal_to<mpl::size<TList>,mpl::int_<2> >,
            mpl::identity<A>,
            mpl::identity<C> >::type::type the_type;
    //BOOST_MPL_ASSERT( ( boost::is_same<typename mpl::size<TList>::type, mpl::long_<2l> > ) );
    //BOOST_MPL_ASSERT( ( boost::is_same<the_type, A> ) );
#endif
};


#if 0
/**
   \class Visitor
  *\ingroup Core
  *\brief The base class of any Acyclic Visitor

   This specialization is not present in the book. It makes it
   easier to define Visitors for multiple types in a shot by using
   a typelist. Example:

   \code
   class SomeVisitor :
   public VisitorBase // required
   public Visitor<TYPELIST_2(RasterBitmap, Paragraph)>,
   public Visitor<Paragraph>
   {
   public:
   void visit(RasterBitmap&); // visit a RasterBitmap
   void visit(Paragraph &);   // visit a Paragraph
   };
   \endcode

   @author Christophe Prud'homme <Christophe.Prudhomme@ann.jussieu.fr>
   @see Gamma, Helm, Johnson, Vlissides,   Design Patterns  Pub: Addison Wesley
   @see Modern C++ Design by Andrei Alexandrescu Pub: Addison Wesley
*/
template <
class Head,
      class Tail,
      typename R
      >
class Visitor< mpl::list<Head, Tail>, R >
    :
public Visitor<Head, R>,
public Visitor<Tail, R>
{
public:
    typedef R return_type;
    using Visitor<Head, R>::visit;
    using Visitor<Tail, R>::visit;
};

template <
class Head,
      typename R
      >
class Visitor< mpl::list<Head>, R>
    :
public Visitor<Head, R>
{
public:
    typedef R return_type;
    using Visitor<Head, R>::visit;
};
#endif
/**
   
  *\ingroup Core
  *\brief
   Implements non-strict visitation (you can implement only part of the Visit
   functions)
*/
template <
class TList,
      typename R = void
      >
class VisitorBaseImpl;


template <
class Head,
      class Tail,
      typename R
      >
class VisitorBaseImpl< mpl::list<Head, Tail>, R >
    :
public Visitor<Head, R>,
public VisitorBaseImpl<Tail, R>
{
public:
    // using BaseVisitorImpl<Tail, R>::Visit;

    virtual R visit( Head* )
    {
        return R();
    }

};

template <
class Head,
      typename R
      >
class VisitorBaseImpl< mpl::list<Head>, R >
    :
public Visitor<Head, R>
{
public:
    virtual R visit( Head* )
    {
        return R();
    }
};


/**
 * \class VisitableCatchAllDefault
 */
template <
typename R,
         typename Visited
         >
struct VisitableCatchAllDefault
{
    static R onUnknownVisitor( Visited&, VisitorBase& )
    {
        return R();
    }

    static R onUnknownVisitor( Visited*, VisitorBase* )
    {
        return R();
    }

};

/**
 * \class VisitableBase
 */
template
<
typename R = void,
         template <class, class> class CatchAll = VisitableCatchAllDefault
         >
class VisitableBase
{
public:

    typedef R return_type;

    virtual ~VisitableBase() {}

    //! accept visitor: use S_DEFINE_VISITABLE() to redefine it
    virtual return_type accept( VisitorBase& ) = 0;

    //! accept visitor: use S_DEFINE_VISITABLE() to redefine it
    virtual return_type accept( VisitorBase* ) = 0;

protected:

    /**
     * give access only to the hierarchy
     *
     * @return the return type object of the visitor
     */
    template <class T>
    static
    return_type
    acceptImpl( T* visited, VisitorBase* guest )
    {
        // Apply the Acyclic Visitor
        if ( Visitor<T>* p = dynamic_cast< Visitor<T>* >( guest ) )
        {
            return p->visit( visited );
        }

        return CatchAll<R, T>::onUnknownVisitor( visited, guest );
    }
};

/**
 * \def FEELPP_DEFINE_VISITABLE()
 *
 * Put it in every class that you want to make visitable (in
 * addition to deriving it from VisitableBase<R>
 */
#define FEELPP_DEFINE_VISITABLE()                          \
    virtual return_type accept( VisitorBase& guest )         \
    {                                                        \
        return this->acceptImpl( this, &guest );             \
    }                                                        \
    virtual return_type accept( VisitorBase* guest )         \
    {                                                        \
        return this->acceptImpl( this, guest );              \
    }

/**
 * \class VisitorCyclic
 *   Put it in every class that you want to make visitable (in addition to
 *   deriving it from VisitableBase<R>
 */
template <
typename R,
         class TList
         >
class VisitorCyclic
    :
public Visitor<TList, R>
{
public:
    typedef R return_type;
    // using Visitor<TList, R>::Visit;

    template <class Visited>
    return_type genericVisit( Visited* host )
    {
        Visitor<Visited, return_type>& subObj = *this;
        return subObj.visit( host );
    }

    template <class Visited>
    return_type genericVisit( Visited& host )
    {
        Visitor<Visited, return_type>& subObj = *this;
        return subObj.visit( host );
    }

};

/**
 * \def FEELPP_DEFINE_CYCLIC_VISITABLE(SomeVisitor)
 *
 * Put it in every class that you want to make visitable by a cyclic visitor
 */
#define FEELPP__DEFINE_CYCLIC_VISITABLE(SomeVisitor)                    \
  virtual SomeVisitor::return_type Accept(SomeVisitor& guest)     \
  {                                                               \
    return guest.genericVisit(*this);                             \
  }                                                               \
  virtual SomeVisitor::return_type Accept(SomeVisitor* guest)     \
  {                                                               \
    return guest->genericVisit(*this);                            \
  }

}
#endif /* __Visitor_H */
