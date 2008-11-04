/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-04-07

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file flags.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-04-07
 */
#ifndef __Flags_H
#define __Flags_H 1

namespace Life
{
/// \cond detail
namespace detail
{

class Flag
{
    int i;
public:
    inline Flag(int i);
    inline operator int() const { return i; }
};

inline Flag::Flag(int ai) : i(ai) {}


template<typename Enum>
class Flags
{
    typedef void **Zero;
    int i;
public:
    typedef Enum enum_type;
    inline Flags(const Flags &f) : i(f.i) {}
    inline Flags(Enum f) : i(f) {}
    inline Flags(Zero = 0) : i(0) {}
    inline Flags(Flag f) : i(f) {}

    inline Flags &operator=(const Flags &f) { i = f.i; return *this; }
    inline Flags &operator&=(int mask) { i &= mask; return *this; }
    inline Flags &operator&=(uint mask) { i &= mask; return *this; }
    inline Flags &operator|=(Flags f) { i |= f.i; return *this; }
    inline Flags &operator|=(Enum f) { i |= f; return *this; }
    inline Flags &operator^=(Flags f) { i ^= f.i; return *this; }
    inline Flags &operator^=(Enum f) { i ^= f; return *this; }

    inline operator int() const { return i; }

    inline Flags operator|(Flags f) const { Flags g; g.i = i | f.i; return g; }
    inline Flags operator|(Enum f) const { Flags g; g.i = i | f; return g; }
    inline Flags operator^(Flags f) const { Flags g; g.i = i ^ f.i; return g; }
    inline Flags operator^(Enum f) const { Flags g; g.i = i ^ f; return g; }
    inline Flags operator&(int mask) const { Flags g; g.i = i & mask; return g; }
    inline Flags operator&(uint mask) const { Flags g; g.i = i & mask; return g; }
    inline Flags operator&(Enum f) const { Flags g; g.i = i & f; return g; }
    inline Flags operator~() const { Flags g; g.i = ~i; return g; }

    inline bool operator!() const { return !i; }

    inline bool testFlag(Enum f) const { return i & f; }

};

} // detail
/// \endcond


#define LIFE_DECLARE_FLAGS(Flags, Enum)         \
    typedef detail::Flags<Enum> Flags;

} // Life
#endif /* __Flags_H */
