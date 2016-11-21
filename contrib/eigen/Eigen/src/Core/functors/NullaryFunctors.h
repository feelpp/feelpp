// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008-2010 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_NULLARY_FUNCTORS_H
#define EIGEN_NULLARY_FUNCTORS_H

namespace Eigen {

namespace internal {

template<typename Scalar>
struct scalar_constant_op {
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE scalar_constant_op(const scalar_constant_op& other) : m_other(other.m_other) { }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE scalar_constant_op(const Scalar& other) : m_other(other) { }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar operator() () const { return m_other; }
  template<typename PacketType>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const PacketType packetOp() const { return internal::pset1<PacketType>(m_other); }
  const Scalar m_other;
};
template<typename Scalar>
struct functor_traits<scalar_constant_op<Scalar> >
{ enum { Cost = 0 /* as the constant value should be loaded in register only once for the whole expression */,
         PacketAccess = packet_traits<Scalar>::Vectorizable, IsRepeatable = true }; };

template<typename Scalar> struct scalar_identity_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_identity_op)
  template<typename IndexType>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar operator() (IndexType row, IndexType col) const { return row==col ? Scalar(1) : Scalar(0); }
};
template<typename Scalar>
struct functor_traits<scalar_identity_op<Scalar> >
{ enum { Cost = NumTraits<Scalar>::AddCost, PacketAccess = false, IsRepeatable = true }; };

template <typename Scalar, typename Packet, bool RandomAccess, bool IsInteger> struct linspaced_op_impl;

// linear access for packet ops:
// 1) initialization
//   base = [low, ..., low] + ([step, ..., step] * [-size, ..., 0])
// 2) each step (where size is 1 for coeff access or PacketSize for packet access)
//   base += [size*step, ..., size*step]
//
// TODO: Perhaps it's better to initialize lazily (so not in the constructor but in packetOp)
//       in order to avoid the padd() in operator() ?
template <typename Scalar, typename Packet>
struct linspaced_op_impl<Scalar,Packet,/*RandomAccess*/false,/*IsInteger*/false>
{
  linspaced_op_impl(const Scalar& low, const Scalar& high, Index num_steps) :
    m_low(low), m_step(num_steps==1 ? Scalar() : (high-low)/Scalar(num_steps-1)),
    m_packetStep(pset1<Packet>(unpacket_traits<Packet>::size*m_step)),
    m_base(padd(pset1<Packet>(low), pmul(pset1<Packet>(m_step),plset<Packet>(-unpacket_traits<Packet>::size)))) {}

  template<typename IndexType>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar operator() (IndexType i) const
  { 
    m_base = padd(m_base, pset1<Packet>(m_step));
    return m_low+Scalar(i)*m_step; 
  }

  template<typename IndexType>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Packet packetOp(IndexType) const { return m_base = padd(m_base,m_packetStep); }

  const Scalar m_low;
  const Scalar m_step;
  const Packet m_packetStep;
  mutable Packet m_base;
};

// random access for packet ops:
// 1) each step
//   [low, ..., low] + ( [step, ..., step] * ( [i, ..., i] + [0, ..., size] ) )
template <typename Scalar, typename Packet>
struct linspaced_op_impl<Scalar,Packet,/*RandomAccess*/true,/*IsInteger*/false>
{
  linspaced_op_impl(const Scalar& low, const Scalar& high, Index num_steps) :
    m_low(low), m_step(num_steps==1 ? Scalar() : (high-low)/Scalar(num_steps-1)),
    m_lowPacket(pset1<Packet>(m_low)), m_stepPacket(pset1<Packet>(m_step)), m_interPacket(plset<Packet>(0)) {}

  template<typename IndexType>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar operator() (IndexType i) const { return m_low+i*m_step; }

  template<typename IndexType>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Packet packetOp(IndexType i) const
  { return internal::padd(m_lowPacket, pmul(m_stepPacket, padd(pset1<Packet>(Scalar(i)),m_interPacket))); }

  const Scalar m_low;
  const Scalar m_step;
  const Packet m_lowPacket;
  const Packet m_stepPacket;
  const Packet m_interPacket;
};

template <typename Scalar, typename Packet>
struct linspaced_op_impl<Scalar,Packet,/*RandomAccess*/true,/*IsInteger*/true>
{
  linspaced_op_impl(const Scalar& low, const Scalar& high, Index num_steps) :
    m_low(low), m_length(high-low), m_divisor(convert_index<Scalar>(num_steps==1?1:num_steps-1)), m_interPacket(plset<Packet>(0))
  {}

  template<typename IndexType>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  const Scalar operator() (IndexType i) const {
    return m_low + (m_length*Scalar(i))/m_divisor;
  }

  template<typename IndexType>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  const Packet packetOp(IndexType i) const {
    return internal::padd(pset1<Packet>(m_low), pdiv(pmul(pset1<Packet>(m_length), padd(pset1<Packet>(Scalar(i)),m_interPacket)),
                                                     pset1<Packet>(m_divisor))); }

  const Scalar m_low;
  const Scalar m_length;
  const Scalar  m_divisor;
  const Packet m_interPacket;
};

// ----- Linspace functor ----------------------------------------------------------------

// Forward declaration (we default to random access which does not really give
// us a speed gain when using packet access but it allows to use the functor in
// nested expressions).
template <typename Scalar, typename PacketType, bool RandomAccess = true> struct linspaced_op;
template <typename Scalar, typename PacketType, bool RandomAccess> struct functor_traits< linspaced_op<Scalar,PacketType,RandomAccess> >
{
  enum
  {
    Cost = 1,
    PacketAccess =   packet_traits<Scalar>::HasSetLinear
                  && ((!NumTraits<Scalar>::IsInteger) || packet_traits<Scalar>::HasDiv),
    IsRepeatable = true
  };
};
template <typename Scalar, typename PacketType, bool RandomAccess> struct linspaced_op
{
  linspaced_op(const Scalar& low, const Scalar& high, Index num_steps)
    : impl((num_steps==1 ? high : low),high,num_steps)
  {}

  template<typename IndexType>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar operator() (IndexType i) const { return impl(i); }

  template<typename Packet,typename IndexType>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Packet packetOp(IndexType i) const { return impl.packetOp(i); }

  // This proxy object handles the actual required temporaries, the different
  // implementations (random vs. sequential access) as well as the
  // correct piping to size 2/4 packet operations.
  // As long as we don't have a Bresenham-like implementation for linear-access and integer types,
  // we have to by-pass RandomAccess for integer types. See bug 698.
  const linspaced_op_impl<Scalar,PacketType,(NumTraits<Scalar>::IsInteger?true:RandomAccess),NumTraits<Scalar>::IsInteger> impl;
};

// Linear access is automatically determined from the operator() prototypes available for the given functor.
// If it exposes an operator()(i,j), then we assume the i and j coefficients are required independently
// and linear access is not possible. In all other cases, linear access is enabled.
// Users should not have to deal with this struture.
template<typename Functor> struct functor_has_linear_access { enum { ret = !has_binary_operator<Functor>::value }; };

} // end namespace internal

} // end namespace Eigen

#endif // EIGEN_NULLARY_FUNCTORS_H
