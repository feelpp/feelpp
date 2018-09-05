// Digit sequence 2-adic arithmetic

#ifndef _CL_2DS_H
#define _CL_2DS_H

namespace cln {

// div2adic(a_len,a_LSDptr,b_len,b_LSDptr,dest_LSDptr);
// dividiert die UDS a_LSDptr[-a_len..-1] mod 2^(intDsize*b_len)
// durch die ungerade UDS b_LSDptr[-b_len..-1] mod 2^(intDsize*b_len)
// (a_len >= b_len > 0) und liefert
// den Quotienten q als UDS dest_LSDptr[-b_len..-1] mod 2^(intDsize*b_len) und
// den "Rest" (a-b*q)/2^(intDsize*b_len) als UDS dest_LSDptr[-a_len..-b_len-1].
// Falls a_len > b_len, wird b implizit als durch Nullen fortgesetzt angenommen.
  extern void div2adic (uintC a_len, const uintD* a_LSDptr, uintC b_len, const uintD* b_LSDptr, uintD* dest_LSDptr);

// div2adic(len,a_LSDptr,b_LSDptr,dest_LSDptr);
// dividiert die UDS a_LSDptr[-len..-1] mod 2^(intDsize*len)
// durch die ungerade UDS b_LSDptr[-len..-1] mod 2^(intDsize*len) (len>0) und
// liefert den Quotienten als UDS dest_LSDptr[-len..-1] mod 2^(intDsize*len).
  inline void div2adic (uintC len, const uintD* a_LSDptr, const uintD* b_LSDptr, uintD* dest_LSDptr)
  { div2adic(len,a_LSDptr,len,b_LSDptr,dest_LSDptr); }

// recip2adic(len,a_LSDptr,dest_LSDptr);
// bildet den Kehrwert der ungeraden UDS a_LSDptr[-len..-1] mod 2^(intDsize*len)
// (len>0) und liefert sie als UDS dest_LSDptr[-len..-1] mod 2^(intDsize*len).
  extern void recip2adic (uintC len, const uintD* a_LSDptr, uintD* dest_LSDptr);

}  // namespace cln

#endif /* _CL_2DS_H */
