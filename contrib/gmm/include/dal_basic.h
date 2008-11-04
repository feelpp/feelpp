/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_basic.h : Basic container, addition of space by mean of  */
/*                          elements pack allocation.                      */
/*                                                                         */
/* Date :  June 01, 1995.                                                  */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1995-2001  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */


#ifndef DAL_BASIC_H__
#define DAL_BASIC_H__

/* *********************************************************************** */
/* Remarks for future improvements:                                        */
/*									   */
/* - For the moment, the allocation is made by new and delete[].           */
/*									   */
/* *********************************************************************** */
#include <limits.h>

#include <dal_std.h>
#include <vector>
#include <dal_algobase.h>

namespace dal
{
  template<class T, unsigned char pks = 5> class dynamic_array;

  /// Iterator class for dynamic array.
  template<class T, unsigned char pks> struct dna_iterator {
    typedef T             value_type;
    typedef value_type*   pointer;
    typedef value_type&   reference;
    typedef size_t        size_type;
    typedef ptrdiff_t     difference_type;
    typedef std::random_access_iterator_tag iterator_category;

#   define DNAMPKS__ ((size_type(1) << pks) - 1)
    dynamic_array<T,pks> *p;
    size_type in;
    pointer pT;

    dna_iterator(void) {}
    dna_iterator(dynamic_array<T,pks> &da, size_type ii)
      { p = &da; in = ii; pT = p->pt_to(in); }

    inline size_type index(void) const { return in; }
    /// next element.
    dna_iterator operator ++(int) {
      dna_iterator tmp = *this;
      if ((++in)&DNAMPKS__) pT++; else pT=p->pt_to(in); return tmp;
    }
    /// previous element.
    dna_iterator operator --(int) {
      dna_iterator tmp = *this;
      if ((in--)&DNAMPKS__) pT--; else pT=p->pt_to(in); return tmp;
    }
    /// next element.
    dna_iterator &operator ++()
      { if ((++in)&DNAMPKS__) pT++; else pT=p->pt_to(in); return *this; }
    /// previous element.
    dna_iterator &operator --()
      { if ((in--)&DNAMPKS__) pT--; else pT=p->pt_to(in); return *this; }
    /// go i elements forward.
    dna_iterator &operator +=(difference_type i)
      { in += i; pT=p->pt_to(in); return *this; }
    /// go i elements backward.
    dna_iterator &operator -=(difference_type i)
      { in -= i; pT=p->pt_to(in); return *this; }
    /// gives an iterator pointing i elements forward.
    dna_iterator operator +(difference_type i) const
      { dna_iterator it = *this; return (it += i); }
    /// gives an iterator pointing i elements backward.
    dna_iterator operator -(difference_type i) const
      { dna_iterator it = *this; return (it -= i); }
    /// Gives the difference, in term of elements between two iterators.
    difference_type operator -(const dna_iterator &i) const
      { return difference_type(in - i.in); }

    reference operator  *() const { return (*pT); }
    pointer   operator ->() const { return pT;    }
    reference operator [](size_type ii) const { return (*p)[in+ii]; }

    bool operator ==(const dna_iterator &i) const { return (i.in==in);}
    bool operator !=(const dna_iterator &i) const { return (i.in!=in);}
    bool operator < (const dna_iterator &i) const { return ( in<i.in);}
  };

  /// Constant iterator class for dynamic array.
  template<class T, unsigned char pks> struct dna_const_iterator {
    typedef T                  value_type;
    typedef const value_type*  pointer;
    typedef const value_type&  reference;
    typedef size_t             size_type;
    typedef ptrdiff_t          difference_type;
    typedef std::random_access_iterator_tag iterator_category;

#   define DNAMPKS__ ((size_type(1) << pks) - 1)
    const dynamic_array<T,pks> *p;
    size_type in;
    pointer pT;

    dna_const_iterator(void) {}
    dna_const_iterator(const dynamic_array<T,pks> &da, size_type ii)
      { p = &da; in = ii; pT = p->pt_to(in); }
    dna_const_iterator(const dna_iterator<T, pks> &it)
      : p(it.p), in(it.in), pT(it.pT) {}

    inline size_type index(void) const { return in; }
    dna_const_iterator operator ++(int) {
      dna_const_iterator tmp = *this;
      if ((++in)&DNAMPKS__) pT++; else pT=p->pt_to(in); return tmp;
    }
    dna_const_iterator operator --(int) {
      dna_const_iterator tmp = *this;
      if ((in--)&DNAMPKS__) pT--; else pT=p->pt_to(in); return tmp;
    }
    dna_const_iterator &operator ++()
      { if ((++in)&DNAMPKS__) pT++; else pT=p->pt_to(in); return *this; }
    dna_const_iterator &operator --()
      { if ((in--)&DNAMPKS__) pT--; else pT=p->pt_to(in); return *this; }
    dna_const_iterator &operator +=(difference_type i)
      { in += i; pT=p->pt_to(in); return *this; }
    dna_const_iterator &operator -=(difference_type i)
      { in -= i; pT=p->pt_to(in); return *this; }
    dna_const_iterator operator +(difference_type i) const
      { dna_const_iterator it = *this; return (it += i); }
    dna_const_iterator operator -(difference_type i) const
      { dna_const_iterator it = *this; return (it -= i); }
    difference_type operator -(const dna_const_iterator &i) const
      { return difference_type(in - i.in); }

    reference operator  *() const { return (*pT); }
    pointer   operator ->() const { return pT;    }
    reference operator [](size_type ii) const { return (*p)[in+ii]; }

    bool operator ==(const dna_const_iterator &i) const
      { return (i.in == in); }
    bool operator !=(const dna_const_iterator &i) const
      { return (i.in != in); }
    bool operator < (const dna_const_iterator &i) const
      { return (in < i.in); }
  };

  /**  Dynamic Array. Defines the basic container of the library which is
   *  dal::dynamic\_array$<$T, pks$>$. This container is virtually an
   * infinite array of element of type T. When a random acces tab[i] is
   * called, a control is made on i and an allocation is made if
   * needed. The allocation is made by blocks of n elements, where
   * $n = 2^{pks}$. $pks$ is an optional parameter assumed to be 5.
   * \subsubsection*{Example of code}
   *  If T is any type (with or without trivial constructor/destructor,
   *  and with constructor T(0) and T(1)), the
   *  following code is valid: \\ \\
   * { \tt
   *  \#include<dal_basic.h> \\ \\
   *  dal::dynamic_array<T> tab; \\
   *  tab[50] = T(0); // to be sure to have at least 50 elements \\
   *  std::fill(tab.begin(), tab.end(), T(0)); // at least 50 elements
   *  are initialized \\ \\
   *  dal::dynamic_array<T>::iterator it = tab.begin();
   *  dal::dynamic_array<T>::iterator ite = it + 50;  \\
   *  for( ; it != ite; ++it)
   *    { *it = T(1); } // only the 50 first elements are changed.
   * }
   */
  template<class T, unsigned char pks> class dynamic_array {
  public :

    typedef T                    value_type;
    typedef value_type*          pointer;
    typedef const value_type*    const_pointer;
    typedef value_type&          reference;
    typedef const value_type&    const_reference;
    typedef size_t               size_type;
    typedef ptrdiff_t            difference_type;
    typedef unsigned char        pack_size_type;
    typedef std::vector<pointer> pointer_array;
    typedef dna_iterator<T, pks> iterator;
    typedef dna_const_iterator<T, pks> const_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;

  protected :

#   define DNAMPKS__ ((size_type(1) << pks) - 1)
    pointer_array array;
    pack_size_type ppks;   /* size of pointer packs (2^ppks).            */
    size_type m_ppks;      /* = (2^ppks) - 1.                            */
    size_type last_ind;    /* allocated = 0 .. last_ind-1.               */
    size_type last_accessed; /* valid = 0 .. last_accessed-1.            */

  public :

    /// Number of allocated elements.
    size_type size(void) const { return last_accessed; }
    size_type capacity(void) const { return last_ind; }
    size_type max_size(void) const { return (size_type(-1)) / 2; }
    /// True if no space is allocated.
    bool empty(void) const { return last_accessed == 0; }
    /// Iterator on the first element.
    iterator begin(void) { return iterator(*this, 0); }
    /// Constant iterator on the first element.
    const_iterator begin(void) const { return const_iterator(*this, 0); }
    /// Iterator on the last + 1 element.
    iterator end(void) { return iterator(*this, size()); }
    /// Constant iterator on the last + 1 element.
    const_iterator end(void) const
      { return const_iterator(*this, size()); }
    reverse_iterator rbegin(void) { return reverse_iterator(end()); }
    const_reverse_iterator rbegin(void) const
      { return const_reverse_iterator(end()); }
    reverse_iterator rend(void) { return reverse_iterator(begin()); }
    const_reverse_iterator rend(void) const
      { return const_reverse_iterator(begin()); }

    reference front(void) { return *begin(); }
    const_reference front(void) const { return *begin(); }
    reference back(void) { return *(end() - 1); }
    const_reference back(void) const { return *(end() - 1); }

    void swap(dynamic_array<T,pks> &da);
    /// Clear and desallocate all the elements.
    void clear(void);
    dynamic_array<T,pks> &operator =(const dynamic_array<T,pks> &da);

  protected:

    void init(void)
      { last_accessed = last_ind = 0; array.resize(8); ppks = 3; m_ppks = 7; }


  public:

    dynamic_array(const dynamic_array<T,pks> &da) { init(); *this = da; }
    dynamic_array(void) { init(); }
    ~dynamic_array(void) { clear(); }
    inline pointer pt_to(size_type ii) /* used by iterators.             */
      { return (ii >=last_ind) ? NULL : &((array[ii>>pks])[ii&DNAMPKS__]); }
    inline const_pointer pt_to(size_type ii) const
      { return (ii >=last_ind) ? NULL : &((array[ii>>pks])[ii&DNAMPKS__]); }

    /// Gives a constant reference on element ii.
    const_reference operator [](size_type ii) const;
    /// Gives a reference on element ii.
    reference operator [](size_type ii);
    void resize(size_type i) { (*this)[i-1]; }

    /// Gives the total memory occupied by the array.
    size_type memsize(void) const {
      return sizeof(pointer) * array.capacity()
	+ last_ind*sizeof(T) + sizeof(dynamic_array<T,pks>);
    }

    /// Swap element i1 and i2.
    void swap(size_type i1, size_type i2)
      { std::swap((*this)[i1], (*this)[i2]); }
  };


  /* ********************************************************************* */
  /* Menbers functions							   */
  /* ********************************************************************* */


  template<class T, unsigned char pks>
  void dynamic_array<T,pks>::swap(dynamic_array<T,pks> &da) {
    array.swap(da.array); std::swap(last_ind, da.last_ind);
    std::swap(ppks, da.ppks); std::swap(m_ppks, da.m_ppks);
  }

  template<class T, unsigned char pks>
  void dynamic_array<T,pks>::clear(void) {
    typename pointer_array::iterator it  = array.begin();
    typename pointer_array::iterator ite = it+ ((last_ind + DNAMPKS__) >> pks);
    while (it != ite) delete[] *it++;
    array.clear(); init();
  }

  template<class T, unsigned char pks> dynamic_array<T,pks>
  &dynamic_array<T,pks>::operator = (const dynamic_array<T,pks> &da) {
    clear(); /* evitable ... ? */
    array.resize(da.array.size());
    last_ind = da.last_ind;
    last_accessed = da.last_accessed;
    ppks = da.ppks; m_ppks = da.m_ppks;
    typename pointer_array::iterator it = array.begin();
    typename pointer_array::const_iterator ita = da.array.begin();
    typename pointer_array::iterator ite = it+ ((last_ind + DNAMPKS__) >> pks);
    while (it != ite) {
      register pointer p = *it++ = new T[DNAMPKS__+1];
      register pointer pe = p + (DNAMPKS__+1);
      register const_pointer pa = *ita++;
      while (p != pe) *p++ = *pa++;
    }
    return *this;
  }

  template<class T, unsigned char pks>
    typename dynamic_array<T,pks>::const_reference
      dynamic_array<T,pks>::operator [](size_type ii) const {
    static T *f = NULL;
    if (f == NULL) { f = new T(); }
    return (ii<last_ind) ? (array[ii>>pks])[ii&DNAMPKS__] : *f;
  }

  template<class T, unsigned char pks> typename dynamic_array<T,pks>::reference
    dynamic_array<T,pks>::operator [](size_type ii) {
    if (ii >= last_accessed) {
      if (ii >= INT_MAX) {
	DAL_THROW(std::out_of_range, "index" << long(ii) << " out of range.");
      }

      last_accessed = ii + 1;
      if (ii >= last_ind) {
	if ((ii >> (pks+ppks)) > 0) {
	  while ((ii >> (pks+ppks)) > 0) ppks++;
	  array.resize(m_ppks = (size_type(1) << ppks)); m_ppks--;
	}
	for (size_type jj = (last_ind >> pks); ii >= last_ind;
	     jj++, last_ind += (DNAMPKS__ + 1))
	  { array[jj] = new T[DNAMPKS__ + 1]; }
      }
    }
    return (array[ii >> pks])[ii & DNAMPKS__];
  }


  /* ********************************************************************* */
  /* Templates functions					      	   */
  /* ********************************************************************* */

  template<class T, unsigned char pks>
    bool operator==(const dynamic_array<T,pks> &x,
		    const dynamic_array<T,pks> &y) {
    typename dynamic_array<T,pks>::const_iterator itxb=x.begin(), itxe=x.end();
    typename dynamic_array<T,pks>::const_iterator ityb=y.begin(), itye=y.end();
    typename dynamic_array<T,pks>::size_type d = std::min(itxe-itxb,itye-ityb);
    typename dynamic_array<T,pks>::const_iterator itxc = itxb+d, ityc = ityb+d;

    if  (!std::equal(itxb, itxc, ityb)) return false;
    for (; itxc != itxe; itxc++) if (*itxc != T()) return false;
    for (; ityc != itye; ityc++) if (*ityc != T()) return false;
    return true;
  }

  template<class T, unsigned char pks>
    bool operator < (const dynamic_array<T,pks> &x,
		     const dynamic_array<T,pks> &y)
  { return std::lexicographical_compare(x.begin(),x.end(),y.begin(),y.end()); }

  template<class T, unsigned char pks> inline
    void swap(const dynamic_array<T,pks> &x, const dynamic_array<T,pks> &y)
  { x.swap(y); }


  /**
      a very simple pointer collection
      which destroys the content of its pointers
  */
  template<class T> class ptr_collection : public std::vector<T*> {
  public:
    typedef typename std::vector<T*>::size_type size_type;
    typedef typename std::vector<T*>::iterator iterator;
    typedef typename std::vector<T*>::const_iterator const_iterator;
    ptr_collection() : std::vector<T*>() {}
    ptr_collection(size_type sz) : std::vector<T*>(sz) { std::fill(this->begin(),this->end(),0); }
    ~ptr_collection() { for (iterator i=this->begin(); i != this->end(); ++i) { if (*i) delete *i; *i = 0; } }
  };

}
#endif /* DAL_BASIC_H__  */
