/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_ref.h : Structures which refere to containers.           */
/*     									   */
/*                                                                         */
/* Date : August 26, 2000.                                                 */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2000-2002  Yves Renard.                                   */
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


#ifndef DAL_REF_H__
#define DAL_REF_H__

/* *********************************************************************** */
/* WARNING : modifiying the container infirm the validity of references.   */
/* *********************************************************************** */


#include <iterator>
#include <dal_basic.h>

namespace dal
{

  /* ********************************************************************* */
  /* Simple reference.                                                     */
  /* ********************************************************************* */

  template<typename ITER> class tab_ref
  {
    protected :

      ITER begin_, end_;

    public :

      typedef typename std::iterator_traits<ITER>::value_type  value_type;
      typedef typename std::iterator_traits<ITER>::pointer     pointer;
      typedef typename std::iterator_traits<ITER>::pointer     const_pointer;
      typedef typename std::iterator_traits<ITER>::reference   reference;
      typedef typename std::iterator_traits<ITER>::reference   const_reference;
      typedef typename std::iterator_traits<ITER>::difference_type
	                                                       difference_type;
      typedef ITER                            iterator;
      typedef ITER                            const_iterator;
      typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
      typedef std::reverse_iterator<iterator> reverse_iterator;
      typedef size_t size_type;
    
      bool empty(void) const { return begin_ == end_; }
      size_type size(void) const { return end_ - begin_; }

      const iterator &begin(void) { return begin_; }
      const const_iterator &begin(void) const { return begin_; }
      const iterator &end(void) { return end_; }
      const const_iterator &end(void) const { return end_; }
      reverse_iterator rbegin(void) { return reverse_iterator(end()); }
      const_reverse_iterator rbegin(void) const
      { return const_reverse_iterator(end()); }
      reverse_iterator rend(void) { return reverse_iterator(begin()); }
      const_reverse_iterator rend(void) const
      { return const_reverse_iterator(begin()); }

      reference front(void) { return *begin(); }
      const_reference front(void) const { return *begin(); }
      reference back(void) { return *(--(end())); }
      const_reference back(void) const { return *(--(end())); }
      void pop_front(void) { ++begin_; }

      const_reference operator [](size_type ii) const { return *(begin_ + ii);}
      reference operator [](size_type ii) { return *(begin_ + ii); }

      tab_ref(void) {}
      tab_ref(const ITER &b, const ITER &e) : begin_(b), end_(e) {}
  };


  /* ********************************************************************* */
  /* Reference with index.                                                 */
  /* ********************************************************************* */

//   template<typename ITER> struct tab_ref_index_iterator_
//     : public dynamic_array<size_t>::const_iterator
//   {
//     typedef typename std::iterator_traits<ITER>::value_type  value_type;
//     typedef typename std::iterator_traits<ITER>::pointer     pointer;
//     typedef typename std::iterator_traits<ITER>::reference   reference;
//     typedef typename std::iterator_traits<ITER>::difference_type  
//     difference_type;
//     typedef std::random_access_iterator_tag iterator_category;
//     typedef size_t size_type;
//     typedef dynamic_array<size_type>::const_iterator dnas_iterator_;
//     typedef tab_ref_index_iterator_<ITER> iterator;
    

//     ITER piter;
    
//     iterator operator ++(int)
//     { iterator tmp = *this; ++(*((dnas_iterator_ *)(this))); return tmp; }
//     iterator operator --(int)
//     { iterator tmp = *this; --(*((dnas_iterator_ *)(this))); return tmp; }
//     iterator &operator ++()
//     { ++(*((dnas_iterator_ *)(this))); return *this; }
//     iterator &operator --()
//     { --(*((dnas_iterator_ *)(this))); return *this; }
//     iterator &operator +=(difference_type i)
//     { (*((dnas_iterator_ *)(this))) += i; return *this; }
//     iterator &operator -=(difference_type i)
//     { (*((dnas_iterator_ *)(this))) -= i; return *this; }
//     iterator operator +(difference_type i) const
//     { iterator it = *this; return (it += i); }
//     iterator operator -(difference_type i) const
//     { iterator it = *this; return (it -= i); }
//     difference_type operator -(const iterator &i) const
//     { return *((dnas_iterator_ *)(this)) - *((dnas_iterator_ *)(&i)); }
	
//     reference operator *() const
//     { return *(piter + *((*((dnas_iterator_ *)(this))))); }
//     reference operator [](int ii)
//     { return *(piter + *((*((dnas_iterator_ *)(this+ii))))); }
    
//     bool operator ==(const iterator &i) const
//     { 
//       return ((piter) == ((i.piter))
//        && *((dnas_iterator_ *)(this)) == *((*((dnas_iterator_ *)(this)))));
//     }
//     bool operator !=(const iterator &i) const
//     { return !(i == *this); }
//     bool operator < (const iterator &i) const
//     { 
//       return ((piter) == ((i.piter))
// 	 && *((dnas_iterator_ *)(this)) < *((*((dnas_iterator_ *)(this)))));
//     }

//     tab_ref_index_iterator_(void) {}
//     tab_ref_index_iterator_(const ITER &iter, const dnas_iterator_ &dnas_iter)
//       : dnas_iterator_(dnas_iter), piter(iter) {}
//   };


//   template<typename ITER> class tab_ref_index
//   {
//     public :

//       typedef typename std::iterator_traits<ITER>::value_type value_type;
//       typedef typename std::iterator_traits<ITER>::pointer    pointer;
//       typedef typename std::iterator_traits<ITER>::pointer    const_pointer;
//       typedef typename std::iterator_traits<ITER>::reference  reference;
//       typedef typename std::iterator_traits<ITER>::reference  const_reference;
//       typedef typename std::iterator_traits<ITER>::difference_type
// 	                                                       difference_type;
//       typedef size_t size_type; 
//       typedef tab_ref_index_iterator_<ITER> iterator;
//       typedef iterator                          const_iterator;
//       typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
//       typedef std::reverse_iterator<iterator> reverse_iterator;
    
//     protected :

//       ITER begin_;
//       dynamic_array<size_type> index_;

//     public :

//       bool empty(void) const { return index_.empty(); }
//       size_type size(void) const { return index_.size(); }


//       iterator begin(void) { return iterator(begin_, index_.begin()); }
//       const_iterator begin(void) const
//       { return iterator(begin_, index_.begin()); }
//       iterator end(void) { return iterator(begin_, index_.end()); }
//       const_iterator end(void) const { return iterator(begin_, index_.end()); }
//       reverse_iterator rbegin(void) { return reverse_iterator(end()); }
//       const_reverse_iterator rbegin(void) const
//       { return const_reverse_iterator(end()); }
//       reverse_iterator rend(void) { return reverse_iterator(begin()); }
//       const_reverse_iterator rend(void) const
//       { return const_reverse_iterator(begin()); }


//       reference front(void) { return *(begin_ +index_[0]); }
//       const_reference front(void) const { return *(begin_ +index_[0]); }
//       reference back(void) { return *(--(end())); }
//       const_reference back(void) const { return *(--(end())); }
   
//       tab_ref_index(void) {}
//       tab_ref_index(const ITER &b, const dynamic_array<size_type> &ind)
//       { begin_ = b; index_ = ind; }

//     // to be changed in a const_reference ?
//       value_type operator [](size_type ii) const
//       { return *(begin_ + index_[ii]);}
//       reference operator [](size_type ii) { return *(begin_ + index_[ii]); }

//   };


  /* ********************************************************************* */
  /* Reference with reference on index.                                    */
  /* ********************************************************************* */

  template<typename ITER, typename ITER_INDEX>
    struct tab_ref_index_ref_iterator_
    {
      typedef typename std::iterator_traits<ITER>::value_type value_type;
      typedef typename std::iterator_traits<ITER>::pointer    pointer;
      typedef typename std::iterator_traits<ITER>::reference  reference;
      typedef typename std::iterator_traits<ITER>::difference_type
                                                              difference_type;
      typedef std::random_access_iterator_tag iterator_category;
      typedef tab_ref_index_ref_iterator_<ITER, ITER_INDEX> iterator;
      typedef size_t size_type;

      ITER piter;
      ITER_INDEX iter_index;
      
      iterator operator ++(int)
      { iterator tmp = *this; ++iter_index; return tmp; }
      iterator operator --(int)
      { iterator tmp = *this; --iter_index; return tmp; }
      iterator &operator ++() { ++iter_index; return *this; }
      iterator &operator --() { --iter_index; return *this; }
      iterator &operator +=(difference_type i)
      { iter_index += i; return *this; }
      iterator &operator -=(difference_type i)
      { iter_index -= i; return *this; }
      iterator operator +(difference_type i) const
      { iterator it = *this; return (it += i); }
      iterator operator -(difference_type i) const
      { iterator it = *this; return (it -= i); }
      difference_type operator -(const iterator &i) const
      { return iter_index - i.iter_index; }
	
      reference operator *() const
      { return *(piter + *iter_index); }
      reference operator [](int ii) const
      { return *(piter + *(iter_index+ii)); }
      
      bool operator ==(const iterator &i) const
      { return ((piter) == ((i.piter)) && iter_index == i.iter_index); }
      bool operator !=(const iterator &i) const { return !(i == *this); }
      bool operator < (const iterator &i) const
      { return ((piter) == ((i.piter)) && iter_index < i.iter_index); }

      tab_ref_index_ref_iterator_(void) {}
      tab_ref_index_ref_iterator_(const ITER &iter, 
				  const ITER_INDEX &dnas_iter)
	: piter(iter), iter_index(dnas_iter) {}
      
    };

  /** 
      convenience template function for quick obtention of a indexed iterator
      without having to specify its (long) typename
  */
  template<typename ITER, typename ITER_INDEX>
  tab_ref_index_ref_iterator_<ITER,ITER_INDEX>
  index_ref_iterator(ITER it, ITER_INDEX it_i) {
    return tab_ref_index_ref_iterator_<ITER,ITER_INDEX>(it, it_i);
  }

  template<typename ITER, typename ITER_INDEX> class tab_ref_index_ref {
  public :
    
    typedef std::iterator_traits<ITER>            traits_type;
    typedef typename traits_type::value_type      value_type;
    typedef typename traits_type::pointer         pointer;
    typedef typename traits_type::pointer         const_pointer;
    typedef typename traits_type::reference       reference;
    typedef typename traits_type::reference       const_reference;
    typedef typename traits_type::difference_type difference_type;
    typedef size_t                                size_type;
    typedef tab_ref_index_ref_iterator_<ITER, ITER_INDEX>   iterator;
    typedef iterator                              const_iterator;
    typedef std::reverse_iterator<const_iterator>     const_reverse_iterator;
    typedef std::reverse_iterator<iterator>           reverse_iterator;
    
  protected :

    ITER begin_;
    ITER_INDEX index_begin_, index_end_;

  public :
    
    bool empty(void) const { return index_begin_ == index_end_; }
    size_type size(void) const { return index_end_ - index_begin_; }
    
    iterator begin(void) { return iterator(begin_, index_begin_); }
    const_iterator begin(void) const
    { return iterator(begin_, index_begin_); }
    iterator end(void) { return iterator(begin_, index_end_); }
    const_iterator end(void) const { return iterator(begin_, index_end_); }
    reverse_iterator rbegin(void) { return reverse_iterator(end()); }
    const_reverse_iterator rbegin(void) const
    { return const_reverse_iterator(end()); }
    reverse_iterator rend(void) { return reverse_iterator(begin()); }
    const_reverse_iterator rend(void) const
    { return const_reverse_iterator(begin()); }
    
    reference front(void) { return *(begin_ + *index_begin_); }
    const_reference front(void) const { return *(begin_ + *index_begin_); }
    reference back(void) { return *(--(end())); }
    const_reference back(void) const { return *(--(end())); }
    void pop_front(void) { ++index_begin_; }
    
    tab_ref_index_ref(void) {}
    tab_ref_index_ref(const ITER &b, const ITER_INDEX &bi,
		      const ITER_INDEX &ei)
      : begin_(b), index_begin_(bi), index_end_(ei) {}
    
    // to be changed in a const_reference ?
    const_reference operator [](size_type ii) const
    { return *(begin_ + index_begin_[ii]);}
    reference operator [](size_type ii)
    { return *(begin_ + index_begin_[ii]); }

  };


  /* ********************************************************************* */
  /* Reference on regularly spaced elements.                               */
  /* ********************************************************************* */

  template<typename ITER> struct tab_ref_reg_spaced_iterator_ {
    
    typedef typename std::iterator_traits<ITER>::value_type value_type;
    typedef typename std::iterator_traits<ITER>::pointer    pointer;
    typedef typename std::iterator_traits<ITER>::reference  reference;
    typedef typename std::iterator_traits<ITER>::difference_type
                                                            difference_type;
    typedef typename std::iterator_traits<ITER>::iterator_category
                                                            iterator_category;
    typedef size_t size_type;
    typedef tab_ref_reg_spaced_iterator_<ITER> iterator;
    
    ITER it;
    size_type N;
    
    iterator operator ++(int) { iterator tmp = *this; it += N; return tmp; }
    iterator operator --(int) { iterator tmp = *this; it -= N; return tmp; }
    iterator &operator ++()   { it += N; return *this; }
    iterator &operator --()   { it -= N; return *this; }
    iterator &operator +=(difference_type i) { it += i * N; return *this; }
    iterator &operator -=(difference_type i) { it -= i * N; return *this; }
    iterator operator +(difference_type i) const 
    { iterator itt = *this; return (itt += i); }
    iterator operator -(difference_type i) const
    { iterator itt = *this; return (itt -= i); }
    difference_type operator -(const iterator &i) const
    { return (it - i.it) / N; }

    reference operator *() const { return *it; }
    reference operator [](int ii) const { return *(it + ii * N); }

    bool operator ==(const iterator &i) const { return (it == i.it); }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const { return (it < i.it); }

    tab_ref_reg_spaced_iterator_(void) {}
    tab_ref_reg_spaced_iterator_(const ITER &iter, size_type n)
      : it(iter), N(n) { }
    
  };

  /** 
      convenience template function for quick obtention of a strided iterator
      without having to specify its (long) typename
  */
  template<typename ITER> tab_ref_reg_spaced_iterator_<ITER> 
  reg_spaced_iterator(ITER it, size_t stride) {
    return tab_ref_reg_spaced_iterator_<ITER>(it, stride);
  }

  template<typename ITER> class tab_ref_reg_spaced
  {
    public :

      typedef typename std::iterator_traits<ITER>::value_type value_type;
      typedef typename std::iterator_traits<ITER>::pointer    pointer;
      typedef typename std::iterator_traits<ITER>::pointer    const_pointer;
      typedef typename std::iterator_traits<ITER>::reference  reference;
      typedef typename std::iterator_traits<ITER>::reference  const_reference;
      typedef typename std::iterator_traits<ITER>::difference_type
	                                                       difference_type;
      typedef size_t size_type;
      typedef tab_ref_reg_spaced_iterator_<ITER> iterator;
      typedef iterator                          const_iterator;
      typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
      typedef std::reverse_iterator<iterator> reverse_iterator;
    
    protected :

      ITER begin_, end_;
      size_type N;

    public :

      bool empty(void) const { return begin_ == end_; }
      size_type size(void) const { return (end_ - begin_) / N; }

      iterator begin(void) { return iterator(begin_, N); }
      const_iterator begin(void) const { return iterator(begin_, N); }
      iterator end(void) { return iterator(end_, N); }
      const_iterator end(void) const { return iterator(end_, N); }
      reverse_iterator rbegin(void) { return reverse_iterator(end()); }
      const_reverse_iterator rbegin(void) const
      { return const_reverse_iterator(end()); }
      reverse_iterator rend(void) { return reverse_iterator(begin()); }
      const_reverse_iterator rend(void) const
      { return const_reverse_iterator(begin()); }

      reference front(void) { return *begin_; }
      const_reference front(void) const { return *begin_; }
      reference back(void) { return *(--(end())); }
      const_reference back(void) const { return *(--(end())); }
      void pop_front(void) { begin_ += N; }

      tab_ref_reg_spaced(void) {}
      tab_ref_reg_spaced(const ITER &b, const ITER &e, size_type n)
	: begin_(b), end_(e), N(n) {}


      const_reference operator [](size_type ii) const
      { return *(begin_ + ii * N);}
      reference operator [](size_type ii) { return *(begin_ + ii * N); }

  };

  /* ********************************************************************* */
  /* Reference elements selected with a condition.                         */
  /* ********************************************************************* */

  template<typename ITER, typename COND> 
    struct tab_ref_with_selection_iterator_ : public ITER
  {
    typedef typename std::iterator_traits<ITER>::value_type value_type;
    typedef typename std::iterator_traits<ITER>::pointer    pointer;
    typedef typename std::iterator_traits<ITER>::reference  reference;
    typedef typename std::iterator_traits<ITER>::difference_type
                                                              difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef tab_ref_with_selection_iterator_<ITER, COND> iterator;
    const COND cond;
    
    void forward(void) { while (!(cond)(*this)) ITER::operator ++(); }
    iterator &operator ++()
    { ITER::operator ++(); forward(); return *this; }
    iterator operator ++(int)
    { iterator tmp = *this; ++(*this); return tmp; }
    
    tab_ref_with_selection_iterator_(void) {}
    tab_ref_with_selection_iterator_(const ITER &iter, const COND c)
      : ITER(iter), cond(c) {}
    
  };

  template<typename ITER, typename COND> class tab_ref_with_selection
  {
    public :

      typedef typename std::iterator_traits<ITER>::value_type value_type;
      typedef typename std::iterator_traits<ITER>::pointer    pointer;
      typedef typename std::iterator_traits<ITER>::pointer    const_pointer;
      typedef typename std::iterator_traits<ITER>::reference  reference;
      typedef typename std::iterator_traits<ITER>::reference  const_reference;
      typedef size_t  size_type;
      typedef tab_ref_with_selection_iterator_<ITER, COND> iterator;
      typedef iterator   const_iterator;
    
    protected :

      ITER begin_, end_;
      COND cond;

    public :

      iterator begin(void) const
      { iterator it(begin_, cond); it.forward(); return it; }
      iterator end(void) const { return iterator(end_, cond); }
      bool empty(void) const { return begin_ == end_; }

      value_type front(void) const { return *begin(); }
      void pop_front(void) { ++begin_; begin_ = begin(); }

      COND &condition(void) { return cond; }
      const COND &condition(void) const { return cond; }
   
      tab_ref_with_selection(void) {}
      tab_ref_with_selection(const ITER &b, const ITER &e, const COND &c)
	: begin_(b), end_(e), cond(c) { begin_ = begin(); }

  };

}

#endif /* DAL_REF_H__  */
