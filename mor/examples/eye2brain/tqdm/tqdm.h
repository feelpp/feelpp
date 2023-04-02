#pragma once

// From https://github.com/tqdm/tqdm.cpp

/**
Customisable progressbar decorator for iterators.
Includes a default range iterator printing to stderr.

TODO:
* handle s(tep) with operator+/-(=)
* chrono and delay printing
* class status printer
* iterator-less: fake range iterable, update() increments value

Usage:
  # include "tqdm/tqdm.h"
  for(int i : tqdm::range(4))
  // same as:
  //   for (int i : tqdm::tqdm({0, 1, 2, 3}))
  // or:
  //   std::vector<int> v{0, 1, 2, 3};
  //   for (int &i : tqdm::tqdm(v.begin(), v.end())
    ...

@author Casper dC-L <github.com/casperdcl>
*/

#include <cassert>      // assert
#include <cinttypes>    // PRIu64
#include <cstddef>      // ptrdiff_t, size_t
#include <cstdint>      // int64_t
#include <cstdio>       // printf
#include <iterator>     // iterator
#include <limits>       // numeric_limits
#include <stdexcept>    // throw
#include <string>       // string
#include <type_traits>  // is_pointer, ...
#include <utility>      // swap
#include "utils.h"

#ifndef SIZE_T_MAX
constexpr size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();
#endif

namespace tqdm {

struct Params {
  std::string desc;
  size_t total = -1;
  bool leave = true;
  FILE *f = stderr;
  int ncols = -1;
  float mininterval = 0.1f, maxinterval = 10.0f;
  unsigned miniters = -1;
  std::string ascii = " 123456789#";
  bool disable = false;
  std::string unit = "it";
  bool unit_scale = false;
  bool dynamic_ncols = false;
  float smoothing = 0.3f;
  std::string bar_format;
  size_t initial = 0;
  int position = -1;
  bool gui = false;
};

template <typename _Iterator>
class Tqdm : public MyIteratorWrapper<_Iterator> {
private:
  using TQDM_IT = MyIteratorWrapper<_Iterator>;
  _Iterator e;  // end
  Params self;  // ha, ha

public:
  /**
   containter-like methods
   */
  // actually current value
  // virtual _Iterator begin() { return this->get(); }
  Tqdm &begin() { return *this; }
  const Tqdm &begin() const { return *this; }
  // virtual _Iterator end() { return e; }
  Tqdm end() const { return Tqdm(e, e); }

  explicit operator _Iterator() { return this->get(); }

  /** constructors
   */
  explicit Tqdm(_Iterator begin, _Iterator end)
      : TQDM_IT(begin), e(end), self() {
    self.total = size_t(end - begin);
  }

  explicit Tqdm(_Iterator begin, size_t total)
      : TQDM_IT(begin), e(begin + total), self() {
    self.total = total;
  }

  // Tqdm(const Tqdm& other)
  //     : TQDM_IT(other.get()),
  //       e(other.end().get()),
  //       self(other.self),
  // {
  //   // std::memcpy(this, &other, sizeof(Tqdm));
  // }

  template <typename _Container,
            typename = typename std::enable_if<
                !std::is_same<_Container, Tqdm>::value>::type>
  Tqdm(_Container &v) : TQDM_IT(std::begin(v)), e(std::end(v)), self() {
    self.total = e - this->get();
  }

  explicit operator bool() const { return this->get() != e; }

  /** TODO: magic methods */
  virtual void _incr() const override {
    if (this->get() == e)
      throw std::out_of_range(
          "exhausted");  // TODO: don't throw, just double total

    TQDM_IT::_incr();
    if (this->get() == e) {
      printf("\nfinished: %" PRIu64 "/%" PRIu64 "\n",
        static_cast<std::uint64_t>(self.total),
        static_cast<std::uint64_t>(self.total));
    } else
      printf("\r%" PRIi64 " left", (int64_t)(e - this->get()));
  }
  virtual void _incr() override { ((Tqdm const &)*this)._incr(); }
};

template <typename _Iterator, typename _Tqdm = Tqdm<_Iterator>>
_Tqdm tqdm(_Iterator begin, _Iterator end) {
  return _Tqdm(begin, end);
}

template <typename _Iterator, typename _Tqdm = Tqdm<_Iterator>>
_Tqdm tqdm(_Iterator begin, size_t total) {
  return _Tqdm(begin, total);
}

template <typename _Container,
          typename _Tqdm = Tqdm<typename _Container::iterator>>
_Tqdm tqdm(_Container &v) {
  return _Tqdm(v);
}

template <size_t N, typename T, typename _Tqdm = Tqdm<T *>>
_Tqdm tqdm(T (&tab)[N]) {
  return _Tqdm(tab, N);
}

template <typename SizeType = int>
using RangeTqdm = Tqdm<RangeIterator<SizeType>>;
template <typename SizeType> RangeTqdm<SizeType> range(SizeType n) {
  return RangeTqdm<SizeType>(RangeIterator<SizeType>(n),
                             RangeIterator<SizeType>(n));
}
template <typename SizeType>
RangeTqdm<SizeType> range(SizeType start, SizeType end) {
  return RangeTqdm<SizeType>(RangeIterator<SizeType>(start, end),
                             RangeIterator<SizeType>(start, end));
}
template <typename SizeType>
RangeTqdm<SizeType> range(SizeType start, SizeType end, SizeType step) {
  return RangeTqdm<SizeType>(RangeIterator<SizeType>(start, end, step),
                             RangeIterator<SizeType>(start, end, step));
}

}  // tqdm

