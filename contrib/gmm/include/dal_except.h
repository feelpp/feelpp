/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_except.h : Exceptions.                                   */
/*     									   */
/*                                                                         */
/* Date : September 01, 2002.                                              */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*          Julien Pommier, Julien.pommier@gmm.insa-tlse.fr                */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
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

#ifndef DAL_EXCEPT_H__
#define DAL_EXCEPT_H__

#include <dal_std.h>

namespace dal {

/* *********************************************************************** */
/*	Getfem++ generic errors.                     			   */
/* *********************************************************************** */

  /* errors definitions  */

  using std::invalid_argument;

  class dimension_error : public std::logic_error {
  public:
    dimension_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class file_not_found_error : public std::logic_error {
  public:
    file_not_found_error(const std::string& what_arg)
      : std::logic_error (what_arg) { }
  };

  class internal_error : public std::logic_error {
  public:
    internal_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class failure_error : public std::logic_error {
  public:
    failure_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class not_linear_error : public std::logic_error {
  public:
    not_linear_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class to_be_done_error : public std::logic_error {
  public:
    to_be_done_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  #define DAL_STANDARD_CATCH_ERROR   catch(std::logic_error e) \
  { \
    cerr << "============================================\n";\
    cerr << "|      An error has been detected !!!      |\n";\
    cerr << "============================================\n";\
    cerr << e.what() << endl << endl;\
    exit(1);\
  }\
  catch(std::runtime_error e)\
  {\
    cerr << "============================================\n";\
    cerr << "|      An error has been detected !!!      |\n";\
    cerr << "============================================\n";\
    cerr << e.what() << endl << endl;\
    exit(1);\
  }\
  catch(std::bad_alloc) {\
    cerr << "============================================\n";\
    cerr << "|  A bad allocation has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  }\
  catch(std::bad_typeid) { \
    cerr << "============================================\n";\
    cerr << "|  A bad typeid     has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  } \
  catch(std::bad_exception) { \
    cerr << "============================================\n";\
    cerr << "|  A bad exception  has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  } \
  catch(std::bad_cast) { \
    cerr << "============================================\n";\
    cerr << "|    A bad cast  has been detected !!!     |\n";\
    cerr << "============================================\n";\
    exit(1);\
  } \
  catch(...) {\
    cerr << "============================================\n";\
    cerr << "|  An unknown error has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  }
//   catch(ios_base::failure) {
//     cerr << "============================================\n";
//     cerr << "| A ios_base::failure has been detected !!!|\n";
//     cerr << "============================================\n";
//     exit(1);
//   }

  /* callback handler for gmm/getfem exceptions */
  struct exception_callback {
    virtual ~exception_callback() {}
    virtual void callback(const std::string&) = 0; //{};

    static exception_callback *which_except(exception_callback *p = 0) {
      static exception_callback *exc_cback = 0;
      if (p != 0) exc_cback = p;
      return exc_cback;
    }

    static void do_exception_callback(const std::string &msg)
      { if (which_except()) which_except()->callback(msg); }

    static void set_exception_callback(exception_callback *e)
      { which_except(e); }

  };

  /* crashing callback for debug mode */
  struct exception_callback_debug : public dal::exception_callback  {
    virtual void callback(const std::string& msg)
    { cerr << msg << endl; *(int *)(0) = 0; }
  };

  /** user function for changing the default exception callback */
  inline void set_exception_callback(exception_callback *e)
  { exception_callback::which_except(e); }

#ifdef GETFEM_HAVE_PRETTY_FUNCTION
#  define DAL_PRETTY_FUNCTION __PRETTY_FUNCTION__
#else
#  define DAL_PRETTY_FUNCTION ""
#endif

#define DAL_THROW(type, thestr) {                                    \
    std::stringstream msg;                                           \
    msg << "Error in "__FILE__ << ", line "                          \
        << __LINE__ << " " << DAL_PRETTY_FUNCTION << ": \n" << thestr << ends; \
    dal::exception_callback::do_exception_callback(msg.str());       \
    throw (type)(msg.str());                                         \
  }

#ifdef DEBUG_MODE
#  define DAL_INTERNAL_ERROR(thestr) { \
  cerr << "Internal error: " << DAL_PRETTY_FUNCTION << " " << thestr << endl; \
   ::abort(); \
   }
#else
#  define DAL_INTERNAL_ERROR(thestr) DAL_THROW(dal::internal_error, "Internal error: " << thestr)
#endif


  struct warning_level {
    static int level(int l = -2)
    { static int level_ = 3; return (l != -2) ? (level_ = l) : level_; }
  };

  inline void set_warning_level(int l) { warning_level::level(std::max(0,l)); }

#define DAL_WARNING(level_, thestr) {                                 \
    std::stringstream msg;                                            \
    msg << "Level " << level_ << " Warning in "__FILE__ << ", line "  \
        << __LINE__ << " " << DAL_PRETTY_FUNCTION << ": " << thestr << ends; \
    if ((level_) <= dal::warning_level::level())                      \
       std::cerr << msg.str() << std::endl;                           \
  }

  // Warning levels : 0 always printed
  //                  1 very important : specify a possible error in the code.
  //                  2 important : specify a default of optimization for inst.
  //                  3 remark
  //                  4 ignored by default.

} /* end of namespace dal.                                                */



#endif /* DAL_EXCEPT_H__ */
