//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//! 
//! Copyright (C) 2017-present Feel++ Consortium
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <mpi4py/mpi4py.h>
#include <petsc4py/petsc4py.h>

#include "petsc_casters.hpp"

#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/binding.hpp>
//#include <boost/parameter/python.hpp>
#include <boost/mpl/vector.hpp>

#include<feel/feelalg/backend.hpp>
#include<feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/matrixpetsc.hpp>
#include <pybind11/stl_bind.h>


namespace py = pybind11;

using namespace Feel;

class PyVectorDouble : Vector<double,unsigned int>
{
    using super = Vector<double,unsigned int>;
    using super::super;
    using size_type = super::size_type;

    size_type size() const override { PYBIND11_OVERLOAD_PURE( size_type, super, size ); }
    size_type localSize() const override { PYBIND11_OVERLOAD_PURE( size_type, super, localSize ); }
    void close() override { PYBIND11_OVERLOAD_PURE(void, super, close ); }
    void zero() override { PYBIND11_OVERLOAD_PURE( void, super, zero ); }
    void zero(size_type start, size_type stop) override {
        PYBIND11_OVERLOAD_PURE( void, super, zero, start, stop );
    }
    void setConstant(value_type v) override { PYBIND11_OVERLOAD_PURE( void, super, setConstant, v );}
    clone_ptrtype clone() const override { PYBIND11_OVERLOAD_PURE( clone_ptrtype, super, clone ); }
    value_type sum() const override { PYBIND11_OVERLOAD_PURE( value_type, super, sum ); }
    real_type min() const override { PYBIND11_OVERLOAD_PURE( real_type, super, min ); }
    real_type max() const override { PYBIND11_OVERLOAD_PURE( real_type, super, max ); }
    real_type l1Norm() const override { PYBIND11_OVERLOAD_PURE( real_type, super, l1Norm ); }
    real_type l2Norm() const override { PYBIND11_OVERLOAD_PURE( real_type, super, l2Norm ); }
    real_type linftyNorm() const override { PYBIND11_OVERLOAD_PURE( real_type, super, linftyNorm );}
    double operator()(const size_type i) const override {
        PYBIND11_OVERLOAD_PURE( double, super, operator(), i ); }
    double& operator()(const size_type i) override {
        PYBIND11_OVERLOAD_PURE( double&, super, operator(), i ); }
    super& operator+=(const super& V) override {
        PYBIND11_OVERLOAD_PURE( super&, super, operator+=, V ); }
    super& operator-=(const super& V) override {
        PYBIND11_OVERLOAD_PURE( super&, super, operator-=, V ); }
    void set(const size_type i, const value_type& value) override {
        PYBIND11_OVERLOAD_PURE( void, super, set, i, value ); }
    void setVector(int* i, int n, value_type* v) override {
        PYBIND11_OVERLOAD_PURE( void, super, setVector, i, n, v ); }
    void add(const size_type i, const value_type& value) override {
        PYBIND11_OVERLOAD_PURE( void, super, add, i, value ); }
    void addVector(int* i, int n, value_type* v, size_type K, size_type K2) override {
        PYBIND11_OVERLOAD_PURE( void, super, addVector, i, n, v, K, K2 ); }
    void add(const value_type& s) override {
        PYBIND11_OVERLOAD_PURE( void, super, add, s ); }
    void add(const super& V) override {
        PYBIND11_OVERLOAD_PURE( void, super, add, V ); }
    void add(const value_type& a, const super& V) override {
        PYBIND11_OVERLOAD_PURE( void, super, add, a, V ); }
    void addVector(const std::vector<double>& v, const std::vector<size_type>& dof_indices) override {
        PYBIND11_OVERLOAD_PURE( void, super, addVector, v, dof_indices ); }
    void addVector(const super& v, const std::vector<size_type>& dof_indices) override {
        PYBIND11_OVERLOAD_PURE( void, super, addVector, v, dof_indices ); }
    void addVector(const super& V_in, const MatrixSparse<double>& A_in) override {
        PYBIND11_OVERLOAD_PURE( void, super, addVector, V_in, A_in ); }
    value_type dot( super const& v ) const override {
        PYBIND11_OVERLOAD_PURE( value_type, super, dot, v ); }
    void insert(const std::vector<double>& v, const std::vector<size_type>& dof_indices) override {
        PYBIND11_OVERLOAD_PURE( void, super, insert, v, dof_indices ); }
    void insert(const super& V, const std::vector<size_type>& dof_indices) override {
        PYBIND11_OVERLOAD_PURE( void, super, insert, V, dof_indices ); }
    void insert(const ublas::vector<double>& V, const std::vector<size_type>& dof_indices) override {
        PYBIND11_OVERLOAD_PURE( void, super, insert, V, dof_indices ); }
    void scale(const double factor) override {
        PYBIND11_OVERLOAD_PURE( void, super, scale, factor ); }
    void printMatlab(const std::string name, bool renumber) const override {
        PYBIND11_OVERLOAD_PURE( void, super, printMatlab, name, renumber ); }
};

class PyMatrixSparseDouble : MatrixSparse<double>
{
    using super = MatrixSparse<double>;
    using self_t = PyMatrixSparseDouble;
    using super::super;
    void init ( const size_type m,
                const size_type n,
                const size_type m_l,
                const size_type n_l,
                const size_type nnz=30,
                const size_type noz=10 ) override {
        PYBIND11_OVERLOAD_PURE(void, super, init, m, n, m_l, n_l, nnz, noz); }
    void init ( const size_type m,
                const size_type n,
                const size_type m_l,
                const size_type n_l,
                graph_ptrtype const& graph ) override 
    {
        PYBIND11_OVERLOAD_PURE(void, super, init, m, n, m_l, n_l); 
    }
    clone_ptrtype clone() const override { PYBIND11_OVERLOAD_PURE( clone_ptrtype, super, clone ); }

    size_type nnz() const override 
    {
        PYBIND11_OVERLOAD_PURE(size_type, super, nnz); 
    }
    void clear () override {
        PYBIND11_OVERLOAD_PURE(void, super, clear); }
    void zero () override {
        PYBIND11_OVERLOAD_PURE(void, super, zero); }
    void zero ( size_type start1, size_type size1,
                size_type start2, size_type size2 ) override {
        PYBIND11_OVERLOAD_PURE(void, super, zero, start1, size1, start2, size2); }
    void close () const override {
        PYBIND11_OVERLOAD_PURE(void, super, close); }
    size_type size1 () const override {
        PYBIND11_OVERLOAD_PURE(size_type, super, size1); }
    size_type size2 () const override {
        PYBIND11_OVERLOAD_PURE(size_type, super, size2); }
    size_type rowStart () const override {
        PYBIND11_OVERLOAD_PURE(size_type, super, rowStart); }
    size_type rowStop () const override {
        PYBIND11_OVERLOAD_PURE(size_type, super, rowStop); }
    void set ( const size_type i,
               const size_type j,
               const value_type& value ) override {
        PYBIND11_OVERLOAD_PURE(void, super, set, i, j, value); }
    void add ( const size_type i,
               const size_type j,
               const value_type& value ) override {
        PYBIND11_OVERLOAD_PURE(void, super, add, i, j, value); }
    void addMatrix ( const ublas::matrix<value_type> &dm,
                     const std::vector<size_type> &rows,
                     const std::vector<size_type> &cols ) override {
        PYBIND11_OVERLOAD_PURE(void, super, addMatrix, dm, rows, cols); }
    void addMatrix ( int* rows, int nrows,
                     int* cols, int ncols,
                     value_type* data, size_type K, size_type K2 ) override {
        PYBIND11_OVERLOAD_PURE(void, super, addMatrix, rows, nrows, cols, ncols, data, K, K2); }
    void addMatrix ( const ublas::matrix<value_type> &dm,
                     const std::vector<size_type> &dof_indices ) override {
        PYBIND11_OVERLOAD_PURE(void, super, addMatrix, dm, dof_indices); }
    void addMatrix ( const double f, super const& m, Feel::MatrixStructure matStruc = Feel::SAME_NONZERO_PATTERN ) override {
        PYBIND11_OVERLOAD_PURE(void, super, addMatrix, f, m, matStruc); }
    void scale ( const double f ) override {
        PYBIND11_OVERLOAD_PURE(void, super, scale, f); }
    double operator () ( const size_type i,
                         const size_type j ) const override {
        PYBIND11_OVERLOAD_PURE(double, super, operator (), i, j); }
    super& operator = ( super const& M ) override {
        PYBIND11_OVERLOAD_PURE(super&, super,  operator=); }
    void diagonal ( Vector<double>& dest ) const override {
        PYBIND11_OVERLOAD_PURE(void, super, diagonal, dest); }
    void transpose( MatrixSparse<value_type>& Mt, size_type options = MATRIX_TRANSPOSE_ASSEMBLED ) const override {
        PYBIND11_OVERLOAD_PURE(void, super, transpose, Mt, options); }
    real_type energy ( vector_type const& v,
                       vector_type const& u,
                       bool transpose = false ) const override {
        PYBIND11_OVERLOAD_PURE(real_type, super, energy, u, v, transpose); }
    real_type l1Norm () const override {
        PYBIND11_OVERLOAD_PURE(real_type, super, l1Norm); }
    real_type linftyNorm () const override {
        PYBIND11_OVERLOAD_PURE(real_type, super, linftyNorm); }
    void zeroRows( std::vector<int> const& rows,
                   Vector<value_type> const& values,
                   Vector<value_type>& rhs,
                   Context const& on_context,
                   value_type value_on_diagonal ) override {
        PYBIND11_OVERLOAD_PURE(void, super, zeroRows, rows, values, rhs, on_context, value_on_diagonal); }
    void updateBlockMat( std::shared_ptr<super > const& m, std::vector<size_type> const& start_i, std::vector<size_type> const& start_j ) override {
        PYBIND11_OVERLOAD_PURE(void, super, updateBlockMat, m, start_i, start_j); }
};

PYBIND11_MODULE(_alg, m )
{
    using namespace Feel;

    if (import_mpi4py()<0 && import_petsc4py()<0) return ;

    py::class_<datamap_t<uint32_type>, datamap_ptr_t<uint32_type>>( m, "DataMap" )
        .def( py::init<std::shared_ptr<WorldComm>>() )
        .def( py::init<uint32_type, uint32_type, std::shared_ptr<WorldComm>>() )
//        .def( py::init<uint32_type, std::vector<int>, std::vector<int>>() )
//        .def( py::init<std::vector<std::shared_ptr<DataMap<uint32_type>>>, std::shared_ptr<WorldComm>>() )
        .def( "nDof", &DataMap<uint32_type>::nDof, "the total number of degrees of freedom in the problem" )
        .def( "nLocalDof", &DataMap<uint32_type>::nLocalDof, "the total number of degrees of freedom in the problem" )
        .def( "nLocalDofWithoutGhost", static_cast<uint32_type (DataMap<uint32_type>::*)() const>(&DataMap<uint32_type>::nLocalDofWithoutGhost), "the total number of degrees of freedom in the problem" )

        ;
    py::class_<Vector<double, uint32_type>, PyVectorDouble, std::shared_ptr<Vector<double, uint32_type>> >(m, "VectorDouble")
        .def(py::init<>())
        ;
    py::class_<VectorPetsc<double>, Vector<double, uint32_type>, std::shared_ptr<VectorPetsc<double>>>( m, "VectorPetscDouble" )
        .def( py::init<>() )
        .def( py::init<int, std::shared_ptr<WorldComm>>() )
        .def( py::init<int, int, std::shared_ptr<WorldComm>>() )
        .def( py::init<datamap_ptr_t<uint32_type>, bool>() )
        .def( "size", &VectorPetsc<double>::size, "return  PETSc Vector size" )
        .def( "clear", &VectorPetsc<double>::clear, "clear PETSc vector" )
        .def( "zero", static_cast<void ( VectorPetsc<double>::* )()>(&VectorPetsc<double>::zero), "zero PETSc vector" )

        .def( "vec", static_cast<Vec ( VectorPetsc<double>::* )() const>( &VectorPetsc<double>::vec ), "return a PETSc Vector" )
        ;
    py::class_<MatrixSparse<double>, PyMatrixSparseDouble, std::shared_ptr<MatrixSparse<double>> >(m, "MatrixSparseDouble")
        .def(py::init<>())
        ;
    py::class_<MatrixPetsc<double>, MatrixSparse<double>, std::shared_ptr<MatrixPetsc<double>>>( m, "MatrixPetscDouble" )
        .def( py::init<worldcomm_ptr_t>() )
        .def( py::init<datamap_ptr_t<uint32_type>,datamap_ptr_t<uint32_type>>() )
        .def( py::init<datamap_ptr_t<uint32_type>,datamap_ptr_t<uint32_type>,worldcomm_ptr_t>() )
        .def( "size1", &MatrixPetsc<double>::size1, "return  PETSc Matrix row size" )
        .def( "size2", &MatrixPetsc<double>::size2, "return  PETSc Matrix column size" )
        .def( "rowStart", &MatrixPetsc<double>::rowStart, "return  PETSc Matrix row start" )
        .def( "rowStop", &MatrixPetsc<double>::rowStop, "return  PETSc Matrix row stop " )
        .def( "clear", &MatrixPetsc<double>::clear, "clear PETSc matrix" )
        .def( "zero", static_cast<void ( MatrixPetsc<double>::* )()>(&MatrixPetsc<double>::zero), "zero PETSc matrix" )
        .def( "mat", static_cast<Mat ( MatrixPetsc<double>::* )() const>( &MatrixPetsc<double>::mat ), "return a PETSc sparse matrix" )
        ;

    py::class_<VectorUblas<double>, Vector<double,uint32_type>, std::shared_ptr<VectorUblas<double>>> vublas(m,"VectorUBlas<double,ublas::vector<double>>");
    vublas.def(py::init<>())
        .def(py::self + py::self )
        .def(py::self - py::self)
        .def(double() + py::self)
        .def(double() - py::self)
        .def(py::self + double())
        .def(py::self - double())
        .def(py::self * double())
        .def(double() * py::self)
        //.def(double() / py::self)
        .def(py::self += py::self)
        .def(py::self -= py::self)
        .def(py::self *= double())
        .def(-py::self)
//        .def(py::self /= double())
//        .def(py::self *= py::self)
//        .def(py::self /= py::self)
        ;

   m.def( "backend_options", &Feel::backend_options, py::arg("prefix"), "create a backend options descriptions with prefix" );
}


#if 0
PYBIND11_MODULE(_alg, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;

    m.def( "backend_options", &Feel::backend_options, py::arg("prefix"), "create a backend options descriptions with prefix" );
    py::class_<VectorUblas<double>,std::shared_ptr<VectorUblas<double>>> vublas(m,"VectorUBlas<double,ublas::vector<double>>");
    vublas.def(py::init<>())
        .def(py::self + py::self )
        .def(py::self - py::self)
        .def(double() + py::self)
        .def(double() - py::self)
        .def(py::self + double())
        .def(py::self - double())
        .def(py::self * double())
        .def(double() * py::self)
        //.def(double() / py::self)
        .def(py::self += py::self)
        .def(py::self -= py::self)
        .def(py::self *= double())
        .def(-py::self)
//        .def(py::self /= double())
//        .def(py::self *= py::self)
//        .def(py::self /= py::self)
        ;
}
#endif