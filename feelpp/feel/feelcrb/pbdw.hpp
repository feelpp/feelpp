/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Romain Hild <romain.hild@cemosis.fr>
       Date: 2020-09-30

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
#ifndef _FEELPP_GEIM_HPP
#define _FEELPP_GEIM_HPP 1

#include <feel/feelalg/backend.hpp>
#include <feel/feelcrb/crbdb.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcore/unwrapptr.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>

#include <limits>

namespace Feel
{

template<typename RBSpace>
class PBDW
{
public:
    using reducedspace_type = RBSpace;
    using reducedspace_ptrtype = std::shared_ptr<reducedspace_type>;
    using space_type = typename reducedspace_type::fespace_type;
    using element_type = typename space_type::element_type;
    using vectorN_type = Eigen::VectorXd;
    using matrixN_type = Eigen::MatrixXd;

    PBDW(std::string const& name,
         reducedspace_ptrtype const& XR,
         std::vector<vector_ptrtype> sigmas);
    int dimensionN() const { return M_XR->size(); }
    int dimensionM() const { return M_sigmas.size(); }
    int dimension() const { return this->dimensionN()+this->dimensionM(); }
    matrixN_type matrix() const { return M_matrix; }
    void offline();
    vectorN_type online(vectorN_type const& yobs) const;
    element_type solution(vectorN_type const& yobs) const;

private:
    std::string M_name;
    reducedspace_ptrtype M_XR;
    std::vector<vector_ptrtype> M_sigmas;
    matrixN_type M_matrix;
};

template<typename RBSpace>
PBDW<RBSpace>::PBDW(std::string const& name,
                    reducedspace_ptrtype const& XR,
                    std::vector<vector_ptrtype> sigmas):
    M_name(name),
    M_XR(XR),
    M_sigmas(sigmas)
{}

template<typename RBSpace>
void
PBDW<RBSpace>::offline()
{
    int MN = this->dimension();
    int M = this->dimensionM();
    int N = this->dimensionN();
    M_matrix = matrixN_type::Zero(MN, MN);
    auto Xh = M_XR->functionSpace();
    auto u = Xh->element();
    auto mf = form2(_test=Xh, _trial=Xh);
    mf = integrate(elements(M_XR->mesh()), inner(id(u),idt(u)));
    auto m = mf.matrixPtr();
    for(int i = 0; i < M; ++i )
    {
        for(int j = 0; j < M; ++j )
            M_matrix(i, j) = m->energy(M_sigmas[i],M_sigmas[j]);
        for(int j = 0; j < N; ++j )
            M_matrix(i, M+j) = m->energy(unwrap_ptr(M_sigmas[i]), M_XR->primalBasisElement(j));
    }
    M_matrix.bottomLeftCorner(N, M) = M_matrix.topRightCorner(M, N).transpose();
}

template<typename RBSpace>
typename PBDW<RBSpace>::vectorN_type
PBDW<RBSpace>::online(vectorN_type const& yobs) const
{
    vectorN_type yobs2 = vectorN_type::Zero(this->dimension());
    yobs2.head(this->dimensionM()) = yobs;
    vectorN_type vn = M_matrix.colPivHouseholderQr().solve(yobs2);
    return vn;
}

template<typename RBSpace>
typename PBDW<RBSpace>::element_type
PBDW<RBSpace>::solution(vectorN_type const& yobs) const
{
    int M = this->dimensionM();
    int N = this->dimensionN();
    auto vn = this->online(yobs);
    auto wn = M_XR->primalRB();
    auto Xh = M_XR->functionSpace();
    auto I = Xh->element();
    I.zero();
    for(int i = 0; i < M; ++i )
        I.add(vn(i), Xh->element(M_sigmas[i]));
    for(int i = 0; i < N; ++i )
        I.add(vn(M+i), unwrap_ptr(wn[i]));
    return I;
}

} // namespace Feel

#endif
