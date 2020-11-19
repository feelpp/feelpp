/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-07-15

  Copyright (C) 2014-2016 Feel++ Consortium

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
#ifndef FEELPP_GLPK_HPP
#define FEELPP_GLPK_HPP 1

#if defined(FEELPP_HAS_GLPK_H)
#include <glpk.h>
#include <iostream>
#include <vector>
#include <exception>

namespace Feel
{
namespace opt
{

class OptimizationLinearProgramming
{
public:
    OptimizationLinearProgramming(int direction = GLP_MIN, std::string name = "");
    ~OptimizationLinearProgramming();
    int solve();
    void addRow(std::string name, int type, double lb, double ub);
    void addColumn(std::string name, double coef, int type, double lb, double ub);
    void setMatrix(std::vector<std::vector<double> > matrix);
    double getObjectiveValue();
    double getColumnPrimalValue(int i); // index start at 1 !
    static std::string const& check(int e);
    glp_prob* problem() const { return M_pb; }
    glp_smcp* params() const { return M_params; }

private:
    glp_prob* M_pb;
    glp_smcp* M_params;
    int M_scaling;

}; // class OptimizationLinearProgramming

struct FeelGlpkException : std::exception
{
    explicit FeelGlpkException(int e);
    const char* what() const noexcept;

    int M_e;
};

} // namespace opt
} // namespace Feel

#endif /* FEELPP_HAS_GLPK_H */
#endif
