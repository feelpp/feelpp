/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
  This file is part of the Feel library

  Copyright (C) 2023 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#pragma once
#include <string>
#include <utility>

namespace Feel
{

/**
 * @brief Check if any material in the model has the given property.
 *
 * @param model The model object
 * @param propName The property name to check for
 * @return true if at least one material has the property, false otherwise
 */
template<typename ModelType>
bool hasAnyMaterialWithProperty(ModelType const& model, std::string const& propName)
{
    for (auto const& [physicName, physicData] : model.physicsFromCurrentType())
    {
        for (std::string const& matName : model.materialsProperties()->physicToMaterials(physicName))
        {
            if (model.materialsProperties()->hasProperty(matName, propName))
                return true;
        }
    }
    return false;
}

template<typename ModelType>
bool hasMaterialWithProperty(ModelType const& model, std::string const& matName, std::string const& propName)
{
    return model.materialsProperties()->hasProperty(matName, propName);
}

/**
 * @brief Apply a function to each material in the model.
 *
 * This function iterates over all materials defined in the model and applies
 * the provided function to each material along with its associated mesh element range.
 *
 * @param model The model object
 * @param f The function to apply, which takes a material name and a range of mesh elements
 */
template<typename ModelType, typename Func>
void forEachMaterial(const ModelType& model, Func&& f)
{
    for ( auto const& [physicName, physicData] : model.physicsFromCurrentType() )
    {
        for ( std::string const& matName : model.materialsProperties()->physicToMaterials(physicName) )
        {
            auto const& range = model.materialsProperties()->rangeMeshElementsByMaterial(model.mesh(), matName);
            f(matName, range);
        }
    }
}

/**
 * Apply a function to each material in the model, but only if a predicate is satisfied.
 * 
 * This function iterates over all materials defined in the model and applies
 * the provided function to each material whose name satisfies the given predicate.
 * @param model The model object
 * @param pred The predicate function that takes a material name and returns true if the material should be processed
 * @param f The function to apply, which takes a material name and a range of mesh elements
 */
template<typename ModelType, typename Predicate, typename Func>
void forEachMaterialIf(const ModelType& model, Predicate&& pred, Func&& f)
{
    forEachMaterial(model,
        [&](std::string const& matName, auto const& range)
        {
            if (pred(matName))
                f(matName, range);
        });
}

/**
 * Apply a function to all materials that have a specific property.
 * This function iterates over all materials defined in the model and applies
 * the provided function to each material that has the specified property.
 * @param model The model object
 * @param propName The name of the property to check for
 * @param f The function to apply, which takes a material name, a range of mesh elements,
 *          and the material property associated with that name.
 */
template<typename ModelType, typename Func>
void forEachMaterialWithProperty(const ModelType& model, std::string const& propName, Func&& f)
{
    forEachMaterial(model,
        [&](std::string const& matName, auto const& range)
        {
            if (model.materialsProperties()->hasProperty(matName, propName))
            {
                auto const& prop = model.materialsProperties()->materialProperty(matName, propName);
                f(matName, range, prop);
            }
        });
}

/**
 * Apply a function to materials that have a property, with rank-aware expr<M,N>.
 * This function iterates over all materials defined in the model and applies
 * the provided function to each material that has the specified property,
 * using an expression of rank M x N.
 * @param model The model object
 * @param propName The name of the property to check for
 * @param symbolsExpr The symbols expression to use in the coefficient expression
 * @param f The function to apply, which takes a material name, a range of mesh elements,
 *          and the coefficient expression associated with that property.
 *          The coefficient expression is constructed using expr<M,N>(prop.expr(), symbolsExpr).
 *          The function should accept the coefficient expression as a parameter.
 *          The coefficient expression is of type expr<M,N>.
 *          The function should accept the coefficient expression as a parameter.
 * @tparam M The rank of the coefficient expression (default is 1)
 * @tparam N The rank of the coefficient expression (default is 1)
 * @tparam ModelType The type of the model object
 * @tparam Func The type of the function to apply
 * @tparam SymbolsExpr The type of the symbols expression
 */
template<int M = 1, int N = 1, typename ModelType, typename Func, typename SymbolsExpr>
void forEachMaterialWithCoefficientExpr(const ModelType& model,
                                        std::string const& propName,
                                        SymbolsExpr const& symbolsExpr,
                                        Func&& f)
{
    forEachMaterialWithProperty(model, propName,
        [&](std::string const& matName, auto const& range, auto const& prop)
        {
            auto coeffExpr = expr(prop.template expr<M, N>(), symbolsExpr);
            f(matName, range, coeffExpr);
        });
}

/**
 * @brief Apply a function to each material with a single coefficient property.
 *
 * This function iterates over all materials defined in the model and applies
 * the provided function to each material that has the specified property.
 * The function receives the material name, the range of mesh elements, and the coefficient expression.
 *
 * @param propName The name of the property to check for
 * @param symbolsExpr The symbols expression to use in the coefficient expression
 * @param functor The function to apply, which takes a material name, a range of mesh elements,
 *                and the coefficient expression associated with that property.
 */
template<typename ModelType, typename SymbolsExpr, typename Functor>
void forEachMaterialWith2Coefficients(ModelType const& model,
                                       std::string const& propName1,
                                       std::string const& propName2,
                                       SymbolsExpr const& symbolsExpr,
                                       Functor&& functor)
{
    forEachMaterial(
        [&](std::string const& matName, auto const& range)
        {
            constexpr int nDim = ModelType::nDim;
            auto props = model.materialsProperties();
            if (props->hasProperty(matName, propName1) && props->hasProperty(matName, propName2))
            {
                auto prop1 = props->materialProperty(matName, propName1);
                auto prop2 = props->materialProperty(matName, propName2);

                auto expr1 = expr(prop1.template expr<nDim,1>(), symbolsExpr); // convection alpha
                auto expr2 = expr(prop2.template expr<nDim,nDim>(), symbolsExpr); // diffusion c

                functor(matName, range, expr1, expr2);
            }
        });
}
} // namespace Feel
