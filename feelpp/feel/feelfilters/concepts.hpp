/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.
*/
#ifndef FEELPP_FEELFILTERS_CONCEPTS_HPP
#define FEELPP_FEELFILTERS_CONCEPTS_HPP 1

#include <concepts>
#include <type_traits>

namespace Feel
{

template <typename T>
concept FiltersMeshConcept = requires {
    { T::nDim } -> std::convertible_to<int>;
    { T::nOrder } -> std::convertible_to<int>;
    typename T::element_type;
    typename T::face_type;
    typename T::point_type;
};

template <typename T>
concept FiltersExporterConcept = requires( T t, double time ) {
    t.save();
    t.step( time );
};

template <typename T>
concept FiltersImporterConcept = requires( T t, typename T::mesh_type* mesh ) {
    t.visit( mesh );
};

} // namespace Feel

#endif // FEELPP_FEELFILTERS_CONCEPTS_HPP
