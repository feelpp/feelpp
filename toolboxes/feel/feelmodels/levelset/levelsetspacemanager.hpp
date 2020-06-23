/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 Date: 2018-05-18

 Copyright (C) 2018 Université de Strasbourg

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
/**
 \file levelsetspacemanager.hpp
 \author Thibaut Metivet <metivet@math.unistra.fr>
 \date 2018-05-18
 */
#ifndef _LEVELSETSPACEMANAGER_HPP
#define _LEVELSETSPACEMANAGER_HPP 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelmodels/modelcore/utils.hpp>

namespace Feel {
namespace FeelModels {

template<
    typename ConvexType, typename BasisType, typename PeriodicityType = NoPeriodicity, 
    typename BasisPnType = BasisType
    >
class LevelSetSpaceManager
{
    typedef LevelSetSpaceManager<ConvexType, BasisType, PeriodicityType, BasisPnType> self_type;

public:
    static const uint16_type Order = BasisType::nOrder;

    //--------------------------------------------------------------------//
    // Mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef elements_reference_wrapper_t<mesh_type> range_elements_type;

    //--------------------------------------------------------------------//
    // Periodicity
    typedef PeriodicityType periodicity_type;

    //--------------------------------------------------------------------//
    // Worldcomm
    using worldcomm_type = WorldComm;
    using worldcomm_ptrtype = std::shared_ptr<WorldComm>;
    using worldscomm_ptrtype = std::vector<std::shared_ptr<WorldComm>>;

    //--------------------------------------------------------------------//
    // Default scalar and vectorial spaces
    typedef BasisType basis_scalar_type;
    //typedef FunctionSpace< mesh_type, bases<basis_scalar_type>, Periodicity<periodicity_type> > space_scalar_type;
    //typedef std::shared_ptr<space_scalar_type> space_scalar_ptrtype;
    //typedef typename space_scalar_type::element_type element_scalar_type;
    //typedef std::shared_ptr<element_scalar_type> element_scalar_ptrtype;

    typedef typename detail::ChangeBasisPolySet<Vectorial, basis_scalar_type>::type basis_vectorial_type;
    typedef FunctionSpace<mesh_type, bases<basis_vectorial_type>/*, Periodicity<periodicity_type>*/ > space_vectorial_type;
    typedef std::shared_ptr<space_vectorial_type> space_vectorial_ptrtype;
    typedef typename space_vectorial_type::element_type element_vectorial_type;
    typedef std::shared_ptr<element_vectorial_type> element_vectorial_ptrtype;

    typedef typename space_vectorial_type::component_functionspace_type space_scalar_type;
    typedef std::shared_ptr<space_scalar_type> space_scalar_ptrtype;
    typedef typename space_scalar_type::element_type element_scalar_type;
    typedef std::shared_ptr<element_scalar_type> element_scalar_ptrtype;

    typedef typename space_scalar_type::value_type value_type;

    //--------------------------------------------------------------------//
    // PN (iso) spaces
    typedef BasisPnType basis_scalar_PN_type;
    typedef FunctionSpace< mesh_type, bases<basis_scalar_PN_type>, Periodicity<periodicity_type> > space_scalar_PN_type;
    typedef std::shared_ptr<space_scalar_PN_type> space_scalar_PN_ptrtype;
    typedef typename space_scalar_PN_type::element_type element_scalar_PN_type;
    typedef std::shared_ptr<element_scalar_PN_type> element_scalar_PN_ptrtype;

    typedef typename detail::ChangeBasisPolySet<Vectorial, basis_scalar_PN_type>::type basis_vectorial_PN_type;
    typedef FunctionSpace<mesh_type, bases<basis_vectorial_PN_type>, Periodicity<periodicity_type> > space_vectorial_PN_type;
    typedef std::shared_ptr<space_vectorial_PN_type> space_vectorial_PN_ptrtype;
    typedef typename space_vectorial_PN_type::element_type element_vectorial_PN_type;
    typedef std::shared_ptr<element_vectorial_PN_type> element_vectorial_PN_ptrtype;

    //--------------------------------------------------------------------//
    // High-order visu spaces
    typedef Mesh< Simplex< nDim, 1, nDim > > mesh_hovisu_type;
    typedef std::shared_ptr<mesh_hovisu_type> mesh_hovisu_ptrtype;
    typedef Lagrange<1, Vectorial, Continuous, PointSetFekete> basis_vectorial_hovisu_type;
    typedef FunctionSpace<mesh_hovisu_type, bases<basis_vectorial_hovisu_type> > space_vectorial_hovisu_type;
    typedef std::shared_ptr<space_vectorial_hovisu_type> space_vectorial_hovisu_ptrtype;
    typedef typename space_vectorial_hovisu_type::element_type element_vectorial_hovisu_type;
    typedef std::shared_ptr<element_vectorial_hovisu_type> element_vectorial_hovisu_ptrtype;

    //typedef Lagrange<1, Scalar, Continuous, PointSetFekete> basis_scalar_hovisu_type;
    //typedef FunctionSpace< mesh_type, bases<basis_scalar_hovisu_type> > space_scalar_hovisu_type;
    typedef typename space_vectorial_hovisu_type::component_functionspace_type space_scalar_hovisu_type;
    typedef std::shared_ptr<space_scalar_hovisu_type> space_scalar_hovisu_ptrtype;
    typedef typename space_scalar_hovisu_type::element_type element_scalar_hovisu_type;
    typedef std::shared_ptr<element_scalar_hovisu_type> element_scalar_hovisu_ptrtype;


    //--------------------------------------------------------------------//
    // Lagrange P1isoPn operators
    typedef OperatorLagrangeP1<space_scalar_PN_type> op_lagrangeP1_type;
    typedef std::shared_ptr<op_lagrangeP1_type> op_lagrangeP1_ptrtype;

    //--------------------------------------------------------------------//
    // Lagrange P1hovisu operators
    typedef OperatorLagrangeP1<space_scalar_type> op_lagrangeP1_hovisu_type;
    typedef std::shared_ptr<op_lagrangeP1_hovisu_type> op_lagrangeP1_hovisu_ptrtype;

    //--------------------------------------------------------------------//
    // isoPN interpolation operators
    typedef OperatorInterpolation<
        space_scalar_type, // from space
        space_scalar_PN_type // to space
        > op_interpolation_scalar_to_PN_type;
    typedef std::shared_ptr<op_interpolation_scalar_to_PN_type> op_interpolation_scalar_to_PN_ptrtype;
    typedef OperatorInterpolation<
        space_scalar_PN_type, // from space
        space_scalar_type // to space
        > op_interpolation_scalar_from_PN_type;
    typedef std::shared_ptr<op_interpolation_scalar_from_PN_type> op_interpolation_scalar_from_PN_ptrtype;

    typedef OperatorInterpolation<
        space_vectorial_type, // from space
        space_vectorial_PN_type // to space
        > op_interpolation_vectorial_to_PN_type;
    typedef std::shared_ptr<op_interpolation_vectorial_to_PN_type> op_interpolation_vectorial_to_PN_ptrtype;
    typedef OperatorInterpolation<
        space_vectorial_PN_type, // from space
        space_vectorial_type // to space
        > op_interpolation_vectorial_from_PN_type;
    typedef std::shared_ptr<op_interpolation_vectorial_from_PN_type> op_interpolation_vectorial_from_PN_ptrtype;

    //--------------------------------------------------------------------//
    // Hovisu interpolation operators
    typedef OperatorInterpolation<
        space_scalar_type, // from space
        space_scalar_hovisu_type // to space
        > op_interpolation_scalar_to_hovisu_type;
    typedef std::shared_ptr<op_interpolation_scalar_to_hovisu_type> op_interpolation_scalar_to_hovisu_ptrtype;

    typedef OperatorInterpolation<
        space_vectorial_type, // from space
        space_vectorial_hovisu_type // to space
        > op_interpolation_vectorial_to_hovisu_type;
    typedef std::shared_ptr<op_interpolation_vectorial_to_hovisu_type> op_interpolation_vectorial_to_hovisu_ptrtype;

    //--------------------------------------------------------------------//
    // Space markers P0
    typedef Lagrange<0, Scalar, Discontinuous> basis_markers_type;
    typedef FunctionSpace<mesh_type, bases<basis_markers_type>, value_type, Periodicity<NoPeriodicity> > space_markers_type;
    typedef std::shared_ptr<space_markers_type> space_markers_ptrtype;
    typedef typename space_markers_type::element_type element_markers_type;
    typedef std::shared_ptr<element_markers_type> element_markers_ptrtype;

    //--------------------------------------------------------------------//
    // Tensor2 symmetric function space
    //typedef Lagrange<Order, Tensor2Symm> basis_tensor2symm_type;
    typedef typename detail::ChangeBasisPolySet<Tensor2Symm, basis_scalar_type>::type basis_tensor2symm_type;
    typedef FunctionSpace<mesh_type, bases<basis_tensor2symm_type>, Periodicity<periodicity_type> > space_tensor2symm_type;
    typedef std::shared_ptr<space_tensor2symm_type> space_tensor2symm_ptrtype;
    typedef typename space_tensor2symm_type::element_type element_tensor2symm_type;
    typedef std::shared_ptr<element_tensor2symm_type> element_tensor2symm_ptrtype;

public:
    LevelSetSpaceManager( 
            mesh_ptrtype const& mesh,
            std::string const& prefix = "",
            std::string const& rootRepository = soption( _name="exporter.directory" )
            );

    // Extended doftable
    bool buildExtendedDofTable() const { return M_buildExtendedDofTable; }
    void setBuildExtendedDofTable( bool buildExtDT ) { M_buildExtendedDofTable = buildExtDT; }
    // Periodicity
    void setPeriodicity( periodicity_type const& p );
    periodicity_type const& periodicity() const { return M_periodicity; }

    // Default minimal function spaces
    void createFunctionSpaceDefault();
    // Iso PN minimal function spaces
    void createFunctionSpaceIsoPN();

    // Hovisu function spaces
    void createFunctionSpaceHovisu();
    
    // Tensor2Symm function space
    void createFunctionSpaceTensor2Symm();

    mesh_ptrtype const& mesh() const { return M_mesh; }
    mesh_ptrtype const& meshIsoPN() const { return M_meshIsoPN; }
    mesh_hovisu_ptrtype const& meshHovisu() const { return M_meshHovisu; }

    range_elements_type const& rangeMeshElements() { return M_rangeMeshElements; }
    range_elements_type const& rangeMeshPNElements() { return M_rangeMeshElements; }
    range_elements_type const& rangeMeshIsoPNElements() { return M_rangeMeshIsoPNElements; }
    range_elements_type const& rangeMeshHovisu() { return M_rangeMeshHovisuElements; }

    space_scalar_ptrtype const& functionSpaceScalar() const { return M_spaceScalar; }
    space_vectorial_ptrtype const& functionSpaceVectorial() const { return M_spaceVectorial; }
    space_markers_ptrtype const& functionSpaceMarkers() const { return M_spaceMarkers; }
    space_scalar_PN_ptrtype const& functionSpaceScalarPN() const { return M_spaceScalarPN; }
    space_vectorial_PN_ptrtype const& functionSpaceVectorialPN() const { return M_spaceVectorialPN; }
    space_scalar_ptrtype const& functionSpaceScalarIsoPN() const { return M_spaceScalarIsoPN; }
    space_vectorial_ptrtype const& functionSpaceVectorialIsoPN() const { return M_spaceVectorialIsoPN; }
    space_scalar_hovisu_ptrtype const& functionSpaceScalarHovisu() const { return M_spaceScalarHovisu; }
    space_vectorial_hovisu_ptrtype const& functionSpaceVectorialHovisu() const { return M_spaceVectorialHovisu; }
    space_markers_ptrtype const& functionSpaceMarkersIsoPN() const { return M_spaceMarkersIsoPN; }
    space_tensor2symm_ptrtype const& functionSpaceTensor2Symm() const { return M_spaceTensor2Symm; }

    op_interpolation_scalar_to_PN_ptrtype const& opInterpolationScalarToPN() const { return M_opInterpolationScalarToPN; }
    op_interpolation_scalar_from_PN_ptrtype const& opInterpolationScalarFromPN() const { return M_opInterpolationScalarFromPN; }
    op_interpolation_vectorial_to_PN_ptrtype const& opInterpolationVectorialToPN() const { return M_opInterpolationVectorialToPN; }
    op_interpolation_vectorial_from_PN_ptrtype const& opInterpolationVectorialFromPN() const { return M_opInterpolationVectorialFromPN; }
    op_interpolation_scalar_to_hovisu_ptrtype const& opInterpolationScalarToHovisu() const { return M_opInterpolationScalarToHovisu; }
    op_interpolation_vectorial_to_hovisu_ptrtype const& opInterpolationVectorialToHovisu() const { return M_opInterpolationVectorialToHovisu; }

private:
    worldscomm_ptrtype const& worldsComm() const { return M_worldsComm; }

private:
    std::string M_prefix;
    std::string M_rootRepository;
    //--------------------------------------------------------------------//
    // Meshes
    mesh_ptrtype M_mesh;
    mesh_ptrtype M_meshIsoPN;
    mesh_hovisu_ptrtype M_meshHovisu;

    // Ranges
    range_elements_type M_rangeMeshElements;
    range_elements_type M_rangeMeshIsoPNElements;
    range_elements_type M_rangeMeshHovisuElements;
    //--------------------------------------------------------------------//
    // WorldsComm
    worldscomm_ptrtype M_worldsComm;
    // Extended doftable
    bool M_buildExtendedDofTable;
    // Periodicity
    periodicity_type M_periodicity;
    bool M_functionSpaceCreated;
    //--------------------------------------------------------------------//
    // Default function spaces
    space_scalar_ptrtype M_spaceScalar;
    space_vectorial_ptrtype M_spaceVectorial;
    space_markers_ptrtype M_spaceMarkers;
    // PN function spaces
    space_scalar_PN_ptrtype M_spaceScalarPN;
    space_vectorial_PN_ptrtype M_spaceVectorialPN;
    // Hovisu function spaces
    space_scalar_hovisu_ptrtype M_spaceScalarHovisu;
    space_vectorial_hovisu_ptrtype M_spaceVectorialHovisu;
    // IsoPN function spaces
    space_scalar_ptrtype M_spaceScalarIsoPN;
    space_vectorial_ptrtype M_spaceVectorialIsoPN;
    space_markers_ptrtype M_spaceMarkersIsoPN;
    // Tensor2Symm function space
    space_tensor2symm_ptrtype M_spaceTensor2Symm;
    //--------------------------------------------------------------------//
    // Lagrange P1 operators
    op_lagrangeP1_ptrtype M_opLagrangeP1isoPN;
    op_lagrangeP1_hovisu_ptrtype M_opLagrangeP1Hovisu;
    //--------------------------------------------------------------------//
    // Interpolation operators
    op_interpolation_scalar_to_PN_ptrtype M_opInterpolationScalarToPN;
    op_interpolation_scalar_from_PN_ptrtype M_opInterpolationScalarFromPN;
    op_interpolation_vectorial_to_PN_ptrtype M_opInterpolationVectorialToPN;
    op_interpolation_vectorial_from_PN_ptrtype M_opInterpolationVectorialFromPN;
    op_interpolation_scalar_to_hovisu_ptrtype M_opInterpolationScalarToHovisu;
    op_interpolation_vectorial_to_hovisu_ptrtype M_opInterpolationVectorialToHovisu;
};

} // namespace FeelModels
} // namespace Feel

#endif // _LEVELSETSPACEMANAGER_HPP

