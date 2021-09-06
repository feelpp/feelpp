/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 This file is part of the Feel library

 Copyright (C) 2010 University of Coimbra

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
 \file ale.hpp
 \author Goncalo Pena <gpena@mat.uc.pt>
 \date 2010-10-12
 */


#ifndef FEELPP_MODELS_ALE_IMPL_H
#define FEELPP_MODELS_ALE_IMPL_H 1

#define ALE_WITH_BOUNDARYELEMENT 0


#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
//#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/interpolate.hpp>

#include <feel/feelmodels/modelcore/feelmodelscoreconstconfig.hpp>
#include <feel/feelmodels/modelcore/modelbase.hpp>
#include <feel/feelmodels/modelmesh/ale.hpp>
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
#include <feel/feelmodels/modelmesh/harmonicextension.hpp>
#endif
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
#include <feel/feelmodels/modelmesh/winslow.hpp>
#endif

/*
 * Class to construct ALE maps
 */

namespace Feel
{
namespace FeelModels
{
namespace ALE_IMPL
{

template< class Convex, int Order = 1 >
class ALE : public Feel::FeelModels::ALE<Convex,Order>
{
public :
    typedef Feel::FeelModels::ALE<Convex,Order> super_type;

    typedef ALE< Convex,Order> self_type;
    /*
     * Reference mesh typedefs
     */
    typedef Convex convex_type;

    static const uint16_type Dim = convex_type::nDim;
    static const uint16_type Order_low = convex_type::nOrder;
    using size_type = typename super_type::size_type;
    typedef Mesh< convex_type > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    using range_elements_type = elements_reference_wrapper_t<mesh_type>;

    static const bool isEqualOrderAndOrderLow = boost::is_same<mpl::int_<Order>,mpl::int_<Order_low> >::type::value;

    typedef typename super_type::template MyReferenceFunctionSpace<Order_low>::type space_low_type;
    typedef typename super_type::template MyReferenceFunctionSpace<Order>::type space_high_type;
    typedef typename super_type::template MyReferenceFunctionSpace<Order_low>::ptrtype space_low_ptrtype;
    typedef typename super_type::template MyReferenceFunctionSpace<Order>::ptrtype space_high_ptrtype;

    typedef typename  super_type::template MyReferenceFunctionSpace<Order_low>::elt_type element_low_type;
    typedef typename  super_type::template MyReferenceFunctionSpace<Order>::elt_type element_high_type;
    typedef typename  super_type::template MyReferenceFunctionSpace<Order_low>::elt_ptrtype element_low_ptrtype;
    typedef typename  super_type::template MyReferenceFunctionSpace<Order>::elt_ptrtype element_high_ptrtype;

    typedef typename super_type::ale_map_basis_type ale_map_basis_type;
    typedef typename super_type::ale_map_functionspace_type ale_map_functionspace_type;
    typedef typename super_type::ale_map_functionspace_ptrtype ale_map_functionspace_ptrtype;
    typedef typename super_type::ale_map_element_type ale_map_element_type;

    // Backend typedefs
    typedef Backend<double> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef Preconditioner<double> preconditioner_type;
    typedef std::shared_ptr<preconditioner_type> preconditioner_ptrtype;


#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
    typedef HarmonicExtension<mesh_type,Order_low> harmonicextension_type;
    typedef std::shared_ptr<harmonicextension_type> harmonicextension_ptrtype;
#endif
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
    typedef Winslow<mesh_type,Order_low/*+1*/ > winslow_type;
    typedef std::shared_ptr<winslow_type> winslow_ptrtype;
#endif
    /**
     * constructor
     *
     */
    ALE( std::string const& prefix="", worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(),
         ModelBaseRepository const& modelRep = ModelBaseRepository() );
    ALE( mesh_ptrtype mesh, std::string const& prefix="", worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(),
         ModelBaseRepository const& modelRep = ModelBaseRepository() );
    ALE( mesh_ptrtype mesh, range_elements_type const& rangeElt, std::string const& prefix="", worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(),
         ModelBaseRepository const& modelRep = ModelBaseRepository() );

    /**
     * copy constructor
     */
    ALE( ALE const& ) = default;

    /**
     * desctructor
     */
    ~ALE() {}

    void init() override;

    std::shared_ptr<std::ostringstream> getInfo() const override;

    /**
     * verbose
     */
    bool verboseSolverTimer() const { return M_verboseSolverTimer; }
    bool verboseSolverTimerAllProc() const { return M_verboseSolverTimerAllProc; }

    /**
     * \return the reference mesh
     */
    mesh_ptrtype referenceMesh();

    /**
     * \calculates the high order ALE map given the boundary's displacement
     */
    void generateMap( ale_map_element_type const & dispOnBoundary,
                      ale_map_element_type const & oldDisp ) override;

    /**
     * \return the low order displacement
     */
    element_low_type const& displacementLowOrder() const { return *M_displacementLow; }

    /**
     * \return the high order ALE map
     */
    ale_map_element_type const& map() const override { return map( mpl::bool_< (Order>Order_low) >() ); }
    element_high_type const& map( mpl::bool_<true> ) const { return *M_aleHigh; }
    element_low_type const& map( mpl::bool_<false> ) const { return *M_aleLow; }

    /**
     * \return the high order displacement
     */
    ale_map_element_type const& displacement() const override { return displacement( mpl::bool_< (Order>Order_low) >() ); }
    element_high_type const& displacement( mpl::bool_<true> ) const { return *M_displacementHigh; }
    element_low_type const& displacement( mpl::bool_<false> ) const { return *M_displacementLow; }

    /**
     * \returns an element containing the position of the points in the reference mesh
     */
    ale_map_element_type const& identity() const { return displacement( mpl::bool_< (Order>Order_low) >() ); }
    element_high_type const& identity( mpl::bool_<true> ) const { return *M_identityHigh; }
    element_low_type const& identity( mpl::bool_<false> ) const { return *M_identityLow; }

    /**
     * \return the functionspace that contains the high order ALE map
     */
    ale_map_functionspace_ptrtype const& functionSpace() const override { return functionSpace( mpl::bool_< (Order>Order_low) >() ); }
    space_high_ptrtype const& functionSpace( mpl::bool_<true> ) const { return M_fspaceHigh; }
    space_low_ptrtype const& functionSpace( mpl::bool_<false> ) const { return M_fspaceLow; }

    /**
     * reset all the data from the class
     */
    //void restart( mesh_ptrtype mesh );

    void initMetricMeshAdaptation() override;
    void updateMetricMeshAdaptation( Expr<GinacExVF<2>> const& e ) override;
    void updateMetricMeshAdaptation( typename super_type::metricmeshadaptation_type::element_scalar_type const& u ) override;

private:

    void createALE( std::optional<range_elements_type> const& rangeElt = std::nullopt );

#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
    void createHarmonicExtension();
#endif
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
    void createWinslow();
#endif

    /**
     * Precomputes the matrix associated with the modified harmonic extension operator
     * This operator has a penalization term that allows to reduce the distortion of the elements
     * in the new mesh, with respect to the reference mesh
     */
    void preCompute();
    void preComputeHO( mpl::true_ );
    void preComputeHO( mpl::false_ );

    /**
     * Creates the low order ALE map, given the boundary's displacement
     */
    void generateLowOrderMap_HARMONIC( ale_map_element_type const & dispOnBoundary );
    void generateLowOrderMap_WINSLOW( ale_map_element_type const & dispOnBoundary,
                                      ale_map_element_type const & oldDisp );

    /**
     * Pass information from low order maps to high order ones
     */
    void interpolateLow2High( mpl::true_ );
    void interpolateLow2High( mpl::false_ );

    //-----------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------//
    // compute ho correction
    void updateBoundaryElements( ale_map_element_type const & dispOnBoundary,
                                 mpl::bool_<true> /**/ );
    //-----------------------------------------------------------------------------------------//
    // compute ho correction : method 0 (problem on whole mesh)
    void updateBoundaryElements( ale_map_element_type const & dispOnBoundary,
                                 mpl::bool_<true> /**/, mpl::int_<0> /**/ );
    // with degree of in volume(3d) or face(2d) : order > 2
    void updateBoundaryElements( ale_map_element_type const & dispOnBoundary,
                                 mpl::bool_<true> /**/, mpl::int_<0> /**/, mpl::bool_<true> /**/ );
    // with no degree of in volume(3d) or face(2d) : order == 2
    void updateBoundaryElements( ale_map_element_type const & dispOnBoundary,
                                 mpl::bool_<true> /**/, mpl::int_<0> /**/, mpl::bool_<false> /**/ );
    //-----------------------------------------------------------------------------------------//
    // compute ho correction : method 1 (problem only on boundary elements)
    void updateBoundaryElements( ale_map_element_type const & dispOnBoundary,
                                 mpl::bool_<true> /**/, mpl::int_<1> /**/ );
    //-----------------------------------------------------------------------------------------//
    // no ho correction
    void updateBoundaryElements( ale_map_element_type const & dispOnBoundary,
                                 mpl::bool_<false> /**/ );

    // metric mesh adaptation
    void updateMetricMeshAdaptationForUse();

private :
    bool M_verboseSolverTimer,M_verboseSolverTimerAllProc;


    mesh_ptrtype M_reference_mesh;

    space_low_ptrtype M_fspaceLow;
    space_high_ptrtype M_fspaceHigh;
    space_high_ptrtype M_fspaceHighLocal;

    element_low_ptrtype M_aleLow, M_displacementLow, M_identityLow;
    element_high_ptrtype M_aleHigh, M_displacementHigh, M_identityHigh;

    // matrix/vector for ho correction
    backend_ptrtype M_bHigh;
    sparse_matrix_ptrtype M_harmonicHigh;
    vector_ptrtype M_rhsHigh;
    preconditioner_ptrtype M_preconditionerHigh;

    std::string M_alemeshTypeName;
    bool M_doHoCorrection;

#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
    harmonicextension_ptrtype M_harmonicextensionFactory;
#endif
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
    winslow_ptrtype M_winslowFactory;
#endif
    bool M_isInitHarmonicExtension, M_isInitWinslow;

    std::map<std::string,std::set<size_type> > M_dofsHighOnBoundary;
};

} // namespace ALE_IMPL
} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_ALE_IMPL_H
