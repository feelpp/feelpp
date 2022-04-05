//! -*- mode: c++; coding: utf-9; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
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
//! @file levelsetredistanciation_hj.hpp
//! @author Thibaut Metivet <thibaut.metivet@gmail.com>
//! @date 06 Sep 2019
//! @copyright 2019 Feel++ Consortium
//!
#ifndef _LEVELSET_REDISTANCIATION_HJ_HPP
#define _LEVELSET_REDISTANCIATION_HJ_HPP 1

#include <feel/feelmodels/levelset/levelsetredistanciation.hpp>
#include <feel/feeldiscr/projector.hpp>
//#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>
#include <feel/feelmodels/levelset/levelsetdeltaexpr.hpp>

namespace Feel
{

template<typename FunctionSpaceType>
class LevelSetRedistanciationHJ
: public LevelSetRedistanciation<FunctionSpaceType>
{
public:
    //--------------------------------------------------------------------//
    // Typedefs
    typedef LevelSetRedistanciation<FunctionSpaceType> super_type;
    typedef LevelSetRedistanciationHJ<FunctionSpaceType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef FunctionSpaceType functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename functionspace_type::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;

    typedef typename functionspace_type::periodicity_0_type periodicity_type;
    static const bool is_periodic = functionspace_type::is_periodic;

    typedef FunctionSpace<mesh_type, bases<Lagrange<1,Vectorial,Continuous>>> functionspace_P1v_type;
    typedef std::shared_ptr<functionspace_P1v_type> functionspace_P1v_ptrtype;

    typedef Projector<functionspace_P1v_type, functionspace_P1v_type> projectorL2_vectorial_type;
    typedef std::shared_ptr<projectorL2_vectorial_type> projectorL2_vectorial_ptrtype;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Scalar,Discontinuous> > > functionspace_P0_type;
    typedef std::shared_ptr<functionspace_P0_type> functionspace_P0_ptrtype;

    //--------------------------------------------------------------------//
    // Hamilton-Jacobi advection
#if 0
    template<typename SpaceType>
    class AdvectionHJ
        : public Feel::FeelModels::AdvDiffReac<SpaceType/*,TODO*/>
        , public std::enable_shared_from_this< AdvectionHJ<SpaceType> >
    {
    public:
        typedef Feel::FeelModels::AdvDiffReac<SpaceType> super_type;

        typedef AdvectionHJ<SpaceType> self_type;
        typedef std::shared_ptr<self_type> self_ptrtype;

        typedef SpaceType functionspace_type;
        typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;

        // Constructor
        AdvectionHJ(
                std::string const& prefix = ""
                );
        
        static self_ptrtype New(
                std::string const& prefix = ""
                );

        // Initialization
        void init( functionspace_ptrtype const& space, bool buildModelAlgebraicFactory = true );

        // BC and source term assembly
        void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const;
        void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;
        bool hasSourceTerm() const { return false; }
    };

    //typedef AdvectionHJ<FunctionSpaceType> advectionhj_type;
    typedef Feel::FeelModels::AdvDiffReac<FunctionSpaceType> advectionhj_type;
    typedef std::shared_ptr<advectionhj_type> advectionhj_ptrtype;
#endif

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Constructor
    LevelSetRedistanciationHJ( 
            functionspace_ptrtype const& space,
            std::string const& prefix = "" );
    //--------------------------------------------------------------------//
    // Parameters
    void loadParametersFromOptionsVm();

    double tolerance() const { return M_tolerance; }
    void setTolerance( double tol ) { M_tolerance = tol; }
    
    double timeStep() const { return 1./*M_advectionHJ->timeStep()*/; }
    void setTimeStep( double dt ) { /*M_advectionHJ->setTimeStep( dt );*/ }

    int maxIterations() const { return M_maxIterations; }
    void setMaxIterations( int max ) { M_maxIterations = max; }

    double thicknessHeaviside() const { return M_thicknessHeaviside; }
    void setThicknessHeaviside( double eps ) { M_thicknessHeaviside = eps; }

    functionspace_P0_ptrtype functionSpaceP0() const { return M_functionSpaceP0; }
    //--------------------------------------------------------------------//
    // Run redistanciation
    element_type run( element_type const& phi ) const;

private:
    //advectionhj_ptrtype M_advectionHJ;

    double M_tolerance;
    double M_timeStep;
    int M_maxIterations;
    double M_thicknessHeaviside;

    mutable int M_nGlobalIter;

    functionspace_P1v_ptrtype M_functionSpaceP1Vec;
    projectorL2_vectorial_ptrtype M_projectorL2Vec;

    bool M_useVolumeConstraint;
    functionspace_P0_ptrtype M_functionSpaceP0;
};


} // namespace Feel

#endif
