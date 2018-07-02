/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 Date: 2018-07-02

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
 \file levelsettoolmanager.hpp
 \author Thibaut Metivet <metivet@math.unistra.fr>
 \date 2018-07-02
 */
#ifndef _LEVELSETTOOLMANAGER_HPP
#define _LEVELSETTOOLMANAGER_HPP 1

namespace Feel {
namespace FeelModels {

template<
    typename ConvexType, typename BasisType, typename PeriodicityType = NoPeriodicity, 
    typename BasisPnType = BasisType
    >
class LevelSetToolManager
{
    typedef LevelSetToolManager<ConvexType, BasisType, PeriodicityType, BasisPnType> self_type;

public:
    typedef typename BasisType::value_type value_type;
    //--------------------------------------------------------------------//
    // Function space manager
    typedef LevelSetSpaceManager<ConvexType, BasisType, PeriodicityType, BasisPnType> levelset_space_manager_type;
    typedef boost::shared_ptr<levelset_space_manager_type> levelset_space_manager_ptrtype;
    // Default scalar and vectorial spaces
    typedef typename levelset_space_manager_type::space_scalar_type space_scalar_type;
    typedef typename levelset_space_manager_type::space_scalar_ptrtype space_scalar_ptrtype;
    typedef typename levelset_space_manager_type::space_vectorial_type space_vectorial_type;
    typedef typename levelset_space_manager_type::space_vectorial_ptrtype space_vectorial_ptrtype;
    // Tensor2 symmetric function space
    typedef typename levelset_space_manager_type::space_tensor2symm_type space_tensor2symm_type;
    typedef typename levelset_space_manager_type::space_tensor2symm_ptrtype space_tensor2symm_ptrtype;
    //--------------------------------------------------------------------//
    // Backend
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    //--------------------------------------------------------------------//
    // Projectors: scalar
    typedef Projector<space_scalar_type, space_scalar_type> projector_scalar_type;
    typedef boost::shared_ptr<projector_scalar_type> projector_scalar_ptrtype;
    // Projectors: vectorial
    typedef Projector<space_vectorial_type, space_vectorial_type> projector_vectorial_type;
    typedef boost::shared_ptr<projector_vectorial_type> projector_vectorial_ptrtype;
    // Projectors: tensor2symm
    typedef Projector<space_tensor2symm_type, space_tensor2symm_type> projector_tensor2symm_type;
    typedef boost::shared_ptr<projector_tensor2symm_type> projector_tensor2symm_ptrtype;

public:
    LevelSetToolManager( 
            levelset_space_manager_ptrtype const& spaceManager,
            std::string const& prefix
            );

    void createProjectorL2Default();
    void createProjectorSMDefault( double smoothCoeff, bool rebuild = false );
    void createProjectorL2Tensor2Symm();

    levelset_space_manager_ptrtype const& functionSpaceManager() const { return M_spaceManager; }

    projector_scalar_ptrtype const& projectorL2Scalar() const { return M_projectorL2Scalar; }
    projector_vectorial_ptrtype const& projectorL2Vectorial() const { return M_projectorL2Vectorial; }
    projector_tensor2symm_ptrtype const& projectorL2Tensor2Symm() const { return M_projectorL2Tensor2Symm; }

    projector_scalar_ptrtype const& projectorSMScalar() const { return M_projectorSMScalar; }
    projector_vectorial_ptrtype const& projectorSMVectorial() const { return M_projectorSMVectorial; }

private:
    std::string M_prefix;
    //--------------------------------------------------------------------//
    // Function space manager
    levelset_space_manager_ptrtype M_spaceManager;
    //--------------------------------------------------------------------//
    // Backends
    backend_ptrtype M_backendProjectorL2Scalar;
    backend_ptrtype M_backendProjectorL2Vectorial;
    backend_ptrtype M_backendProjectorSMScalar;
    backend_ptrtype M_backendProjectorSMVectorial;
    backend_ptrtype M_backendProjectorL2Tensor2Symm;
    //--------------------------------------------------------------------//
    // Projectors L2
    projector_scalar_ptrtype M_projectorL2Scalar;
    projector_vectorial_ptrtype M_projectorL2Vectorial;
    projector_tensor2symm_ptrtype M_projectorL2Tensor2Symm;
    // Projectors SMOOTH
    projector_scalar_ptrtype M_projectorSMScalar;
    projector_vectorial_ptrtype M_projectorSMVectorial;

};

#define LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS \
    template< typename ConvexType, typename BasisType, typename PeriodicityType, typename BasisPnType > \
        /**/
#define LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE \
    LevelSetToolManager<ConvexType, BasisType, PeriodicityType, BasisPnType> \
        /**/

LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS
LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::LevelSetToolManager( 
        levelset_space_manager_ptrtype const& spaceManager
        std::string const& prefix
        ) :
    M_spaceManager( spaceManager ),
    M_prefix( prefix )
{
}

LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::createProjectorL2Default()
{
    if( !M_projectorL2Scalar )
    {
        auto backendName = prefixvm( this->M_prefix, "projector-l2-scalar" );
        M_backendProjectorL2Scalar = backend_type::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpaceManager()->functionSpaceScalar()->worldComm()
                );
        M_projectorL2Scalar = projector(
                this->functionSpaceManager()->functionSpaceScalar(),
                this->functionSpaceManager()->functionSpaceScalar(),
                this->M_backendProjectorL2Scalar
                );
    }
    if( !M_projectorL2Vectorial )
    {
        auto backendName = prefixvm( this->M_prefix, "projector-l2-vectorial" );
        M_backendProjectorL2Vectorial = backend_type::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpaceManager()->functionSpaceVectorial()->worldComm()
                );
        M_projectorL2Vectorial = projector(
                this->functionSpaceManager()->functionSpaceVectorial(),
                this->functionSpaceManager()->functionSpaceVectorial(),
                this->M_backendProjectorL2Vectorial
                );
    }
}

LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::createProjectorSMDefault( double smoothCoeff, bool rebuild )
{
    if( !M_projectorSMScalar || rebuild )
    {
        auto backendName = prefixvm( this->M_prefix, "smoother-scalar" );
        M_backendProjectorSMScalar = backend_type::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpaceManager()->functionSpaceScalar()->worldComm()
                );
        M_projectorSMScalar = projector(
                this->functionSpaceManager()->functionSpaceScalar(),
                this->functionSpaceManager()->functionSpaceScalar(),
                this->M_backendProjectorSMScalar,
                DIFF, smoothCoeff, 30
                );
    }
    if( !M_projectorSMVectorial || rebuild )
    {
        auto backendName = prefixvm( this->M_prefix, "smoother-vectorial" );
        M_backendProjectorSMVectorial = backend_type::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpaceManager()->functionSpaceVectorial()->worldComm()
                );
        M_projectorL2Vectorial = projector(
                this->functionSpaceManager()->functionSpaceVectorial(),
                this->functionSpaceManager()->functionSpaceVectorial(),
                this->M_backendProjectorL2Vectorial,
                DIFF, smoothCoeff, 30
                );
    }
}

LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::createProjectorL2Tensor2Symm()
{
    if( !M_projectorL2Tensor2Symm )
    {
        auto backendName = prefixvm( this->M_prefix, "projector-l2-tensor2symm" );
        M_backendProjectorL2Tensor2Symm = backend_type::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpaceManager()->functionSpaceTensor2Symm()->worldComm()
                );
        M_projectorL2Tensor2Symm = projector(
                this->functionSpaceManager()->functionSpaceTensor2Symm(),
                this->functionSpaceManager()->functionSpaceTensor2Symm(),
                this->M_backendProjectorL2Tensor2Symm
                );
    }
}


}
}

#endif

