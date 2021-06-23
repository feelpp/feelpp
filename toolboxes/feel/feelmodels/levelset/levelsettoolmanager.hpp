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

#include <feel/feelmodels/levelset/levelsetcurvaturediffusion.hpp>

namespace Feel {
namespace FeelModels {

//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
// Derivation methods
enum class LevelSetDerivationMethod { 
    NODAL_PROJECTION, L2_PROJECTION, SMOOTH_PROJECTION/*, PN_NODAL_PROJECTION*/
};
typedef boost::bimap<std::string, LevelSetDerivationMethod> levelsetderivationmethod_maptype;
static const levelsetderivationmethod_maptype LevelSetDerivationMethodMap = boost::assign::list_of< levelsetderivationmethod_maptype::relation >
    ( "nodal-projection", LevelSetDerivationMethod::NODAL_PROJECTION )
    ( "l2-projection", LevelSetDerivationMethod::L2_PROJECTION )
    ( "smooth-projection", LevelSetDerivationMethod::SMOOTH_PROJECTION )
    //( "pn-nodal-projection", LevelSetDerivationMethod::PN_NODAL_PROJECTION )
;

enum class LevelSetCurvatureMethod { 
    NODAL_PROJECTION, L2_PROJECTION, SMOOTH_PROJECTION, /*PN_NODAL_PROJECTION,*/
    DIFFUSION_ORDER1, DIFFUSION_ORDER2
};
typedef boost::bimap<std::string, LevelSetCurvatureMethod> levelsetcurvaturemethod_maptype;
static const levelsetcurvaturemethod_maptype LevelSetCurvatureMethodMap = boost::assign::list_of< levelsetcurvaturemethod_maptype::relation >
    ( "nodal-projection", LevelSetCurvatureMethod::NODAL_PROJECTION )
    ( "l2-projection", LevelSetCurvatureMethod::L2_PROJECTION )
    ( "smooth-projection", LevelSetCurvatureMethod::SMOOTH_PROJECTION )
    //( "pn-nodal-projection", LevelSetCurvatureMethod::PN_NODAL_PROJECTION )
    ( "diffusion-order1", LevelSetCurvatureMethod::DIFFUSION_ORDER1 )
    ( "diffusion-order2", LevelSetCurvatureMethod::DIFFUSION_ORDER2 )
;

//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
// LevelSetToolManager class
template<
    typename ConvexType, typename BasisType, typename PeriodicityType = NoPeriodicity, 
    typename BasisPnType = BasisType
    >
class LevelSetToolManager
{
    typedef LevelSetToolManager<ConvexType, BasisType, PeriodicityType, BasisPnType> self_type;

public:
    //--------------------------------------------------------------------//
    // Function space manager
    typedef LevelSetSpaceManager<ConvexType, BasisType, PeriodicityType, BasisPnType> levelset_space_manager_type;
    typedef std::shared_ptr<levelset_space_manager_type> levelset_space_manager_ptrtype;

    typedef typename levelset_space_manager_type::value_type value_type;
    // Default scalar and vectorial spaces and elements
    typedef typename levelset_space_manager_type::space_scalar_type space_scalar_type;
    typedef typename levelset_space_manager_type::space_scalar_ptrtype space_scalar_ptrtype;
    typedef typename space_scalar_type::element_type element_scalar_type;
    typedef typename space_scalar_type::element_ptrtype element_scalar_ptrtype;
    typedef typename levelset_space_manager_type::space_vectorial_type space_vectorial_type;
    typedef typename levelset_space_manager_type::space_vectorial_ptrtype space_vectorial_ptrtype;
    typedef typename space_vectorial_type::element_type element_vectorial_type;
    typedef typename space_vectorial_type::element_ptrtype element_vectorial_ptrtype;
    // Tensor2 symmetric function space and element
    typedef typename levelset_space_manager_type::space_tensor2symm_type space_tensor2symm_type;
    typedef typename levelset_space_manager_type::space_tensor2symm_ptrtype space_tensor2symm_ptrtype;
    typedef typename space_tensor2symm_type::element_type element_tensor2symm_type;
    typedef typename space_tensor2symm_type::element_ptrtype element_tensor2symm_ptrtype;
    //--------------------------------------------------------------------//
    // Backend
    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;
    //--------------------------------------------------------------------//
    // Projectors: scalar
    typedef Projector<space_scalar_type, space_scalar_type> projector_scalar_type;
    typedef std::shared_ptr<projector_scalar_type> projector_scalar_ptrtype;
    // Projectors: vectorial
    typedef Projector<space_vectorial_type, space_vectorial_type> projector_vectorial_type;
    typedef std::shared_ptr<projector_vectorial_type> projector_vectorial_ptrtype;
    // Projectors: tensor2symm
    typedef Projector<space_tensor2symm_type, space_tensor2symm_type> projector_tensor2symm_type;
    typedef std::shared_ptr<projector_tensor2symm_type> projector_tensor2symm_ptrtype;
    //--------------------------------------------------------------------//
    // Curvature diffusion method
    typedef LevelSetCurvatureDiffusion<space_scalar_type> levelset_curvaturediffusion_type;
    typedef std::shared_ptr<levelset_curvaturediffusion_type> levelset_curvaturediffusion_ptrtype;

public:
    LevelSetToolManager( 
            levelset_space_manager_ptrtype const& spaceManager,
            std::string const& prefix
            );

    void createProjectorL2Default();
    void createProjectorSMDefault();
    void createProjectorL2Tensor2Symm();
    void createProjectorL2IsoPN();
    void createProjectorSMIsoPN();

    void createCurvatureDiffusion();

    levelset_space_manager_ptrtype const& functionSpaceManager() const { return M_spaceManager; }

    projector_scalar_ptrtype const& projectorL2Scalar() const { return M_projectorL2Scalar; }
    projector_vectorial_ptrtype const& projectorL2Vectorial() const { return M_projectorL2Vectorial; }
    projector_tensor2symm_ptrtype const& projectorL2Tensor2Symm() const { return M_projectorL2Tensor2Symm; }

    projector_scalar_ptrtype const& projectorSMScalar() const { return M_projectorSMScalar; }
    projector_vectorial_ptrtype const& projectorSMVectorial() const { return M_projectorSMVectorial; }

    projector_scalar_ptrtype const& projectorL2ScalarIsoPN() const { return M_projectorL2ScalarIsoPN; }
    projector_vectorial_ptrtype const& projectorL2VectorialIsoPN() const { return M_projectorL2VectorialIsoPN; }
    projector_scalar_ptrtype const& projectorSMScalarIsoPN() const { return M_projectorSMScalarIsoPN; }
    projector_vectorial_ptrtype const& projectorSMVectorialIsoPN() const { return M_projectorSMVectorialIsoPN; }

    element_vectorial_type grad( element_scalar_type const& phi, LevelSetDerivationMethod method ) const;
    element_vectorial_type grad( element_scalar_ptrtype const& phi, LevelSetDerivationMethod method ) const { return this->grad( *phi, method ); }
    element_scalar_type modGrad( element_scalar_type const& phi, LevelSetDerivationMethod method ) const;
    element_scalar_type modGrad( element_scalar_ptrtype const& phi, LevelSetDerivationMethod method ) const { return this->modGrad( *phi, method ); }

    levelset_curvaturediffusion_ptrtype const& curvatureDiffusion() const { return M_curvatureDiffusion; }

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

    backend_ptrtype M_backendProjectorL2ScalarIsoPN;
    backend_ptrtype M_backendProjectorL2VectorialIsoPN;
    backend_ptrtype M_backendProjectorSMScalarIsoPN;
    backend_ptrtype M_backendProjectorSMVectorialIsoPN;
    //--------------------------------------------------------------------//
    // Projectors L2
    projector_scalar_ptrtype M_projectorL2Scalar;
    projector_vectorial_ptrtype M_projectorL2Vectorial;
    projector_tensor2symm_ptrtype M_projectorL2Tensor2Symm;

    projector_scalar_ptrtype M_projectorL2ScalarIsoPN;
    projector_vectorial_ptrtype M_projectorL2VectorialIsoPN;
    projector_tensor2symm_ptrtype M_projectorL2Tensor2SymmIsoPN;
    // Projectors SMOOTH
    double M_projectorSMScalarCoeff;
    projector_scalar_ptrtype M_projectorSMScalar;
    double M_projectorSMVectorialCoeff;
    projector_vectorial_ptrtype M_projectorSMVectorial;

    double M_projectorSMScalarIsoPNCoeff;
    projector_scalar_ptrtype M_projectorSMScalarIsoPN;
    double M_projectorSMVectorialIsoPNCoeff;
    projector_vectorial_ptrtype M_projectorSMVectorialIsoPN;
    //--------------------------------------------------------------------//
    // Curvature diffusion method
    levelset_curvaturediffusion_ptrtype M_curvatureDiffusion;

};

#define LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS \
    template< typename ConvexType, typename BasisType, typename PeriodicityType, typename BasisPnType > \
        /**/
#define LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE \
    LevelSetToolManager<ConvexType, BasisType, PeriodicityType, BasisPnType> \
        /**/

LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS
LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::LevelSetToolManager( 
        levelset_space_manager_ptrtype const& spaceManager,
        std::string const& prefix
        ) :
    M_prefix( prefix ),
    M_spaceManager( spaceManager )
{
    double h = M_spaceManager->mesh()->hAverage();
    double O = BasisType::nOrder;
    double On = BasisPnType::nOrder;
    M_projectorSMScalarCoeff = h / O * doption( _name="smooth-coeff", _prefix=prefixvm(M_prefix,"projector-sm-scalar") );
    M_projectorSMVectorialCoeff = h / O * doption( _name="smooth-coeff", _prefix=prefixvm(M_prefix,"projector-sm-vectorial") );
    M_projectorSMScalarIsoPNCoeff = h / On * doption( _name="smooth-coeff", _prefix=prefixvm(M_prefix,"projector-sm-scalar-isopn") );
    M_projectorSMVectorialIsoPNCoeff = h / On * doption( _name="smooth-coeff", _prefix=prefixvm(M_prefix,"projector-sm-vectorial-isopn") );
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
                this->functionSpaceManager()->functionSpaceScalar()->worldCommPtr()
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
                this->functionSpaceManager()->functionSpaceVectorial()->worldCommPtr()
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
LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::createProjectorSMDefault()
{
    if( !M_projectorSMScalar )
    {
        auto backendName = prefixvm( this->M_prefix, "projector-sm-scalar" );
        M_backendProjectorSMScalar = backend_type::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpaceManager()->functionSpaceScalar()->worldCommPtr()
                );
        M_projectorSMScalar = projector(
                this->functionSpaceManager()->functionSpaceScalar(),
                this->functionSpaceManager()->functionSpaceScalar(),
                this->M_backendProjectorSMScalar,
                DIFF, M_projectorSMScalarCoeff, 30
                );
    }
    if( !M_projectorSMVectorial )
    {
        auto backendName = prefixvm( this->M_prefix, "projector-sm-vectorial" );
        M_backendProjectorSMVectorial = backend_type::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpaceManager()->functionSpaceVectorial()->worldCommPtr()
                );
        M_projectorSMVectorial = projector(
                this->functionSpaceManager()->functionSpaceVectorial(),
                this->functionSpaceManager()->functionSpaceVectorial(),
                this->M_backendProjectorSMVectorial,
                DIFF, M_projectorSMVectorialCoeff, 30
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
                this->functionSpaceManager()->functionSpaceTensor2Symm()->worldCommPtr()
                );
        M_projectorL2Tensor2Symm = projector(
                this->functionSpaceManager()->functionSpaceTensor2Symm(),
                this->functionSpaceManager()->functionSpaceTensor2Symm(),
                this->M_backendProjectorL2Tensor2Symm
                );
    }
}

LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::createProjectorL2IsoPN()
{
    if( !M_projectorL2ScalarIsoPN )
    {
        auto backendName = prefixvm( this->M_prefix, "projector-l2-scalar-isopn" );
        M_backendProjectorL2ScalarIsoPN = backend_type::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpaceManager()->functionSpaceScalarIsoPN()->worldCommPtr()
                );
        M_projectorL2ScalarIsoPN = projector(
                this->functionSpaceManager()->functionSpaceScalarIsoPN(),
                this->functionSpaceManager()->functionSpaceScalarIsoPN(),
                this->M_backendProjectorL2ScalarIsoPN
                );
    }
    if( !M_projectorL2VectorialIsoPN )
    {
        auto backendName = prefixvm( this->M_prefix, "projector-l2-vectorial-isopn" );
        M_backendProjectorL2VectorialIsoPN = backend_type::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpaceManager()->functionSpaceVectorialIsoPN()->worldCommPtr()
                );
        M_projectorL2VectorialIsoPN = projector(
                this->functionSpaceManager()->functionSpaceVectorialIsoPN(),
                this->functionSpaceManager()->functionSpaceVectorialIsoPN(),
                this->M_backendProjectorL2VectorialIsoPN
                );
    }
}

LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::createProjectorSMIsoPN()
{
    if( !M_projectorSMScalarIsoPN )
    {
        auto backendName = prefixvm( this->M_prefix, "projector-sm-scalar-isopn" );
        M_backendProjectorSMScalarIsoPN = backend_type::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpaceManager()->functionSpaceScalarIsoPN()->worldCommPtr()
                );
        M_projectorSMScalarIsoPN = projector(
                this->functionSpaceManager()->functionSpaceScalarIsoPN(),
                this->functionSpaceManager()->functionSpaceScalarIsoPN(),
                this->M_backendProjectorSMScalarIsoPN,
                DIFF, M_projectorSMScalarIsoPNCoeff, 30
                );
    }
    if( !M_projectorSMVectorialIsoPN )
    {
        auto backendName = prefixvm( this->M_prefix, "projector-sm-vectorial-isopn" );
        M_backendProjectorSMVectorialIsoPN = backend_type::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpaceManager()->functionSpaceVectorialIsoPN()->worldCommPtr()
                );
        M_projectorSMVectorialIsoPN = projector(
                this->functionSpaceManager()->functionSpaceVectorialIsoPN(),
                this->functionSpaceManager()->functionSpaceVectorialIsoPN(),
                this->M_backendProjectorSMVectorialIsoPN,
                DIFF, M_projectorSMVectorialIsoPNCoeff, 30
                );
    }
}

LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::createCurvatureDiffusion()
{
    if( !M_curvatureDiffusion )
    {
        M_curvatureDiffusion = std::make_shared<levelset_curvaturediffusion_type>( 
                this->functionSpaceManager()->functionSpaceScalar(), 
                prefixvm( this->M_prefix, "curvature-diffusion" ) 
                );
    }
}

LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::element_vectorial_type
LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::grad( element_scalar_type const& phi, LevelSetDerivationMethod method ) const
{
    switch( method )
    {
        case LevelSetDerivationMethod::NODAL_PROJECTION:
            return vf::project( 
                    _space=this->functionSpaceManager()->functionSpaceVectorial(),
                    _range=this->functionSpaceManager()->rangeMeshElements(),
                    _expr=trans(gradv(phi))
                    );
        case LevelSetDerivationMethod::L2_PROJECTION:
            //return this->projectorL2Vectorial()->project( _expr=trans(gradv(phi)) );
            return this->projectorL2Vectorial()->derivate( idv(phi) );
        case LevelSetDerivationMethod::SMOOTH_PROJECTION:
            return this->projectorSMVectorial()->project( trans(gradv(phi)) );
        //case LevelSetDerivationMethod::PN_NODAL_PROJECTION:
            //CHECK( M_useSpaceIsoPN ) << "use-space-iso-pn must be enabled to use PN_NODAL_PROJECTION \n";
            //auto phiPN = this->functionSpaceManager()->opInterpolationScalarToPN()->operator()( phi );
            //auto gradPhiPN = vf::project(
                    //_space=this->functionSpaceManager()->functionSpaceVectorialPN(),
                    //_range=this->functionSpaceManager()->rangeMeshPNElements(),
                    //_expr=trans(gradv(phiPN))
                    //);
            //return this->functionSpaceManager()->opInterpolationVectorialFromPN()->operator()( gradPhiPN );
    }
}

LEVELSETTOOLMANAGER_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::element_scalar_type
LEVELSETTOOLMANAGER_CLASS_TEMPLATE_TYPE::modGrad( element_scalar_type const& phi, LevelSetDerivationMethod method ) const
{
    switch( method )
    {
        case LevelSetDerivationMethod::NODAL_PROJECTION:
            return vf::project( 
                    _space=this->functionSpaceManager()->functionSpaceScalar(),
                    _range=this->functionSpaceManager()->rangeMeshElements(),
                    _expr=sqrt( gradv(phi)*trans(gradv(phi)) )
                    );
        case LevelSetDerivationMethod::L2_PROJECTION:
            return this->projectorL2Scalar()->project( sqrt( gradv(phi)*trans(gradv(phi)) ) );
        case LevelSetDerivationMethod::SMOOTH_PROJECTION:
            return this->projectorSMScalar()->project( sqrt( gradv(phi)*trans(gradv(phi)) ) );
        //case LevelSetDerivationMethod::PN_NODAL_PROJECTION:
            //CHECK( M_useSpaceIsoPN ) << "use-space-iso-pn must be enabled to use PN_NODAL_PROJECTION \n";
            //auto phiPN = this->functionSpaceManager()->opInterpolationScalarToPN()->operator()( phi );
            //auto modGradPhiPN = vf::project(
                    //_space=this->functionSpaceManager()->functionSpaceScalarPN(),
                    //_range=this->functionSpaceManager()->rangeMeshPNElements(),
                    //_expr=sqrt( gradv(phiPN)*trans(gradv(phiPN)) )
                    //);
            //return this->functionSpaceManager()->opInterpolationScalarFromPN()->operator()( modGradPhiPN );
    }
}


}
}

#endif

