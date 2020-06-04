/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2011-07-05

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
 \file aitkenrelaxationfsi.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-05
 */

#ifndef FEELPP_MODELS_AITKENRELAXATION_FSI_H
#define FEELPP_MODELS_AITKENRELAXATION_FSI_H 1


#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/aitken.hpp>


namespace Feel
{
namespace FeelModels
{


template <class SolidType >
class AitkenRelaxationFSI
{
public :
    typedef AitkenRelaxationFSI<SolidType> self_type;
    typedef SolidType solid_type;
    typedef std::shared_ptr<solid_type> solid_ptrtype;

    typedef typename solid_type::space_displacement_type space_disp_type;
    typedef typename solid_type::element_displacement_type element_disp_type;
    typedef typename solid_type::element_displacement_ptrtype element_disp_ptrtype;
    typedef Aitken<space_disp_type> aitken_type;
    typedef std::shared_ptr<aitken_type> aitken_ptrtype;

    typedef typename solid_type::solid_1dreduced_type::space_displacement_component_type space_disp_1dreduced_type;
    typedef typename solid_type::solid_1dreduced_type::element_displacement_component_type element_disp_1dreduced_type;
    typedef typename solid_type::solid_1dreduced_type::element_displacement_component_ptrtype element_disp_1dreduced_ptrtype;
    typedef Aitken<space_disp_1dreduced_type> aitken_1dreduced_type;
    typedef std::shared_ptr<aitken_1dreduced_type> aitken_1dreduced_ptrtype;

    //-----------------------------------------------------------------------------------//

    AitkenRelaxationFSI( solid_ptrtype solid,
                         std::string aitkenType = "method1",
                         double initialTheta = 1.0,
                         double tolPtFixe =1.0e-6,
                         double minTheta =1e-4 );
    AitkenRelaxationFSI( self_type const & M ) = default;

    //-----------------------------------------------------------------------------------//

    void saveOldSolution();

    bool isFinished();

    void restart();

    void printInfo();

    void applyRelaxation();

    void shiftRight();

    uint16_type nIterations();

    double residualNorm();

    void setTheta(double v);

private :

    solid_ptrtype M_solid;

    aitken_ptrtype M_aitken;
    element_disp_ptrtype M_oldSol;
    element_disp_ptrtype M_residual;

    aitken_1dreduced_ptrtype M_aitken1dReduced;
    element_disp_1dreduced_ptrtype M_oldSol1dReduced;
    element_disp_1dreduced_ptrtype M_residual1dReduced;

    //double M_areaFSIinterface;
};

template <class SolidType >
class FixPointConvergenceFSI
{
public :
    typedef FixPointConvergenceFSI<SolidType> self_type;
    typedef SolidType solid_type;
    typedef std::shared_ptr<solid_type> solid_ptrtype;

    typedef typename solid_type::space_displacement_type space_disp_type;
    typedef typename solid_type::element_displacement_type element_disp_type;
    typedef typename solid_type::element_displacement_ptrtype element_disp_ptrtype;

    typedef typename solid_type::solid_1dreduced_type::space_displacement_component_type space_disp_1dreduced_type;
    typedef typename solid_type::solid_1dreduced_type::element_displacement_component_type element_disp_1dreduced_type;
    typedef typename solid_type::solid_1dreduced_type::element_displacement_component_ptrtype element_disp_1dreduced_ptrtype;

    FixPointConvergenceFSI( solid_ptrtype solid );
    FixPointConvergenceFSI( self_type const & M ) = default;

    void saveOldSolution();
    double computeConvergence();

private :
    solid_ptrtype M_solid;
    element_disp_ptrtype M_oldSol;
    element_disp_ptrtype M_residual;
    element_disp_1dreduced_ptrtype M_oldSol1dReduced;
    element_disp_1dreduced_ptrtype M_residual1dReduced;
};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_AITKENRELAXATION_FSI_H
