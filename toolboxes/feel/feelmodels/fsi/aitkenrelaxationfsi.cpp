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
 \file aitkenrelaxationfsi.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-05
 */

#include <feel/feelmodels/fsi/aitkenrelaxationfsi.hpp>


namespace Feel
{
namespace FeelModels
{

template <class SolidType >
AitkenRelaxationFSI<SolidType>::AitkenRelaxationFSI( solid_ptrtype solid,
                                                     std::string aitkenType,
                                                     double initialTheta,
                                                     double tolPtFixe,
                                                     double minTheta )
    :
    M_solid(solid)
{

    if (M_solid->isStandardModel())
    {
        M_oldSol.reset( new element_disp_type( M_solid->functionSpaceDisplacement() ) );
        M_residual.reset( new element_disp_type( M_solid->functionSpaceDisplacement() ) );
        M_aitken = aitkenPtr(_space=M_solid->functionSpaceDisplacement(),
                             _initial_theta=initialTheta,
                             _type=aitkenType,
                             _tolerance=tolPtFixe,
                             _min_theta=minTheta);

        // redifine here! not in initialize funcion because not compile with incompressibiliy model
        /*auto currentEltSolid = vf::project(_space=M_solid->functionSpaceDisplacement(),
         _range=elements(M_solid->mesh()), 
         _expr=vf::idv(M_solid->getDisplacement()));*/
        M_aitken->initialize( _residual=*M_residual,
                              _currentElt=/**M_residual*//*currentEltSolid*/M_solid->fieldDisplacement() );

        //M_areaFSIinterface = integrate(_range=markedfaces(M_oldSol->mesh(),M_solid->getMarkerNameFSI().front()),
        //_expr=cst(1.)).evaluate()(0,0);

    }
    else if ( M_solid->is1dReducedModel() )
    {
        M_oldSol1dReduced.reset( new element_disp_1dreduced_type( M_solid->solid1dReduced()->fieldDisplacementScal1dReduced().functionSpace() ) );
        M_residual1dReduced.reset( new element_disp_1dreduced_type( M_solid->solid1dReduced()->fieldDisplacementScal1dReduced().functionSpace() ) );
        M_aitken1dReduced = aitkenPtr(_space=M_solid->solid1dReduced()->fieldDisplacementScal1dReduced().functionSpace(),
                                      _initial_theta=initialTheta,
                                      _type=aitkenType,
                                      _tolerance=tolPtFixe,
                                      _min_theta=minTheta);
        M_aitken1dReduced->initialize( _residual=*M_residual1dReduced,
                                       _currentElt=M_solid->solid1dReduced()->fieldDisplacementScal1dReduced() );
    }
}

//-----------------------------------------------------------------------------------//

template <class SolidType >
void
AitkenRelaxationFSI<SolidType>::saveOldSolution()
{
    if (M_solid->isStandardModel())
    {
        M_oldSol->on(_range=elements(M_oldSol->mesh()),
                     //_range=markedfaces(M_oldSol->mesh(),"ParoiFSI"),
                     _expr=vf::idv(M_solid->fieldDisplacement()) );
    }
    else if ( M_solid->is1dReducedModel() )
    {
        M_oldSol1dReduced->on(_range=elements(M_oldSol1dReduced->mesh()),
                              _expr=vf::idv(M_solid->solid1dReduced()->fieldDisplacementScal1dReduced()) );
    }
}

//-----------------------------------------------------------------------------------//

template <class SolidType >
bool
AitkenRelaxationFSI<SolidType>::isFinished()
{
    bool res=false;
    if (M_solid->isStandardModel())
        res=M_aitken->isFinished();
    else if ( M_solid->is1dReducedModel() )
        res=M_aitken1dReduced->isFinished();
    return res;
}

//-----------------------------------------------------------------------------------//

template <class SolidType >
void
AitkenRelaxationFSI<SolidType>::restart()
{
    if (M_solid->isStandardModel())
        M_aitken->restart();
    else if ( M_solid->is1dReducedModel() )
        M_aitken1dReduced->restart();
}

//-----------------------------------------------------------------------------------//

template <class SolidType >
void
AitkenRelaxationFSI<SolidType>::printInfo()
{
    std::cout << std::scientific;
    if (M_solid->isStandardModel())
        M_aitken->printInfo();
    else if ( M_solid->is1dReducedModel() )
        M_aitken1dReduced->printInfo();
}

//-----------------------------------------------------------------------------------//

template <class SolidType >
void
AitkenRelaxationFSI<SolidType>::applyRelaxation()
{
    if (M_solid->isStandardModel())
    {
        M_residual->on(_range=elements(M_residual->mesh()),
                       //_range=markedfaces(M_oldSol->mesh(),"ParoiFSI"),
                       _expr=vf::idv(M_solid->fieldDisplacement() ));
        *M_residual -= *M_oldSol;

        M_aitken->apply2(_newElt=M_solid->fieldDisplacement(),
                         _residual=*M_residual,
                         _currentElt=M_solid->fieldDisplacement() );
    }
    else if ( M_solid->is1dReducedModel() )
    {
        M_residual1dReduced->on(_range=elements(M_residual1dReduced->mesh()),
                                _expr=vf::idv(M_solid->solid1dReduced()->fieldDisplacementScal1dReduced() ));
        *M_residual1dReduced -= *M_oldSol1dReduced;

        M_aitken1dReduced->apply2(_newElt=M_solid->solid1dReduced()->fieldDisplacementScal1dReduced(),
                                  _residual=*M_residual1dReduced,
                                  _currentElt=M_solid->solid1dReduced()->fieldDisplacementScal1dReduced() );
    }
}

//-----------------------------------------------------------------------------------//

template <class SolidType >
void
AitkenRelaxationFSI<SolidType>::shiftRight()
{
    if (M_solid->isStandardModel())
        M_aitken->shiftRight();
    else if ( M_solid->is1dReducedModel() )
        M_aitken1dReduced->shiftRight();
}

//-----------------------------------------------------------------------------------//

template <class SolidType >
uint16_type
AitkenRelaxationFSI<SolidType>::nIterations()
{
    uint16_type res=0;
    if (M_solid->isStandardModel())
        res = M_aitken->nIterations();
    else if ( M_solid->is1dReducedModel() )
        res = M_aitken1dReduced->nIterations();

    return res;
}

//-----------------------------------------------------------------------------------//

template <class SolidType >
double
AitkenRelaxationFSI<SolidType>::residualNorm()
{
    double res=0;
    if (M_solid->isStandardModel())
        res=M_aitken->residualNorm();
    else if ( M_solid->is1dReducedModel() )
        res=M_aitken1dReduced->residualNorm();

    return res;
}

//-----------------------------------------------------------------------------------//

template <class SolidType >
void
AitkenRelaxationFSI<SolidType>::setTheta(double v)
{
    if (M_solid->isStandardModel())
        M_aitken->setTheta(v);
    else if (M_solid->is1dReducedModel())
        M_aitken1dReduced->setTheta(v);
}

//-----------------------------------------------------------------------------------//

template <class SolidType >
FixPointConvergenceFSI<SolidType>::FixPointConvergenceFSI( solid_ptrtype solid )
    :
    M_solid(solid)
{
    if (M_solid->isStandardModel())
    {
        M_oldSol.reset( new element_disp_type( M_solid->functionSpaceDisplacement() ) );
        M_residual.reset( new element_disp_type( M_solid->functionSpaceDisplacement() ) );
    }
    else if ( M_solid->is1dReducedModel() )
    {
        M_oldSol1dReduced.reset( new element_disp_1dreduced_type( M_solid->solid1dReduced()->fieldDisplacementScal1dReduced().functionSpace() ) );
        M_residual1dReduced.reset( new element_disp_1dreduced_type( M_solid->solid1dReduced()->fieldDisplacementScal1dReduced().functionSpace() ) );
    }
}

template <class SolidType >
void
FixPointConvergenceFSI<SolidType>::saveOldSolution()
{
    if (M_solid->isStandardModel())
    {
        M_oldSol->on(_range=elements(M_oldSol->mesh()),
                     //_range=markedfaces(M_oldSol->mesh(),"ParoiFSI"),
                     _expr=vf::idv(M_solid->fieldDisplacement()) );
    }
    else if ( M_solid->is1dReducedModel() )
    {
        M_oldSol1dReduced->on(_range=elements(M_oldSol1dReduced->mesh()),
                              _expr=vf::idv(M_solid->solid1dReduced()->fieldDisplacementScal1dReduced()) );
    }
}

template <class SolidType >
double
FixPointConvergenceFSI<SolidType>::computeConvergence()
{
    double residualConvergence=1;

    if (M_solid->isStandardModel())
    {
        M_residual->on(_range=elements(M_residual->mesh()),
                       _expr=vf::idv(M_solid->fieldDisplacement() ));
        *M_residual -= *M_oldSol;
        double oldEltL2Norm = M_oldSol->l2Norm();
        if ( oldEltL2Norm > 1e-13 )
            residualConvergence = M_residual->l2Norm()/oldEltL2Norm;
        else
            residualConvergence = M_residual->l2Norm();
    }
    else if ( M_solid->is1dReducedModel() )
    {
        M_residual1dReduced->on(_range=elements(M_residual1dReduced->mesh()),
                                _expr=vf::idv(M_solid->solid1dReduced()->fieldDisplacementScal1dReduced() ));
        *M_residual1dReduced -= *M_oldSol1dReduced;
        double oldEltL2Norm = M_oldSol1dReduced->l2Norm();
        if ( oldEltL2Norm > 1e-13 )
            residualConvergence = M_residual1dReduced->l2Norm()/oldEltL2Norm;
        else
            residualConvergence = M_residual1dReduced->l2Norm();
    }

    return residualConvergence;
}


} // namespace FeelModels
} // namespace Feel
