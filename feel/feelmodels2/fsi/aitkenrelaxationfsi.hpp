/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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


template <class  SolidMech >
class AitkenRelaxationFSI

{
public :

    typedef SolidMech solid_type;
    typedef boost::shared_ptr<solid_type> solid_ptrtype;

    typedef typename solid_type::space_displacement_type space_disp_type;
    //typedef typename solid_type::element_displacement_type element_disp_type;
    typedef typename solid_type::element_vectorial_type element_disp_type;
    typedef boost::shared_ptr<element_disp_type> element_disp_ptrtype;

    typedef Aitken<space_disp_type> aitken_type;
    typedef boost::shared_ptr<aitken_type> aitken_ptrtype;

    typedef typename solid_type::space_1dreduced_type space_disp_1dreduced_type;
    typedef typename solid_type::element_1dreduced_type element_disp_1dreduced_type;
    typedef boost::shared_ptr<element_disp_1dreduced_type> element_disp_1dreduced_ptrtype;

    typedef Aitken<space_disp_1dreduced_type> aitken_1dreduced_type;
    typedef boost::shared_ptr<aitken_1dreduced_type> aitken_1dreduced_ptrtype;

    //-----------------------------------------------------------------------------------//

    AitkenRelaxationFSI( solid_ptrtype solid,
                         std::string aitkenType = "method1",
                         double initialTheta = 1.0,
                         double tolPtFixe =1.0e-6,
                         double minTheta =1e-4 )
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
                M_oldSol1dReduced.reset( new element_disp_1dreduced_type( M_solid->fieldDisplacementScal1dReduced().functionSpace() ) );
                M_residual1dReduced.reset( new element_disp_1dreduced_type( M_solid->fieldDisplacementScal1dReduced().functionSpace() ) );
                M_aitken1dReduced = aitkenPtr(_space=M_solid->fieldDisplacementScal1dReduced().functionSpace(),
                                              _initial_theta=initialTheta,
                                              _type=aitkenType,
                                              _tolerance=tolPtFixe,
                                              _min_theta=minTheta);
                M_aitken1dReduced->initialize( _residual=*M_residual1dReduced,
                                               _currentElt=M_solid->fieldDisplacementScal1dReduced() );
            }
        }

    //-----------------------------------------------------------------------------------//

    void saveOldSolution()
        {
            if (M_solid->isStandardModel())
            {
                *M_oldSol = Feel::vf::project(M_oldSol->functionSpace(),
                                              elements(M_oldSol->mesh()),
                                              //markedfaces(M_oldSol->mesh(),"ParoiFSI"),
                                              vf::idv(M_solid->fieldDisplacement()) );
            }
            else if ( M_solid->is1dReducedModel() )
            {
                *M_oldSol1dReduced = Feel::vf::project(M_oldSol1dReduced->functionSpace(),
                                                       elements(M_oldSol1dReduced->mesh()),
                                                       vf::idv(M_solid->fieldDisplacementScal1dReduced()) );
            }
        }

    //-----------------------------------------------------------------------------------//

    bool isFinished()
        {
            bool res=false;
            if (M_solid->isStandardModel())
                res=M_aitken->isFinished();
            else if ( M_solid->is1dReducedModel() )
                res=M_aitken1dReduced->isFinished();
            return res;
        }

    //-----------------------------------------------------------------------------------//

    void restart()
        {
            if (M_solid->isStandardModel())
                M_aitken->restart();
            else if ( M_solid->is1dReducedModel() )
                M_aitken1dReduced->restart();
        }

    //-----------------------------------------------------------------------------------//

    void printInfo()
        {
            if (M_solid->isStandardModel())
                M_aitken->printInfo();
            else if ( M_solid->is1dReducedModel() )
                M_aitken1dReduced->printInfo();
        }

    //-----------------------------------------------------------------------------------//

    void applyRelaxation()
        {
            if (M_solid->isStandardModel())
            {
                *M_residual = vf::project(M_residual->functionSpace(),
                                          elements(M_residual->mesh()),
                                          //markedfaces(M_oldSol->mesh(),"ParoiFSI"),
                                          vf::idv(M_solid->fieldDisplacement() ));
                *M_residual -= *M_oldSol;

                M_aitken->apply2(_newElt=M_solid->fieldDisplacement(),
                                 _residual=*M_residual,
                                 _currentElt=M_solid->fieldDisplacement() );
            }
            else if ( M_solid->is1dReducedModel() )
            {
                *M_residual1dReduced = vf::project(M_residual1dReduced->functionSpace(),
                                                   elements(M_residual1dReduced->mesh()),
                                                   vf::idv(M_solid->fieldDisplacementScal1dReduced() ));
                *M_residual1dReduced -= *M_oldSol1dReduced;

                M_aitken1dReduced->apply2(_newElt=M_solid->fieldDisplacementScal1dReduced(),
                                          _residual=*M_residual1dReduced,
                                          _currentElt=M_solid->fieldDisplacementScal1dReduced() );
            }
        }

    //-----------------------------------------------------------------------------------//

    void shiftRight()
        {
            if (M_solid->isStandardModel())
                M_aitken->shiftRight();
            else if ( M_solid->is1dReducedModel() )
                M_aitken1dReduced->shiftRight();
        }

    //-----------------------------------------------------------------------------------//

    uint16_type nIterations()
        {
            uint16_type res=0;
            if (M_solid->isStandardModel())
                res = M_aitken->nIterations();
            else if ( M_solid->is1dReducedModel() )
                res = M_aitken1dReduced->nIterations();

            return res;
        }

    //-----------------------------------------------------------------------------------//

    double residualNorm()
        {
            double res=0;
            if (M_solid->isStandardModel())
                res=M_aitken->residualNorm();
            else if ( M_solid->is1dReducedModel() )
                res=M_aitken1dReduced->residualNorm();

            return res;
        }

    //-----------------------------------------------------------------------------------//

    void setTheta(double v)
        {
            if (M_solid->isStandardModel())
                M_aitken->setTheta(v);
            else if (M_solid->is1dReducedModel())
                M_aitken1dReduced->setTheta(v);
        }

    //-----------------------------------------------------------------------------------//

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

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_AITKENRELAXATION_FSI_H
