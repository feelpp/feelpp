/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-05-25

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file stvenant_kirchhoff_base.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-05-25
 */
#ifndef __StructureBase_H
#define __StructureBase_H 1

#include <life/lifecore/life.hpp>
#include <life/lifecore/typetraits.hpp>
#include <life/lifealg/glas.hpp>
#include <life/lifecore/application.hpp>


namespace Life
{
/**
 * \class StructureBase
 * \brief Base class for structure models
 *
 *  @author Christophe Prud'homme
 *  @see
 */
class StructureBase
{
public:


    /** @name Typedefs
     */
    //@{

    typedef Life::node<double>::type node_type;

    typedef StructureBase structure_type;
    typedef boost::shared_ptr<structure_type> structure_ptrtype;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    static structure_ptrtype New( Life::po::variables_map const& vm );

    StructureBase( int d );
    StructureBase( int d, Life::po::variables_map const& vm );
    StructureBase( StructureBase const & );
    ~StructureBase();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    Life::po::variables_map vm() const { return M_vm; }

    double d() const { return M_dimension; }
    double h() const { return M_h; }

    /* time discretisation data */
    double T0() const { return M_T0; }
    double T() const { return M_T; }
    double dt() const { return M_dt; }

    int spaceOrder() const { return M_sorder; }
    int timeOrder() const { return M_torder; }

    //! penalisation parameter for weak Dirichlet handling
    double gammaBc() const { return M_gammabc; }

    int nSubSteps() const { return M_nsubsteps; }

    std::vector<std::string> const& dirichletMarkers() const { return M_dirichlet; }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    std::string createMesh();

    void print() const;

    static Life::AboutData makeAbout();
    static Life::po::options_description makeOptions();

    virtual void run() = 0;


    //@}



protected:


    Life::po::variables_map M_vm;

    int M_dimension;

    double M_h;
    double M_T0;
    double M_T;
    double M_dt;
    int M_nsubsteps;

    int M_sorder;
    int M_torder;

    double M_gammabc;

    std::vector<std::string> M_dirichlet;
    std::vector<std::string> M_neumann;

};
} // Life
#endif /* __StructureBase_H */
