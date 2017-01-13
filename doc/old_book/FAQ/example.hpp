/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-10

  Copyright (C) 2010-2014 Feel++ Consortium

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
\page FAQDocumentationDoxygenHowto How to generate Doxygen documentation pages
\author Feel++ Consortium
\date 2013-05-06

\tableofcontents

In this section we present some useful commands to edit new pages.
You can create cross references using \ref Layout .

Use section / subsection commands to structure your page.
The following template is prefered

\section Example_Description Description
This is our first section.
You can eventually describe here the physical problem.

\subsection Description_First First sub section


\section Example_Theory Theory
Explain mathematical theory / variational formulation.
You can use Latex formulas style. There are some syntax example :

Inline formulas \f$ \delta=\int_\Omega 1 dx\f$ works and Latex keywords aswell :
\verbatim
\f{equation}
\left\{
\begin{aligned}
   -\Delta u & =  f & \text{on}\;\Omega \;, \\
            u & =  0 & \text{on}\;\partial\Omega \;,\\
\end{aligned}
\right.
\f}
where \f$u\in\Omega\f$ is the unknown "trial" function and \f$\Omega\f$ the domain.
\endverbatim
results in :
\f{equation}
\left\{
\begin{aligned}
   -\Delta u & =  f & \text{on}\;\Omega \;, \\
            u & =  0 & \text{on}\;\partial\Omega \;,\\
\end{aligned}
\right.
\f}
where \f$u\in\Omega\f$ is the unknown "trial" function and \f$\Omega\f$ the domain.

\section Example_Implementation Implementation
There we present the Feel++ code using the following command :
\verbatim
\snippet example.cpp marker1
\endverbatim
Where example.cpp contains following line :

\verbatim
/// [marker1]
int main(){
  int example = 0;
  return example;
}
/// [marker1]
\endverbatim

which results in
\snippet example.cpp marker1


\section Example_Results Results
Present your results in this section.
You can use html syntax to integrate figures.
There is an example :

<center>
<table border=0px>
<tr>
  <td>\image html mode-0.png</td>
  <td>\image html mode-1.png</td>
  <td>\image html mode-2.png</td>
</tr>
<tr>
  <td><center>first mode</center></td>
  <td><center>second mode</center></td>
  <td><center>third mode</center></td>
</tr>
<tr>
<th align="center" colspan="3">Three eigenmodes</th>
</tr>
</table>
</center>

Note that image files have to be in one of the following directory :
\li "@FEELPP_HOME_DIR@/doc/manual/pngs/"
\li "@CMAKE_CURRENT_SOURCE_DIR@/pngs/"
\li "@FEELPP_HOME_DIR@/doc/figures/"

*/


namespace Feel
{

/**
 * \brief Brief description of the class Example
 *
 * there is the detailed description
 * which can be written on several lines
 * The detailed description ends there.
 *
 * \author Feel++ Consortium
 * \see http://www/feelpp.org
 * \date 2013-05-02
 */
class Example
{
public:

    /** \name Constructors, destructor
     */
    //@{

    Example( type1 param1, type2 param2=value2 );
    Example();
    virtual ~Example();

    //@}


    /** \name Operator overloads
     */
    //@{

    /**
     * \brief assignment operator
     */
    Example& operator=( Example const& __g );

    //@}


    /** \name Accessors
     */
    //@{

    /**
     * \return data1 value
     */
    type1 accessor1() const
        {
            return M_data1;
        }

    /**
     * \brief get the value of the second data
     * \return the data2 value
     */
    type2 acessor2() const
        {
            return M_data2;
        }

    //@}

    /** \name  Mutators
     */
    //@{

    /**
     * \brief set the data1
     */
    type1 setData1( int dim )
        {
            M_data1 = dim;
        }

    //@}


    /** \name  Methods
     */
    //@{

    /**
     * \brief a brief description of methode1
     *
     * there is a detailed description of
     * the methode1
     * \param parameter1 this a description of the first parameter
     * \param parameter2 this a description of the second parameter
     * \return this is a description of the output
     */
     type1 methode1( ptype1 parameter1, ptype2 parameter2);

    //@}

private:

   // Private data wont be documented

protected:
    /// this is a short description of object1
    object1;

};

}


/// \cond DETAIL

// All the code between '\cond DETAIL' and '\endcond' wont be documented.
// Use this command to hide obsolete or useless parts.

/// \endcond

// The following boost function wont be documented
BOOST_PARAMETER_FUNCTION(
    ( typename exampletype ), // return type
    examplefun,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( param, * )
        )
    ( optional
         )
{ }
                         )
