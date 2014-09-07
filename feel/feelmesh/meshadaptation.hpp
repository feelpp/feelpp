/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
       Date: 2011-16-12

       Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)
       Copyright (C) CNRS 2012

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file meshadaptation.hpp
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2012-08-04
 */

#ifndef __MESHADAPTATION_HPP
#define __MESHADAPTATION_HPP 1

/** include linear algebra backend */
#include <feel/feelalg/backend.hpp>

/** include function space class */
#include <feel/feeldiscr/functionspace.hpp>

/** include gmsh mesh importer */
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>

#include <feel/feelfilters/exporter.hpp>

/** include  the header for the variational formulation language (vf) aka FEEL++ */
#include <feel/feelvf/vf.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feeldiscr/elementdiv.hpp>

#include <boost/assert.hpp>
#include <boost/regex.hpp>

#include <fstream>

#if defined( FEELPP_HAS_GMSH_H )
#include <Gmsh.h>
#include <GModel.h>
#include <Context.h>
#include <MElement.h>
#include <OpenFile.h>
#include <StringUtils.h>
#include <Field.h>
#include <PView.h>
#include <PViewData.h>
#endif

namespace Feel
{

    BOOST_PARAMETER_NAME(initMesh)
    BOOST_PARAMETER_NAME(geofile)
    BOOST_PARAMETER_NAME(adaptType)
    BOOST_PARAMETER_NAME(var)
    BOOST_PARAMETER_NAME(metric)
    BOOST_PARAMETER_NAME(hMin)
    BOOST_PARAMETER_NAME(hMax)
    BOOST_PARAMETER_NAME(tol)

    template<int Dim,
             int Order,
             int OrderGeo>
    class MeshAdaptation
    {
    public:
        static const bool isP1 = (Order == 1);

        //typedef double value_type;
        typedef Eigen::Matrix<double, Dim, 1> vectorN_type;
        typedef Eigen::Matrix<double, Dim, Dim> matrixN_type;

        //! linear algebra backend factory
        typedef Backend<double> backend_type;
        //! linear algebra backend factory shared_ptr<> type
        typedef boost::shared_ptr<backend_type> backend_ptrtype;

        //! sparse matrix type associated with backend
        typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
        //! sparse matrix type associated with backend (shared_ptr<> type)
        typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
        //! vector type associated with backend
        typedef typename backend_type::vector_type vector_type;
        //! vector type associated with backend (shared_ptr<> type)
        typedef typename backend_type::vector_ptrtype vector_ptrtype;

        typedef Simplex<Dim,OrderGeo> convex_type;
        //! mesh type
        typedef Mesh<convex_type> mesh_type;
        //! mesh shared_ptr<> type
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

        //! the exporter factory type
        typedef Exporter<mesh_type> export_type;
        //! the exporter factory (shared_ptr<> type)
        typedef boost::shared_ptr<export_type> export_ptrtype;

        typedef bases<Lagrange<Order, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type> space_type;
        //! the approximation function space type (shared_ptr<> type)
        typedef boost::shared_ptr<space_type> space_ptrtype;
        //! an element type of the approximation function space
        typedef typename space_type::element_type element_type;

        //! Scalar P0 space
        typedef bases<Lagrange<0,Scalar, Discontinuous > > p0_basis_type;
        typedef FunctionSpace<mesh_type, p0_basis_type> p0_space_type;
        typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
        typedef typename p0_space_type::element_type p0_element_type;

        //! Scalar P1 space
        typedef bases<Lagrange<1, Scalar> > p1_basis_type;
        typedef FunctionSpace<mesh_type, p1_basis_type> p1_space_type;
        typedef boost::shared_ptr<p1_space_type> p1_space_ptrtype;
        typedef typename p1_space_type::element_type p1_element_type;

        //! Vectorial P1 space
        typedef bases<Lagrange<1, Vectorial> > p1vec_basis_type;
        typedef FunctionSpace<mesh_type, p1vec_basis_type> p1vec_space_type;
        typedef boost::shared_ptr<p1vec_space_type> p1vec_space_ptrtype;
        typedef typename p1vec_space_type::element_type p1vec_element_type;

        typedef typename mesh_type::point_type point_type;

        //! Constructor
        MeshAdaptation()
        {
            // build default value for _var parameter of interface
            element_type initVar;
            std::pair<element_type, std::string> initPairVar = std::make_pair( initVar, "defaultVar");
            M_defaultVar.push_back(initPairVar);

            // build default value for _metric parameter of interface
            p1_element_type initMetric1;
            std::vector<p1_element_type> initMetric;
            initMetric.push_back(initMetric1);
            std::pair<std::vector<p1_element_type>, std::string> initPairMetric = std::make_pair(initMetric, "defaultMetric");
            M_defaultMetric.push_back(initPairMetric);
        }

        //! Set of measures on a dof
        std::vector<vectorN_type> const measures()
        {
            return M_measures;
        }

        //! Set of directions on a dof
        std::vector<matrixN_type> const directions()
        {
            return M_directions;
        }

        //! Interface
        mesh_ptrtype adaptMeshImpl( const mesh_ptrtype& initMesh, std::string geofile, std::string adaptType,
                                         std::list< std::pair<element_type, std::string> > var,
                                         std::list< std::pair< std::vector<p1_element_type>, std::string> > metric,
                                         double hMin, double hMax, double tol);

        //! Generates a GMSH post processing file (for mesh adaptation)
        std::string createPosfile(std::string name_var, const p1_element_type& bbNewMap, const mesh_ptrtype& mesh);
        std::string createPosfileAnisotropic(std::string nameVar, const std::vector<p1_element_type>& bbNewMap, const mesh_ptrtype& mesh);

        //! Adapt geofile for mesh adaptation
        std::string createAdaptedGeo(std::string geofile, std::string name, std::vector<std::string> posfiles, bool aniso);

        //! Build adapted mesh
        std::string buildAdaptedMesh(std::string geofile, std::string name, std::vector<std::string> posfiles , bool aniso);

        //! Compute the metric for mesh adaptation using hessian matrix
        void computeMetric(const double tol, const double h_min, const double h_max, const matrixN_type & hessian_matrix,
                           matrixN_type & M, double & max_eigenvalue, int dofId);

        //! Mesh adaptation from Hessian matrix
        //1 : proj u Pk -> P1 and Hessian P1
        //2 : Hessian P(k-2) and proj hessian P(k-2) -> P1 (possibility to choose L2 or H1 projection)
        std::string adaptMeshHess1(element_type& U, const mesh_ptrtype& mesh, double hMin, double hMax,
                                   std::string name, std::string geofile, double tol, bool aniso);

        std::string adaptMeshHess2(element_type& U, const mesh_ptrtype& mesh, double hMin, double hMax,
                                   std::string name, std::string geofile, double tol, bool aniso);
        std::string adaptMeshHess2(element_type& U, const mesh_ptrtype& mesh, double hMin, double hMax,
                                   std::string name, std::string geofile, double tol, bool aniso,
                                   mpl::bool_<true>);
        std::string adaptMeshHess2(element_type& U, const mesh_ptrtype& mesh, double hMin, double hMax,
                                   std::string name, std::string geofile, double tol, bool aniso,
                                   mpl::bool_<false>);

        //! Mesh adaptation interface
        // _initMesh = initial mesh
        // _geofile = complete path to access geofile
        // _adaptType = isotropic | anisotropic
        // _var = list of (element_type, string) = list of (variable, name of the variable)
        // _metric = list of (vector<p1_element_type>, string) = list of (hsize(isotropic) | matrix(anisotropic), name for metric)
        // _hMin, _hMax = limits for hsize
        // _tol = geometric tolerance (for adaptation from hessian)

        BOOST_PARAMETER_MEMBER_FUNCTION(
                                        (mesh_ptrtype),
                                        adaptMesh,
                                        tag,
                                        (required
                                         (initMesh,(mesh_ptrtype)) // initial mesh
                                         (geofile, (std::string)) // geometry
                                         (adaptType, (std::string))) //type of adaptation : isotropic | anisotropic
                                        (optional
                                         (var, (std::list< std::pair<element_type, std::string> >), M_defaultVar)
                                         (metric, (std::list< std::pair< std::vector<p1_element_type>, std::string> >), M_defaultMetric)
                                         (hMin, (double), 1.0e-3)
                                         (hMax, (double), 100.0)
                                         (tol, (double), 1.0))
                                        )
        {
            mesh_ptrtype newMesh = this-> adaptMeshImpl( initMesh, geofile, adaptType, var, metric, hMin, hMax, tol);
            return newMesh;
        }

    private :
        std::list< std::pair<element_type, std::string> > M_defaultVar;
        std::list< std::pair<std::vector<p1_element_type>, std::string> > M_defaultMetric;

        std::vector<vectorN_type> M_measures;
        std::vector<matrixN_type> M_directions;

    };

    template<int Dim,
             int Order,
             int OrderGeo>
    const bool
    MeshAdaptation<Dim, Order, OrderGeo>::isP1;

    template<int Dim,
             int Order,
             int OrderGeo>
    typename MeshAdaptation<Dim, Order, OrderGeo>::mesh_ptrtype
    MeshAdaptation<Dim,
                   Order,
                   OrderGeo>::adaptMeshImpl( const mesh_ptrtype& initMesh, std::string geofile, std::string adaptType,
                                         std::list< std::pair<element_type, std::string> > var,
                                         std::list< std::pair< std::vector<p1_element_type>, std::string> > metric,
                                         double hMin, double hMax, double tol)
    {
        std::vector<std::string> posfiles;
        std::string finalName = "";
        if ( !(var.front().second == "defaultVar") ) // adptation from a variable
            {
                boost::for_each( var, [&initMesh, &hMin, &hMax, &geofile, &tol, &adaptType, &posfiles, &finalName, this]
                                 (std::pair<element_type, std::string> _var)
                                 {
                                     std::string newPosfile = this->adaptMeshHess2( _var.first, initMesh, hMin, hMax, _var.second, geofile,
                                                                              tol, adaptType=="anisotropic");
                                     posfiles.push_back( newPosfile );
                                     finalName += "-" + _var.second;
                                 });
            }

        else if ( !(metric.front().second == "defaultMetric") ) //adaptation from a metric
            {
                boost::for_each(metric, [&initMesh, &adaptType, &posfiles, &finalName, this]
                                (std::pair<std::vector<p1_element_type>, std::string> _metric)
                                {
                                    std::string newPosfile;
                                    if (adaptType == "isotropic")
                                        {
                                            BOOST_ASSERT(_metric.first.size()==1);
                                            newPosfile = this->createPosfile(_metric.second, _metric.first[0], initMesh);
                                        }
                                    else if( adaptType == "anisotropic")
                                        {
                                            BOOST_ASSERT(_metric.first.size()==Dim*Dim);
                                            newPosfile = this->createPosfileAnisotropic(_metric.second, _metric.first, initMesh);
                                        }
                                    posfiles.push_back( newPosfile );
                                    finalName += "-" + _metric.second;
                                });
            }

        std::string newMeshName = this->buildAdaptedMesh( geofile, finalName, posfiles, adaptType=="anisotropic");

        // Update the mesh
        LOG(INFO) << "load the new mesh...\n";
        mesh_ptrtype newMesh = mesh_type::New();
        newMesh = loadGMSHMesh( _mesh = new mesh_type,
                                _filename = newMeshName,
                                _update=MESH_UPDATE_FACES | MESH_UPDATE_EDGES);
        LOG(INFO) << "new mesh loaded \n";

        return newMesh;
    }

    template<int Dim,
             int Order,
             int OrderGeo>
    std::string
    MeshAdaptation<Dim,
                   Order,
                   OrderGeo>::createPosfile(std::string nameVar, const p1_element_type& bbNewMap, const mesh_ptrtype& mesh)
    {
        using namespace Feel::vf;

        std::string posFormat = "pos";
        std::string mshFormat = "msh";

        p1_space_ptrtype P1h = p1_space_type::New( mesh );

        std::ofstream newPosFile( (boost::format( "%1%.%2%" ) % nameVar %posFormat ).str() );

        newPosFile << "View \" background mesh \" { \n";
        auto eltIt = mesh->beginElement();
        auto eltEnd = mesh->endElement();
        for ( ; eltIt != eltEnd; eltIt++)
            {
                std::vector<point_type> eltPoints( eltIt->nPoints() );
                if (Dim == 2)
                    newPosFile << "ST(";
                if (Dim == 3)
                    newPosFile << "SS(";
                for (int i=0; i < eltIt->nPoints(); i++)
                    {
                        eltPoints[i] = eltIt->point(i);
                        for (int j=0; j< Dim; j++)
                            {
                                newPosFile << eltPoints[i](j);
                                if ( Dim == 2 || (Dim == 3 &&  j!=Dim-1) )
                                    {
                                        newPosFile << ", ";
                                    }
                            }

                        if (Dim == 2)
                            newPosFile << "0"; //2D case => z coordinate = 0
                        if (i!= eltIt->nPoints() - 1)
                            newPosFile << ", ";
                    }

                newPosFile << "){";

                for (size_t k=0; k<eltPoints.size(); k++)
                    {
                        auto dofIndex = boost::get<0>( P1h->dof()->localToGlobal( eltIt->id(), k) );
                        newPosFile << bbNewMap[ dofIndex ];
                        if ( k!= eltPoints.size() - 1)
                            newPosFile << ", ";
                    }

                newPosFile << "};\n";
                eltPoints.clear();
            }

        newPosFile << "};";
        newPosFile.close();

        return (boost::format( "%1%.%2%" ) % nameVar %posFormat ).str();
    }

    template<int Dim,
             int Order,
             int OrderGeo>
    std::string
    MeshAdaptation<Dim,
                   Order,
                   OrderGeo>::createPosfileAnisotropic(std::string nameVar, const std::vector<p1_element_type>& bbNewMap,
                                                       const mesh_ptrtype& mesh)
    {
        using namespace Feel::vf;

        std::string posFormat = "pos";
        std::string mshFormat = "msh";

        p1_space_ptrtype P1h = p1_space_type::New( mesh );

        std::ofstream newPosFile( (boost::format( "%1%.%2%" ) % nameVar %posFormat ).str() );

        newPosFile << "View \" background mesh \" { \n";
        auto eltIt = mesh->beginElement();
        auto eltEnd = mesh->endElement();
        for ( ; eltIt != eltEnd; eltIt++)
            {
                std::vector<point_type> eltPoints( eltIt->nPoints() );
                newPosFile << "T";
                if (Dim == 2)
                    newPosFile << "T";
                else
                    newPosFile << "S";

                newPosFile << "(";
                for (int i=0; i < eltIt->nPoints(); i++)
                    {
                        eltPoints[i] = eltIt->point(i);
                        for (int j=0; j< Dim; j++)
                            {
                                newPosFile << eltPoints[i](j);
                                if ( Dim == 2 || (Dim == 3 &&  j!=Dim-1) )
                                    newPosFile << ", ";
                            }
                        if (Dim == 2)
                            newPosFile << "0"; //2D case => z coordinate = 0
                        if (i!= eltIt->nPoints() - 1)
                            newPosFile << ", ";
                    }

                newPosFile << "){";

                for (size_t k=0; k<eltPoints.size(); k++)
                    {
                        auto dofIndex = boost::get<0>( P1h->dof()->localToGlobal( eltIt->id(), k) );

                        int num = 0;
                        newPosFile << (bbNewMap[num++])[ dofIndex ];
                        newPosFile << ", ";

                        newPosFile << (bbNewMap[num++])[ dofIndex ] << ", ";

                        if (Dim == 3)
                            newPosFile << (bbNewMap[num++])[ dofIndex ] << ", ";
                        else
                            newPosFile << "0, ";

                        newPosFile << (bbNewMap[num++])[ dofIndex ] << ", ";
                        newPosFile << (bbNewMap[num++])[ dofIndex ] << ", ";

                        if (Dim == 3)
                            {
                                newPosFile << (bbNewMap[num++])[ dofIndex ] << ", ";
                                newPosFile << (bbNewMap[num++])[ dofIndex ] << ", ";
                                newPosFile << (bbNewMap[num++])[ dofIndex ] << ", ";
                                newPosFile << (bbNewMap[num++])[ dofIndex ];
                            }
                        else
                            {
                                newPosFile << "0, ";
                                newPosFile << "0, ";
                                newPosFile << "0, ";
                                newPosFile << "1";
                            }

                        if (k!= eltPoints.size() - 1)
                            newPosFile << ", ";
                    }
                newPosFile << "};\n";
                eltPoints.clear();
            }

        newPosFile << "};";
        newPosFile.close();

        return (boost::format( "%1%.%2%" ) % nameVar %posFormat ).str();
    }

    template<int Dim,
             int Order,
             int OrderGeo>
    std::string
    MeshAdaptation<Dim,
                   Order,
                   OrderGeo>::createAdaptedGeo(std::string geofile, std::string name, std::vector<std::string> posfiles,
                                               bool aniso)
    {
        ///// Transform inupt geofile into string
        std::ifstream stringToGeo;
        std::string geofileString;
        stringToGeo.open(geofile);
        stringToGeo.seekg(0, std::ios::end);
        geofileString.reserve(stringToGeo.tellg());
        stringToGeo.seekg(0, std::ios::beg);
        geofileString.assign((std::istreambuf_iterator<char>(stringToGeo)),
                             std::istreambuf_iterator<char>());
        stringToGeo.close();

        int nbPosfiles = posfiles.size();

        //// Detect lines corresponding to algorithm used, fields list, and background field
        std::ostringstream __exprAlgo, __exprFlist, __exprBgm;
        __exprAlgo << "Field.([[:blank:]]*)([[:digit:]]*)([[:blank:]]*).([[:blank:]]*)=([[:blank:]]*)(Min|MinAniso);";
        __exprFlist << "Field.([[:blank:]]*)([[:digit:]]*)([[:blank:]]*)..FieldsList([[:blank:]]*)=([[:blank:]]*).(([[:digit:]].?)*).;";
        __exprBgm << "Background[[:blank:]]Field([[:blank:]]*)=([[:blank:]]*)([[:digit:]]*);";

        boost::regex exprAlgo(__exprAlgo.str().c_str());
        boost::regex exprFlist(__exprFlist.str().c_str());
        boost::regex exprBgm(__exprBgm.str().c_str());
        std::list<boost::regex> endExpr = boost::assign::list_of(exprAlgo)(exprFlist)(exprBgm);

        bool isInserted = false;

        /// Check each field appears in geo file (add missing ones)
        int i=1;
        for (; i <= nbPosfiles ; i++)
            {
                std::ostringstream __expr1, __expr2;
                __expr1 << "Field.([[:blank:]]*)"<< i <<"([[:blank:]]*).([[:blank:]]*)=([[:blank:]]*)PostView;" ;
                __expr2 << "Field.([[:blank:]]*)" << i << "([[:blank:]]*)..IView([[:blank:]]*)=([[:blank:]]*)" << i-1 << ";" ;

                boost::regex expression1( __expr1.str().c_str() );
                boost::regex expression2( __expr2.str().c_str() );

                boost::match_results<std::string::iterator> posFields;
                bool fieldIsFoundL1 = boost::regex_search(geofileString.begin(), geofileString.end(), posFields, expression1);
                bool fieldIsFoundL2 = boost::regex_search(geofileString.begin(), geofileString.end(), posFields, expression2);

                /// If field linked with i-eme posfile doesn't appear => add it
                if (!(fieldIsFoundL1 && fieldIsFoundL2))
                    {
                        std::ostringstream __newExpr;
                        __newExpr << "Field[" << i << "] = PostView; \n";
                        __newExpr << "Field[" << i << "].IView = " << i-1 << ";\n";
                        std::string newExpr = __newExpr.str().c_str();

                        //// Place to add new field (always before algo, field_list and background lines)
                        boost::for_each( endExpr, [&isInserted, &geofileString, &newExpr]( boost::regex _expr)
                                         {
                                             boost::match_results<std::string::iterator> itExpr;
                                             bool exprIsFound = boost::regex_search(geofileString.begin(), geofileString.end(), itExpr, _expr);
                                             if (exprIsFound)
                                                 {
                                                     geofileString.insert(itExpr[0].first, newExpr.begin(), newExpr.end() );
                                                     isInserted = true;
                                                     return;
                                                 }
                                         });
                        if (!isInserted)
                            geofileString.insert(geofileString.end(), newExpr.begin(), newExpr.end() );
                    }
            }

        ///// Check if algo, field list and background are defined
        std::ostringstream __newAlgo, __newList, __newBgm;

        if (aniso)
            __newAlgo << "Field[" << i << "] = MinAniso;\n";
        else
            __newAlgo << "Field[" << i << "] = Min;\n";

        __newList << "Field[" << i << "].FieldsList = {";

        size_t j=1;
        for (; j< posfiles.size(); j++)
            __newList << j << ", ";
        __newList << j <<"};\n";

        if (posfiles.size() > 1)
            __newBgm << "Background Field = " << i << ";\n";
        else
            __newBgm << "Background Field = " << i-1 << ";\n";

        std::string newAlgo = __newAlgo.str();
        std::string newList = __newList.str();
        std::string newBgm = __newBgm.str();

        std::list<std::string> newEndExpr;
        if (posfiles.size() > 1)
            {
                newEndExpr.push_back(newAlgo);
                newEndExpr.push_back(newList);
            }

        newEndExpr.push_back(newBgm);

        std::list<std::string>::iterator itNew = newEndExpr.begin();

        //// Update place of matches (since fields have eventually been added)
        //// If fields have changed => algo, field list and background have to be replaced
        isInserted = false;
        for (std::list<boost::regex>::iterator itExpr = endExpr.begin(); itExpr != endExpr.end(); itExpr++)
            {
                boost::match_results<std::string::iterator> itPlace;
                bool exprIsFound = boost::regex_search(geofileString.begin(), geofileString.end(), itPlace, *itExpr);
                if (exprIsFound)
                    geofileString = boost::regex_replace(geofileString, *itExpr, *itNew);

                else
                    {
                        std::list<boost::regex>::iterator itNextExpr = itExpr;
                        itNextExpr++;
                        for ( ; itNextExpr != endExpr.end(); itNextExpr++)
                            {
                                bool nextIsFound = boost::regex_search(geofileString.begin(), geofileString.end(), itPlace, *itNextExpr);
                                if (nextIsFound)
                                    {
                                        geofileString.insert(itPlace[0].first, (*itNew).begin(), (*itNew).end() );
                                        isInserted = true;
                                        break;
                                    }
                            }
                        if (!isInserted)
                            geofileString.insert(geofileString.end(), (*itNew).begin(), (*itNew).end() );
                    }
                itNew++;
                if( itNew == newEndExpr.end() )
                    break;
            }

        ///// Geofile_s is complete => transform it into new file
        std::ofstream newGeofile;

        // Find name of the geofile (without extension), from the complete path
        size_t geoNameBegin = geofile.find_last_of("/");
        size_t extensionBegin = geofile.find(".geo") - 1;
        std::string geofileName = geofile.substr( geoNameBegin+1, extensionBegin - geoNameBegin);

        if( boost::filesystem::exists((boost::format( "./new-%1%.geo" ) % geofileName ).str()) )
            remove((boost::format( "./new-%1%.geo" ) % geofileName ).str().c_str() );

        newGeofile.open( (boost::format( "./new-%1%.geo" ) % geofileName ).str().c_str() );
        newGeofile << geofileString;
        newGeofile.close();

        std::string newAccessGeofile = (boost::format( "./new-%1%.geo" ) % geofileName ).str();

        return newAccessGeofile;
    }

    template<int Dim,
             int Order,
             int OrderGeo>
    std::string
    MeshAdaptation<Dim,
                   Order,
                   OrderGeo>::buildAdaptedMesh(std::string geofile, std::string name, std::vector<std::string> posfiles,
                                          bool aniso)
    {
        // Find name of the geofile (without extension), from the complete path
        size_t geoNameBegin = geofile.find_last_of("/");
        size_t extensionBegin = geofile.find(".geo") - 1;
        std::string geofileName = geofile.substr( geoNameBegin+1, extensionBegin - geoNameBegin);

        std::string prefix = (boost::format( "./%1%" ) % geofileName ).str();
        std::string mshFormat = "msh";
        std::string newMeshName =  (boost::format( "%1%%2%.%3%" ) % prefix % name % mshFormat ).str();

        int nbPosfiles = posfiles.size();

        if ( !mpi::environment::initialized() || ( mpi::environment::initialized()  &&  Environment::worldComm().globalRank() == Environment::worldComm().masterRank() ) )
            {
                /// *********** GMSH call to build new mesh ********** ///
#if FEELPP_HAS_GMSH
#if defined(FEELPP_HAS_GMSH_H)
                //#if 0

                std::string geofileNameWE = (boost::format( "%1%.geo" ) % geofileName).str();
                std::string exeName = "";

                // Load geofile (.geo) and post processing (.pos) files
                int argcGmsh = nbPosfiles + 2;
                char **argvGmsh = new char * [argcGmsh]; // geofile + posfiles

                /// argv[0] is the executable name
                char *argExe = new char[exeName.size() + 1];
                copy(exeName.begin(), exeName.end(), argExe);
                argExe[exeName.size()] = '\0';
                argvGmsh[0] = argExe;


                /// argv[1] => geofile
                char *argGeo = new char[geofileNameWE.size() + 1];
                copy(geofileNameWE.begin(), geofileNameWE.end(), argGeo);
                argGeo[geofileNameWE.size()] = '\0';
                argvGmsh[1] = argGeo;

                /// argv[2,...,n] => posfiles
                for (int i=0; i<nbPosfiles; i++)
                    {
                        int j = i + 2; // Shift (0 = exe, 1 = geo)
                        char *argPos = new char[posfiles[i].size() + 1];
                        copy(posfiles[i].begin(), posfiles[i].end(), argPos);
                        argPos[posfiles[i].size()] = '\0';
                        argvGmsh[j] = argPos;
                    }

                // Initializing
                GmshInitialize(argcGmsh, argvGmsh);

                GmshSetOption("Mesh", "Algorithm", 5.);
                ::GModel *newGmshModel = new GModel();

                // Read the input geometry file
                newGmshModel->readGEO(geofileNameWE);

                // Retreive list of files
                for (unsigned int i = 0; i < CTX::instance()->files.size(); i++)
                    {
                        LOG(INFO) << "[MeshAdaptation] loaded files : " << CTX::instance()->files[i] << "\n";
                        MergeFile(CTX::instance()->files[i]);
                    }

                // /* Create Fields from PView list */
                ::FieldManager* myFieldManager = newGmshModel->getFields();
                std::list<int> idList;

                for (unsigned int i = 0; i < ::PView::list.size(); i++)
                    {
                        ::PView *v = ::PView::list[i];
                        if (v->getData()->hasModel(::GModel::current()))
                            {
                                Msg::Error("Cannot use view based on current mesh for background mesh: you might"
                                           " want to use a list-based view (.pos file) instead");
                                return 0;
                            }

                        LOG(INFO) << "[MeshAdaptation] PostView : " << v->getData()->getFileName() << "\n";

                        /// Add new PView as post processing file
                        int id = myFieldManager->newId();
                        myFieldManager->newField(id, "PostView");

                        ::Field *f = myFieldManager->get(id);
                        f->options["IView"]->numericalValue(i);
                        idList.push_back(id);
                    }

                // /* Create minAniso field, with all the posfiles (intersection) */
                int id = myFieldManager->newId();

                if (aniso)
                    myFieldManager->newField(id, "MinAniso");
                else
                    myFieldManager->newField(id, "Min");
                ::Field *f = myFieldManager->get(id);

                /// Check if list of field for MinAniso is empty
                assert( f->options["FieldsList"]->list().size() == 0);
                /// Copy idlist vector into algorithm fieldlist
                //std::copy(idList.begin(), idList.end(), std::back_inserter(f->options["FieldsList"]->list()) );
                f->options["FieldsList"]->list(idList);

                // /* Now create the adapted mesh */

                /// Define background field from Fields
                myFieldManager->setBackgroundFieldId(id);
                f = myFieldManager->get(myFieldManager->getBackgroundField());

                /// Algo for remeshing
                CTX::instance()->mesh.algo2d = ALGO_2D_BAMG;
                if (Dim == 3)
                    CTX::instance()->mesh.algo3d = ALGO_3D_MMG3D;

                newGmshModel->deleteMesh(); //Delete current mesh
                newGmshModel->mesh(Dim);

                LOG(INFO) << "[MeshAdaptation] New mesh built : " << newMeshName << "\n";
                newGmshModel->writeMSH(newMeshName);
                LOG(INFO) << "[MeshAdaptation] vertices : " << newGmshModel->getNumMeshVertices() << "\n";
                LOG(INFO) << "[MeshAdaptation] elements : " << newGmshModel->getNumMeshElements() << "\n";

                //// Delete list of loaded files
                //CTX::instance()->files.erase(CTX::instance()->files.begin(), CTX::instance()->files.end());

                //// Clear PView list to update with new posfile
                ::PView::list.erase(::PView::list.begin(), ::PView::list.end() );

                // Cleanup memory
                newGmshModel->destroy();
                delete newGmshModel;
                delete[] argGeo;
                delete[] argExe;
                delete[] argvGmsh;
                GmshFinalize();

#else
                std::string newAccessGeofile = createAdaptedGeo(geofile, name, posfiles, aniso);

                // Execute gmsh command to generate new mesh from new geofile built
                std::ostringstream __str;
                __str << "gmsh" << " " << newAccessGeofile << " "
                      << "-o " << newMeshName << " ";
                for (int p=0; p<nbPosfiles; p++)
                    {
                        __str << posfiles[p] << " ";
                    }
                if (aniso)
                    {
                        __str << "-algo bamg" << " ";
                        if( Dim == 3)
                            __str << "-algo mmg3d" << " ";
                    }
                __str << "-" << Dim;

                ::system(__str.str().c_str());
#endif
#else
                throw std::invalid_argument("Gmsh is not available on this system");
#endif
            }

        if ( mpi::environment::initialized() )
            {
                Environment::worldComm().barrier();
                LOG(INFO) << "Broadcast adapted mesh file " << newMeshName << " to all other mpi processes\n";
                mpi::broadcast( Environment::worldComm().globalComm(), newMeshName, 0 );
            }

        return newMeshName;
    }

    template<int Dim,
             int Order,
             int OrderGeo>
    void
    MeshAdaptation<Dim,
                   Order,
                   OrderGeo>::computeMetric(const double tol, const double hMin, const double hMax,
                                            const matrixN_type & hessianMatrix, matrixN_type & metric, double & maxEigenvalue, int dofId)
    {
        using namespace Feel::vf;

        Eigen::EigenSolver< matrixN_type > eigenSolver;
        eigenSolver.compute(hessianMatrix);
        vectorN_type eigenvalues = eigenSolver.eigenvalues().array().abs();

        /***** Use of Frey formula - bamg user manual
               eigenvalues *= 1./(tol*norm);
        ******/

        /***** Use of Alauzet formula - p35 Metric-Based Anisotropic Mesh Adaptation ****/
        //eigenvalues *= 1./(tol*norm);
        eigenvalues *= 1./tol;
        eigenvalues *= 0.5*(Dim/(double) (Dim+1))*(Dim/(double) (Dim+1));
        for (int i=0; i<Dim; i++)
            {
                eigenvalues(i)=std::min(std::max(eigenvalues(i),1/(hMax*hMax)), 1/(hMin*hMin));
            }
        /******/

        matrixN_type S = matrixN_type::Zero();
        for (int i=0; i<Dim; i++)
            S(i,i)=eigenvalues(i); //1/math::sqrt(eigenvalues(i));

        matrixN_type R = matrixN_type::Zero();
        for (int j=0; j<Dim; j++)
            for (int i=0; i<Dim; i++)
                R(i,j)=real((eigenSolver.eigenvectors())(i,j)); // because eigenvectors are complex

        metric = (R * S) * R.transpose();
        maxEigenvalue = eigenvalues.maxCoeff();

        for(int i=0; i<Dim; i++)
            M_measures[dofId](i) = 1.0/math::sqrt(eigenvalues(i));

        M_directions[dofId] = R;
    }

    template<int Dim,
             int Order,
             int OrderGeo>
    std::string
    MeshAdaptation<Dim,
                   Order,
                   OrderGeo>::adaptMeshHess1(element_type& var, const mesh_ptrtype& mesh, double hMin, double hMax, std::string name,
                                             std::string geofile, double tol, bool aniso)
    {
        using namespace Feel::vf;

        p0_space_ptrtype P0h = p0_space_type::New( mesh ); //P0 space
        p1_space_ptrtype P1h = p1_space_type::New( mesh ); //P1 space

        // Initialize M_measures and M_directions size
        M_measures.resize(P1h->nLocalDof());
        M_directions.resize(P1h->nLocalDof());

        // Store measure on each point
        int bbItemSize;
        if (aniso)
            bbItemSize = Dim*Dim;
        else
            bbItemSize = 1;

        //cout << "n_bb_item=" << bbItemSize << endl;
        std::vector<p1_element_type> bbNewMap(bbItemSize, P1h->element()); // map of adapted measures

        // // Store max and min of solution U
        // double max_val = U.max();
        // double min_val = U.min();
        // double cutoff = 1.e-5;
        // double norm = std::max(std::max(fabs(max_val),fabs(min_val)), cutoff);

        // ******* From U in P1 space -> approximate components of hessian matrix on P1 space ******** //
        // Project u on P1 space
        auto varP1 = P1h->element();
        varP1 = vf::project( P1h, elements(mesh), idv(var) );

        // First derivatives of U (P1) are P0
        std::vector<p0_element_type> dvard1P0(Dim, P0h->element());
        std::vector<p1_element_type> dvard1P1(Dim, P1h->element());
        for (int i=0; i<Dim; i++)
            {
                dvard1P0[i] = vf::project(P0h, elements(mesh), gradv(varP1)(0,i) );
                dvard1P1[i] = div( sum(P1h, idv(dvard1P0[i])*meas()), sum(P1h, meas() ) ); //L2 proj -> P1
            }

        // Second derivatives of proj_P1(U) are P0
        std::vector<p0_element_type> hessP0(Dim*Dim, P0h->element());
        std::vector<p1_element_type> hessP1(Dim*Dim, P1h->element());
        for (int i=0; i<Dim; i++)
            {
                for (int j=0; j<Dim; j++)
                    {
                        // Hessian components are entered column by column
                        hessP0[i+j*Dim] = vf::project(P0h, elements(mesh), gradv( dvard1P1[j] )(0,i) );
                        hessP1[i+j*Dim] = div( sum(P1h, idv(hessP0[i+j*Dim])*meas()), sum(P1h, meas() ) ); //L2 proj -> P1
                    }
            }
        // ******************************************************************************************* //

        auto dofptItP1 = P1h->dof()->dofPointBegin();
        auto dofptEnP1 = P1h->dof()->dofPointEnd();

        // double hsizeMin = (1.0/1000)*meshSize;
        // std::cout << "hsize min = " << hsizeMin << std::endl;
        // double hsizeMax = 1000*meshSize;
        // std::cout << "hsize max = " << hsizeMax << std::endl;

        for ( ; dofptItP1 != dofptEnP1; dofptItP1++)
            {
                auto dofptCoordP1 = boost::get<0>(*dofptItP1); //dofptItP1->get<0>();
                auto dofptIdP1 = boost::get<1>(*dofptItP1); //dofptItP1->get<1>();

                matrixN_type hessianMatrix;
                Eigen::EigenSolver< matrixN_type > eigenSolver;

                // Associate each dof with the P1 hessian matrix projection
                for (int i=0; i<Dim; i++)
                    {
                        for (int j=0; j<Dim; j++)
                            {
                                hessianMatrix(i,j) = (hessP1[i+j*Dim])[dofptIdP1];
                            }
                    }

                double maxEigenvalue;
                matrixN_type metrics;
                computeMetric(tol, hMin, hMax, hessianMatrix, metrics, maxEigenvalue, dofptIdP1);

                if (aniso)
                    {
                        for (int j=0; j<Dim; j++)
                            {
                                for (int i=0; i<Dim; i++)
                                    {
                                        (bbNewMap[i+j*Dim])[dofptIdP1] = metrics(i,j);
                                    }
                            }
                    }
                else
                    {
                        auto newHsize = 1.0/math::sqrt( maxEigenvalue);
                        (bbNewMap[0])[dofptIdP1] = newHsize;
                    }
            }

        /// Build posfile
        std::string posfileName;
        if (aniso)
            posfileName = createPosfileAnisotropic(name, bbNewMap, mesh);
        else
            posfileName = createPosfile(name, bbNewMap[0], mesh);

        return posfileName;

    }

    template<int Dim,
             int Order,
             int OrderGeo>
    std::string
    MeshAdaptation<Dim,
                   Order,
                   OrderGeo>::adaptMeshHess2(element_type& var, const mesh_ptrtype& mesh, double hMin, double hMax,
                                             std::string name, std::string geofile, double tol, bool aniso)
    {
        return adaptMeshHess2(var, mesh, hMin, hMax, name, geofile, tol, aniso, Feel::mpl::bool_< isP1 >() );
    }

    template<int Dim,
             int Order,
             int OrderGeo>
    std::string
    MeshAdaptation<Dim,
                   Order,
                   OrderGeo>::adaptMeshHess2(element_type& var, const mesh_ptrtype& mesh, double hMin, double hMax,
                                             std::string name, std::string geofile, double tol, bool aniso,
                                             mpl::bool_<true>)
    {
        return adaptMeshHess1(var, mesh, hMin, hMax, name, geofile, tol, aniso);
    }

    template<int Dim,
             int Order,
             int OrderGeo>
    std::string
    MeshAdaptation<Dim,
                   Order,
                   OrderGeo>::adaptMeshHess2(element_type& var, const mesh_ptrtype& mesh, double hMin, double hMax,
                                             std::string name, std::string geofile, double tol, bool aniso,
                                             mpl::bool_<false>)
    {
        p1_space_ptrtype P1h = p1_space_type::New( mesh ); //P1 space
        p1vec_space_ptrtype P1hvec = p1vec_space_type::New( mesh ); //P1 (vectorial) space

        // Initialize M_measures and M_directions size
        M_measures.resize(P1h->nLocalDof());
        M_directions.resize(P1h->nLocalDof());


        // Define P(k-2) space from Order parameter
        typedef bases<Lagrange<Order-2,Scalar, Discontinuous> > p_km2_basis_type;
        typedef FunctionSpace<mesh_type, p_km2_basis_type> p_km2_space_type;
        typedef boost::shared_ptr<p_km2_space_type> p_km2_space_ptrtype;
        typedef typename p_km2_space_type::element_type p_km2_element_type;

        p_km2_space_ptrtype Pkm2H = p_km2_space_type::New( mesh );

        // Store measure on each point
        int bbItemSize; // = (aniso) ? Dim*Dim : 1;
        if (aniso)
            bbItemSize = Dim*Dim;
        else
            bbItemSize = 1;

        //cout << "n_bb_item=" << bbItemSize << endl;
        std::vector<p1_element_type> bbNewMap(bbItemSize, P1h->element()); // map of adapted measures

        // Store max and min of solution U
        // double max_val = U.max();
        // double min_val = U.min();
        // double cutoff = 1.e-5;
        // double norm = std::max(std::max(fabs(max_val),fabs(min_val)), cutoff);

        // Second derivatives of U are P(k-2)
        std::vector<p_km2_element_type> dvard2Pkm2(Dim*Dim, Pkm2H->element() );
        for (int i=0; i<Dim; i++)
            {
                for (int j=0; j<Dim; j++)
                    {
                        dvard2Pkm2[i+j*Dim] = vf::project( Pkm2H, elements(mesh),hessv(var)(i,j));
                    }
            }

        p1_element_type V(P1h, "V");

        /// Proj -> P1 for each hessian matrix components :
        /// Find \int proj_P1(U)V = \int UV \forall V
        std::vector<p1_element_type> dvard2P1(Dim*Dim, P1h->element());
        //auto l2proj = opProjection( _domainSpace=P1h, _imageSpace=P1h );

        for(int i=0; i<Dim; i++)
            {
                for (int j=0; j<Dim; j++)
                    {
                        //dvard2P1[i+j*Dim] = l2proj->project( _expr=idv( dvard2Pkm2[i+j*Dim]) );
                        dvard2P1[i+j*Dim] = vf::project( P1h, elements(mesh), idv( dvard2Pkm2[i+j*Dim]) );
                    }
            }
        // ******************************************************************************************* //

        auto dofptItP1 = P1h->dof()->dofPointBegin();
        auto dofptEnP1 = P1h->dof()->dofPointEnd();

        // double hsizeMin = (1.0/10000)*meshSize;
        // std::cout << "hsize min = " << hsizeMin << std::endl;
        // double hsizeMax = 10000*meshSize;
        // std::cout << "hsize max = " << hsizeMax << std::endl;

        std::vector<p1_element_type> hessianComponents(4);

        for ( ; dofptItP1 != dofptEnP1; dofptItP1++)
            {
                auto dofptCoordP1 = boost::get<0>(*dofptItP1); //dofptItP1->get<0>();
                auto dofptIdP1 = boost::get<1>(*dofptItP1); //dofptItP1->get<1>();

                matrixN_type hessianMatrix;

                // Associate each dof with the P1 hessian matrix projection
                for (int i=0; i<Dim; i++)
                    {
                        for (int j=0; j<Dim; j++)
                            {
                                hessianMatrix(i,j) = (dvard2P1[i+j*Dim])[dofptIdP1];
                            }
                    }

                double maxEigenvalue;
                matrixN_type metrics;
                computeMetric(tol, hMin, hMax, hessianMatrix, metrics, maxEigenvalue, dofptIdP1);

                if (aniso)
                    {
                        for (int j=0; j<Dim; j++)
                            {
                                for (int i=0; i<Dim; i++)
                                    {
                                        (bbNewMap[i+j*Dim])[dofptIdP1] = metrics(i,j);
                                    }
                            }
                    }
                else
                    {
                        auto newHsize = 1.0/math::sqrt( maxEigenvalue);
                        (bbNewMap[0])[dofptIdP1] = newHsize;

                    }
            }

        /// Build posfile
        std::string posfileName;
        if (aniso)
            posfileName = createPosfileAnisotropic(name, bbNewMap, mesh);
        else
            posfileName = createPosfile(name, bbNewMap[0], mesh);

        return posfileName;

    }

}
#endif // __MESHADAPTATION_HPP
