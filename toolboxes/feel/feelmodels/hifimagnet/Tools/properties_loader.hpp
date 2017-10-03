/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
  Author(s): Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
       Date: 2011-16-12

       Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)
       Copyright (C) CNRS

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
   \file properties_loader
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2013-12-03
 */

#ifndef __PROP_LOADER_HPP
#define __PROP_LOADER_HPP 1

#include <applications/Tools/material.hpp>
#include <applications/Tools/boundary_condition.hpp>
#include <boost/tokenizer.hpp>

template<int Dim>
class properties_loader
{
public :

    typedef boost::tokenizer<boost::char_separator<char> > c_tokenizer;

    // Material pointer type, element_type specified
    typedef Hmaterial<Dim> material_elt;
    typedef boost::shared_ptr<material_elt> material_ptrelt;
    typedef std::map<std::string, material_ptrelt> material_map_type;
    // Condition pointer type, element_type specified
    typedef condition<Dim> condition_elt;
    typedef boost::shared_ptr<condition_elt> condition_ptrelt;
    typedef std::map<std::string, std::vector<condition_ptrelt> > condition_map_type;

    material_map_type load_materials(std::string materials_filename_access, std::vector<std::string> params=std::vector<std::string>() );
#if defined(CRB_MODEL_linear) || defined(CRB_MODEL_nonlinear)
    condition_map_type load_conditions(std::string conditions_filename_access,
                                       std::vector<std::string> params=std::vector<std::string>());
#else
    condition_map_type load_conditions(std::string conditions_filename_access, material_map_type materials,
                                       std::vector<std::string> params=std::vector<std::string>());
#endif
    void write_json(material_map_type mat_map, condition_map_type bc_map, std::string variable, std::string filename="mymodel");
};

template<int Dim>
typename properties_loader<Dim>::material_map_type
properties_loader<Dim>::load_materials(std::string materials_filename_access, std::vector<std::string> params)
{
    material_map_type materials;
    using boost::spirit::ascii::space;
    typedef std::string::const_iterator iterator_type;

    int proc_rank = Environment::worldComm().globalRank();

    // Transform materials files content into string
    std::ifstream materials_file(materials_filename_access);
    if(!materials_file)
        {
            std::ostringstream error_msg;
            error_msg << "[properties_loader] load_materials : " << materials_filename_access << " no such file\n";
            throw std::logic_error( error_msg.str() );
        }


    std::string materials_file_content, tmp_content;
    while(materials_file)
        {
            std::getline(materials_file, tmp_content);
            if(tmp_content.compare(0,1,"#") != 0) // Ignore comments line (#)
                materials_file_content += tmp_content;
        }

    // Split materials files into vector of string (one string per material/bc)
    std::vector<std::string> list_materials;
    boost::char_separator<char> sep("[");
    c_tokenizer materials_tokens(materials_file_content, sep);
    for (c_tokenizer::iterator tok_iter_mat = materials_tokens.begin(); tok_iter_mat != materials_tokens.end(); ++tok_iter_mat)
        list_materials.push_back(*tok_iter_mat);

    int idx = 0;
    for(std::string _current_mat: list_materials)
        {
            material_ptrelt mat = material_ptrelt( new material_elt() );
            mat->idx_ = idx;

            std::string::const_iterator iter = _current_mat.begin();
            std::string::const_iterator end = _current_mat.end();

            bool parse_ok = parse_material(iter, end, *mat, params);
            if( parse_ok )
                {
                    Feel::cout << "Material " << mat->name() << " correctly parsed" << std::endl << std::flush;
                }
            else
                {
                    throw std::logic_error( "[properties_loader] Error loading Material " + mat->name() );
                }
            materials.insert(std::pair<std::string, material_ptrelt>(mat->name(), mat));
            idx++;
        };

    LOG(INFO) << "Materials properties collected \n";
    return materials;

}

#if defined(CRB_MODEL_linear) || defined(CRB_MODEL_nonlinear)

// Load boundary conditions without taking materials into account (for CRB models)
template<int Dim>
typename properties_loader<Dim>::condition_map_type
properties_loader<Dim>::load_conditions(std::string conditions_filename_access, std::vector<std::string> params)
{
    condition_map_type conditions;

    using boost::spirit::ascii::space;
    typedef std::string::const_iterator iterator_type;

    int proc_rank = Environment::worldComm().globalRank();

    // Transform bc files content into string
    std::ifstream conditions_file(conditions_filename_access);
    if(!conditions_file)
        {
            std::ostringstream error_msg;
            error_msg << "[properties_loader] load_conditions : " << conditions_filename_access << " no such file\n";
            throw std::logic_error( error_msg.str() );
        }

    std::string conditions_file_content, tmp_content;
    while(conditions_file)
        {
            std::getline(conditions_file, tmp_content);
            if(tmp_content.compare(0,1,"#") != 0) // Ignore comments line (#)
                conditions_file_content += tmp_content;
        }

    // Split bc files into vector of string (one string per bc)
    std::vector<std::string> list_conditions;
    boost::char_separator<char> sep("[");
    c_tokenizer conditions_tokens(conditions_file_content, sep);
    for (c_tokenizer::iterator tok_iter_bc = conditions_tokens.begin(); tok_iter_bc != conditions_tokens.end(); ++tok_iter_bc)
        list_conditions.push_back(*tok_iter_bc);

    int idx = 0;
    for(std::string _current_c: list_conditions)
        {
            condition_ptrelt cond = condition_ptrelt( new condition_elt() );
            cond->idx_ = idx;
            std::string::const_iterator iter = _current_c.begin();
            std::string::const_iterator end = _current_c.end();

            bool parse_ok = parse_condition(iter, end, *cond, params);
            if( parse_ok )
                {
                    Feel::cout << "Condition " << cond->name() << " correctly parsed" << std::endl << std::flush;
                }
            else
                {
                    std::ostringstream error_msg;
                    error_msg << "[properties_loader] Error loading Condition" << cond->name() << "\n";
                    throw std::logic_error( error_msg.str() );
                }

            conditions[cond->variable()].push_back(cond);
            idx++;
        };

    LOG(INFO) << "Boundary conditions collected \n";
    return conditions;
}

#else

// Load boundary conditions using materials
template<int Dim>
typename properties_loader<Dim>::condition_map_type
properties_loader<Dim>::load_conditions(std::string conditions_filename_access, material_map_type materials, std::vector<std::string> params)
{
    condition_map_type conditions;

    using boost::spirit::ascii::space;
    typedef std::string::const_iterator iterator_type;

    int proc_rank = Environment::worldComm().globalRank();

    // Transform bc files content into string
    std::ifstream conditions_file(conditions_filename_access);
    if(!conditions_file )
        Feel::cout << "Boundary conditions file cannot be read" << std::endl;

    std::string conditions_file_content, tmp_content;
    while(conditions_file)
        {
            std::getline(conditions_file, tmp_content);
            if(tmp_content.compare(0,1,"#") != 0) // Ignore comments line (#)
                conditions_file_content += tmp_content;
        }

    // Split bc files into vector of string (one string per bc)
    std::vector<std::string> list_conditions;
    boost::char_separator<char> sep("[");
    c_tokenizer conditions_tokens(conditions_file_content, sep);
    for (c_tokenizer::iterator tok_iter_bc = conditions_tokens.begin(); tok_iter_bc != conditions_tokens.end(); ++tok_iter_bc)
        list_conditions.push_back(*tok_iter_bc);

    int idx = 0;
    for(std::string _current_c: list_conditions)
        {
            condition_ptrelt cond = condition_ptrelt( new condition_elt() );
            cond->idx_ = idx;
            std::string::const_iterator iter = _current_c.begin();
            std::string::const_iterator end = _current_c.end();

            bool parse_ok = parse_condition(iter, end, *cond, params);
            if(parse_ok)
                Feel::cout << "Condition " << cond->name() << " correctly parsed" << std::endl;

            auto mat_depend = cond->material_name();
            if( mat_depend )
                {
                    if ( materials.find(*mat_depend) == materials.end() )
                        throw std::logic_error( "[" + *mat_depend + "] : no such material declared in Materials in " + cond->name() + " BC def" );

                    cond->elec_cond(materials[*mat_depend]->elec_cond());
                    cond->thermic_cond(materials[*mat_depend]->thermic_cond());
                    cond->lorentz(materials[*mat_depend]->lorentz());
                    cond->alpha(materials[*mat_depend]->alpha());
                    cond->magnetic_permeability(materials[*mat_depend]->magnetic_permeability());
                    cond->young_mod(materials[*mat_depend]->young_mod());
                    cond->poisson_coeff(materials[*mat_depend]->poisson_coeff());
                    cond->dilatation_coeff(materials[*mat_depend]->dilatation_coeff());
                    if( materials[*mat_depend]->intensity() )
                        cond->intensity(*materials[*mat_depend]->intensity());
                }

            conditions[cond->variable()].push_back(cond);
            idx++;
        };

    LOG(INFO) << "Boundary conditions collected \n";
    return conditions;
}
#endif

template<int Dim>
void
properties_loader<Dim>::write_json(material_map_type mat_map, condition_map_type bc_map, std::string variable, std::string filename)
{
    std::ofstream jsonfile;
    std::string name = filename + ".json";
    jsonfile.open(name, std::ios::trunc);

    jsonfile << "// -*- mode: javascript -*- \n{\n";

    jsonfile << "\t \"Name\": \"" << filename << "\",\n";
    jsonfile << "\t \"Model\": \"" << variable << "\",\n"; //TEMP

    jsonfile << "\t \"BoundaryConditions\": \n";
    jsonfile << "\t { \n";

    // // Loop on map (variables)
    // for( auto bcmap_it=bc_map.begin(); bcmap_it!=bc_map.end; bcmap_it++)
    //     {
    // jsonfile << "\"" << bcmap_it->first << "\":{ \n";

    //jsonfile << "\t \t \"" << variable << "\":{ \n";

    jsonfile << "\t \t \"u\":{ \n"; //TEMPORARY

    //for( condition_ptrelt cond : bcmap_it->second )
    std::string conditions_dir_str, conditions_neu_str;
    std::string conditions_phi_str;

    int count_dirichlet = 0;
    int count_neumann = 0;

    for( condition_ptrelt cond : bc_map[variable] )
        {
            Feel::cout << "[" << cond->physical_name() << "]: " << cond->type() << std::endl << std::flush;
            if (cond->type() == "Dirichlet")
                {
                    count_dirichlet +=1;
                    conditions_dir_str += "\n \t \t \t \t \"" + cond->physical_name() + "\": {";
#ifdef SADDLEPOINT
                    conditions_phi_str += "\n \t \t \t \t \"" + cond->physical_name() + "\": {";
#endif
                    if( cond->is_vectorial() )
                        {
                            conditions_dir_str += "\"expr\": \"{";
                            for(int i=0; i<Dim-1; i++)
                                conditions_dir_str += cond->vectorial_constant()[i]->func_str() + ",";
                            conditions_dir_str += cond->vectorial_constant()[Dim-1]->func_str() + "}:x:y:z\"},";
                        }
                    else
                        conditions_dir_str += "\"expr\": \"" + cond->constant()->func_str() + ":x:y:z\"},";
#ifdef SADDLEPOINT
                    conditions_phi_str += "\"expr\": \"0:x:y:z\" },";
#endif
                }

            if (cond->type() == "Neumann")
                {
                    count_neumann +=1;
                    conditions_neu_str += "\n \t \t \t \t \"" + cond->physical_name() + "\": { ";
#ifdef SADDLEPOINT
                    conditions_phi_str += "\n \t \t \t \t \"" + cond->physical_name() + "\": { ";
#endif

                    if( cond->is_vectorial() )
                        {
                            conditions_neu_str += "\"expr\": \"{";
                            for(int i=0; i<Dim-1; i++)
                                conditions_neu_str += cond->vectorial_constant()[i]->func_str() + ",";
                            conditions_neu_str += cond->vectorial_constant()[Dim-1]->func_str() + "}:x:y:z\"},";
                        }
                    else
                        conditions_neu_str += "\"expr\": \"" + cond->constant()->func_str() + ":x:y:z\"},";
#ifdef SADDLEPOINT
                    conditions_phi_str += "\"expr\": \"0:x:y:z\" },";
#endif
                }
        }
    //Remove comma
    conditions_dir_str = conditions_dir_str.substr(0, conditions_dir_str.length()-1);
    conditions_neu_str = conditions_neu_str.substr(0, conditions_neu_str.length()-1);
#ifdef SADDLEPOINT
    conditions_phi_str = conditions_phi_str.substr(0, conditions_phi_str.length()-1);
#endif

    //Dirichlet u
    if ( count_dirichlet )
        {
            jsonfile << "\t \t \t \"Dirichlet\":{";
            jsonfile << conditions_dir_str + "\n \t \t \t }";
        }
    //Neumann u
    if ( count_neumann )
        {
            jsonfile << ",\n";
            jsonfile << "\t \t \t \"Neumann\":{";
            jsonfile << conditions_neu_str + "\n \t \t \t }\n";
        }
    else
        {
            jsonfile << "\n";
        }
    jsonfile << "\t \t  }";

#ifdef SADDLEPOINT
    jsonfile << ",\n";
    jsonfile << "\t \t \"phi\":{ \n"; //TEMPORARY
    if ( count_dirichlet )
        {
    //Dirichlet phi
            jsonfile << "\t \t \t \"Dirichlet\":{";
            jsonfile << conditions_phi_str + "\n \t \t \t }\n";
        }
    jsonfile << "\t \t }\n";
#else
    jsonfile << "\n";
#endif

    //} //end for loop on variables
    jsonfile << "\t }\n";
    jsonfile << "}\n"; //end file

    jsonfile.close();
}

#endif //__PROP_LOADER_HPP
