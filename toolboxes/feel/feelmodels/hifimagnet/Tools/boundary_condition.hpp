/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \file thermoelectric.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2011-16-12
 */

#ifndef __BOUNDARY_CONDITION_HPP
#define __BOUNDARY_CONDITION_HPP 1

#include <boost/spirit/include/qi.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/spirit/include/qi_eol.hpp>
#include <boost/spirit/include/qi_rule.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

#include <applications/Tools/parsedfunction_ginac.hpp>

namespace phx = boost::phoenix;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

template<int Dim>
class condition
{
public :

#if defined(CRB_MODEL_linear) || defined(CRB_MODEL_nonlinear)
    //CRB : conditions need to be evaluate in mu (cst_ref())
    typedef vf::Cst<boost::reference_wrapper<double> > cst_double;
    typedef parsedfunc_ginac<Dim, cst_double> parsedfunc;
#else
    typedef parsedfunc_ginac<Dim> parsedfunc;
#endif
    typedef boost::shared_ptr<parsedfunc> parsedfunc_ptr;

    condition()
        :
        elec_cond_(parsedfunc_ptr( new parsedfunc("0.0")) ),
        thermic_cond_(parsedfunc_ptr( new parsedfunc("0.0")) ),
        magnetic_permeability_(parsedfunc_ptr( new parsedfunc("0.0")) ),
        coeff_var_(parsedfunc_ptr( new parsedfunc("0.0")) ),
        constant_(parsedfunc_ptr( new parsedfunc("0.0")) ),
        vectorial_constant_( std::vector<parsedfunc_ptr>(Dim, parsedfunc_ptr( new parsedfunc("0.0"))) ),
        is_vectorial_( false )
    {
    }

    // Access to members
    std::string name() const {return name_;}
    std::string variable() const {return variable_;}
    std::string physical_name() const {return physical_name_;}
    std::string type() const {return type_;}
    parsedfunc_ptr elec_cond() const {return elec_cond_;}
    parsedfunc_ptr thermic_cond() const {return thermic_cond_;}
    double lorentz() const {return lorentz_;}
    double alpha() const {return alpha_;}
    double young_mod() const {return young_mod_;}
    double poisson_coeff() const {return poisson_coeff_;}
    double dilatation_coeff() const {return dilatation_coeff_;}
    parsedfunc_ptr coeff_var() const {return coeff_var_;}
    parsedfunc_ptr constant() const {return constant_;}
    parsedfunc_ptr magnetic_permeability() const {return magnetic_permeability_;}
    //double magnetic_permeability() const {return magnetic_permeability_;}
    std::vector<parsedfunc_ptr> vectorial_constant() const {return vectorial_constant_;}
    boost::optional<double> intensity() const {return intensity_;}
    boost::optional<std::string> material_name() const {return mat_name_;}
    //std::list<double>& params() {return params_;}
    std::map<std::string, double>& params() {return params_;}
    bool is_vectorial() {return is_vectorial_;}

    // Change members value
    void name(std::string name) {name_ = name;}
    void setName(std::string name) {name_ = name;}
    void variable(std::string var) {variable_ = var;}
    void setVariable(std::string var) {variable_ = var;}
    void physical_name(std::string physical_name) {physical_name_ = physical_name;}
    void setPhysical_name(std::string physical_name) {physical_name_ = physical_name;}
    void type(std::string type) {type_ = type;}
    void setType(std::string type) {type_ = type;}
    void elec_cond(parsedfunc_ptr elec_cond) {elec_cond_ = elec_cond;}
    void setElec_cond(parsedfunc_ptr elec_cond) {elec_cond_ = elec_cond;}
    void thermic_cond(parsedfunc_ptr thermic_cond) {thermic_cond_ = thermic_cond;}
    void setThermic_cond(parsedfunc_ptr thermic_cond) {thermic_cond_ = thermic_cond;}
    void magnetic_permeability(parsedfunc_ptr magnetic_permeability) {magnetic_permeability_ = magnetic_permeability;}
    void setMagnetic_permeability(parsedfunc_ptr magnetic_permeability) {magnetic_permeability_ = magnetic_permeability;}
    //void magnetic_permeability(double magnetic_permeability) {magnetic_permeability_ = magnetic_permeability;}
    void lorentz(double L){lorentz_ = L;}
    void setLorentz(double L){lorentz_ = L;}
    void alpha( double alpha){alpha_ = alpha;}
    void setAlpha( double alpha){alpha_ = alpha;}
    void young_mod(double ym) {young_mod_ = ym;}
    void setYoung_mod(double ym) {young_mod_ = ym;}
    void poisson_coeff(double p) {poisson_coeff_ = p;}
    void setPoisson_coeff(double p) {poisson_coeff_ = p;}
    void dilatation_coeff(double d) {dilatation_coeff_ = d;}
    void setDilatation_coeff(double d) {dilatation_coeff_ = d;}
    void coeff_var(parsedfunc_ptr cv) {coeff_var_ = cv;}
    void setCoeff_var(parsedfunc_ptr cv) {coeff_var_ = cv;}
    void constant(parsedfunc_ptr cs) {constant_ = cs;}
    void setConstant(parsedfunc_ptr cs) {constant_ = cs;}
    void vectorial_constant(int idx, parsedfunc_ptr f) { if(idx < Dim){vectorial_constant_[idx] = f;} }
    void setVectorial_constant(int idx, parsedfunc_ptr f) { if(idx < Dim){vectorial_constant_[idx] = f;} }
    void intensity(double i) {intensity_ = i;}
    void setIntensity(double i) {intensity_ = i;}
    void material_name(std::string mat_name) {mat_name_ = mat_name;}
    void setMaterial_name(std::string mat_name) {mat_name_ = mat_name;}

    void elec_cond_str(std::string elec_cond, std::vector<std::string> params)
    {
        elec_cond_ = parsedfunc_ptr( new parsedfunc(elec_cond, (boost::format("ElecCondBC%1%") % idx_ ).str(), params) );
    }
    void thermic_cond_str(std::string thermic_cond, std::vector<std::string> params)
    {
        thermic_cond_ = parsedfunc_ptr( new parsedfunc(thermic_cond, (boost::format("ThermCondBC%1%") % idx_ ).str(), params) );
    }
    void magnetic_permeability_str(std::string permeability, std::vector<std::string> params)
    {
        magnetic_permeability_ = parsedfunc_ptr(new parsedfunc(permeability, (boost::format("MagneticPermeability%1%") % idx_ ).str(), params ) );
    }
    void coeff_var_str(std::string cv, std::vector<std::string> params)
    {
        coeff_var_ = parsedfunc_ptr( new parsedfunc(cv, (boost::format("CoeffVar%1%") % idx_ ).str(), params) );
    }
    void constant_str(std::string cs, std::vector<std::string> params)
    {
        constant_ = parsedfunc_ptr( new parsedfunc(cs, (boost::format("Cst%1%") % idx_ ).str(), params) );
    }
    void vectorial_constant_str(int idx, std::string f, std::vector<std::string> params)
    {
        is_vectorial_=true;
        if(idx < Dim)
            vectorial_constant_[idx] = parsedfunc_ptr( new parsedfunc(f,(boost::format("vCst%1%comp%2%") % idx_ % idx ).str(), params) );
    }

    int idx_;
    std::string name_;
    std::string variable_;
    std::string physical_name_;
    std::string type_;
    parsedfunc_ptr elec_cond_;
    parsedfunc_ptr thermic_cond_;
    parsedfunc_ptr magnetic_permeability_;
    //double magnetic_permeability_;
    double lorentz_;
    double alpha_;
    double young_mod_;
    double poisson_coeff_;
    double dilatation_coeff_;

    parsedfunc_ptr coeff_var_;
    parsedfunc_ptr constant_;
    std::vector<parsedfunc_ptr> vectorial_constant_;
    bool is_vectorial_;
    std::map<std::string, double> params_;

    boost::optional<double> intensity_;
    boost::optional<std::string> mat_name_;
};

template <typename Iterator, int Dim>
bool parse_condition(Iterator first, Iterator last, condition<Dim>& bc, std::vector<std::string> params = std::vector<std::string>() )
    {
        using qi::int_;
        using qi::lit;
        using qi::double_;
        using qi::lexeme;
        using qi::eol;
        using ascii::char_;
        using ascii::space;
        using qi::_1;
        using qi::phrase_parse;
        using boost::optional;
        using qi::rule;
        using qi::locals;

        qi::rule<Iterator, std::string(), ascii::space_type> quoted_string;
        quoted_string %= lexeme['"' >> +(char_ - '"') >> '"'];

        qi::rule<Iterator, void(), ascii::space_type> name_rule;
        name_rule = quoted_string[phx::bind(&condition<Dim>::setName, &bc, _1)] >> *(char_(']')) >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> variable_rule;
        variable_rule = char_('_')
            >> (lit("variable") | lit("var") )
            >> char_('=') >> quoted_string[phx::bind(&condition<Dim>::setVariable, &bc, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> physical_name_rule;
        physical_name_rule = char_('_')
            >> (lit("physical_name") | lit("physical") )
            >> char_('=') >> quoted_string[phx::bind(&condition<Dim>::setPhysical_name, &bc, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> type_rule;
        type_rule = char_('_')
            >> (lit("condition_type") | lit("type") )
            >> char_('=') >> quoted_string[phx::bind(&condition<Dim>::setType, &bc, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> elec_cond_rule;
        elec_cond_rule = char_('_')
            >> (lit("electrical_conductivity") | lit("elec_cond") )
            >> char_('=') >> quoted_string[phx::bind(&condition<Dim>::elec_cond_str, &bc, _1, params)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> thermal_cond_rule;
        thermal_cond_rule = char_('_')
            >> (lit("thermal_conductivity") | lit("thermal_cond") )
            >> char_('=') >> quoted_string[phx::bind(&condition<Dim>::thermic_cond_str, &bc, _1, params)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> magnetic_permeability_rule;
        magnetic_permeability_rule = char_('_')
            >> (lit("magnetic_permeability")  | lit("mu") )
            >> char_('=') >> quoted_string[phx::bind(&condition<Dim>::magnetic_permeability_str, &bc, _1, params)] >> *(eol);
        // magnetic_permeability_rule = char_('_')
        //     >> (lit("magnetic_permeability")  | lit("mu") )
        //     >> char_('=') >> double_[phx::bind(&condition<Dim>::magnetic_permeability, &bc, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> lorentz_rule;
        lorentz_rule = char_('_')
            >> (lit("lorentz") | lit("L") )
            >> char_('=') >> double_[phx::bind(&condition<Dim>::setLorentz, &bc, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> alpha_rule;
        alpha_rule = char_('_')
            >> lit("alpha")
            >> char_('=') >> double_[phx::bind(&condition<Dim>::setAlpha, &bc, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> young_rule;
        young_rule = char_('_')
            >> (lit("young") | lit("E") | lit("young_modulus"))
            >> char_('=') >> double_[phx::bind(&condition<Dim>::setYoung_mod, &bc, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> poisson_rule;
        poisson_rule = char_('_')
            >> (lit("poisson") | lit("nu") )
            >> char_('=') >> double_[phx::bind(&condition<Dim>::setPoisson_coeff, &bc, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> dilatation_rule;
        dilatation_rule = char_('_')
            >> (lit("dilatation") | lit("alpha") )
            >> char_('=') >> double_[phx::bind(&condition<Dim>::setDilatation_coeff, &bc, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> coeff_var_rule;
        coeff_var_rule = char_('_')
            >> (lit("robin_variable_coeff") | lit("variable_coeff") | lit("var_coeff") )
            >> char_('=') >> quoted_string[phx::bind(&condition<Dim>::coeff_var_str, &bc, _1, params)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> constant_rule;
        constant_rule = char_('_')
            >> (lit("constant") | lit("cst") )
            >> char_('=') >> quoted_string[phx::bind(&condition<Dim>::constant_str, &bc, _1, params)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> vectorial_constant_rule;
        vectorial_constant_rule = char_('_')
            >> (lit("constant") | lit("cst") )
            >> char_('=') >> char_("(")
            >> quoted_string[phx::bind(&condition<Dim>::vectorial_constant_str, &bc, 0, _1, params)]
            >> char_(',')
            >> quoted_string[phx::bind(&condition<Dim>::vectorial_constant_str, &bc, 1, _1, params)]
            >> -( char_(',') )
            >> -( quoted_string[phx::bind(&condition<Dim>::vectorial_constant_str, &bc, 2, _1, params)] )
            >> char_(")")
            >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> intensity_rule;
        intensity_rule = char_('_')
            >> (lit("intensity") | lit("I") )
            >> char_('=') >> double_[phx::bind(&condition<Dim>::setIntensity, &bc, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> material_rule;
        material_rule = char_('_')
            >> (lit("material") | lit("material_name") )
            >> char_('=') >> quoted_string[phx::bind(&condition<Dim>::setMaterial_name, &bc, _1)] >> *(eol);

        bool r = phrase_parse(first, last,

                              //  Begin grammar
                              (
                               name_rule
                               >> variable_rule
                               >> type_rule
                               >> physical_name_rule
                               >> -(elec_cond_rule)
                               >> -(thermal_cond_rule)
                               >> -(lorentz_rule)
                               >> -(alpha_rule)
                               >> -(magnetic_permeability_rule)
                               >> -(young_rule)
                               >> -(poisson_rule)
                               >> -(dilatation_rule)
                               >> -(coeff_var_rule)
                               >> -(constant_rule | vectorial_constant_rule)
                               >> -(intensity_rule)
                               >> -(material_rule)
                               ),
                              //  End grammar

                              space);

        if (!r || first != last) // fail if we did not get a full match
            return false;

        return r;
    }

#endif //__BOUNDARY_CONDITION_HPP
