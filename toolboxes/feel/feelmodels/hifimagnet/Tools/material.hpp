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
   \file material
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2013-12-03
 */

#ifndef __MATERIAL_HPP
#define __MATERIAL_HPP 1

#include <boost/spirit/include/qi.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/spirit/include/qi_eol.hpp>
#include <boost/spirit/include/qi_rule.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

#include <applications/Tools/parsedfunction_ginac.hpp>

namespace phx = boost::phoenix;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

template<int Dim>
class Hmaterial
{
public :

    typedef parsedfunc_ginac<Dim> parsedfunc;

    typedef boost::shared_ptr<parsedfunc> parsedfunc_ptr;

    Hmaterial()
        :
        magnetic_permeability_( parsedfunc_ptr( new parsedfunc("1")) ),
        saturation_induction_(0.0),
        js_( std::vector<parsedfunc_ptr>(Dim, parsedfunc_ptr( new parsedfunc("0.0"))) ),
        density_( parsedfunc_ptr( new parsedfunc("0.0")) ),
        specific_heat_( parsedfunc_ptr( new parsedfunc("0.0")) ),
        volumic_forces_( std::vector<parsedfunc_ptr>(Dim, parsedfunc_ptr( new parsedfunc("0.0"))) )
    {
    }

    // Acces to members
    std::string name() const {return name_;}
    std::string physical_name() const {return physical_name_;}
    parsedfunc_ptr elec_cond() const {return elec_cond_;}
    parsedfunc_ptr thermic_cond() const {return thermic_cond_;}
    double lorentz() const {return lorentz_;}
    double alpha() const {return alpha_;}
    double young_mod() const {return young_mod_;}
    double poisson_coeff() const {return poisson_coeff_;}
    double dilatation_coeff() const {return dilatation_coeff_;}
    boost::optional<double> ddp() const {return ddp_;}
    boost::optional<int> nbr_turns() const {return nbr_turns_;}
    std::vector<parsedfunc_ptr> js() const {return js_;}
    boost::optional<double> rhos() const {return rhos_;}
    boost::optional<double> intensity() const {return intensity_;}
    parsedfunc_ptr density() const {return density_;}
    parsedfunc_ptr specific_heat() const {return specific_heat_;}
    parsedfunc_ptr magnetic_permeability() const {return magnetic_permeability_;}
    double saturation_induction() const {return saturation_induction_;}
    double relative_permeability() const {return relative_permeability_;}
    std::vector<parsedfunc_ptr> volumic_forces() const {return volumic_forces_;}

    // Change members value
    void name(std::string name) {name_ = name;}
    void physical_name(std::string physical_name) {physical_name_ = physical_name;}
    void elec_cond(parsedfunc_ptr elec_cond) {elec_cond_ = elec_cond;}
    void thermic_cond(parsedfunc_ptr thermic_cond) {thermic_cond_ = thermic_cond;}
    void magnetic_permeability(parsedfunc_ptr mag_permeability) {magnetic_permeability_ = mag_permeability;}
    void saturation_induction(double sat) {saturation_induction_ = sat;}
    void relative_permeability(double rel_permeability) {relative_permeability_ = rel_permeability;}
    void lorentz(double L){lorentz_ = L;}
    void alpha(double alpha){alpha_ = alpha;}
    void ddp(double U){ddp_ = U;}
    void nbr_turns(int nt) {nbr_turns_ = nt;}
    void js(int idx, parsedfunc_ptr j) {if(idx < Dim){js_[idx] = j;} }
    void rhos(double rhos) {rhos_ = rhos;}
    void intensity(double i) {intensity_ = i;}
    void density(parsedfunc_ptr density) {density_ = density;}
    void specific_heat(parsedfunc_ptr specific_heat) {specific_heat_ = specific_heat;}
    void young_mod(double ym) {young_mod_ = ym;}
    void poisson_coeff(double p) {poisson_coeff_ = p;}
    void dilatation_coeff(double d) {dilatation_coeff_ = d;}
    void volumic_forces(int idx, parsedfunc_ptr f) { if(idx < Dim){volumic_forces_[idx] = f;} }

    void elec_cond_str(std::string elec_cond, std::vector<std::string> params)
    {
        elec_cond_ = parsedfunc_ptr( new parsedfunc(elec_cond, (boost::format("ElecCond%1%") % idx_ ).str(), params ) );
    }
    void thermic_cond_str(std::string thermic_cond, std::vector<std::string> params)
    {
        thermic_cond_ = parsedfunc_ptr( new parsedfunc(thermic_cond,(boost::format("ThermCond%1%") % idx_ ).str(), params ) );
    }
    void magnetic_permeability_str(std::string permeability, std::vector<std::string> params)
    {
        magnetic_permeability_ = parsedfunc_ptr(new parsedfunc(permeability, (boost::format("MagneticPermeability%1%") % idx_ ).str(), params ) );
    }
    void density_str(std::string density, std::vector<std::string> params)
    {
        density_ = parsedfunc_ptr( new parsedfunc(density, (boost::format("Density%1%") % idx_ ).str(), params ) );
    }
    void specific_heat_str(std::string specific_heat, std::vector<std::string> params)
    {
        specific_heat_ = parsedfunc_ptr( new parsedfunc(specific_heat, (boost::format("SpeHeat%1%") % idx_ ).str(), params ) );
    }
    void js_str(int idx, std::string f, std::vector<std::string> params)
    {
        if(idx < Dim)
            js_[idx] = parsedfunc_ptr( new parsedfunc(f, (boost::format("Js%1%comp%2%") % idx_ % idx ).str(), params ) );
    }
    void volumic_forces_str(int idx, std::string f, std::vector<std::string> params)
    {
        if(idx < Dim)
            volumic_forces_[idx] = parsedfunc_ptr( new parsedfunc(f, (boost::format("vForces%1%comp%2%") % idx_ % idx ).str(), params ) );
    }

    int idx_;
    std::string name_;
    std::string physical_name_;
    parsedfunc_ptr elec_cond_;
    parsedfunc_ptr thermic_cond_;
    parsedfunc_ptr magnetic_permeability_;
    double saturation_induction_;
    double relative_permeability_;

    double lorentz_;
    double alpha_;
    double young_mod_;
    double poisson_coeff_;
    double dilatation_coeff_;

    boost::optional<double> ddp_;
    boost::optional<int> nbr_turns_;
    std::vector<parsedfunc_ptr> js_;
    boost::optional<double> rhos_;
    boost::optional<double> intensity_;
    parsedfunc_ptr density_;
    parsedfunc_ptr specific_heat_;

    std::vector<parsedfunc_ptr> volumic_forces_;
};

template <typename Iterator, int Dim>
bool parse_material(Iterator first, Iterator last, Hmaterial<Dim>& mat, std::vector<std::string> params=std::vector<std::string>() )
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
        name_rule = quoted_string[phx::bind(&Hmaterial<Dim>::name, &mat, _1)] >> *(char_(']')) >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> physical_name_rule;
        physical_name_rule = char_('_')
            >> (lit("physical_name") | lit("physical") )
            >> char_('=') >> quoted_string[phx::bind(&Hmaterial<Dim>::physical_name, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> elec_cond_rule;
        elec_cond_rule = char_('_')
            >> (lit("electrical_conductivity") | lit("elec_cond") )
            >> char_('=') >> quoted_string[phx::bind(&Hmaterial<Dim>::elec_cond_str, &mat, _1, params)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> thermal_cond_rule;
        thermal_cond_rule = char_('_')
            >> (lit("thermal_conductivity") | lit("thermal_cond") )
            >> char_('=') >> quoted_string[phx::bind(&Hmaterial<Dim>::thermic_cond_str, &mat, _1, params)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> lorentz_rule;
        lorentz_rule = char_('_')
            >> (lit("lorentz") | lit("L") )
            >> char_('=') >> double_[phx::bind(&Hmaterial<Dim>::lorentz, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> alpha_rule;
        alpha_rule = char_('_')
            >> lit("alpha")
            >> char_('=') >> double_[phx::bind(&Hmaterial<Dim>::alpha, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> magnetic_permeability_rule;
        magnetic_permeability_rule = char_('_')
            >> (lit("magnetic_permeability")  | lit("mu") )
            >> char_('=') >> quoted_string[phx::bind(&Hmaterial<Dim>::magnetic_permeability_str, &mat, _1, params)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> saturation_induction_rule;
        saturation_induction_rule = char_('_')
            >> (lit("saturation_induction")  | lit("Bs") )
            >> char_('=') >> double_[phx::bind(&Hmaterial<Dim>::saturation_induction, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> relative_permeability_rule;
        relative_permeability_rule = char_('_')
            >> (lit("relative_permeability")  | lit("uri") )
            >> char_('=') >> double_[phx::bind(&Hmaterial<Dim>::relative_permeability, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> ddp_rule;
        ddp_rule = char_('_')
            >> (lit("potential_difference") | lit("ddp") | lit("voltage") )
            >> char_('=') >> double_[phx::bind(&Hmaterial<Dim>::ddp, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> nbr_turns_rule;
        nbr_turns_rule = char_('_')
            >> (lit("nb_turns") | lit("turns") )
            >> char_('=') >> int_[phx::bind(&Hmaterial<Dim>::nbr_turns, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> js_rule;
        js_rule = char_('_')
            >> (lit("current_density_source") | lit("density_source") | lit("js") )
            >> char_('=') >> char_('(')
            >> quoted_string[phx::bind(&Hmaterial<Dim>::js_str, &mat, 0, _1, params)]
            >> char_(',')
            >> quoted_string[phx::bind(&Hmaterial<Dim>::js_str, &mat, 1, _1, params)]
            >> -( char_(',') )
            >> -( quoted_string[phx::bind(&Hmaterial<Dim>::js_str, &mat, 2, _1, params)] )
            >> char_(')')
            >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> rhos_rule;
        rhos_rule = char_('_')
            >> (lit("resistivity_source") | lit("rhos") )
            >> char_('=') >> double_[phx::bind(&Hmaterial<Dim>::rhos, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> intensity_rule;
        intensity_rule = char_('_')
            >> (lit("intensity") | lit("I") )
            >> char_('=') >> double_[phx::bind(&Hmaterial<Dim>::intensity, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> density_rule;
        density_rule = char_('_')
            >> (lit("density") | lit("rho") )
            >> char_('=') >> quoted_string[phx::bind(&Hmaterial<Dim>::density_str, &mat, _1, params)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> specific_heat_rule;
        specific_heat_rule = char_('_')
            >> (lit("specific_heat") | lit("cp") )
            >> char_('=') >> quoted_string[phx::bind(&Hmaterial<Dim>::specific_heat_str, &mat, _1, params)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> young_rule;
        young_rule = char_('_')
            >> (lit("young") | lit("E") | lit("young_modulus"))
            >> char_('=') >> double_[phx::bind(&Hmaterial<Dim>::young_mod, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> poisson_rule;
        poisson_rule = char_('_')
            >> (lit("poisson") | lit("nu") )
            >> char_('=') >> double_[phx::bind(&Hmaterial<Dim>::poisson_coeff, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> dilatation_rule;
        dilatation_rule = char_('_')
            >> (lit("dilatation") | lit("alpha") )
            >> char_('=') >> double_[phx::bind(&Hmaterial<Dim>::dilatation_coeff, &mat, _1)] >> *(eol);

        qi::rule<Iterator, void(), ascii::space_type> volumic_forces_rule;
        volumic_forces_rule = char_('_')
            >> (lit("vforces") | lit("volumic_forces") )
            >> char_('=') >> char_('(')
            >> quoted_string[phx::bind(&Hmaterial<Dim>::volumic_forces_str, &mat, 0, _1, params)]
            >> char_(',')
            >> quoted_string[phx::bind(&Hmaterial<Dim>::volumic_forces_str, &mat, 1, _1, params)]
            >> -( char_(',') )
            >> -( quoted_string[phx::bind(&Hmaterial<Dim>::volumic_forces_str, &mat, 2, _1, params)] )
            >> char_(')')
            >> *(eol);

        bool r = phrase_parse(first, last,

                              //  Begin grammar
                              (
                               name_rule
                               >> physical_name_rule
                               >> elec_cond_rule
                               >> thermal_cond_rule
                               >> lorentz_rule
                               >> alpha_rule
                               >> -(magnetic_permeability_rule)
                               >> -(saturation_induction_rule)
                               >> -(relative_permeability_rule)
                               >> -(ddp_rule)
                               >> -(nbr_turns_rule)
                               >> -(js_rule)
                               >> -(rhos_rule)
                               >> -(intensity_rule)
                               >> -(density_rule)
                               >> -(specific_heat_rule)
                               >> -(young_rule)
                               >> -(poisson_rule)
                               >> -(dilatation_rule)
                               >> -(volumic_forces_rule)
                               ),
                              //  End grammar

                              space);

        if (!r || first != last) // fail if we did not get a full match
            return false;

        return r;
    }

#endif //__MATERIAL_HPP
