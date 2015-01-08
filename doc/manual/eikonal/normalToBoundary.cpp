/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2008-02-07

 Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
 \file dist2walls.cpp
 \author Vincent Doyeux <vincent.doyeux@ujf-grenoble.fr>
 \date 2014-01-21
 */


#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelpde/reinit_fms.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feelvf/vf.hpp>

#define DIM 2

using namespace Feel;

inline
Feel::po::options_description
makeOptions()
{
	Feel::po::options_description stokesoptions( "Normal options" );
	stokesoptions.add_options()
		( "proj_type", Feel::po::value<std::string>()->default_value( "L2" ), "Projector Type" )
		( "epsilon", Feel::po::value<double>()->default_value( 0.01 ), "epsilon" )
		( "gamma", Feel::po::value<double>()->default_value( 0.01 ), "gamma" )
		;
	return stokesoptions.add( Feel::feel_options()) ;
}

namespace Feel
{
	class normal : public Application{

		public:
			typedef Application super;
			typedef Mesh< Simplex<DIM> > mesh_type;

			typedef Lagrange<1, Scalar> basis_1_type;
			typedef Lagrange<2, Scalar> basis_2_type;
			typedef Lagrange<1,Vectorial> basis_v_type;

			typedef FunctionSpace<mesh_type, basis_1_type> space_P1_type;
			typedef boost::shared_ptr<space_P1_type> space_P1_ptrtype;

			typedef FunctionSpace<mesh_type, basis_2_type> space_P2_type;
			typedef boost::shared_ptr<space_P2_type> space_P2_ptrtype;

			typedef FunctionSpace<mesh_type, basis_v_type> space_v_type;
			typedef boost::shared_ptr<space_v_type> space_v_ptrtype;

			typedef boost::tuple<boost::mpl::size_t<MESH_ELEMENTS>,
							typename MeshTraits<mesh_type>::element_const_iterator,
							typename MeshTraits<mesh_type>::element_const_iterator> range_visu_ho_type;

			typedef OperatorInterpolation<space_P1_type, //espace depart
							space_P2_type, //espace arrivee
							range_visu_ho_type> op_inte_P1_to_P2_type;

			typedef boost::shared_ptr<backend_type> backend_ptrtype;
			backend_ptrtype M_proj_backend;

			normal();
			void run();
	};

	normal::normal()
		: super(),
		M_proj_backend(backend_type::build("petsc"))
	{
	}

	void normal::run()
	{

		auto mesh = loadMesh( _mesh=new mesh_type );

		//space_P1_ptrtype XhP1 = space_P1_type::New(mesh);
		//space_P2_ptrtype XhP2 = space_P2_type::New(mesh);
		//space_v_ptrtype XhV = space_v_type::New(mesh);
		auto XhP1 = Pch<1>(mesh);
		auto XhP2 = Pch<2>(mesh);
		auto XhV = Pchv<2>(mesh);

		/*
		 * Compute the distance
		 */
		auto thefms = fms( XhP1 );
		auto phio = XhP1->element();
		phio = vf::project(XhP1, elements(mesh), h() );
		phio +=vf::project(XhP1, boundaryfaces(mesh), -idv(phio) - h()/100. );
		auto phi_P1 = thefms->march(phio);

		/*
		 * Project it on P2
		 */
		auto op_inte_P1_to_P2 = opInterpolation(_domainSpace=XhP1, _imageSpace=XhP2, _backend=Feel::backend(_rebuild=true));
		auto phi_P2 = XhP2->element();
		op_inte_P1_to_P2->apply(phi_P1,phi_P2);	

		/*
		 * Compute the normal
		 */
		std::map<std::string, ProjectorType> m_projectorType;
		m_projectorType["NODAL"] = NODAL;
		m_projectorType["L2"] = L2;
		m_projectorType["H1"] = H1;
		m_projectorType["DIFF"] = DIFF;
		m_projectorType["HDIV"] = HDIV;
		m_projectorType["HCURL"] = HCURL;
		m_projectorType["LIFT"] = LIFT;
		m_projectorType["CIP"] = CIP;
		CHECK(m_projectorType.count(soption("proj_type"))) << soption("proj_type") <<" is not in the list of possible projector methods\n";
		auto l2pv = projector( XhV ,
				XhV,
				M_proj_backend,
				m_projectorType.find(soption("proj_type"))->second,
				doption("epsilon"),
				doption("gamma"));
		auto n_phi_p   = vf::project(_space=XhV, _range=boundaryfaces(mesh), _expr=trans(gradv(phi_P2))/sqrt( gradv(phi_P2) * trans(gradv(phi_P2))));
		auto n_phi_s = l2pv->project(_space=XhV, _range=boundaryfaces(mesh), _expr=trans(gradv(phi_P2)));
		n_phi_s = vf::project(XhV, boundaryfaces(mesh),idv(n_phi_s)/sqrt(idv(n_phi_s.comp(X))*idv(n_phi_s.comp(X))+idv(n_phi_s.comp(Y))*idv(n_phi_s.comp(Y))));

		auto Nmesh = vf::project(XhV, boundaryfaces(mesh),-N());
		auto Nexact = vf::project( _space=XhV, _range=boundaryfaces(mesh), _expr=vec(Px(),Py())/sqrt(Px()*Px()+Py()*Py()));

		auto err_phi_p = normL2( _range=markedfaces(mesh,"surface"), _expr=idv(n_phi_p)-idv(Nexact))/normL2( _range=markedfaces(mesh,"surface"), _expr=idv(Nexact));
		auto err_phi_s = normL2( _range=markedfaces(mesh,"surface"), _expr=idv(n_phi_s)-idv(Nexact))/normL2( _range=markedfaces(mesh,"surface"), _expr=idv(Nexact));
		auto err_N     = normL2( _range=markedfaces(mesh,"surface"), _expr=idv(Nmesh  )-idv(Nexact))/normL2( _range=markedfaces(mesh,"surface"), _expr=idv(Nexact));

		std::cout << 
			err_phi_p<< "\t" << 
			err_phi_s<< "\t" << 
			err_N    << std::endl; 


		auto exp = exporter(_mesh=mesh, _name="NormalToBoundary");
		exp->step(0)->add("phi_P2", phi_P2);
		exp->step(0)->add("phi_P1", phi_P1);
		exp->step(0)->add("n_phi_p", n_phi_p);
		exp->step(0)->add("n_phi_s", n_phi_s);
		exp->step(0)->add("Nmesh", Nmesh);
		exp->step(0)->add("Nexact", Nexact);
		exp->save();
	}

}

int main( int argc, char** argv )
{
	Feel::Environment env( _argc=argc, _argv=argv,
			_desc=makeOptions(),
			_about=about(_name="NormalToBoundary",
				_author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));
	normal n;
	n.run();
	return 0;
}
