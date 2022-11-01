/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 07 juil. 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_EXTENDERFROMINTERFACE_HPP
#define FEELPP_EXTENDERFROMINTERFACE_HPP 1


#include <feel/feells/reinit_fms.hpp>

namespace Feel {

template<int Dim,int GeoOrder=1, template<uint16_type,int,uint16_type> class Convex = Simplex>
class ExtenderFromInterface
{
public:
    
    enum state_type {TODO=0, DONE=1};

    using convex_type = Convex<Dim,GeoOrder,Dim>;
    using mesh_ptrtype = ls_mesh_ptrtype<Dim,GeoOrder,Convex>;
    using space_ptrtype = ls_space_ptrtype<Dim,GeoOrder,Convex>;
    using element_type = ls_element_type<Dim,GeoOrder,Convex>;
    template<typename Storage>
    using element_s_type = typename ls_space_type<Dim,GeoOrder,Convex>::template Element<double,Storage>;

    ExtenderFromInterface( space_ptrtype Xh, element_type const& phi ) 
        : 
        M_Xh( Xh ), 
        M_phi( Xh->element() ),
        M_marker( Xh->element() ),
        M_states_init(M_Xh->nDof(), TODO),
        M_interf_id()
        {
            this->update(phi);
        }

    template<typename Storage>
    void extendFromInterface( element_s_type<Storage>& field );

    double dist(double a, double b ) const
        {
            return (a*a + b*b);
        }
    element_type const& marker() const { return M_marker; }

    void update( element_type const& phi );

private:
    
    space_ptrtype M_Xh;
    element_type M_phi,M_marker;
    std::vector<state_type> M_states_init;
    std::vector<size_type> M_interf_id;
}; // ExtenderFromInterface

template<int Dim,int GeoOrder, template<uint16_type,int,uint16_type> class Convex>
void
ExtenderFromInterface<Dim,GeoOrder,Convex>::update( element_type const& phi )
{
    M_phi = phi;
    M_phi.close();
    
    auto it_elt = M_Xh->mesh()->beginElement();
    auto en_elt = M_Xh->mesh()->endElement();

    int nb_elt_crossed = 0;
    for (;it_elt!=en_elt; it_elt++)
    {
        auto const& elt = it_elt->second;
        int nbplus=0;
        int nbminus=0;
        std::vector<size_type> indices_nodes( convex_type::numPoints );
        for (int j=0; j<convex_type::numPoints; j++)
        {
            size_type index = M_phi.start() + M_Xh->dof()->localToGlobal( elt.id(), j, 0 ).index();
            //double index = velocX.localToGlobal(it_elt->id(), j, 0);
            indices_nodes[j]=index;

            if ( (M_phi)[index] > 0.)
                ++nbplus;
            else if ((M_phi)[index] < 0.)
                ++nbminus;
        }

        //if elt crossed by interface -> store its informations
        if ( (nbminus != convex_type::numPoints) && (nbplus!=convex_type::numPoints) )
        {
            LOG(INFO) << "element crossed " << elt.id();
            nb_elt_crossed++;
            for (int j=0; j<convex_type::numPoints; j++)
            {
                if (M_states_init[indices_nodes[j]]!=DONE)
                {
                    M_interf_id.push_back(indices_nodes[j]);
                    M_marker[indices_nodes[j]]=1; 
                    M_states_init[indices_nodes[j]] = DONE;
                }
            }
        }
    }
    LOG(INFO) << "interf " << M_interf_id;
    LOG(INFO) << "states_init " << M_states_init;

    if (M_interf_id.size() == 0 )
    {
        LOG(WARNING) <<"no element in the interface\n";
    }
    // **************************************************************

}

template<int Dim,int GeoOrder, template<uint16_type,int,uint16_type> class Convex>
template<typename Storage>
void
ExtenderFromInterface<Dim,GeoOrder,Convex>::extendFromInterface( element_s_type<Storage>& field )
{
    std::vector<state_type> states( M_states_init );
    LOG(INFO) << "states " << M_states_init;

    auto posX = M_Xh->element(),posY=M_Xh->element();
    posX.on(_range=elements(M_Xh->mesh()), _expr=Px());
    posY.on(_range=elements(M_Xh->mesh()), _expr=Py());
    
    auto dofpt_it = M_Xh->dof()->dofPointBegin();
    auto dofpt_en = M_Xh->dof()->dofPointEnd();
    
    for (int compt=0; dofpt_it != dofpt_en; dofpt_it++)
    {
        
        auto dofpt_coord = dofpt_it->second.template get<0>();
        auto dofpt_id = dofpt_it->second.template get<1>();

        if (states[dofpt_id]==TODO)
        {

            LOG(INFO) << "dof " << dofpt_coord << " " << dofpt_id << "  TOBEDONE";
            compt++;

            double mindist = 10000000;
            int ind;

            for (int i=0; i<M_interf_id.size(); i++)
            {
                //point of the interface to track
                auto __px = posX[M_interf_id[i]];
                auto __py = posY[M_interf_id[i]];

                double t_dist = this->dist((dofpt_coord[0] - __px), (dofpt_coord[1] - __py));

                if (t_dist < mindist)
                {
                    mindist = t_dist;
                    ind = i;
                }
            }
            
            field[dofpt_id] = field[M_interf_id[ind]];
            LOG(INFO) << "set " << dofpt_id << " to " << field[M_interf_id[ind]];
            states[dofpt_id]=DONE;
        }
    }

}

}
#endif

