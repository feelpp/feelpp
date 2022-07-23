/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): JB Wahl <wahl@math.unistra.fr>
 Date: 18 April 2017

 Copyright (C) 2016 Feel++ Consortium

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

#ifndef FEELPP_SER_HPP
#define FEELPP_SER_HPP 1

#include <boost/shared_ptr.hpp>

namespace Feel{

template <typename CRBType>
class SER
{
public :
    typedef CRBType crb_type;
    typedef std::shared_ptr<crb_type> crb_ptrtype;

    typedef typename crb_type::model_type crbmodel_type;
    typedef std::shared_ptr<crbmodel_type> crbmodel_ptrtype;

    typedef typename crbmodel_type::model_type model_type;
    typedef std::shared_ptr<model_type> model_ptrtype;


    SER( crb_ptrtype crb, crbmodel_ptrtype crbmodel );

    void run();

    crb_ptrtype crb( int const& level )
    {
        return M_crbs[level];
    }

private :
    std::vector<crbmodel_ptrtype> M_models;
    std::vector<crb_ptrtype> M_crbs;
    std::string M_prefix;
    int M_nbLevels;
    bool M_useRbInEimMuSelect;
    bool M_useRbInEimBasis;
    int M_eimFreq;
}; //class SER


template <typename CRBType>
SER<CRBType>::SER( crb_ptrtype crb, crbmodel_ptrtype crbmodel )
    : M_prefix(crbmodel->prefix()),
      M_nbLevels(ioption( _prefix=M_prefix,_name="ser.nb-levels" )),
      M_useRbInEimMuSelect(boption(_prefix=M_prefix,_name="ser.use-rb-in-eim-mu-selection")),
      M_useRbInEimBasis(boption(_prefix=M_prefix,_name="ser.use-rb-in-eim-basis-build")),
      M_eimFreq(ioption(_prefix=M_prefix,_name="ser.eim-frequency"))
{
    M_models.push_back( crbmodel );
    M_crbs.push_back( crb );
} // SER constructor


template <typename CRBType>
void
SER<CRBType>::run()
{
    bool do_offline_eim = false;

    for( int ser_level=0; ser_level < M_nbLevels; ++ser_level )
    {
        if ( ser_level > 0 ) // create new crb and model
        {
            auto model = std::make_shared<crbmodel_type>( M_crbs.front()->name(), crb::stage::offline, ser_level );
            auto crb = crb_type::New( M_crbs.front()->name(), model,
                                      crb::stage::offline, (boost::format("ser%1%")%ser_level).str() );
            M_models.push_back( model );
            M_crbs.push_back( crb );
        }

        auto model = M_models.back();
        auto crb = M_crbs.back();
        auto crbInEIM = ( ser_level > 0 )? M_crbs[ser_level-1] : crb;

        auto eim_sc_vector = model->scalarContinuousEim();
        auto eim_sd_vector = model->scalarDiscontinuousEim();
        for( auto eim_sc : eim_sc_vector )
            do_offline_eim = do_offline_eim || eim_sc->offlineStep();
        for( auto eim_sd : eim_sd_vector )
            do_offline_eim = do_offline_eim || eim_sd->offlineStep();

        auto deim_vector = model->deimVector();
        auto mdeim_vector = model->mdeimVector();
        for ( auto deim : deim_vector )
            do_offline_eim = do_offline_eim || deim->offlineStep();
        for ( auto mdeim : mdeim_vector )
            do_offline_eim = do_offline_eim || mdeim->offlineStep();

        do
        {
            //Begin with rb since first eim has already been built in initModel
            //this->loadDB(); // update AffineDecomposition and enrich RB database
            tic();
            crb->setOfflineStep( true );
            do  // SER r-adaptation for RB
            {
                crb->setAdaptationSER( false ); //re-init to false
                crb->offline();
            }
            while( crb->adaptationSER() );
            toc("SER - crb offline", FLAGS_v>0);

            crb->setRebuild( false ); //do not rebuild since co-build is not finished
            int use_rb = M_useRbInEimMuSelect || M_useRbInEimBasis;

            tic();
            if( do_offline_eim && crb->offlineStep() ) //Continue to enrich EIM functionspace only is RB is not complete
            {
                // Be careful if you have multiple EIMs to not overstep the bounds of the vectors when computing beqtQm (see #1130)
                do_offline_eim = false; //re-init
                for( auto const& eim_sc : eim_sc_vector )
                {
                    //if ( ser_level > 0 )
                    eim_sc->setDBSubDirectory( (boost::format("eim_ser%1%")%ser_level).str() );
                    eim_sc->setRestart( false ); //do not restart since co-build is not finished

                    if( use_rb )
                    {
                        eim_sc->setRB( crbInEIM ); // update rb model member to be used in eim offline
                    }
                    do //r-adaptation for EIM
                    {
                        eim_sc->setAdaptationSER( false ); //re-init to false
                        eim_sc->offline();
                    } while( eim_sc->adaptationSER() );

                    do_offline_eim = do_offline_eim || eim_sc->offlineStep();
                }
                for( auto const& eim_sd : eim_sd_vector )
                {
                    //if ( ser_level > 0 )
                    eim_sd->setDBSubDirectory( (boost::format("eim_ser%1%")%ser_level).str() );
                    eim_sd->setRestart( false ); //do not restart since co-build is not finished

                    if( use_rb )
                    {
                        eim_sd->setRB( crbInEIM ); // update rb model member to be used in eim offline
                    }
                    do //r-adaptation for EIM
                    {
                        eim_sd->setAdaptationSER( false ); //re-init to false
                        eim_sd->offline();
                    } while( eim_sd->adaptationSER() );

                    do_offline_eim = do_offline_eim || eim_sd->offlineStep();
                }

                for ( auto const& deim : deim_vector )
                {
                    deim->setRestart( false );
                    deim->setSerFrequency( M_eimFreq );
                    deim->setSerUseRB( use_rb );

                    if ( use_rb )
                    {
                        deim->setRB( crbInEIM ); //update rb model member to be used in eim offline
                    }

                    if ( deim->offlineStep() )
                        deim->offline();

                    do_offline_eim = do_offline_eim || deim->offlineStep();
                }

                for ( auto const& mdeim : mdeim_vector )
                {
                    mdeim->setRestart( false );
                    mdeim->setSerFrequency( M_eimFreq );
                    mdeim->setSerUseRB( use_rb );

                    if ( use_rb )
                    {
                        mdeim->setRB( crbInEIM ); //update rb model member to be used in eim offline
                    }

                    if ( mdeim->offlineStep() )
                        mdeim->offline();

                    do_offline_eim = do_offline_eim || mdeim->offlineStep();
                }

                model->assemble(); //Affine decomposition has changed since eim has changed
            }
            toc("SER - eim offline + re-assemble", FLAGS_v>0);

        } while( crb->offlineStep() );

    } // for( int ser_level=0; ser_level < nbLevels; ++ser_level )


} // run


} //namespace Feel

#endif // FEELPP_SER_HPP
