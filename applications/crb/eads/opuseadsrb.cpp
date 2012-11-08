// -*- c++ -*-
/**
 *  @file  opuseadsrb.c
 *  @brief The opuseadsrb declaration file
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>

#include <opusdata.hpp>
#include <opusmodelbase.hpp>
#include <opusmodelfactory.hpp>

#include "opuseadsrb.hpp"

namespace Feel
{
class OpusApp   : public Application
{
    typedef Application super;
public:

    typedef OpusModelBase opus_type;
    typedef boost::shared_ptr<opus_type> opus_ptrtype;

    OpusApp( AboutData const& ad )
        :
        super( ad )
    {
        using namespace Feel;
        typedef OpusModelBase opus_type;
        typedef boost::shared_ptr<opus_type> opus_ptrtype;
        this->changeRepository( boost::format( "%1%" )
                                % this->about().appName()
                              );
        M_opus = OpusModelFactory::New( 2 );



    }

    void run()
    {

        M_opus->run();
    }
    void run( const double * X, unsigned long N,
              double * Y, unsigned long P )
    {
        LOG(INFO) << "run from OT\n";

        for ( int i = 0; i < N; ++i )
            LOG(INFO)<< "[opuseadsrb::run] X[" << i << "]="<< X[i] << "\n";

        if ( M_opus )
            M_opus->run( X, N, Y, P );

        for ( int i = 0; i < P; ++i )
            LOG(INFO)<< "[opuseadsrb::run] Y[" << i << "]="<< Y[i] << "\n";

        LOG(INFO) << "done run from OT\n";
    }
private:

    opus_ptrtype M_opus;

}; // Opus

/**
 * \fn makeAbout()
 * \brief Create the About data of the OpusApp
 *
 */
AboutData
makeAbout()
{
    Feel::AboutData about( "opuseadsrb" ,
                           "opuseadsrb" ,
                           "0.1",
                           "2D OPUS/EADS Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

} // Feel


/* myCFuntion returns a non-null value when it fails */
int opuseadsrb( const double * X, unsigned long N,
                double * Y, unsigned long P )
{
    Feel::OpusApp app( Feel::makeAbout() );
    app.run( X, N, Y, P );
    return 0;
}

