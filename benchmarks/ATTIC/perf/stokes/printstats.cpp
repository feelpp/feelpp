
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>

#include <boost/any.hpp>
#include <boost/utility.hpp>
#include <feel/feelcore/feel.hpp>

void
printStats( std::ostream& out, std::vector<std::map<std::string,boost::any> > & stats )
{

    for ( auto it = stats.begin(), en =stats.end(); it!=en; ++it )
    {
        std::for_each( stats.begin(), stats.end(), [] ( std::pair<std::string,boost::any> const& p )
        {
            std::cout << p.first << "\n";
        } );
    }

    out << std::setw( 10 ) << std::right << "levels"
        << std::setw( 10 ) << std::right << "h"
        << std::setw( 15 ) << std::right << "e_u"
        << std::setw( 15 ) << std::right << "ROC"
        << std::setw( 15 ) << std::right << "e_p"
        << std::setw( 15 ) << std::right << "ROC"
        << "\n" ;
    int l=1;

    for ( auto it = stats.begin(), en =stats.end(); it!=en; ++it,++l )
    {
        //std::for_each( it->begin(),it->end(), []( std::pair<std::string,boost::any> const& o ) { std::cout << o.first << "\n"; } );
        std::map<std::string,boost::any> data = *it;
        std::map<std::string,boost::any> datap;
        double rocu = 1, rocp=1;
        double h  = boost::any_cast<double>( data.find( "mesh.h" )->second ) ;
        double u  = boost::any_cast<double>( data.find( "e.l2.u" )->second );
        double p  =  boost::any_cast<double>( data.find( "e.l2.p" )->second );

        if ( l > 1 )
        {
            datap = *boost::prior( it );

            double hp  = boost::any_cast<double>( datap.find( "mesh.h" )->second );
            double up  = boost::any_cast<double>( datap.find( "e.l2.u" )->second );
            double pp  = boost::any_cast<double>( datap.find( "e.l2.p" )->second );
            rocu = std::log10( up/u )/std::log10( hp/h );
            rocp = std::log10( pp/p )/std::log10( hp/h );
        }

        out << std::right << std::setw( 10 ) << l
            << std::right << std::setw( 10 ) << std::fixed  << std::setprecision( 4 ) << h
            << std::right << std::setw( 15 ) << std::scientific << std::setprecision( 2 ) << u
            << std::right << std::setw( 15 ) << std::fixed << std::setprecision( 2 ) << rocu
            << std::right << std::setw( 15 ) << std::scientific << std::setprecision( 2 ) << p
            << std::right << std::setw( 15 ) << std::fixed << std::setprecision( 2 ) << rocp
            << "\n";
    }

    out << std::setw( 10 ) << std::right << "levels"
        << std::setw( 10 ) << std::right << "h"
        << std::setw( 10 ) << std::right << "nElts"
        << std::setw( 10 ) << std::right << "nDof"
        << std::setw( 10 ) << std::right << "nDofu"
        << std::setw( 10 ) << std::right << "nDofp"
        << std::setw( 10 ) << std::right << "nNz"
        << std::setw( 10 ) << std::right << "T_space"
        << std::setw( 10 ) << std::right << "T_matrix"
        << std::setw( 15 ) << std::right << "T_m_assembly"
        << std::setw( 15 ) << std::right << "T_v_assembly"
        << std::setw( 15 ) << std::right << "T_assembly"
        << std::setw( 10 ) << std::right << "T_solve"
        << "\n" ;
    l=1;

    for ( auto it = stats.begin(), en =stats.end(); it!=en; ++it,++l )
    {
        //std::for_each( it->begin(),it->end(), []( std::pair<std::string,boost::any> const& o ) { std::cout << o.first << "\n"; } );
        std::map<std::string,boost::any> data = *it;
        std::map<std::string,boost::any> datap;
        double rocu = 1, rocp=1;
        double h  = boost::any_cast<double>( data.find( "h" )->second ) ;
        using namespace Feel;
        out << std::right << std::setw( 10 ) << l
            << std::right << std::setw( 10 ) << std::fixed  << std::setprecision( 4 ) << h
            << std::right << std::setw( 10 ) << boost::any_cast<size_type>( data.find( "n.mesh.elts" )->second )
            << std::right << std::setw( 10 ) << boost::any_cast<size_type>( data.find( "n.space.dof" )->second )
            << std::right << std::setw( 10 ) << boost::any_cast<size_type>( data.find( "n.space.dof.u" )->second )
            << std::right << std::setw( 10 ) << boost::any_cast<size_type>( data.find( "n.space.dof.p" )->second )
            << std::right << std::setw( 10 ) << boost::any_cast<size_type>( data.find( "n.matrix.nnz" )->second )
            << std::right << std::setw( 10 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "space.time" )->second )
            << std::right << std::setw( 10 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "matrix.init" )->second )
            << std::right << std::setw( 15 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "matrix.assembly" )->second )
            << std::right << std::setw( 15 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "vector.assembly" )->second )
            << std::right << std::setw( 15 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "vector.assembly" )->second )+boost::any_cast<double>( data.find( "matrix.assembly" )->second )
            << std::right << std::setw( 10 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "solver.time" )->second )
            << "\n";
    }
}
