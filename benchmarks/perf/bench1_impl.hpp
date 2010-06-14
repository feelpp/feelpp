/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-06-14

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
/**
   \file bench1_impl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-06-14
 */
#include <bench1.hpp>

namespace Life
{

template<typename FSType>
void
Bench1::R( boost::shared_ptr<FSType> const& Xh  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );


    typedef fusion::vector<Lagrange<0, Scalar> > dp0_basis_type;
    typedef FunctionSpace<typename FSType::mesh_type, dp0_basis_type, Discontinuous> dp0_space_type;
    typename dp0_space_type::pointer_type P0h = dp0_space_type::New( Xh->mesh() );
    typename dp0_space_type::element_type w( P0h );
    w = vf::project( P0h, elements(Xh->mesh()), Px() );
    form2(Xh,Xh,M,_init=true);
    boost::timer timer;

#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    //
    // R
    //
    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  idt(u)*id(v) );
    Log() << " o- R<const> time : " << timer.elapsed() << " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  idv(w)*idt(u)*id(v) );
    Log() << " o- R<p0 const> time : " << timer.elapsed() << " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  ((Px()^(3))+(Py()^(2))*Pz())*idt(u)*id(v) );
    Log() << " o-   R<xyz> time : " << timer.elapsed() << " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  idv(w)*((Px()^(3))+(Py()^(2))*Pz())*idt(u)*id(v) );
    Log() << " o-   R<dp0 xyz> time : " << timer.elapsed() << " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  val((Px()^(3))+(Py()^(2))*Pz())*idt(u)*id(v) );
    Log() << " o-   R<val xyz> time : " << timer.elapsed() << " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  idv(w)*val((Px()^(3))+(Py()^(2))*Pz())*idt(u)*id(v) );
    Log() << " o-   R<dp0 val xyz> time : " << timer.elapsed() << " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";
#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}

template<typename FSType>
void
Bench1::D( boost::shared_ptr<FSType> const& Xh  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    form2(Xh,Xh,M,_init=true);
    boost::timer timer;

#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    //
    // D
    //
    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  gradt(u)*trans(grad(v)) );
    Log() << " o- D<const> time : " << timer.elapsed()  << " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  dxt(u)*dx(v)+dyt(u)*dy(v) );
    Log() << " o- D<const2> time : " << timer.elapsed()  << " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  dxt(u)*dx(v)+dyt(u)*dy(v)+dzt(u)*dz(v) );
    Log() << " o- D<const3> time : " << timer.elapsed()  << " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  val((Px()^(3))+(Py()^(2))*Pz())*gradt(u)*trans(grad(v)) );
    Log() << " o-   D<xyz> time : " << timer.elapsed()  << " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}
template<typename FSType>
void
Bench1::DR( boost::shared_ptr<FSType> const& Xh )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    form2(Xh,Xh,M,_init=true);
    boost::timer timer;

#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    //
    // DR
    //
    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  gradt(u)*trans(grad(v))+idt( u )*id( v ));
    Log() << " o- DR<const> time : " << timer.elapsed() <<  " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  val((Px()^(3))+(Py()^(2))*Pz())*(gradt(u)*trans(grad(v))+idt( u )*id( v ) ));
    Log() << " o-   DR<xyz> time : " << timer.elapsed() <<  " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";
#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}
template<typename FSType>
void
Bench1::ADR( boost::shared_ptr<FSType> const& Xh, mpl::int_<1>  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    form2(Xh,Xh,M,_init=true);
    boost::timer timer;

#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    //
    // ADR
    //
    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),
                                 gradt(u)*trans(grad(v))+idt( u )*id( v ) +
                                 (gradt(u)*vec(constant(1.0)))*id(v));
    Log() << " o- ADR<const> time : " << timer.elapsed() <<  " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),
                                 val((Px()^(3))+(Py()^(2))*Pz())*(gradt(u)*trans(grad(v))+idt( u )*id( v )) +
                                 (gradt(u)*vec(val((Px()^(3))+(Py()^(2))*Pz())))*id(v));
    Log() << " o-   ADR<xyz> time : " << timer.elapsed() <<  " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";
#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}

template<typename FSType>
void
Bench1::ADR( boost::shared_ptr<FSType> const& Xh, mpl::int_<2>  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    form2(Xh,Xh,M,_init=true);
    boost::timer timer;

#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    //
    // ADR
    //
    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),
                                           gradt(u)*trans(grad(v))+idt( u )*id( v ) +
                                           (gradt(u)*vec(constant(1.0),constant(1.0)))*id(v));
    Log() << " o- ADR<const> time : " << timer.elapsed() <<  " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),
                                           val((Px()^(3))+(Py()^(2))*Pz())*(gradt(u)*trans(grad(v))+idt( u )*id( v )) +
                                 (gradt(u)*vec(val((Px()^(3))+(Py()^(2))*Pz()),val((Px()^(3))+(Py()^(2)))))*id(v));
    Log() << " o-   ADR<xyz> time : " << timer.elapsed() <<  " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";
#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}
template<typename FSType>
void
Bench1::ADR( boost::shared_ptr<FSType> const& Xh, mpl::int_<3>  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype  M( M_backend->newMatrix( Xh, Xh ) );

    form2(Xh,Xh,M,_init=true);
    boost::timer timer;

#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    //
    // ADR
    //
    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),
                                           gradt(u)*trans(grad(v))+idt( u )*id( v ) +
                                           (gradt(u)*vec(constant(1.0),constant(1.0),constant(1.0)))*id(v));
    Log() << " o- ADR<const> time : " << timer.elapsed() <<  " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),
                                 val((Px()^(3))+(Py()^(2))*Pz())*(gradt(u)*trans(grad(v))+idt( u )*id( v )) +
                                 (gradt(u)*vec(val((Px()^(3)+(Py()^(2)))*Pz()),val((Px()^(3))+(Py()^(2))),val(Px()^(3))))*id(v));
    Log() << " o-   ADR<xyz> time : " << timer.elapsed() <<  " " << timer.elapsed()*1e6/Xh->mesh()->numElements() << "\n";
#if defined(HAVE_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}
template<typename MeshType, int Order>
void
Bench1::bench1( boost::shared_ptr<MeshType> & mesh )
{


    const int nDim = MeshType::nDim;
    std::ostringstream femstr;
    femstr << "FEM_PK(" << nDim << "," << Order << ")";
    std::ostringstream geostr;
    geostr << "GT_PK(" << nDim << "," << 1 << ")";

    Log() << "------------------------------------------------------------" << "\n";
    Log() << "dimension : " << nDim << "\n";
    Log() << "      fem : " << femstr.str() << "\n";
    Log() << "      geo : " << geostr.str() << "\n";
    Log() << "++++++++++++++++++++++++++++++" << "\n";

    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;
    typedef FunctionSpace<MeshType, basis_type> space_type;
    typename space_type::pointer_type Xh = space_type::New( mesh );



    //IM<nDim, 2*Order, double, Simplex> im;
    //IM<nDim, 2*Order,double> im;
    //IM_PK<nDim, 4> im;
    //IMSimplex<2,2*Order> im;

    //v.space()->gm()->setCacheInformation( QDR, mesh->numElements(), QDR_NPTS);
    //v.space()->fe()->setCacheInformation( QDR, mesh->numElements(), QDR_NPTS);

    Log() << "dof : " << Xh->nDof() << "\n"
                     << "elt : " << Xh->mesh()->numElements() << "\n";


    boost::timer timer;

    R( Xh );
    D( Xh );
    DR( Xh );
    ADR( Xh, mpl::int_<nDim>() );
    Log() << "------------------------------------------------------------" << "\n";

}



}
