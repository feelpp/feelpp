AUTO( def,  0.5*( grad( v )  + trans( grad( v ) ) ) );
AUTO( deft, 0.5*( gradt( u ) + trans( gradt( u ) ) ) );
form2( Xh, Xh, M ) =
    integrate(
        elements( Xh->mesh() ), IM,
        alpha*trans( idt( u ) )*id( v )
        + 2.0*nu*trace( trans( deft )*def )
        + trans( gradt( u )*idv( beta ) )*id( v )
        - div( v )*idt( p ) + divt( u )*id( q )
    );
