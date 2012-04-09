AUTO( coeff, $\gamma_ {\beta} + \module {\beta\SCAL n}$ );
form2( Xh, Xh, M ) +=
    integrate(
        internalfaces( Xh->mesh() ), IM,
        coeff * $h^2/N^ {3.5}$ * (
            trans( jumpt( gradt( u ) ) )
            * jump( grad( v ) )
        )
    );
