#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/mpfr.hpp>



namespace Feel
{
namespace mpfr
{
MpfrClass::PrecisionType &MpfrClass::CurrPrecision = __gmpfr_default_fp_bit_precision;
//MpfrClass::PrecisionType &MpfrClass::CurrPrecision = ( uint )mpfr_get_default_prec();
// ??? MpfrClass::PrecisionType &MpfrClass::CurrPrecision = mpfr_get_default_prec(); does not do what is wanted?
RoundingMode &MpfrClass::CurrRndMode = __gmpfr_default_rounding_mode;



//--------------------------------------------------------------
//
// Constructors and destructors
//
//--------------------------------------------------------------

MpfrClass::MpfrClass ()
{
    mpfr_init2( mpfr_rep, CurrPrecision );
    nbref = new RefCounter( 0 );
    inexact=new InexactFlag();
}

MpfrClass::MpfrClass ( double d,
                       RoundingMode rnd,
                       PrecisionType prec )
{
    mpfr_init2( mpfr_rep, prec );

    nbref=new RefCounter( 1 );
    inexact=new InexactFlag();
    inexact->refvalue() = mpfr_set_d( mpfr_rep, d, rnd );
}

MpfrClass::MpfrClass ( long double d,
                       RoundingMode rnd,
                       PrecisionType prec )
{
    mpfr_init2( mpfr_rep, prec );

    nbref=new RefCounter( 1 );
    inexact=new InexactFlag();
    inexact->refvalue() = mpfr_set_ld( mpfr_rep, d, rnd );
}

MpfrClass::MpfrClass ( int i,
                       RoundingMode rnd,
                       PrecisionType prec )
{
    mpfr_init2( mpfr_rep, prec );
    nbref=new RefCounter( 1 );
    inexact=new InexactFlag();
    inexact->refvalue() = mpfr_set_si( mpfr_rep, ( long )( i ), rnd );
}

MpfrClass::MpfrClass ( unsigned int i,
                       RoundingMode rnd,
                       PrecisionType prec )
{
    mpfr_init2( mpfr_rep, prec );
    nbref=new RefCounter( 1 );
    inexact=new InexactFlag();
    inexact->refvalue() = mpfr_set_ui( mpfr_rep, ( unsigned long int )i, rnd );
}

MpfrClass::MpfrClass ( long int i,
                       RoundingMode rnd,
                       PrecisionType prec )
{
    mpfr_init2( mpfr_rep, prec );
    nbref=new RefCounter( 1 );
    inexact=new InexactFlag();
    inexact->refvalue() = mpfr_set_si( mpfr_rep, i, rnd );
}

MpfrClass::MpfrClass ( unsigned long int i,
                       RoundingMode rnd,
                       PrecisionType prec )
{
    mpfr_init2( mpfr_rep, prec );
    nbref=new RefCounter( 1 );
    inexact=new InexactFlag();
    inexact->refvalue() = mpfr_set_ui( mpfr_rep, i, rnd );
}

MpfrClass::MpfrClass( string s )
{
    mpfr_init2( mpfr_rep,MpfrClass::CurrPrecision );
    nbref=new RefCounter( 1 );
    inexact=new InexactFlag();

    if ( mpfr_set_str( mpfr_rep, const_cast<char *>( s.c_str() ), 10, MpfrClass::CurrRndMode ) )
    {
        cerr << "Pb while reading " << s << endl ;
        // maybe throw an exception...
    }
}

MpfrClass::MpfrClass ( mpz_srcptr z,
                       RoundingMode rnd,
                       PrecisionType prec )
{
    mpfr_init2( mpfr_rep, prec );
    nbref=new RefCounter( 1 );
    inexact=new InexactFlag();
    inexact->refvalue() = mpfr_set_z( mpfr_rep, z, rnd );
}

MpfrClass::MpfrClass ( mpq_srcptr q,
                       RoundingMode rnd,
                       PrecisionType prec )
{
    mpfr_init2( mpfr_rep, prec );
    nbref=new RefCounter( 1 );
    inexact=new InexactFlag();
    inexact->refvalue() = mpfr_set_q( mpfr_rep, q, rnd );
}

MpfrClass::MpfrClass ( mpfr_t r,
                       RoundingMode rnd,
                       PrecisionType prec )
{
    mpfr_init2( mpfr_rep, prec );
    nbref = new RefCounter( 1 );
    inexact=new InexactFlag();
    inexact->refvalue() = mpfr_set ( mpfr_rep, r, rnd );
}

MpfrClass::MpfrClass ( const MpfrClass& r )
{
    mpfr_rep[0]=r.mpfr_rep[0];
    inexact = r.inexact;
    nbref=r.nbref;
    nbref->incr();
}

MpfrClass::~MpfrClass ()
{
    if ( nbref->decr() <= 0 )
    {
        mpfr_clear( mpfr_rep );
        delete inexact;
        delete nbref;
    }
}


//--------------------------------------------------------------
//
// Assignment and physical copy
//
//--------------------------------------------------------------

// Assignment is a logical copy
MpfrClass& MpfrClass::operator = ( const MpfrClass& r )
{
    //(*this).copy(r);
    if ( this == &r )
        return *this;

    else if ( nbref->decr() <= 0 )
    {
        mpfr_clear( mpfr_rep );
        delete nbref;
        delete inexact;
    }

    mpfr_rep[0]=r.mpfr_rep[0];
    nbref=r.nbref;
    nbref->incr();
    inexact = r.inexact;
    return *this;
}

// Physical copy
MpfrClass& MpfrClass::copy ( const MpfrClass& r, RoundingMode rnd, PrecisionType prec )
{
    if ( nbref->decr() <= 0 )
    {
        mpfr_clear( mpfr_rep );
        delete nbref;
        delete inexact;
    }

    nbref=new RefCounter( 1 );
    inexact = new InexactFlag();

    // If the desired precision is different from the current one
    if ( prec != GetDefaultPrecision() )
    {
        mpfr_init2( mpfr_rep, prec );
        inexact->refvalue() = mpfr_set( mpfr_rep, r.mpfr_rep, rnd );
    }

    else
        inexact->refvalue() = mpfr_init_set( mpfr_rep, r.mpfr_rep, rnd );

    return *this;
}


//--------------------------------------------------------------
//
// Precision and rounding: consultation, modification
//
//--------------------------------------------------------------

void MpfrClass::SetDefaultPrecision ( MpfrClass::PrecisionType newprec )
{
    mpfr_set_default_prec( newprec );
}

// ??? Rien sur les flags inexact ?
void MpfrClass::SetPrecision ( MpfrClass::PrecisionType newprec )
{
    nbref->refvalue() = 0;
    inexact->refvalue() = UNAFFECTED_INEXACT_FLAG;
    mpfr_set_prec( mpfr_rep, newprec );
}

const MpfrClass::PrecisionType MpfrClass::GetDefaultPrecision ()
{
    return CurrPrecision;
}
// idem: { return mpfr_get_default_prec(); }

MpfrClass::PrecisionType MpfrClass::GetPrecision () const
{
    return mpfr_get_prec( mpfr_rep );
}

void MpfrClass::SetDefaultRndMode ( RoundingMode newrndmode )
{
    mpfr_set_default_rounding_mode( newrndmode );
}

const RoundingMode MpfrClass::GetDefaultRndMode ()
{
    return MpfrClass::CurrRndMode;
}

void MpfrClass::ChangePrec ( RoundingMode rnd, PrecisionType prec )
// { inexact->refvalue() = mpfr_round_prec(mpfr_rep, rnd, prec); }
// above: old version for mpfr < 2.0.2
{
    if ( nbref->getvalue() <= 0 )
    {
        inexact->refvalue() = mpfr_prec_round( mpfr_rep, prec, rnd );
    }

    else
    {
        MpfrClass tmp ( *this );
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter();
        inexact = new InexactFlag();
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_set( mpfr_rep, tmp.mpfr_rep, rnd );
    }
}


// Rounding in the direction of rnd2 when the result is computed
// in the direction of rnd1, with a number of wrong bits given as argument.
// Should be in place or not? In place
int MpfrClass::CanRound ( mp_prec_t nb_wrong_bits, RoundingMode rnd1,
                          RoundingMode rnd2, PrecisionType prec )
{
    return mpfr_can_round( mpfr_rep, nb_wrong_bits, rnd1, rnd2, prec );
}


//--------------------------------------------------------------
//
// Constants: Pi, log(2) and Euler constant
//
//--------------------------------------------------------------

MpfrClass MpfrClass::Pi ( RoundingMode rnd, PrecisionType prec )
{
    MpfrClass res( 0,MpfrClass::CurrRndMode, prec );
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_const_pi( res.mpfr_rep, rnd );
    return res;
}

MpfrClass MpfrClass::Log2 ( RoundingMode rnd, PrecisionType prec )
{
    MpfrClass res( 0,MpfrClass::CurrRndMode, prec );
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_const_log2( res.mpfr_rep, rnd );
    return res;
}

MpfrClass MpfrClass::Euler ( RoundingMode rnd, PrecisionType prec )
{
    MpfrClass res( 0,MpfrClass::CurrRndMode, prec );
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_const_euler( res.mpfr_rep, rnd );
    return res;
}


//--------------------------------------------------------------
//
// Input Output
//
//--------------------------------------------------------------

// Outputs the value on the stream o in base 'base' with 'nb_digits'
// with rounding mode rnd
ostream& MpfrClass::put( ostream& o,
                         RoundingMode rnd, PrecisionType /*prec*/,
                         int base, int nb_digits ) const
{
    bool neg = ( sign( *this )<0 );

    if ( isinf( *this ) )
    {
        if ( neg )
            o << "-";

        o << "Inf";
        return o;
    }

    if ( isnan( *this ) )
    {
        o << "NaN";
        return o;
    }

    if ( iszero( *this ) )
    {
        if ( neg )
            o << "-";

        o << "0";
        return o;
    }

    // Based on mpfr_get_str to create the string
    // Output done in another way (digit . fraction x E...).
    char *s, *fraction, first_digit;
    mp_exp_t e;
    int size_s;
    s = mpfr_get_str( NULL, &e, base, nb_digits, mpfr_rep, rnd );
    size_s = strlen( s );

    // to print the floating point
    if ( s[0] == '-' )
    {
        o << s[0];
        first_digit = s[1];
        fraction = &s[2];
        size_s++;
    }

    else
    {
        first_digit = s[0];
        fraction = &s[1];
    }

    o << first_digit << "." << fraction;
    ( *__gmp_free_func )( s, size_s );

    //Recommended way: cf. below
    //free(s, size_s);
    if ( --e )
    {
        o << "E" << e ;
    }

    return o;
}

ostream& operator << ( ostream& o, const MpfrClass& r )
{
    return r.put( o );
}

istream& operator >> ( istream& i, MpfrClass& r )
// Mimicking mpfr_set_str and Integer::>> of Givaro
{
    if ( !i.good() )
    {
        cerr << " Pb reading on the input stream " << endl;
        // maybe throw an exception...
    }

    if ( !i )
        return i;

    // read white spaces
    i >> ws;

    // reading the string,
    // the conversion from a string to a mpfr is done by mpfr_set_str
    bool noend = true;
    int nread=0;
    size_t alloc_size=100;
    char c, *str=( char * )( *__gmp_allocate_func )( alloc_size );

    // read the characters on i until a white space is read
    while ( noend )
    {
        // read until alloc_size char are read, or a white space is encountered
        while ( ( noend ) && ( nread < alloc_size ) )
        {
            i.get( c );

            if ( i.eof() )
            {
                noend=false;
            }

            else if ( isspace( c ) )
            {
                noend=false;
                i.putback( c );
            }

            else
            {
                str[nread++] = c;
            }
        }

        if ( nread >= alloc_size )
        {
            size_t old_alloc_size = alloc_size;
            alloc_size = alloc_size * 3 / 2;
            str = ( char * )( *__gmp_reallocate_func )( str, old_alloc_size, alloc_size );
        }
    }

    str[nread] = '\0';

    if ( mpfr_set_str( r.mpfr_rep, str, 10, MpfrClass::CurrRndMode ) )
    {
        cerr << " Pb reading on the input stream " << endl;
        // maybe throw an exception...
    }

    return i;
}


//--------------------------------------------------------------
//
// Comparisons
//
//--------------------------------------------------------------

int compare ( const MpfrClass& r1, const MpfrClass& r2 )
{
    return mpfr_cmp( r1.mpfr_rep, r2.mpfr_rep );
}

int compare ( const MpfrClass& r1, const double r2 )
{
    return mpfr_cmp_d( r1.mpfr_rep, r2 );
}

int compare ( const MpfrClass& r1, const int r2 )
{
    return mpfr_cmp_si( r1.mpfr_rep, long( r2 ) );
}

int compare ( const MpfrClass& r1, const unsigned int r2 )
{
    return mpfr_cmp_ui( r1.mpfr_rep, ( unsigned long ) r2 ) ;
}

int compare ( const MpfrClass& r1, const long int r2 )
{
    return mpfr_cmp_si( r1.mpfr_rep, r2 );
}

int compare ( const MpfrClass& r1, const unsigned long r2 )
{
    return mpfr_cmp_ui( r1.mpfr_rep, r2 );
}

bool operator == ( const MpfrClass& r1, const MpfrClass& r2 )
{
    return ( compare( r1,r2 ) == 0 );
}

bool operator == ( const double r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) == 0 );
}

bool operator == ( const int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) == 0 );
}

bool operator == ( const unsigned int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) == 0 );
}

bool operator == ( const long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) == 0 );
}

bool operator == ( const unsigned long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) == 0 );
}

bool operator == ( const MpfrClass& r1, const double r2 )
{
    return ( compare( r1,r2 ) == 0 );
}

bool operator == ( const MpfrClass& r1, const int r2 )
{
    return ( compare( r1,r2 ) == 0 );
}

bool operator == ( const MpfrClass& r1, const unsigned int r2 )
{
    return ( compare( r1,r2 ) == 0 );
}

bool operator == ( const MpfrClass& r1, const long int r2 )
{
    return ( compare( r1,r2 ) == 0 );
}

bool operator == ( const MpfrClass& r1, const unsigned long int r2 )
{
    return ( compare( r1,r2 ) == 0 );
}

bool operator != ( const MpfrClass& r1, const MpfrClass& r2 )
{
    return ( compare( r1,r2 ) != 0 );
}

bool operator != ( const double r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) != 0 );
}

bool operator != ( const int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) != 0 );
}

bool operator != ( const unsigned int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) != 0 );
}

bool operator != ( const long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) != 0 );
}

bool operator != ( const unsigned long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) != 0 );
}

bool operator != ( const MpfrClass& r1, const double r2 )
{
    return ( compare( r1,r2 ) != 0 );
}

bool operator != ( const MpfrClass& r1, const int r2 )
{
    return ( compare( r1,r2 ) != 0 );
}

bool operator != ( const MpfrClass& r1, const unsigned int r2 )
{
    return ( compare( r1,r2 ) != 0 );
}

bool operator != ( const MpfrClass& r1, const long int r2 )
{
    return ( compare( r1,r2 ) != 0 );
}

bool operator != ( const MpfrClass& r1, const unsigned long int r2 )
{
    return ( compare( r1,r2 ) != 0 );
}

bool operator <  ( const MpfrClass& r1, const MpfrClass& r2 )
{
    return ( compare( r1,r2 ) < 0 );
}

bool operator <  ( const double r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) >  0 );
}

bool operator <  ( const int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) >  0 );
}

bool operator <  ( const unsigned int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) >  0 );
}

bool operator <  ( const long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) >  0 );
}

bool operator <  ( const unsigned long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) >  0 );
}

bool operator <  ( const MpfrClass& r1, const double r2 )
{
    return ( compare( r1,r2 ) <  0 );
}

bool operator <  ( const MpfrClass& r1, const int r2 )
{
    return ( compare( r1,r2 ) <  0 );
}

bool operator <  ( const MpfrClass& r1, const unsigned int r2 )
{
    return ( compare( r1,r2 ) <  0 );
}

bool operator <  ( const MpfrClass& r1, const long int r2 )
{
    return ( compare( r1,r2 ) <  0 );
}

bool operator <  ( const MpfrClass& r1, const unsigned long int r2 )
{
    return ( compare( r1,r2 ) <  0 );
}

bool operator <= ( const MpfrClass& r1, const MpfrClass& r2 )
{
    return ( compare( r1,r2 ) <= 0 );
}

bool operator <= ( const double r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) >= 0 );
}

bool operator <= ( const int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) >= 0 );
}

bool operator <= ( const unsigned int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) >= 0 );
}

bool operator <= ( const long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) >= 0 );
}

bool operator <= ( const unsigned long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) >= 0 );
}

bool operator <= ( const MpfrClass& r1, const double r2 )
{
    return ( compare( r1,r2 ) <= 0 );
}

bool operator <= ( const MpfrClass& r1, const int r2 )
{
    return ( compare( r1,r2 ) <= 0 );
}

bool operator <= ( const MpfrClass& r1, const unsigned int r2 )
{
    return ( compare( r1,r2 ) <= 0 );
}

bool operator <= ( const MpfrClass& r1, const long int r2 )
{
    return ( compare( r1,r2 ) <= 0 );
}

bool operator <= ( const MpfrClass& r1, const unsigned long int r2 )
{
    return ( compare( r1,r2 ) <= 0 );
}

bool operator >  ( const MpfrClass& r1, const MpfrClass& r2 )
{
    return ( compare( r1,r2 ) > 0 );
}

bool operator >  ( const double r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) < 0 );
}

bool operator >  ( const int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) < 0 );
}

bool operator >  ( const unsigned int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) < 0 );
}

bool operator >  ( const long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) < 0 );
}

bool operator >  ( const unsigned long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) < 0 );
}

bool operator >  ( const MpfrClass& r1, const double r2 )
{
    return ( compare( r1,r2 ) >  0 );
}

bool operator >  ( const MpfrClass& r1, const int r2 )
{
    return ( compare( r1,r2 ) >  0 );
}

bool operator >  ( const MpfrClass& r1, const unsigned int r2 )
{
    return ( compare( r1,r2 ) >  0 );
}

bool operator >  ( const MpfrClass& r1, const long int r2 )
{
    return ( compare( r1,r2 ) >  0 );
}

bool operator >  ( const MpfrClass& r1, const unsigned long int r2 )
{
    return ( compare( r1,r2 ) >  0 );
}

bool operator >= ( const MpfrClass& r1, const MpfrClass& r2 )
{
    return ( compare( r1,r2 ) >= 0 );
}

bool operator >= ( const double r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) <= 0 );
}

bool operator >= ( const int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) <= 0 );
}

bool operator >= ( const unsigned int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) <= 0 );
}

bool operator >= ( const long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) <= 0 );
}

bool operator >= ( const unsigned long int r1, const MpfrClass& r2 )
{
    return ( compare( r2,r1 ) <= 0 );
}

bool operator >= ( const MpfrClass& r1, const double r2 )
{
    return ( compare( r1,r2 ) >= 0 );
}

bool operator >= ( const MpfrClass& r1, const int r2 )
{
    return ( compare( r1,r2 ) >= 0 );
}

bool operator >= ( const MpfrClass& r1, const unsigned int r2 )
{
    return ( compare( r1,r2 ) >= 0 );
}

bool operator >= ( const MpfrClass& r1, const long int r2 )
{
    return ( compare( r1,r2 ) >= 0 );
}

bool operator >= ( const MpfrClass& r1, const unsigned long int r2 )
{
    return ( compare( r1,r2 ) >= 0 );
}

// To be checked: IEEE recommendation when one operand is an exception
MpfrClass min ( const MpfrClass& a, MpfrClass& b )
{
    //MpfrClass res;
    //res.nbref->refvalue() = 1;
    //res.inexact->refvalue() = mpfr_min(res.mpfr_rep, a.mpfr_rep, b.mpfr_rep, rnd);
    //return res;

    if ( a <= b )
        return a;

    else
        return b;
}

MpfrClass max ( const MpfrClass& a, MpfrClass& b )
//const MpfrClass max (const MpfrClass& a, MpfrClass& b, RoundingMode rnd)
{
    //MpfrClass res;
    //res.nbref->refvalue() = 1;
    //res.inexact->refvalue() = mpfr_max(res.mpfr_rep, a.mpfr_rep, b.mpfr_rep, rnd);
    //return res;
    if ( a >= b ) return a;

    else return b;
}



//--------------------------------------------------------------
//
// Arithmetic operations
// Pb with operators: it is impossible to stick to MPFR philosophy
// which states that the result of an arithmetic op. is computed
// with the precision of the result, since the result is not yet known.
//
//--------------------------------------------------------------

// Addition-----------------------------------------------------
MpfrClass MpfrClass::operator + ( const MpfrClass& b ) const
{

    if ( iszero( *this ) )
    {
        b.nbref->incr();
        return b;
    }

    else if ( iszero( b ) )
    {
        nbref->incr();
        return *this;
    }

    else          // no need to check nbref counter: res has just been created
    {
        MpfrClass res;
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_add( res.mpfr_rep, mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
        return res;
    }
}

/*
MpfrClass operator + (const MpfrClass& a, const MpfrClass& b)
  {

  if ( iszero(b) )
    { a.nbref->incr(); return a; }
  else if ( iszero(a) )
    { b.nbref->incr(); return b; }
  else
    {
    MpfrClass res;
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_add(res.mpfr_rep, mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode);
    return res;
    }
  }
*/

MpfrClass operator + ( const MpfrClass& a, const double b )
{
    return a + MpfrClass( b, MpfrClass::CurrRndMode, 53 );
}

MpfrClass operator + ( const MpfrClass& a, const int b )
{
    if ( b == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
        return MpfrClass( b, MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );

    else          //use the efficient mpfr_<add,sub>_ui
    {
        // no need to check nbref counter: res has just been created
        MpfrClass res;
        res.nbref->refvalue() = 1;

        if ( b >= 0 )
            res.inexact->refvalue() = mpfr_add_ui( res.mpfr_rep, a.mpfr_rep,
                                                   ( unsigned long int )( b ), MpfrClass::CurrRndMode );

        else
            res.inexact->refvalue() = mpfr_sub_ui( res.mpfr_rep, a.mpfr_rep,
                                                   ( unsigned long int )( -b ), MpfrClass::CurrRndMode );

        return res;
    }
}

MpfrClass operator + ( const MpfrClass& a, const unsigned int b )
{
    if ( b == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
        return MpfrClass( ( unsigned long int ) b, MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );

    else          // no need to check nbref counter: res has just been created
    {
        MpfrClass res;
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_add_ui( res.mpfr_rep, a.mpfr_rep, ( unsigned long int )( b ), MpfrClass::CurrRndMode );
        return res;
    }
}

MpfrClass operator + ( const MpfrClass& a, const long int b )
{
    if ( b == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
        return MpfrClass( b, MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );

    else          //use the efficient mpfr_<add,sub>_ui
    {
        // no need to check nbref counter: res has just been created
        MpfrClass res;
        res.nbref->refvalue() = 1;

        if ( b >= 0 )
            res.inexact->refvalue() = mpfr_add_ui( res.mpfr_rep, a.mpfr_rep,
                                                   ( unsigned long int )( b ), MpfrClass::CurrRndMode );

        else
            res.inexact->refvalue() = mpfr_sub_ui( res.mpfr_rep, a.mpfr_rep,
                                                   ( unsigned long int )( -b ), MpfrClass::CurrRndMode );

        return res;
    }
}

MpfrClass operator + ( const MpfrClass& a, const mpz_srcptr b )
{
    if ( mpz_cmp_si( b, 0 ) == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
        return MpfrClass( b, MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );

    else          // no need to check nbref counter: res has just been created
    {
        MpfrClass res;
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_add_z( res.mpfr_rep, a.mpfr_rep, b, MpfrClass::CurrRndMode );
        return res;
    }
}

MpfrClass operator + ( const MpfrClass& a, const mpq_srcptr b )
{
    if ( mpq_cmp_ui( b,0,1 ) == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
        return MpfrClass( b, MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );

    else          // no need to check nbref counter: res has just been created
    {
        MpfrClass res;
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_add_q( res.mpfr_rep, a.mpfr_rep, b, MpfrClass::CurrRndMode );
        return res;
    }
}

MpfrClass operator + ( const MpfrClass& a, const unsigned long int b )
{
    if ( b == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
        return MpfrClass( b, MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );

    else          // no need to check nbref counter: res has just been created
    {
        MpfrClass res;
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_add_ui( res.mpfr_rep, a.mpfr_rep, b, MpfrClass::CurrRndMode );
        return res;
    }
}

MpfrClass operator + ( const double a, const MpfrClass& b )
{
    return b+a;
}

MpfrClass operator + ( const int a, const MpfrClass& b )
{
    return b+a;
}

MpfrClass operator + ( const unsigned int a, const MpfrClass& b )
{
    return b+a;
}

MpfrClass operator + ( const long int a, const MpfrClass& b )
{
    return b+a;
}

MpfrClass operator + ( const unsigned long int a, const MpfrClass& b )
{
    return b+a;
}

MpfrClass operator + ( const mpz_srcptr a, const MpfrClass& b )
{
    return b+a;
}

MpfrClass operator + ( const mpq_srcptr a, const MpfrClass& b )
{
    return b+a;
}

MpfrClass& MpfrClass::operator += ( const MpfrClass& b )
{
    if ( iszero( b ) )
        return *this;

    if ( iszero( *this ) )
    {
        if ( GetPrecision() == b.GetPrecision() )
        {
            if ( nbref->decr() <= 0 )		  // no other ref. on the value of *this
            {
                mpfr_clear( mpfr_rep );
                delete nbref;
                delete inexact;
            }

            mpfr_rep[0] = b.mpfr_rep[0];
            nbref = b.nbref;
            nbref->incr();
            inexact = b.inexact;
        }

        else
        {
            if ( nbref->decr() <= 0 )			// no other ref. on the value of *this
            {
                // the memory can be reused
                nbref->refvalue() = 1;
                inexact->refvalue() = mpfr_set( mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
            }

            else						// the value and memory used must be preserved
            {
                MpfrClass tmp ( *this );
                PrecisionType prec = GetPrecision();
                mpfr_init2( mpfr_rep, prec );
                nbref = new RefCounter( 1 );
                inexact = new InexactFlag();
                inexact->refvalue() = mpfr_set( mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
            }
        }

        return *this;
    }

    // none of the operands is zero
    if ( nbref->decr() <= 0 )			// no other ref. on the value of *this
    {
        // the memory can be reused
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_add( mpfr_rep, mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
    }

    else						// the value and memory used must be preserved
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_add( mpfr_rep, tmp.mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator += ( const double b )
{
    return ( *this ) += MpfrClass( b, MpfrClass::CurrRndMode, 53 );
}

MpfrClass& MpfrClass::operator += ( const long int b ) 	// using the more efficient mpfr_add_ui
{
    if ( b == 0 )
        return *this;

    if ( iszero( *this ) )
    {
        *this = MpfrClass( b, MpfrClass::CurrRndMode, GetPrecision() );
        return *this;
    }

    // none of the operands is zero
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;

        if ( b > 0 )
            inexact->refvalue() = mpfr_add_ui( mpfr_rep, mpfr_rep, ( unsigned long int )b, MpfrClass::CurrRndMode );

        else
            inexact->refvalue() = mpfr_sub_ui( mpfr_rep, mpfr_rep, ( unsigned long int )( -b ), MpfrClass::CurrRndMode );
    }

    else						// the previous value and memory must be preserved
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();

        if ( b > 0 )
            inexact->refvalue() = mpfr_add_ui( mpfr_rep, tmp.mpfr_rep, ( unsigned long int )b, MpfrClass::CurrRndMode );

        else
            inexact->refvalue() = mpfr_sub_ui( mpfr_rep, tmp.mpfr_rep, ( unsigned long int )( -b ), MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator += ( const unsigned long int b ) 	// using mpfr_add_ui
{
    if ( b == 0 )
        return *this;

    if ( iszero( *this ) )
    {
        *this = MpfrClass( b, MpfrClass::CurrRndMode, GetPrecision() );
        return *this;
    }

    // none of the operands is zero
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        inexact->refvalue() = mpfr_add_ui( mpfr_rep, mpfr_rep, ( unsigned long int )b, MpfrClass::CurrRndMode );
        nbref->refvalue() = 1;
    }

    else
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_add_ui( mpfr_rep, tmp.mpfr_rep, ( unsigned long int )b, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator += ( const int b ) 	// using the efficient mpfr_add_ui
{
    return ( *this )+= long( b );
}

MpfrClass& MpfrClass::operator += ( const unsigned int b ) // using mpfr_add_ui
{
    return ( *this )+= ( unsigned long int )( b );
}

MpfrClass& MpfrClass::operator += ( const mpz_srcptr b ) 	// using mpfr_add_z
{
    if ( mpz_cmp_si( b,0 ) == 0 )
        return *this;

    if ( iszero( *this ) )
    {
        *this = MpfrClass( b, MpfrClass::CurrRndMode, GetPrecision() );
        return *this;
    }

    // none of the operands is zero
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_add_z( mpfr_rep, mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    else
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_add_z( mpfr_rep, tmp.mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator += ( const mpq_srcptr b ) 	// using mpfr_add_q
{
    if ( mpq_cmp_si( b,0,1 ) == 0 )
        return *this;

    if ( iszero( *this ) )
    {
        *this = MpfrClass( b, MpfrClass::CurrRndMode, GetPrecision() );
        return *this;
    }

    // none of the operands is zero
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_add_q( mpfr_rep, mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    else
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_add_q( mpfr_rep, tmp.mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    return *this;
}

void MpfrClass::add ( MpfrClass& res, const MpfrClass& r1, const MpfrClass& r2, RoundingMode rnd )
// This function is the one where the philosophy of MPFR is preserved:
// in this function the result is computed with the precision of the result.
{
    if ( iszero( r1 ) )
    {
        if ( res.GetPrecision() == r2.GetPrecision() )
        {
            if ( res.nbref->decr() <= 0 )			// the memory can be reused
            {
                mpfr_clear ( res.mpfr_rep );
                delete res.nbref;
                delete res.inexact;
            }

            res.mpfr_rep[0] = r2.mpfr_rep[0];
            res.nbref = r2.nbref;
            res.nbref->incr();
            res.inexact = r2.inexact;
        }

        else 		// The result does not have the same precision as r2
        {
            // it must be converted and possibly rounded.
            if ( res.nbref->decr() <= 0 )			// the memory can be reused
            {
                res.nbref->refvalue() = 1;
                res.inexact->refvalue() = mpfr_set( res.mpfr_rep, r2.mpfr_rep, rnd );
            }

            else						// the previous memory must be preserved
            {
                PrecisionType prec = res.GetPrecision();
                mpfr_init2( res.mpfr_rep,prec );
                res.nbref = new RefCounter( 1 );
                res.inexact = new InexactFlag( 1 );
                res.inexact->refvalue() = mpfr_set( res.mpfr_rep, r2.mpfr_rep, rnd );
            }
        }
    } // r1 == 0

    if ( iszero( r2 ) )
    {
        if ( res.GetPrecision() == r1.GetPrecision() )
        {
            if ( res.nbref->decr() <= 0 )			// the memory can be reused
            {
                mpfr_clear ( res.mpfr_rep );
                delete res.nbref;
                delete res.inexact;
            }

            res.mpfr_rep[0] = r1.mpfr_rep[0];
            res.nbref = r1.nbref;
            res.nbref->incr();
            res.inexact = r1.inexact;
        }

        else 		// The result does not have the same precision as r1
        {
            // it must be converted and possibly rounded.
            if ( res.nbref->decr() <= 0 )			// the memory can be reused
            {
                res.nbref->refvalue() = 1;
                res.inexact->refvalue() = mpfr_set( res.mpfr_rep, r1.mpfr_rep, rnd );
            }

            else						// the previous memory must be preserved
            {
                PrecisionType prec = res.GetPrecision();
                mpfr_init2( res.mpfr_rep,prec );
                res.nbref = new RefCounter( 1 );
                res.inexact = new InexactFlag( 1 );
                res.inexact->refvalue() = mpfr_set( res.mpfr_rep, r1.mpfr_rep, rnd );
            }
        }
    } // r2 == 0

    // none of the operands is zero
    if ( res.nbref->decr() <= 0 )			// the memory can be reused
    {
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_add( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, MpfrClass::CurrRndMode );
    }

    else						// the previous memory must be preserved
    {
        PrecisionType prec = res.GetPrecision();
        mpfr_init2( res.mpfr_rep,prec );
        res.nbref = new RefCounter( 1 );
        res.inexact = new InexactFlag( 1 );
        res.inexact->refvalue() = mpfr_add( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, MpfrClass::CurrRndMode );
    }
}

// Subtraction --------------------------------------------------------
MpfrClass MpfrClass::operator - ( const MpfrClass& b ) const
{
    MpfrClass res;
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_sub( res.mpfr_rep, mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator - ( const MpfrClass& a, const double b )
{
    return a - MpfrClass( b, MpfrClass::CurrRndMode, 53 );
}

MpfrClass operator - ( const MpfrClass &a, const int b )
{
    if ( b == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
        return MpfrClass( ( -b ), MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );

    else          //use the efficient mpfr_<add,sub>_ui
    {
        // no need to check nbref counter: res has just been created
        MpfrClass res;
        res.nbref->refvalue() = 1;

        if ( b >= 0 )
            res.inexact->refvalue() = mpfr_sub_ui( res.mpfr_rep, a.mpfr_rep,
                                                   ( unsigned long int )( b ), MpfrClass::CurrRndMode );

        else
            res.inexact->refvalue() = mpfr_add_ui( res.mpfr_rep, a.mpfr_rep,
                                                   ( unsigned long int )( -b ), MpfrClass::CurrRndMode );

        return res;
    }
}

MpfrClass operator - ( const MpfrClass &a, const unsigned int b )
{
    if ( b == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
    {
        int tmp;
        tmp = b;
        return MpfrClass( ( -tmp ), MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );
    }

    else          //use the efficient mpfr_<add,sub>_ui
    {
        // no need to check nbref counter: res has just been created
        MpfrClass res;
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_sub_ui( res.mpfr_rep, a.mpfr_rep,
                                               ( unsigned long int )( b ), MpfrClass::CurrRndMode );
        return res;
    }
}

MpfrClass operator - ( const MpfrClass &a, const long int b )
{
    if ( b == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
        return MpfrClass( ( -b ), MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );

    else          //use the efficient mpfr_<add,sub>_ui
    {
        // no need to check nbref counter: res has just been created
        MpfrClass res;
        res.nbref->refvalue() = 1;

        if ( b >= 0 )
            res.inexact->refvalue() = mpfr_sub_ui( res.mpfr_rep, a.mpfr_rep,
                                                   ( unsigned long int )( b ), MpfrClass::CurrRndMode );

        else
            res.inexact->refvalue() = mpfr_add_ui( res.mpfr_rep, a.mpfr_rep,
                                                   ( unsigned long int )( -b ), MpfrClass::CurrRndMode );

        return res;
    }
}

MpfrClass operator - ( const MpfrClass &a, const unsigned long int b )
{
    if ( b == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
    {
        long int tmp;
        tmp = b;
        return MpfrClass( ( -tmp ), MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );
    }

    else          //use the efficient mpfr_<add,sub>_ui
    {
        // no need to check nbref counter: res has just been created
        MpfrClass res;
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_sub_ui( res.mpfr_rep, a.mpfr_rep,
                                               ( unsigned long int )( b ), MpfrClass::CurrRndMode );
        return res;
    }
}

MpfrClass operator - ( const MpfrClass& a, const mpz_srcptr b )
{
    if ( mpz_cmp_si( b, 0 ) == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
        return ( -MpfrClass( b, MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() ) );

    else          // no need to check nbref counter: res has just been created
    {
        MpfrClass res;
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_sub_z( res.mpfr_rep, a.mpfr_rep, b, MpfrClass::CurrRndMode );
        return res;
    }
}

MpfrClass operator - ( const MpfrClass& a, const mpq_srcptr b )
{
    if ( mpq_cmp_ui( b,0,1 ) == 0 )
    {
        a.nbref->incr();
        return a;
    }

    else if ( iszero( a ) )
        return ( -MpfrClass( b, MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() ) );

    else          // no need to check nbref counter: res has just been created
    {
        MpfrClass res;
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_sub_q( res.mpfr_rep, a.mpfr_rep, b, MpfrClass::CurrRndMode );
        return res;
    }
}

MpfrClass operator - ( const double a, const MpfrClass&  b )
{
    return MpfrClass( a, MpfrClass::CurrRndMode, 53 ) - b;
}

MpfrClass operator - ( const int a, const MpfrClass&  b )
{
    return ( MpfrClass( a, MpfrClass::CurrRndMode, sizeof( int ) ) - b );
}


MpfrClass operator - ( const unsigned int a, const MpfrClass&  b )
{
    if ( a == 0 )
        return ( -b );

    else if ( iszero( b ) )
    {
        return MpfrClass( ( unsigned long int )( a ), MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );
    }

    else          //use the efficient mpfr_ui_sub
    {
        // no need to check nbref counter: res has just been created
        MpfrClass res;
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_ui_sub( res.mpfr_rep, ( unsigned long int )( a ),
                                               b.mpfr_rep, MpfrClass::CurrRndMode );
        return res;
    }
}

MpfrClass operator - ( const long int a, const MpfrClass&  b )
{
    return ( MpfrClass( a, MpfrClass::CurrRndMode, sizeof( long int ) ) - b );
}

MpfrClass operator - ( const unsigned long int a, const MpfrClass&  b )
{
    if ( a == 0 )
        return ( -b );

    else if ( iszero( b ) )
    {
        return MpfrClass( a, MpfrClass::CurrRndMode, MpfrClass::GetDefaultPrecision() );
    }

    else          //use the efficient mpfr_<add,sub>_ui
    {
        // no need to check nbref counter: res has just been created
        MpfrClass res;
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_ui_sub( res.mpfr_rep, a,
                                               b.mpfr_rep, MpfrClass::CurrRndMode );
        return res;
    }
}

MpfrClass operator - ( const mpz_srcptr a, const MpfrClass&  b )
{
    return ( MpfrClass( a, MpfrClass::CurrRndMode, sizeof( int ) ) - b );
}

MpfrClass operator - ( const mpq_srcptr a, const MpfrClass&  b )
{
    return ( MpfrClass( a, MpfrClass::CurrRndMode, sizeof( int ) ) - b );
}

MpfrClass MpfrClass::operator - () const
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_neg( res.mpfr_rep, mpfr_rep, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass& MpfrClass::operator -= ( const MpfrClass& b )
{
    if ( iszero( b ) )
        return *this;

    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;

        if ( iszero( *this ) )
            inexact->refvalue() = mpfr_neg( mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );

        else
            inexact->refvalue() = mpfr_sub( mpfr_rep, mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
    }

    else						// the previous value must be preserved
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep,prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();

        if ( iszero( tmp ) )
            inexact->refvalue() = mpfr_neg( mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );

        else
            inexact->refvalue() = mpfr_sub( mpfr_rep, tmp.mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator -= ( const double b )
{
    return ( *this ) -= MpfrClass( b, MpfrClass::CurrRndMode, 53 );
}

MpfrClass& MpfrClass::operator -= ( const int b )
{
    ( *this ) += -b;
    return *this;
}

MpfrClass& MpfrClass::operator -= ( const long int b )
{
    ( *this ) += -b;
    return *this;
}

MpfrClass& MpfrClass::operator -= ( const unsigned long int b ) 	// using the more efficient mpfr_sub_ui
{
    if ( b == 0 )
        return *this;

    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_sub_ui( mpfr_rep, mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    else						// the previous value must be preserved
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_sub_ui( mpfr_rep, tmp.mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator -= ( const unsigned int b )
{
    ( *this ) -= ( unsigned long int )b;
    return *this;
}

MpfrClass& MpfrClass::operator -= ( const mpz_srcptr b ) 	// using mpfr_sub_z
{
    if ( mpz_cmp_si( b,0 ) == 0 )
        return *this;

    if ( iszero( *this ) )
    {
        *this = - MpfrClass( b, MpfrClass::CurrRndMode, GetPrecision() );
        return *this;
    }

    // none of the operands is zero
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_sub_z( mpfr_rep, mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    else
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_sub_z( mpfr_rep, tmp.mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator -= ( const mpq_srcptr b ) 	// using mpfr_sub_q
{
    if ( mpq_cmp_si( b,0,1 ) == 0 )
        return *this;

    if ( iszero( *this ) )
    {
        *this = - MpfrClass( b, MpfrClass::CurrRndMode, GetPrecision() );
        return *this;
    }

    // none of the operands is zero
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_sub_q( mpfr_rep, mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    else
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_sub_q( mpfr_rep, tmp.mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    return *this;
}

void MpfrClass::neg ( MpfrClass& res, const MpfrClass& r, RoundingMode /*rnd*/ )
{
    if ( res.nbref->decr() <= 0 )			// the memory can be reused
    {
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_neg( res.mpfr_rep, r.mpfr_rep,  MpfrClass::CurrRndMode );
    }

    else						// the previous value must be preserved
    {
        PrecisionType prec = res.GetPrecision();
        mpfr_init2( res.mpfr_rep, prec );
        res.nbref = new RefCounter( 1 );
        res.inexact = new InexactFlag();
        res.inexact->refvalue() = mpfr_neg( res.mpfr_rep, r.mpfr_rep, MpfrClass::CurrRndMode );
    }
}

void MpfrClass::sub ( MpfrClass& res, const MpfrClass& r1, const MpfrClass& r2, RoundingMode /*rnd*/ )
{
    if ( res.nbref->decr() <= 0 )
    {
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_sub( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, MpfrClass::CurrRndMode );
    }

    else						// the previous value must be preserved
    {
        PrecisionType prec = res.GetPrecision();
        mpfr_init2( res.mpfr_rep, prec );
        res.nbref = new RefCounter( 1 );
        res.inexact = new InexactFlag();
        res.inexact->refvalue() = mpfr_sub( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, MpfrClass::CurrRndMode );
    }
}

// Multiplication-----------------------------------------------
MpfrClass MpfrClass::operator * ( const MpfrClass& b ) const
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_mul( res.mpfr_rep, mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator * ( const MpfrClass& a, const double b )
{
    return a * MpfrClass( b, MpfrClass::CurrRndMode, 53 );
}

MpfrClass operator * ( const MpfrClass& a, const int b )
{
    return a * MpfrClass( b, MpfrClass::CurrRndMode, sizeof( int ) );
}

MpfrClass operator * ( const MpfrClass& a, const unsigned int b )
{
    MpfrClass res ;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_mul_ui( res.mpfr_rep, a.mpfr_rep,
                                           ( unsigned long int )( b ), MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator * ( const MpfrClass& a, const long int b )
{
    return a * MpfrClass( b, MpfrClass::CurrRndMode, sizeof( long int ) );
}

MpfrClass operator * ( const MpfrClass& a, const unsigned long int b )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_mul_ui( res.mpfr_rep, a.mpfr_rep,
                                           b, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator * ( const MpfrClass& a, const mpz_srcptr b )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_mul_z( res.mpfr_rep, a.mpfr_rep, b, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator * ( const MpfrClass& a, const mpq_srcptr b )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_mul_q( res.mpfr_rep, a.mpfr_rep, b, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator * ( const double a, const MpfrClass& b )
{
    return MpfrClass( a, MpfrClass::CurrRndMode, 53 ) * b;
}

MpfrClass operator * ( const int a, const MpfrClass& b )
{
    return MpfrClass( a, MpfrClass::CurrRndMode, sizeof( int ) ) * b;
}

MpfrClass operator * ( const unsigned int a, const MpfrClass& b )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_mul_ui( res.mpfr_rep, b.mpfr_rep,
                                           ( unsigned long int )( a ), MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator * ( const long int a, const MpfrClass& b )
{
    return MpfrClass( a, MpfrClass::CurrRndMode, sizeof( long int ) ) * b;
}

MpfrClass operator * ( const mpz_srcptr a, const MpfrClass& b )
{
    return b*a;
}

MpfrClass operator * ( const mpq_srcptr a, const MpfrClass& b )
{
    return b*a;
}

MpfrClass operator * ( const unsigned long int a, const MpfrClass& b )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_mul_ui( res.mpfr_rep, b.mpfr_rep,
                                           a, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass& MpfrClass::operator *= ( const MpfrClass& b )
{
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_mul( mpfr_rep, mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
    }

    else						// the previous value must be preserved
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_mul( mpfr_rep, tmp.mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator *= ( const double b )
{
    return ( *this ) *= MpfrClass( b, MpfrClass::CurrRndMode, 53 );
}

MpfrClass& MpfrClass::operator *= ( const long int b ) 	// using the more efficient mpfr_mul_ui
{
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;

        if ( b >= 0 )
            inexact->refvalue() = mpfr_mul_ui( mpfr_rep, mpfr_rep, ( unsigned long int )( b ), MpfrClass::CurrRndMode );

        else
        {
            int dummy;
            inexact->refvalue() = mpfr_mul_ui( mpfr_rep, mpfr_rep, ( unsigned long int )( -b ), MpfrClass::CurrRndMode );
            dummy = mpfr_neg( mpfr_rep, mpfr_rep, MpfrClass::CurrRndMode );
            inexact->refvalue() = - inexact->getvalue();
        }
    }

    else						// the previous value must be preserved
    {
        MpfrClass tmp ( *this );
        int dummy;
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();

        if ( b >= 0 )
            inexact->refvalue() = mpfr_mul_ui( mpfr_rep, tmp.mpfr_rep, ( unsigned long int )( b ), MpfrClass::CurrRndMode );

        else
        {
            inexact->refvalue() = mpfr_mul_ui( mpfr_rep, tmp.mpfr_rep, ( unsigned long int )( -b ), MpfrClass::CurrRndMode );
            dummy = mpfr_neg( mpfr_rep, mpfr_rep, MpfrClass::CurrRndMode );
            inexact->refvalue() = - inexact->getvalue();
        }
    }

    return *this;
}

MpfrClass& MpfrClass::operator *= ( const int b )
{
    ( *this )*= ( long int ) ( b );
    return *this;
}

MpfrClass& MpfrClass::operator *= ( const unsigned long int b ) 	// the efficient mpfr_mul_ui is used
{
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_mul_ui( mpfr_rep, mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    else						// the previous value must be preserved
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision ();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_mul_ui( mpfr_rep, tmp.mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator *= ( const unsigned int b )
{
    ( *this )*= ( unsigned long int ) ( b );
    return *this;
}

MpfrClass& MpfrClass::operator *= ( const mpz_srcptr b ) 	// using mpfr_mul_z
{
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_mul_z( mpfr_rep, mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    else
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_mul_z( mpfr_rep, tmp.mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator *= ( const mpq_srcptr b ) 	// using mpfr_mul_q
{
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_mul_q( mpfr_rep, mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    else
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_mul_q( mpfr_rep, tmp.mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    return *this;
}

void MpfrClass::mul ( MpfrClass& res, const MpfrClass& r1, const MpfrClass& r2, RoundingMode rnd )
{
    if ( res.nbref->decr() <= 0 )			// the memory can be reused
    {
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_mul( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, MpfrClass::CurrRndMode );
    }

    else						// the previous value must be preserved
    {
        PrecisionType prec = res.GetPrecision();
        mpfr_init2( res.mpfr_rep, prec );
        res.nbref = new RefCounter( 1 );
        res.inexact = new InexactFlag();
        res.inexact->refvalue() = mpfr_mul( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, MpfrClass::CurrRndMode );
    }
}

// Division-----------------------------------------------------
MpfrClass MpfrClass::operator / ( const MpfrClass& b ) const
{

    MpfrClass res;
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_div( res.mpfr_rep, mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator / ( const MpfrClass& a, const double b )
{
    return a / MpfrClass( b, MpfrClass::CurrRndMode, 53 );
}

MpfrClass operator / ( const MpfrClass& a, const int b )
{
    return a / MpfrClass( b, MpfrClass::CurrRndMode, sizeof( int ) );
}

MpfrClass operator / ( const MpfrClass& a, const unsigned int b )
{
    MpfrClass res;
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_div_ui( res.mpfr_rep, a.mpfr_rep,
                                           ( unsigned long int )( b ), MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator / ( const MpfrClass& a, const long int b )
{
    return a / MpfrClass( b, MpfrClass::CurrRndMode, sizeof( long int ) );
}

MpfrClass operator / ( const MpfrClass& a, const unsigned long int b )
{
    MpfrClass res;
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_div_ui( res.mpfr_rep, a.mpfr_rep,
                                           b, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator / ( const MpfrClass& a, const mpz_srcptr b )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_div_z( res.mpfr_rep, a.mpfr_rep, b, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator / ( const MpfrClass& a, const mpq_srcptr b )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_div_q( res.mpfr_rep, a.mpfr_rep, b, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator / ( const double a, const MpfrClass& b )
{
    return MpfrClass( a, MpfrClass::CurrRndMode, 53 )/b;
}

MpfrClass operator / ( const long int a, const MpfrClass& b )
{
    return MpfrClass( a, MpfrClass::CurrRndMode, sizeof( long int ) ) / b;
}

MpfrClass operator / ( const int a, const MpfrClass& b )
{
    return ( long ) a / b;
}

MpfrClass operator / ( const unsigned long int a, const MpfrClass& b )	// using mpfr_ui_div
{
    MpfrClass res;
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_ui_div( res.mpfr_rep, a, b.mpfr_rep, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator / ( const unsigned int a, const MpfrClass& b )
{
    MpfrClass res;
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_ui_div( res.mpfr_rep, ( unsigned long int )( a ),
                                           b.mpfr_rep, MpfrClass::CurrRndMode );
    return res;
}

MpfrClass operator / ( const mpz_srcptr a, const MpfrClass& b )
{
    return MpfrClass( a ) / b;
}

MpfrClass operator / ( const mpq_srcptr a, const MpfrClass& b )
{
    return MpfrClass( a ) / b;
}

MpfrClass& MpfrClass::operator /= ( const MpfrClass& b )
{
    if ( nbref->decr() <= 0 )				// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_div( mpfr_rep, mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
    }

    else							// the previous value must be preserved
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        mpfr_div( mpfr_rep, tmp.mpfr_rep, b.mpfr_rep, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator /= ( const double b )
{
    return ( *this ) /= MpfrClass( b, MpfrClass::CurrRndMode, GetPrecision() );
}

MpfrClass& MpfrClass::operator /= ( const long int b ) 	// using the more efficient mpfr_div_ui
{
    if ( nbref->decr() <= 0 )				// the memory can be reused
    {
        nbref->refvalue() = 1;

        if ( b > 0 )
            inexact->refvalue() = mpfr_div_ui( mpfr_rep, mpfr_rep, ( unsigned long int )b, MpfrClass::CurrRndMode );

        else
        {
            if ( MpfrClass::CurrRndMode == RoundDown )
                inexact->refvalue() = mpfr_div_ui( mpfr_rep, mpfr_rep, ( unsigned long int )( -b ), RoundUp );

            else if ( MpfrClass::CurrRndMode == RoundUp )
                inexact->refvalue() = mpfr_div_ui( mpfr_rep, mpfr_rep, ( unsigned long int )( -b ), RoundDown );

            else
                inexact->refvalue() = mpfr_div_ui( mpfr_rep, mpfr_rep, ( unsigned long int )( -b ), MpfrClass::CurrRndMode );

            int dummy = mpfr_neg( mpfr_rep, mpfr_rep, MpfrClass::CurrRndMode );
            inexact->refvalue() = -inexact->getvalue();
        }
    }

    else							// the previous value must be preserved
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();

        if ( b > 0 )
            inexact->refvalue() = mpfr_div_ui( mpfr_rep, tmp.mpfr_rep, ( unsigned long int )b, MpfrClass::CurrRndMode );

        else
        {
            if ( MpfrClass::CurrRndMode == RoundDown )
                inexact->refvalue() = mpfr_div_ui( mpfr_rep, mpfr_rep, ( unsigned long int )( -b ), RoundUp );

            else if ( MpfrClass::CurrRndMode == RoundUp )
                inexact->refvalue() = mpfr_div_ui( mpfr_rep, mpfr_rep, ( unsigned long int )( -b ), RoundDown );

            else
                inexact->refvalue() = mpfr_div_ui( mpfr_rep, mpfr_rep, ( unsigned long int )( -b ), MpfrClass::CurrRndMode );

            int dummy = mpfr_neg( mpfr_rep, mpfr_rep, MpfrClass::CurrRndMode );
            inexact->refvalue() = -inexact->getvalue();
        }
    }

    return *this;
}

MpfrClass& MpfrClass::operator /= ( const int b )
{
    ( *this ) /= ( long int ) b;
    return *this;
}

MpfrClass& MpfrClass::operator /= ( const unsigned long int b ) 	// using mpfr_div_ui
{
    if ( nbref->decr() <= 0 )				// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_div_ui( mpfr_rep, mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    else							// the previous value must be preserved
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact->refvalue() = mpfr_div_ui( mpfr_rep, tmp.mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator /= ( const unsigned int b )
{
    ( *this ) /= ( unsigned long int ) b;
    return *this;
}

MpfrClass& MpfrClass::operator /= ( const mpz_srcptr b ) 	// using mpfr_div_z
{
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_div_z( mpfr_rep, mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    else
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_div_z( mpfr_rep, tmp.mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    return *this;
}

MpfrClass& MpfrClass::operator /= ( const mpq_srcptr b ) 	// using mpfr_div_q
{
    if ( nbref->decr() <= 0 )			// the memory can be reused
    {
        nbref->refvalue() = 1;
        inexact->refvalue() = mpfr_div_q( mpfr_rep, mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    else
    {
        MpfrClass tmp ( *this );
        PrecisionType prec = GetPrecision();
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
        inexact->refvalue() = mpfr_div_q( mpfr_rep, tmp.mpfr_rep, b, MpfrClass::CurrRndMode );
    }

    return *this;
}

void MpfrClass::div ( MpfrClass& res, const MpfrClass& r1, const MpfrClass& r2, RoundingMode /*rnd*/ )
{
    if ( res.nbref->decr() <= 0 )				// the memory can be reused
    {
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_div( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, MpfrClass::CurrRndMode );
    }

    else							// the previous value must be preserved
    {
        PrecisionType prec = res.GetPrecision();
        mpfr_init2( res.mpfr_rep, prec );
        res.nbref = new RefCounter( 1 );
        res.inexact = new InexactFlag();
        res.inexact->refvalue() = mpfr_div( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, MpfrClass::CurrRndMode );
    }
}

void MpfrClass::fma ( MpfrClass& res, const MpfrClass& r1, const MpfrClass& r2, const MpfrClass& r3, RoundingMode /*rnd*/ )
{
    if ( res.nbref->decr() <= 0 )				// the memory can be reused
    {
        res.nbref->refvalue() = 1;
        res.inexact->refvalue() = mpfr_fma( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, r3.mpfr_rep, MpfrClass::CurrRndMode );
    }

    else							// the previous value must be preserved
    {
        PrecisionType prec = res.GetPrecision();
        mpfr_init2( res.mpfr_rep, prec );
        res.nbref = new RefCounter( 1 );
        res.inexact = new InexactFlag();
        res.inexact->refvalue() = mpfr_fma( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, r3.mpfr_rep, MpfrClass::CurrRndMode );
    }
}


//--------------------------------------------------------------
//
// Mathematical and miscellaneous functions
//
//--------------------------------------------------------------
// NaN of infinity handled by MPFR
// => no test on the value (sign...) of the operand
void MpfrClass::random ( PrecisionType prec )		// member to avoid conflict with GMP random
{
    if ( nbref->decr() <= 0 )				// the memory can be reused
    {
        nbref->refvalue() = 1;
        mpfr_set_prec( mpfr_rep, prec );
    }

    else							// the previous value must be preserved
    {
        mpfr_init2( mpfr_rep, prec );
        nbref = new RefCounter( 1 );
        inexact = new InexactFlag();
    }

    mpfr_random( mpfr_rep );
    inexact->refvalue() = EXACT_FLAG;
}


MpfrClass abs ( const MpfrClass& r, RoundingMode /*rnd*/ )
{
    if ( r >= 0 )
        return r+0;	// changes the precision of r

    // if CurrPrecision != prec of r
    else
        return -r;
}

MpfrClass agm ( const MpfrClass& r1, const MpfrClass& r2, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_agm( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, rnd );
    return res;
}

MpfrClass sqrt ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    // NaN handled by MPFR in case r < 0
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_sqrt( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass exp ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_exp( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass expm1 ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_expm1( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}




MpfrClass log ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    // NaN handled by MPFR in case r < 0
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_log( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass log2 ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    // NaN handled by MPFR in case r < 0
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_log2( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass log10 ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    // NaN handled by MPFR in case r < 0
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_log10( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass log1p ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    // NaN handled by MPFR in case r < 0
    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_log1p( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass sin ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_sin( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass cos ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_cos( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

void sin_cos ( MpfrClass& res_sin, MpfrClass& res_cos, const MpfrClass& r, RoundingMode rnd )
{
    int inexact;

    if ( res_sin.nbref->decr() <= 0 )
    {
        res_sin.nbref->refvalue() = 1;
    }

    else
    {
        MpfrClass::PrecisionType prec = res_sin.GetPrecision();
        mpfr_init2( res_sin.mpfr_rep, prec );
        res_sin.nbref = new RefCounter( 1 );
    }

    if ( res_cos.nbref->decr() <= 0 )
    {
        res_cos.nbref->refvalue() = 1;
    }

    else
    {
        MpfrClass::PrecisionType prec = res_cos.GetPrecision();
        mpfr_init2( res_cos.mpfr_rep, prec );
        res_cos.nbref = new RefCounter( 1 );
    }

    inexact = mpfr_sin_cos( res_sin.mpfr_rep, res_cos.mpfr_rep, r.mpfr_rep, rnd );

    if ( inexact == 0 )
    {
        res_sin.inexact->refvalue() = EXACT_FLAG;
        res_cos.inexact->refvalue() = EXACT_FLAG;
    }

    else
    {
        res_sin.inexact->refvalue() = UNAFFECTED_INEXACT_FLAG;
        res_cos.inexact->refvalue() = UNAFFECTED_INEXACT_FLAG;
    }
}


// new functions 10/04/03

MpfrClass tan ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_tan( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}


MpfrClass acos ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_acos( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}


MpfrClass asin ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_asin( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}


MpfrClass atan ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_atan( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass cosh ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_cosh( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass sinh ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_sinh( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}


MpfrClass tanh ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_tanh( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}


MpfrClass asinh ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_asinh( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}


MpfrClass acosh ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_acosh( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass atanh ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_atanh( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass MpfrClass::pow ( const unsigned long int e, RoundingMode rnd ) const
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_pow_ui( res.mpfr_rep, mpfr_rep, e, rnd );
    return res;
}

MpfrClass MpfrClass::pow ( const long int e, RoundingMode rnd ) const
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_pow_si( res.mpfr_rep, mpfr_rep, e, rnd );
    return res;
}

MpfrClass MpfrClass::pow ( const MpfrClass& e, RoundingMode rnd ) const
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_pow( res.mpfr_rep, mpfr_rep, e.mpfr_rep, rnd );
    return res;
}


MpfrClass cbrt ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_cbrt( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass exp2 ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_exp2( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass gamma ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_gamma( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass erf ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_erf( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

MpfrClass factorial ( const unsigned long int n, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_fac_ui( res.mpfr_rep, n, rnd );
    return res;
}

MpfrClass hypot ( const MpfrClass& r1, const MpfrClass& r2, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_hypot( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, rnd );
    return res;
}


MpfrClass zeta ( const MpfrClass& r, RoundingMode rnd )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_zeta( res.mpfr_rep, r.mpfr_rep, rnd );
    return res;
}

// end new functions
//
//

MpfrClass round ( const MpfrClass& r )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_round( res.mpfr_rep, r.mpfr_rep );
    return res;
}

MpfrClass floor ( const MpfrClass& r )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_floor( res.mpfr_rep, r.mpfr_rep );
    return res;
}

MpfrClass trunc ( const MpfrClass& r )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_trunc( res.mpfr_rep, r.mpfr_rep );
    return res;
}

MpfrClass ceil ( const MpfrClass& r )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_ceil( res.mpfr_rep, r.mpfr_rep );
    return res;
}

MpfrClass frac ( const MpfrClass& r )
{
    MpfrClass res;

    res.nbref->refvalue() = 1;
    res.inexact->refvalue() = mpfr_frac( res.mpfr_rep, r.mpfr_rep, MpfrClass::CurrRndMode );
    return res;
}

long int to_int ( const MpfrClass& r, RoundingMode rnd )
{
    return mpfr_get_si( r.mpfr_rep,rnd );
}

unsigned long int to_uint ( const MpfrClass& r, RoundingMode rnd )
{
    return mpfr_get_ui( r.mpfr_rep,rnd );
}

double to_double ( const MpfrClass& r, RoundingMode rnd )
{
    return mpfr_get_d( r.mpfr_rep,rnd );
}

long double to_ldouble ( const MpfrClass& r, RoundingMode rnd )
{
    return mpfr_get_ld( r.mpfr_rep,rnd );
}

MpfrClass reldiff ( const MpfrClass& r1, const MpfrClass& r2, RoundingMode rnd )
{
    MpfrClass res;

    mpfr_reldiff( res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, rnd );
    return res;
}

MpfrClass nextabove ( const MpfrClass& r )
{
    MpfrClass res;

    res.copy( r );
    mpfr_nextabove( res.mpfr_rep );
    return res;
}

MpfrClass nextbelow ( const MpfrClass& r )
{
    MpfrClass res;

    res.copy( r );
    mpfr_nextbelow( res.mpfr_rep );
    return res;
}

MpfrClass nexttoward ( const MpfrClass& r, const MpfrClass& dir )
{
    MpfrClass res;

    res.copy( r );
    mpfr_nexttoward( res.mpfr_rep, dir.mpfr_rep );
    return res;
}

}
}
