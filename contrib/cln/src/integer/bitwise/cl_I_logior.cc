// logior().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "integer/bitwise/cl_I_log.h"

namespace cln {

// Logische Operationen auf Integers:
// Methode: aus den Längen der beiden Argumente eine obere Schranke für
// die Länge des Ergebnisses berechnen (das Maximum der beiden Längen und
// FN_maxlength), so daß das MSD für unendlich viele Bits steht.
// Dann beide Argumente in gleichgroße Digit sequences umwandeln, Operation
// mit einer einfachen Schleife durchführen.

const cl_I logior (const cl_I& x, const cl_I& y)
    { if (fixnump(x) && fixnump(y)) // Beides Fixnums -> ganz einfach:
        { // bitweise als Fixnum zurück
          return cl_I_from_word(x.word | y.word);
        }
        else
        { CL_ALLOCA_STACK;
          var uintC n; // Anzahl der Digits
         {var uintC nx = I_to_DS_need(x);
          var uintC ny = I_to_DS_need(y);
          n = (nx>=ny ? nx : ny);
         }
         {var uintD* xptr; I_to_DS_n(x,n,xptr=); // Pointer in DS zu x
          var uintD* yptr; I_to_DS_n(y,n,yptr=); // Pointer in DS zu y
          var uintD* zptr = xptr; // Pointer aufs Ergebnis
          or_loop_msp(xptr,yptr,n); // mit OR verknüpfen
          return DS_to_I(zptr,n); // Ergebnis als Integer
    }   }}

}  // namespace cln
