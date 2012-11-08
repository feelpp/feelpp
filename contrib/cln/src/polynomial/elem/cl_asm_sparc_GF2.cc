// Externe Routinen
// Prozessor: SPARC
// Compiler: GNU-C oder SUN-C
// Parameter-Übergabe: in Registern %o0-%o5.
// Einstellungen: intCsize=32, intDsize=32.

#ifdef ASM_UNDERSCORE /* SunOS 4 */
  #if defined(__STDC__) || defined (__cplusplus)
    #define C(entrypoint) _##entrypoint
  #else
    #define C(entrypoint) _/**/entrypoint
  #endif
#else /* SunOS 5 = Solaris 2 */
  #define C(entrypoint) entrypoint
#endif

// When this file is compiled into a shared library, ELF linkers need to
// know which symbols are functions.
#if defined(__NetBSD__) || defined(__OpenBSD__)
  #define DECLARE_FUNCTION(name) .type C(name),@function
#elif defined(__svr4__) || defined(__ELF__)
  #define DECLARE_FUNCTION(name) .type C(name),#function
#else
  #define DECLARE_FUNCTION(name)
#endif

  // Indikatoren für Anweisungen (Instruktionen) in Delay-Slots
  // (diese werden VOR der vorigen Instruktion ausgeführt):
  #define _             // Instruktion, die stets ausgeführt wird
  #define __            // Instruktion, die nur im Sprung-Fall ausgeführt wird
  // Abkürzungen für Anweisungen:
  #define ret   jmp %i7+8    // return from subroutine
  #define retl  jmp %o7+8    // return from leaf subroutine (no save/restore)

        .seg "text"

        .global C(gf2_mul16),C(gf2_mul32)

// extern uint32 gf2_mul16 (uint16 x, uint16 y);
        DECLARE_FUNCTION(gf2_mul16)
C(gf2_mul16:) // Input in %o0,%o1, Output in %o0
        sll %o0,16,%o0
        sll %o1,16,%o1
        srl %o1,16,%o1
        // 16-bit multiply of x and y
        // input %o1 = factor1, %o0 = 2^16*factor2, output %o0
        addcc %o0,%o0,%o0
        bcs Lb01
       _ addcc %o0,%o0,%o0
La01:   bcs Lb02
       _ addcc %o0,%o0,%o0
La02:   bcs Lb03
       _ addcc %o0,%o0,%o0
La03:   bcs Lb04
       _ addcc %o0,%o0,%o0
La04:   bcs Lb05
       _ addcc %o0,%o0,%o0
La05:   bcs Lb06
       _ addcc %o0,%o0,%o0
La06:   bcs Lb07
       _ addcc %o0,%o0,%o0
La07:   bcs Lb08
       _ addcc %o0,%o0,%o0
La08:   bcs Lb09
       _ addcc %o0,%o0,%o0
La09:   bcs Lb10
       _ addcc %o0,%o0,%o0
La10:   bcs Lb11
       _ addcc %o0,%o0,%o0
La11:   bcs Lb12
       _ addcc %o0,%o0,%o0
La12:   bcs Lb13
       _ addcc %o0,%o0,%o0
La13:   bcs Lb14
       _ addcc %o0,%o0,%o0
La14:   bcs Lb15
       _ addcc %o0,%o0,%o0
La15:   bcs Lb16
       _ add %o0,%o0,%o0
La16:   retl
       _ nop
Lb01:   xor %o0,%o1,%o0
        bcc La02
       _ addcc %o0,%o0,%o0
Lb02:   xor %o0,%o1,%o0
        bcc La03
       _ addcc %o0,%o0,%o0
Lb03:   xor %o0,%o1,%o0
        bcc La04
       _ addcc %o0,%o0,%o0
Lb04:   xor %o0,%o1,%o0
        bcc La05
       _ addcc %o0,%o0,%o0
Lb05:   xor %o0,%o1,%o0
        bcc La06
       _ addcc %o0,%o0,%o0
Lb06:   xor %o0,%o1,%o0
        bcc La07
       _ addcc %o0,%o0,%o0
Lb07:   xor %o0,%o1,%o0
        bcc La08
       _ addcc %o0,%o0,%o0
Lb08:   xor %o0,%o1,%o0
        bcc La09
       _ addcc %o0,%o0,%o0
Lb09:   xor %o0,%o1,%o0
        bcc La10
       _ addcc %o0,%o0,%o0
Lb10:   xor %o0,%o1,%o0
        bcc La11
       _ addcc %o0,%o0,%o0
Lb11:   xor %o0,%o1,%o0
        bcc La12
       _ addcc %o0,%o0,%o0
Lb12:   xor %o0,%o1,%o0
        bcc La13
       _ addcc %o0,%o0,%o0
Lb13:   xor %o0,%o1,%o0
        bcc La14
       _ addcc %o0,%o0,%o0
Lb14:   xor %o0,%o1,%o0
        bcc La15
       _ addcc %o0,%o0,%o0
Lb15:   xor %o0,%o1,%o0
        bcc La16
       _ add %o0,%o0,%o0
Lb16:   retl
       _ xor %o0,%o1,%o0

// extern uint32 gf2_mul32 (uint32 x, uint32 y, uint32* plo);
        DECLARE_FUNCTION(gf2_mul32)
C(gf2_mul32:) // Input in %o0,%o1,%o2, Output in [%o2],%o0
#if 0
        sll %o0,16,%o4
        srl %o4,16,%o4          // %o4 = low16(x)
        sll %o1,16,%o5          // %o5 = 2^16*low16(y)
        srl %o0,16,%o0          // %o0 = high16(x)
        srl %o1,16,%o1
        sll %o1,16,%o1          // %o1 = 2^16*high16(y)
        xor %o1,%o5,%o3         // %o3 = 2^16*(high16(y)+low16(y))
        // 16-bit multiply of low16(x) and low16(y)
        // input %o4 = factor1, %o5 = 2^16*factor2, output %o5
        addcc %o5,%o5,%o5
        bcs Ld01
       _ addcc %o5,%o5,%o5
Lc01:   bcs Ld02
       _ addcc %o5,%o5,%o5
Lc02:   bcs Ld03
       _ addcc %o5,%o5,%o5
Lc03:   bcs Ld04
       _ addcc %o5,%o5,%o5
Lc04:   bcs Ld05
       _ addcc %o5,%o5,%o5
Lc05:   bcs Ld06
       _ addcc %o5,%o5,%o5
Lc06:   bcs Ld07
       _ addcc %o5,%o5,%o5
Lc07:   bcs Ld08
       _ addcc %o5,%o5,%o5
Lc08:   bcs Ld09
       _ addcc %o5,%o5,%o5
Lc09:   bcs Ld10
       _ addcc %o5,%o5,%o5
Lc10:   bcs Ld11
       _ addcc %o5,%o5,%o5
Lc11:   bcs Ld12
       _ addcc %o5,%o5,%o5
Lc12:   bcs Ld13
       _ addcc %o5,%o5,%o5
Lc13:   bcs Ld14
       _ addcc %o5,%o5,%o5
Lc14:   bcs Ld15
       _ addcc %o5,%o5,%o5
Lc15:   bcs Ld16
       _ add %o5,%o5,%o5
Lc16:   b Ld17
       _ nop
Ld01:   xor %o5,%o4,%o5
        bcc Lc02
       _ addcc %o5,%o5,%o5
Ld02:   xor %o5,%o4,%o5
        bcc Lc03
       _ addcc %o5,%o5,%o5
Ld03:   xor %o5,%o4,%o5
        bcc Lc04
       _ addcc %o5,%o5,%o5
Ld04:   xor %o5,%o4,%o5
        bcc Lc05
       _ addcc %o5,%o5,%o5
Ld05:   xor %o5,%o4,%o5
        bcc Lc06
       _ addcc %o5,%o5,%o5
Ld06:   xor %o5,%o4,%o5
        bcc Lc07
       _ addcc %o5,%o5,%o5
Ld07:   xor %o5,%o4,%o5
        bcc Lc08
       _ addcc %o5,%o5,%o5
Ld08:   xor %o5,%o4,%o5
        bcc Lc09
       _ addcc %o5,%o5,%o5
Ld09:   xor %o5,%o4,%o5
        bcc Lc10
       _ addcc %o5,%o5,%o5
Ld10:   xor %o5,%o4,%o5
        bcc Lc11
       _ addcc %o5,%o5,%o5
Ld11:   xor %o5,%o4,%o5
        bcc Lc12
       _ addcc %o5,%o5,%o5
Ld12:   xor %o5,%o4,%o5
        bcc Lc13
       _ addcc %o5,%o5,%o5
Ld13:   xor %o5,%o4,%o5
        bcc Lc14
       _ addcc %o5,%o5,%o5
Ld14:   xor %o5,%o4,%o5
        bcc Lc15
       _ addcc %o5,%o5,%o5
Ld15:   xor %o5,%o4,%o5
        bcc Ld17
       _ add %o5,%o5,%o5
Ld16:   xor %o5,%o4,%o5
Ld17:                           // %o5 = low16(x)*low16(y)
        // 16-bit multiply of high16(x) and high16(y)
        // input %o0 = factor1, %o1 = 2^16*factor2, output %o1
        addcc %o1,%o1,%o1
        bcs Lf01
       _ addcc %o1,%o1,%o1
Le01:   bcs Lf02
       _ addcc %o1,%o1,%o1
Le02:   bcs Lf03
       _ addcc %o1,%o1,%o1
Le03:   bcs Lf04
       _ addcc %o1,%o1,%o1
Le04:   bcs Lf05
       _ addcc %o1,%o1,%o1
Le05:   bcs Lf06
       _ addcc %o1,%o1,%o1
Le06:   bcs Lf07
       _ addcc %o1,%o1,%o1
Le07:   bcs Lf08
       _ addcc %o1,%o1,%o1
Le08:   bcs Lf09
       _ addcc %o1,%o1,%o1
Le09:   bcs Lf10
       _ addcc %o1,%o1,%o1
Le10:   bcs Lf11
       _ addcc %o1,%o1,%o1
Le11:   bcs Lf12
       _ addcc %o1,%o1,%o1
Le12:   bcs Lf13
       _ addcc %o1,%o1,%o1
Le13:   bcs Lf14
       _ addcc %o1,%o1,%o1
Le14:   bcs Lf15
       _ addcc %o1,%o1,%o1
Le15:   bcs Lf16
       _ add %o1,%o1,%o1
Le16:   b Lf17
       _ nop
Lf01:   xor %o1,%o0,%o1
        bcc Le02
       _ addcc %o1,%o1,%o1
Lf02:   xor %o1,%o0,%o1
        bcc Le03
       _ addcc %o1,%o1,%o1
Lf03:   xor %o1,%o0,%o1
        bcc Le04
       _ addcc %o1,%o1,%o1
Lf04:   xor %o1,%o0,%o1
        bcc Le05
       _ addcc %o1,%o1,%o1
Lf05:   xor %o1,%o0,%o1
        bcc Le06
       _ addcc %o1,%o1,%o1
Lf06:   xor %o1,%o0,%o1
        bcc Le07
       _ addcc %o1,%o1,%o1
Lf07:   xor %o1,%o0,%o1
        bcc Le08
       _ addcc %o1,%o1,%o1
Lf08:   xor %o1,%o0,%o1
        bcc Le09
       _ addcc %o1,%o1,%o1
Lf09:   xor %o1,%o0,%o1
        bcc Le10
       _ addcc %o1,%o1,%o1
Lf10:   xor %o1,%o0,%o1
        bcc Le11
       _ addcc %o1,%o1,%o1
Lf11:   xor %o1,%o0,%o1
        bcc Le12
       _ addcc %o1,%o1,%o1
Lf12:   xor %o1,%o0,%o1
        bcc Le13
       _ addcc %o1,%o1,%o1
Lf13:   xor %o1,%o0,%o1
        bcc Le14
       _ addcc %o1,%o1,%o1
Lf14:   xor %o1,%o0,%o1
        bcc Le15
       _ addcc %o1,%o1,%o1
Lf15:   xor %o1,%o0,%o1
        bcc Lf17
       _ add %o1,%o1,%o1
Lf16:   xor %o1,%o0,%o1
Lf17:                           // %o1 = high16(x)*high16(y)
        xor %o0,%o4,%o4         // %o4 = high16(x)+low16(x)
        // 16-bit multiply of high16(x)+low16(x) and high16(y)+low16(y)
        // input %o4 = factor1, %o3 = 2^16*factor2, output %o3
        addcc %o3,%o3,%o3
        bcs Lh01
       _ addcc %o3,%o3,%o3
Lg01:   bcs Lh02
       _ addcc %o3,%o3,%o3
Lg02:   bcs Lh03
       _ addcc %o3,%o3,%o3
Lg03:   bcs Lh04
       _ addcc %o3,%o3,%o3
Lg04:   bcs Lh05
       _ addcc %o3,%o3,%o3
Lg05:   bcs Lh06
       _ addcc %o3,%o3,%o3
Lg06:   bcs Lh07
       _ addcc %o3,%o3,%o3
Lg07:   bcs Lh08
       _ addcc %o3,%o3,%o3
Lg08:   bcs Lh09
       _ addcc %o3,%o3,%o3
Lg09:   bcs Lh10
       _ addcc %o3,%o3,%o3
Lg10:   bcs Lh11
       _ addcc %o3,%o3,%o3
Lg11:   bcs Lh12
       _ addcc %o3,%o3,%o3
Lg12:   bcs Lh13
       _ addcc %o3,%o3,%o3
Lg13:   bcs Lh14
       _ addcc %o3,%o3,%o3
Lg14:   bcs Lh15
       _ addcc %o3,%o3,%o3
Lg15:   bcs Lh16
       _ add %o3,%o3,%o3
Lg16:   b Lh17
       _ nop
Lh01:   xor %o3,%o4,%o3
        bcc Lg02
       _ addcc %o3,%o3,%o3
Lh02:   xor %o3,%o4,%o3
        bcc Lg03
       _ addcc %o3,%o3,%o3
Lh03:   xor %o3,%o4,%o3
        bcc Lg04
       _ addcc %o3,%o3,%o3
Lh04:   xor %o3,%o4,%o3
        bcc Lg05
       _ addcc %o3,%o3,%o3
Lh05:   xor %o3,%o4,%o3
        bcc Lg06
       _ addcc %o3,%o3,%o3
Lh06:   xor %o3,%o4,%o3
        bcc Lg07
       _ addcc %o3,%o3,%o3
Lh07:   xor %o3,%o4,%o3
        bcc Lg08
       _ addcc %o3,%o3,%o3
Lh08:   xor %o3,%o4,%o3
        bcc Lg09
       _ addcc %o3,%o3,%o3
Lh09:   xor %o3,%o4,%o3
        bcc Lg10
       _ addcc %o3,%o3,%o3
Lh10:   xor %o3,%o4,%o3
        bcc Lg11
       _ addcc %o3,%o3,%o3
Lh11:   xor %o3,%o4,%o3
        bcc Lg12
       _ addcc %o3,%o3,%o3
Lh12:   xor %o3,%o4,%o3
        bcc Lg13
       _ addcc %o3,%o3,%o3
Lh13:   xor %o3,%o4,%o3
        bcc Lg14
       _ addcc %o3,%o3,%o3
Lh14:   xor %o3,%o4,%o3
        bcc Lg15
       _ addcc %o3,%o3,%o3
Lh15:   xor %o3,%o4,%o3
        bcc Lh17
       _ add %o3,%o3,%o3
Lh16:   xor %o3,%o4,%o3
Lh17:                           // %o3 = (high16(x)+low16(x))*(high16(y)+low16(y))
        // Now %o5 = low16(x)*low16(y)
        //     %o1 = high16(x)*high16(y)
        //     %o3 = (high16(x)+low16(x))*(high16(y)+low16(y))
        // The result is   x*y = 2^32*%o1 + 2^16*(%o3+%o1+%o5) + %o5
        xor %o3,%o1,%o3
        xor %o3,%o5,%o3
        // The result is   x*y = 2^32*%o1 + 2^16*%o3 + %o5
        srl %o3,16,%o0
        xor %o0,%o1,%o0         // high 32 bits in %o0
        sll %o3,16,%o1
        xor %o1,%o5,%o1         // low 32 bits in %o1
        retl
       _ st %o1,[%o2]
#else
        mov 0,%o3
        // 32-bit multiply of x and y
        // input %o1 = factor1, %o0|%o3 = 2^32*factor2, output %o0|%o3
        addcc %o0,%o0,%o0
        bcs Ld01
       _ addcc %o3,%o3,%o3
Lc01:   addxcc %o0,%o0,%o0
        bcs Ld02
       _ addcc %o3,%o3,%o3
Lc02:   addxcc %o0,%o0,%o0
        bcs Ld03
       _ addcc %o3,%o3,%o3
Lc03:   addxcc %o0,%o0,%o0
        bcs Ld04
       _ addcc %o3,%o3,%o3
Lc04:   addxcc %o0,%o0,%o0
        bcs Ld05
       _ addcc %o3,%o3,%o3
Lc05:   addxcc %o0,%o0,%o0
        bcs Ld06
       _ addcc %o3,%o3,%o3
Lc06:   addxcc %o0,%o0,%o0
        bcs Ld07
       _ addcc %o3,%o3,%o3
Lc07:   addxcc %o0,%o0,%o0
        bcs Ld08
       _ addcc %o3,%o3,%o3
Lc08:   addxcc %o0,%o0,%o0
        bcs Ld09
       _ addcc %o3,%o3,%o3
Lc09:   addxcc %o0,%o0,%o0
        bcs Ld10
       _ addcc %o3,%o3,%o3
Lc10:   addxcc %o0,%o0,%o0
        bcs Ld11
       _ addcc %o3,%o3,%o3
Lc11:   addxcc %o0,%o0,%o0
        bcs Ld12
       _ addcc %o3,%o3,%o3
Lc12:   addxcc %o0,%o0,%o0
        bcs Ld13
       _ addcc %o3,%o3,%o3
Lc13:   addxcc %o0,%o0,%o0
        bcs Ld14
       _ addcc %o3,%o3,%o3
Lc14:   addxcc %o0,%o0,%o0
        bcs Ld15
       _ addcc %o3,%o3,%o3
Lc15:   addxcc %o0,%o0,%o0
        bcs Ld16
       _ addcc %o3,%o3,%o3
Lc16:   addxcc %o0,%o0,%o0
        bcs Ld17
       _ addcc %o3,%o3,%o3
Lc17:   addxcc %o0,%o0,%o0
        bcs Ld18
       _ addcc %o3,%o3,%o3
Lc18:   addxcc %o0,%o0,%o0
        bcs Ld19
       _ addcc %o3,%o3,%o3
Lc19:   addxcc %o0,%o0,%o0
        bcs Ld20
       _ addcc %o3,%o3,%o3
Lc20:   addxcc %o0,%o0,%o0
        bcs Ld21
       _ addcc %o3,%o3,%o3
Lc21:   addxcc %o0,%o0,%o0
        bcs Ld22
       _ addcc %o3,%o3,%o3
Lc22:   addxcc %o0,%o0,%o0
        bcs Ld23
       _ addcc %o3,%o3,%o3
Lc23:   addxcc %o0,%o0,%o0
        bcs Ld24
       _ addcc %o3,%o3,%o3
Lc24:   addxcc %o0,%o0,%o0
        bcs Ld25
       _ addcc %o3,%o3,%o3
Lc25:   addxcc %o0,%o0,%o0
        bcs Ld26
       _ addcc %o3,%o3,%o3
Lc26:   addxcc %o0,%o0,%o0
        bcs Ld27
       _ addcc %o3,%o3,%o3
Lc27:   addxcc %o0,%o0,%o0
        bcs Ld28
       _ addcc %o3,%o3,%o3
Lc28:   addxcc %o0,%o0,%o0
        bcs Ld29
       _ addcc %o3,%o3,%o3
Lc29:   addxcc %o0,%o0,%o0
        bcs Ld30
       _ addcc %o3,%o3,%o3
Lc30:   addxcc %o0,%o0,%o0
        bcs Ld31
       _ addcc %o3,%o3,%o3
Lc31:   addxcc %o0,%o0,%o0
        bcs Ld32
       _ addcc %o3,%o3,%o3
Lc32:   b Ld34
       _ addx %o0,%o0,%o0
Ld01:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc02
       _ addcc %o3,%o3,%o3
Ld02:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc03
       _ addcc %o3,%o3,%o3
Ld03:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc04
       _ addcc %o3,%o3,%o3
Ld04:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc05
       _ addcc %o3,%o3,%o3
Ld05:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc06
       _ addcc %o3,%o3,%o3
Ld06:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc07
       _ addcc %o3,%o3,%o3
Ld07:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc08
       _ addcc %o3,%o3,%o3
Ld08:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc09
       _ addcc %o3,%o3,%o3
Ld09:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc10
       _ addcc %o3,%o3,%o3
Ld10:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc11
       _ addcc %o3,%o3,%o3
Ld11:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc12
       _ addcc %o3,%o3,%o3
Ld12:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc13
       _ addcc %o3,%o3,%o3
Ld13:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc14
       _ addcc %o3,%o3,%o3
Ld14:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc15
       _ addcc %o3,%o3,%o3
Ld15:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc16
       _ addcc %o3,%o3,%o3
Ld16:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc17
       _ addcc %o3,%o3,%o3
Ld17:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc18
       _ addcc %o3,%o3,%o3
Ld18:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc19
       _ addcc %o3,%o3,%o3
Ld19:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc20
       _ addcc %o3,%o3,%o3
Ld20:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc21
       _ addcc %o3,%o3,%o3
Ld21:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc22
       _ addcc %o3,%o3,%o3
Ld22:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc23
       _ addcc %o3,%o3,%o3
Ld23:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc24
       _ addcc %o3,%o3,%o3
Ld24:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc25
       _ addcc %o3,%o3,%o3
Ld25:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc26
       _ addcc %o3,%o3,%o3
Ld26:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc27
       _ addcc %o3,%o3,%o3
Ld27:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc28
       _ addcc %o3,%o3,%o3
Ld28:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc29
       _ addcc %o3,%o3,%o3
Ld29:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc30
       _ addcc %o3,%o3,%o3
Ld30:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Lc31
       _ addcc %o3,%o3,%o3
Ld31:   addxcc %o0,%o0,%o0
        xor %o3,%o1,%o3
        bcc Ld33
       _ addcc %o3,%o3,%o3
Ld32:   xor %o3,%o1,%o3
Ld33:   addx %o0,%o0,%o0
Ld34:   // Now x*y = 2^32*%o0+%o3
        retl
       _ st %o3,[%o2]
#endif

