// Externe Routinen zu ARILEV1.D
// Prozessor: SPARC 64-bit
// Compiler: GNU-C oder ...
// Parameter-Übergabe: in Registern %o0-%o5.
// Parameter-Übergabe: in Registern %o0-%o5.
//   Argumente vom Typ uint8, uint16, uint32 sind bereits vom Aufrufer zu
//   uint64 umgewandelt worden (zero-extend, "srl reg,0,reg").
//   Argumente vom Typ sint8, sint16, sint32 sind bereits vom Aufrufer zu
//   sint64 umgewandelt worden (sign-extend, "sra reg,0,reg").
//   Ergebnisse vom Typ uint8, uint16, uint32 müssen vor Rückgabe zu uint64
//   umgewandelt werden (zero-extend, "srl reg,0,reg").
//   Ergebnisse vom Typ sint8, sint16, sint32 müssen vor Rückgabe zu sint64
//   umgewandelt werden (sign-extend, "sra reg,0,reg").
// Einstellungen: intCsize=32, intDsize=32.

#ifdef ASM_UNDERSCORE
  #define C(entrypoint) _##entrypoint
#else
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

        .register %g2,#scratch

        .global C(mulu16_),C(mulu32_),C(mulu32_unchecked),C(mulu64_)
        .global C(divu_6432_3232_),C(divu_3216_1616_)
        .global C(copy_loop_up),C(copy_loop_down),C(fill_loop_up),C(fill_loop_down)
        .global C(clear_loop_up),C(clear_loop_down)
        .global C(test_loop_up),C(test_loop_down)
        .global C(xor_loop_up),C(compare_loop_up),C(shiftleftcopy_loop_up),C(shiftxor_loop_up)
#if CL_DS_BIG_ENDIAN_P
        .global C(or_loop_up),C(and_loop_up),C(eqv_loop_up)
        .global C(nand_loop_up),C(nor_loop_up),C(andc2_loop_up),C(orc2_loop_up)
        .global C(not_loop_up)
        .global C(and_test_loop_up)
        .global C(add_loop_down),C(addto_loop_down),C(inc_loop_down)
        .global C(sub_loop_down),C(subx_loop_down),C(subfrom_loop_down),C(dec_loop_down)
        .global C(neg_loop_down)
        .global C(shift1left_loop_down),C(shiftleft_loop_down),C(shiftleftcopy_loop_down)
        .global C(shift1right_loop_up),C(shiftright_loop_up),C(shiftrightsigned_loop_up),C(shiftrightcopy_loop_up)
        .global C(mulusmall_loop_down),C(mulu_loop_down),C(muluadd_loop_down),C(mulusub_loop_down)
#else
        .global C(or_loop_down),C(xor_loop_down),C(and_loop_down),C(eqv_loop_down)
        .global C(nand_loop_down),C(nor_loop_down),C(andc2_loop_down),C(orc2_loop_down)
        .global C(not_loop_down)
        .global C(and_test_loop_down),C(compare_loop_down)
        .global C(add_loop_up),C(addto_loop_up),C(inc_loop_up)
        .global C(sub_loop_up),C(subx_loop_up),C(subfrom_loop_up),C(dec_loop_up)
        .global C(neg_loop_up)
        .global C(shift1left_loop_up),C(shiftleft_loop_up)
        .global C(shift1right_loop_down),C(shiftright_loop_down),C(shiftrightsigned_loop_down),C(shiftrightcopy_loop_down)
        .global C(mulusmall_loop_up),C(mulu_loop_up),C(muluadd_loop_up),C(mulusub_loop_up)
#endif

#define LOOP_TYPE  1    // 1: Standard-Schleifen
                        // 2: Schleifen ohne Pointer, nur mit Zähler
#define STANDARD_LOOPS  (LOOP_TYPE==1)
#define COUNTER_LOOPS  (LOOP_TYPE==2)

// extern uint32 mulu16_ (uint16 arg1, uint16 arg2);
// ergebnis := arg1*arg2.
        DECLARE_FUNCTION(mulu16_)
C(mulu16_:) // Input in %o0,%o1, Output in %o0
        umul %o0,%o1,%o2
        retl
       _ srl %o2,0,%o0

// extern struct { uint32 lo; uint32 hi; } mulu32_ (uint32 arg1, uint32 arg2);
// 2^32*hi+lo := arg1*arg2.
        DECLARE_FUNCTION(mulu32_)
C(mulu32_:) // Input in %o0,%o1, Output in %o0,%g1
        umul %o0,%o1,%o2
        rd %y,%g1
        retl
       _ srl %o2,0,%o0

// extern uint32 mulu32_unchecked (uint32 x, uint32 y);
// ergebnis := arg1*arg2 < 2^32.
        DECLARE_FUNCTION(mulu32_unchecked)
C(mulu32_unchecked:) // Input in %o0,%o1, Output in %o0
        umul %o0,%o1,%o2
        retl
       _ srl %o2,0,%o0

// extern struct { uint64 lo; uint64 hi; } mulu64_ (uint64 arg1, uint64 arg2);
// 2^64*hi+lo := arg1*arg2.
        DECLARE_FUNCTION(mulu64_)
C(mulu64_:) // Input in %o0,%o1, Output in %o0,%g2
        srlx %o0,32,%o2         // %o2 = high32(arg1)
        srl %o0,0,%o0           // %o0 = low32(arg1)
        srlx %o1,32,%o3         // %o3 = high32(arg2)
        srl %o1,0,%o1           // %o1 = low32(arg2)
        mulx %o2,%o3,%g2        // high part
        mulx %o2,%o1,%o2        // first mid part
        mulx %o0,%o3,%o3        // second mid part
        addcc %o2,%o3,%o2       // sum of mid parts
        mov 0,%o3
        movcs %xcc,1,%o3        // carry from sum of mid parts
        sllx %o3,32,%o3
        add %g2,%o3,%g2         // add to high part
        srlx %o2,32,%o3
        add %g2,%o3,%g2         // add high32(midparts) to high part
        mulx %o0,%o1,%o0        // low part
        sllx %o2,32,%o2
        addcc %o0,%o2,%o0       // add low32(midparts)*2^32 to low part
        add %g2,1,%o3
        retl
       _ movcs %xcc,%o3,%g2     // add carry to high part

// extern struct { uint32 q; uint32 r; } divu_6432_3232_ (uint32 xhi, uint32 xlo, uint32 y);
// x = 2^32*xhi+xlo = q*y+r schreiben. Sei bekannt, daß 0 <= x < 2^32*y .
        DECLARE_FUNCTION(divu_6432_3232_)
C(divu_6432_3232_:) // Input in %o0,%o1,%o2, Output in %o0,%g1
        wr %o0,%g0,%y
        udiv %o1,%o2,%o0        // x durch y dividieren, %o0 := q
        umul %o0,%o2,%g1        // %g1 := (q*y) mod 2^32
        sub %o1,%g1,%g1         // %g1 := (xlo-q*y) mod 2^32 = r
        retl
       _ srl %o0,0,%o0

// extern struct { uint16 q; uint16 r; } divu_3216_1616_ (uint32 x, uint16 y);
// x = q*y+r schreiben. Sei bekannt, daß 0 <= x < 2^16*y .
        DECLARE_FUNCTION(divu_3216_1616_)
C(divu_3216_1616_:) // Input in %o0,%o1, Output in %o0 (Rest und Quotient).
        wr %g0,%g0,%y
        udiv %o0,%o1,%o2        // dividieren, Quotient nach %o2
#if 0 // Who says that %y has some meaningful contents after `udiv' ??
        rd %y,%g1               // Rest aus %y
#else
        umul %o2,%o1,%g1        // %g1 := (q*y) mod 2^32
        sub %o0,%g1,%g1         // %g1 := (x-q*y) mod 2^32 = r
#endif
        sll %g1,16,%g1          // in die oberen 16 Bit schieben
        or %o2,%g1,%o0
        retl
       _ srl %o0,0,%o0

#if !defined(__GNUC__)
        .global C(_get_g1)
// extern uint32 _get_g1 (void);
        DECLARE_FUNCTION(_get_g1)
C(_get_g1:)
        retl
       _ srl %g1,0,%o0
#endif

#if !defined(__GNUC__)
        .global C(_get_g2)
// extern uint64 _get_g2 (void);
        DECLARE_FUNCTION(_get_g2)
C(_get_g2:)
        retl
       _ mov %g2,%o0
#endif

// extern uintD* copy_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
        DECLARE_FUNCTION(copy_loop_up)
C(copy_loop_up:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ nop
1:        ldx [%o0],%o3
          add %o0,8,%o0
          stx %o3,[%o1]
          subcc %o2,1,%o2
          bne,pt %xcc,1b
         _ add %o1,8,%o1
2:      retl
       _ mov %o1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,8,%o1
        sub %g0,%o2,%o2         // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &sourceptr[count]
        sub %o1,%o2,%o1         // %o1 = &destptr[count-1]
1:        ldx [%o0+%o2],%o3     // nächstes Digit holen
          addcc %o2,8,%o2       // Zähler "erniedrigen", Pointer erhöhen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ add %o1,8,%o0
#endif

// extern uintD* copy_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
        DECLARE_FUNCTION(copy_loop_down)
C(copy_loop_down:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o0,8,%o0
1:        ldx [%o0],%o3
          sub %o1,8,%o1
          stx %o3,[%o1]
          subcc %o2,1,%o2
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov %o1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o0,8,%o0
        sllx %o2,3,%o2          // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &sourceptr[-count-1]
        sub %o1,%o2,%o1         // %o1 = &destptr[-count]
1:        ldx [%o0+%o2],%o3     // nächstes Digit holen
          subcc %o2,8,%o2       // Zähler erniedrigen, Pointer erniedrigen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ mov %o1,%o0
#endif

// extern uintD* fill_loop_up (uintD* destptr, uintC count, uintD filler);
        DECLARE_FUNCTION(fill_loop_up)
C(fill_loop_up:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ nop
1:        stx %o2,[%o0]
          subcc %o1,1,%o1
          bne,pt %xcc,1b
         _ add %o0,8,%o0
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %o0,8,%o0
        sub %g0,%o1,%o1         // %o1 = -count
        sllx %o1,3,%o1          // %o1 = -8*count
        sub %o0,%o1,%o0         // %o0 = &destptr[count-1]
1:        addcc %o1,8,%o1       // Zähler "erniedrigen", Pointer erhöhen
          bne,pt %xcc,1b
         _ stx %o2,[%o0+%o1]    // Digit ablegen
2:      retl
       _ add %o0,8,%o0
#endif

// extern uintD* fill_loop_down (uintD* destptr, uintC count, uintD filler);
        DECLARE_FUNCTION(fill_loop_down)
C(fill_loop_down:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %o0,8,%o0
1:        stx %o2,[%o0]
          subcc %o1,1,%o1
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ add %o0,8,%o0
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sllx %o1,3,%o1         // %o1 = 8*count
        sub %o0,%o1,%o0         // %o0 = &destptr[-count]
1:        subcc %o1,8,%o1       // Zähler erniedrigen, Pointer erniedrigen
          bne,pt %xcc,1b
         _ stx %o2,[%o0+%o1]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern uintD* clear_loop_up (uintD* destptr, uintC count);
        DECLARE_FUNCTION(clear_loop_up)
C(clear_loop_up:) // Input in %o0,%o1, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ nop
1:        stx %g0,[%o0]
          subcc %o1,1,%o1
          bne,pt %xcc,1b
         _ add %o0,8,%o0
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %o0,8,%o0
        sub %g0,%o1,%o1         // %o1 = -count
        sllx %o1,3,%o1          // %o1 = -8*count
        sub %o0,%o1,%o0         // %o0 = &destptr[count-1]
1:        addcc %o1,8,%o1       // Zähler "erniedrigen", Pointer erhöhen
          bne,pt %xcc,1b
         _ stx %g0,[%o0+%o1]    // Digit 0 ablegen
2:      retl
       _ add %o0,8,%o0
#endif

// extern uintD* clear_loop_down (uintD* destptr, uintC count);
        DECLARE_FUNCTION(clear_loop_down)
C(clear_loop_down:) // Input in %o0,%o1, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %o0,8,%o0
1:        stx %g0,[%o0]
          subcc %o1,1,%o1
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ add %o0,8,%o0
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sllx %o1,3,%o1         // %o1 = 8*count
        sub %o0,%o1,%o0         // %o0 = &destptr[-count]
1:        subcc %o1,8,%o1       // Zähler erniedrigen, Pointer erniedrigen
          bne,pt %xcc,1b
         _ stx %g0,[%o0+%o1]    // Digit 0 ablegen
2:      retl
       _ nop
#endif

// extern boolean test_loop_up (uintD* ptr, uintC count);
        DECLARE_FUNCTION(test_loop_up)
C(test_loop_up:) // Input in %o0,%o1, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ nop
          ldx [%o0],%o2
1:        add %o0,8,%o0
          brnz,pn %o2,3f
         _ subcc %o1,1,%o1
          bne,a,pt %xcc,1b
         __ ldx [%o0],%o2
2:      retl
       _ mov 0,%o0
3:      retl
       _ mov 1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %g0,%o1,%o1        // %o1 = -count
        sllx %o1,3,%o1          // %o1 = -8*count
        sub %o0,%o1,%o0         // %o0 = &ptr[count]
          ldx [%o0+%o1],%o2     // nächstes Digit holen
1:        brnz,pn %o2,3f        // testen
         _ addcc %o1,8,%o1      // Zähler "erniedrigen", Pointer erhöhen
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o1],%o2   // nächstes Digit holen
2:      retl
       _ mov 0,%o0
3:      retl
       _ mov 1,%o0
#endif

// extern boolean test_loop_down (uintD* ptr, uintC count);
        DECLARE_FUNCTION(test_loop_down)
C(test_loop_down:) // Input in %o0,%o1, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %o0,8,%o0
          ldx [%o0],%o2
1:        sub %o0,8,%o0
          brnz,pn %o2,3f
         _ subcc %o1,1,%o1
          bne,a,pt %xcc,1b
         __ ldx [%o0],%o2
2:      retl
       _ mov 0,%o0
3:      retl
       _ mov 1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sllx %o1,3,%o1         // %o1 = 8*count
        sub %o0,%o1,%o0         // %o0 = &ptr[-count]
        sub %o1,8,%o1
          ldx [%o0+%o1],%o2     // nächstes Digit holen
1:        brnz,pn %o2,3f        // testen
         _ subcc %o1,8,%o1      // Zähler erniedrigen, Pointer erniedrigen
          bcc,a,pt %xcc,1b
         __ ldx [%o0+%o1],%o2   // nächstes Digit holen
2:      retl
       _ mov 0,%o0
3:      retl
       _ mov 1,%o0
#endif

#if CL_DS_BIG_ENDIAN_P

// extern void or_loop_up (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(or_loop_up)
C(or_loop_up:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          or %o3,%o4,%o3        // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ add %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o0,8,%o0
        sub %g0,%o2,%o2         // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[count-1]
        sub %o1,%o2,%o1         // %o1 = &yptr[count]
1:        ldx [%o1+%o2],%o3     // nächstes Digit holen
          addcc %o2,8,%o2       // Zähler "erniedrigen", Pointer erhöhen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          or %o4,%o3,%o3        // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

#endif

// extern void xor_loop_up (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(xor_loop_up)
C(xor_loop_up:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          xor %o3,%o4,%o3       // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ add %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o0,8,%o0
        sub %g0,%o2,%o2         // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[count-1]
        sub %o1,%o2,%o1         // %o1 = &yptr[count]
1:        ldx [%o1+%o2],%o3     // nächstes Digit holen
          addcc %o2,8,%o2       // Zähler "erniedrigen", Pointer erhöhen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          xor %o4,%o3,%o3       // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

#if CL_DS_BIG_ENDIAN_P

// extern void and_loop_up (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(and_loop_up)
C(and_loop_up:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          and %o3,%o4,%o3       // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ add %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o0,8,%o0
        sub %g0,%o2,%o2         // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[count-1]
        sub %o1,%o2,%o1         // %o1 = &yptr[count]
1:        ldx [%o1+%o2],%o3     // nächstes Digit holen
          addcc %o2,8,%o2       // Zähler "erniedrigen", Pointer erhöhen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          and %o4,%o3,%o3       // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void eqv_loop_up (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(eqv_loop_up)
C(eqv_loop_up:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          xnor %o3,%o4,%o3      // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ add %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o0,8,%o0
        sub %g0,%o2,%o2         // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[count-1]
        sub %o1,%o2,%o1         // %o1 = &yptr[count]
1:        ldx [%o1+%o2],%o3     // nächstes Digit holen
          addcc %o2,8,%o2       // Zähler "erniedrigen", Pointer erhöhen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          xnor %o4,%o3,%o3      // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void nand_loop_up (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(nand_loop_up)
C(nand_loop_up:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          and %o3,%o4,%o3       // verknüpfen
          xnor %g0,%o3,%o3
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ add %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o0,8,%o0
        sub %g0,%o2,%o2         // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[count-1]
        sub %o1,%o2,%o1         // %o1 = &yptr[count]
1:        ldx [%o1+%o2],%o3     // nächstes Digit holen
          addcc %o2,8,%o2       // Zähler "erniedrigen", Pointer erhöhen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          and %o4,%o3,%o3       // beide verknüpfen
          xnor %g0,%o3,%o3
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void nor_loop_up (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(nor_loop_up)
C(nor_loop_up:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          or %o3,%o4,%o3        // verknüpfen
          xnor %g0,%o3,%o3
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ add %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o0,8,%o0
        sub %g0,%o2,%o2         // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[count-1]
        sub %o1,%o2,%o1         // %o1 = &yptr[count]
1:        ldx [%o1+%o2],%o3     // nächstes Digit holen
          addcc %o2,8,%o2       // Zähler "erniedrigen", Pointer erhöhen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          or %o4,%o3,%o3        // beide verknüpfen
          xnor %g0,%o3,%o3
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void andc2_loop_up (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(andc2_loop_up)
C(andc2_loop_up:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          andn %o3,%o4,%o3      // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ add %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o0,8,%o0
        sub %g0,%o2,%o2         // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[count-1]
        sub %o1,%o2,%o1         // %o1 = &yptr[count]
1:        ldx [%o1+%o2],%o3     // nächstes Digit holen
          addcc %o2,8,%o2       // Zähler "erniedrigen", Pointer erhöhen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          andn %o4,%o3,%o3      // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void orc2_loop_up (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(orc2_loop_up)
C(orc2_loop_up:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          orn %o3,%o4,%o3       // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ add %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o0,8,%o0
        sub %g0,%o2,%o2         // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[count-1]
        sub %o1,%o2,%o1         // %o1 = &yptr[count]
1:        ldx [%o1+%o2],%o3     // nächstes Digit holen
          addcc %o2,8,%o2       // Zähler "erniedrigen", Pointer erhöhen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          orn %o4,%o3,%o3       // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void not_loop_up (uintD* xptr, uintC count);
        DECLARE_FUNCTION(not_loop_up)
C(not_loop_up:) // Input in %o0,%o1
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ nop
1:        ldx [%o0],%o2
          subcc %o1,1,%o1
          xnor %g0,%o2,%o2
          stx %o2,[%o0]
          bne,pt %xcc,1b
         _ add %o0,8,%o0
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %o0,8,%o0
        sub %g0,%o1,%o1         // %o1 = -count
        sllx %o1,3,%o1          // %o1 = -8*count
        sub %o0,%o1,%o0         // %o0 = &destptr[count-1]
1:        addcc %o1,8,%o1       // Zähler "erniedrigen", Pointer erhöhen
          ldx [%o0+%o1],%o2     // nächstes Digit holen
          xnor %g0,%o2,%o2
          bne,pt %xcc,1b
         _ stx %o2,[%o0+%o1]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern boolean and_test_loop_up (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(and_test_loop_up)
C(and_test_loop_up:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ nop
1:        ldx [%o0],%o3
          ldx [%o1],%o4
          add %o0,8,%o0
          andcc %o3,%o4,%g0
          bne,pn %xcc,3f
         _ subcc %o2,1,%o2
          bne,pt %xcc,1b
         _ add %o1,8,%o1
2:      retl
       _ mov 0,%o0
3:      retl
       _ mov 1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %g0,%o2,%o2        // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[count]
        sub %o1,%o2,%o1         // %o1 = &yptr[count]
          ldx [%o0+%o2],%o3     // nächstes Digit holen
1:        ldx [%o1+%o2],%o4     // noch ein Digit holen
          andcc %o3,%o4,%g0     // beide verknüpfen
          bne,pn %xcc,3f
         _ addcc %o2,8,%o2      // Zähler "erniedrigen", Pointer erhöhen
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o2],%o3   // nächstes Digit holen
2:      retl
       _ mov 0,%o0
3:      retl
       _ mov 1,%o0
#endif

#endif

// extern cl_signean compare_loop_up (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(compare_loop_up)
C(compare_loop_up:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ nop
          ldx [%o0],%o3
1:        ldx [%o1],%o4
          add %o0,8,%o0
          subcc %o3,%o4,%g0
          bne,pn %xcc,3f
         _ add %o1,8,%o1
          subcc %o2,1,%o2
          bne,a,pt %xcc,1b
         __ ldx [%o0],%o3
2:      retl
       _ mov 0,%o0
3:      mov 1,%o0
        movlu %xcc,-1,%o0
        retl
       _ sra %o0,0,%o0          // sign-extend %o0
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %g0,%o2,%o2        // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[count]
        sub %o1,%o2,%o1         // %o1 = &yptr[count]
          ldx [%o0+%o2],%o3     // nächstes Digit holen
1:        ldx [%o1+%o2],%o4     // noch ein Digit holen
          subcc %o3,%o4,%g0     // vergleichen
          bne,pn %xcc,3f
         _ addcc %o2,8,%o2      // Zähler "erniedrigen", Pointer erhöhen
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o2],%o3   // nächstes Digit holen
2:      retl
       _ mov 0,%o0
3:      subcc %o3,%o4,%g0       // nochmals vergleichen
        mov 1,%o0
        movlu %xcc,-1,%o0
        retl
       _ sra %o0,0,%o0          // sign-extend %o0
#endif

#if CL_DS_BIG_ENDIAN_P

// extern uintD add_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
        DECLARE_FUNCTION(add_loop_down)
C(add_loop_down:) // Input in %o0,%o1,%o2,%o3, verändert %g1, Output in %o0
#if STANDARD_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %g0,%g1            // Carry := 0
        sub %o0,8,%o0
1:        ldx [%o0],%o4         // source1-digit
          sub %o1,8,%o1
          ldx [%o1],%o5         // source2-digit
          addcc %o4,%g1,%o4
          movcc %xcc,0,%g1      // %g1|%o4 := %o4 + alter Carry %g1
          addcc %o4,%o5,%o4
          movcs %xcc,1,%g1      // %g1|%o4 := %o4 + alter Carry %g1 + %o5
          sub %o2,8,%o2
          stx %o4,[%o2]         // Digit ablegen
          subcc %o3,1,%o3
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov %g1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %g0,%g1            // Carry := 0
        sub %o0,8,%o0
        sub %o1,8,%o1
        sllx %o3,3,%o3          // %o3 = 8*count
        sub %o0,%o3,%o0         // %o0 = &sourceptr1[-count-1]
        sub %o1,%o3,%o1         // %o1 = &sourceptr2[-count-1]
        sub %o2,%o3,%o2         // %o2 = &destptr[-count]
1:        ldx [%o0+%o3],%o4     // source1-digit
          ldx [%o1+%o3],%o5     // source2-digit
          addcc %o4,%g1,%o4
          movcc %xcc,0,%g1      // %g1|%o4 := %o4 + alter Carry %g1
          addcc %o4,%o5,%o4
          movcs %xcc,1,%g1      // %g1|%o4 := %o4 + alter Carry %g1 + %o5
          subcc %o3,8,%o3
          bne,pt %xcc,1b
         _ stx %o4,[%o2+%o3]    // Digit ablegen
2:      retl
       _ mov %g1,%o0
#endif

// extern uintD addto_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
        DECLARE_FUNCTION(addto_loop_down)
C(addto_loop_down:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ mov %g0,%o5            // Carry := 0
        sub %o0,8,%o0
1:        ldx [%o0],%o3         // source-digit
          sub %o1,8,%o1
          ldx [%o1],%o4         // dest-digit
          addcc %o3,%o5,%o3
          movcc %xcc,0,%o5      // %o5|%o3 := %o3 + alter Carry %o5
          addcc %o3,%o4,%o4
          movcs %xcc,1,%o5      // %o5|%o4 := %o3 + alter Carry %o5 + %o4
          stx %o4,[%o1]         // Digit ablegen
          subcc %o2,1,%o2
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov %o5,%o0
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ mov %g0,%o5            // Carry := 0
        sub %o0,8,%o0
        sub %o1,8,%o1
        sllx %o2,3,%o2          // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &sourceptr[-count-1]
        sub %o1,%o2,%o1         // %o1 = &destptr[-count-1]
          ldx [%o0+%o2],%o3     // source-digit
1:        ldx [%o1+%o2],%o4     // dest-digit
          addcc %o3,%o5,%o3
          movcc %xcc,0,%o5      // %o5|%o3 := %o3 + alter Carry %o5
          addcc %o3,%o4,%o4
          movcs %xcc,1,%o5      // %o5|%o4 := %o3 + alter Carry %o5 + %o4
          stx %o4,[%o1+%o2]     // Digit ablegen
          subcc %o2,8,%o2
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o2],%o3   // source-digit
2:      retl
       _ mov %o5,%o0
#endif

// extern uintD inc_loop_down (uintD* ptr, uintC count);
        DECLARE_FUNCTION(inc_loop_down)
C(inc_loop_down:) // Input in %o0,%o1, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %o0,8,%o0
1:        ldx [%o0],%o2
          addcc %o2,1,%o2
          bne,pn %xcc,3f
         _ stx %o2,[%o0]
          subcc %o1,1,%o1
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov 1,%o0
3:      retl
       _ mov 0,%o0
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %o0,8,%o0
        sllx %o1,3,%o1          // %o1 = 8*count
        sub %o0,%o1,%o0         // %o0 = &ptr[-count-1]
          ldx [%o0+%o1],%o2     // digit holen
1:        addcc %o2,1,%o2       // incrementieren
          bne,pn %xcc,3f
         _ stx %o2,[%o0+%o1]    // ablegen
          subcc %o1,8,%o1       // Zähler erniedrigen, Pointer erniedrigen
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o1],%o2
2:      retl
       _ mov 1,%o0
3:      retl
       _ mov 0,%o0
#endif

// extern uintD sub_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
        DECLARE_FUNCTION(sub_loop_down)
C(sub_loop_down:) // Input in %o0,%o1,%o2,%o3, verändert %g1, Output in %o0
#if STANDARD_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %g0,%g1            // Carry := 0
        sub %o1,8,%o1
1:        ldx [%o1],%o5         // source2-digit
          sub %o0,8,%o0
          ldx [%o0],%o4         // source1-digit
          addcc %o5,%g1,%o5
          movcc %xcc,0,%g1      // %g1|%o5 := %o5 + alter Carry %g1
          subcc %o4,%o5,%o4
          movcs %xcc,1,%g1      // %o4-2^64*%g1 := %o4 - %o5 - alter Carry %g1
          sub %o2,8,%o2
          stx %o4,[%o2]         // Digit ablegen
          subcc %o3,1,%o3
          bne,pt %xcc,1b
         _ sub %o1,8,%o1
2:      retl
       _ mov %g1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %g0,%g1            // Carry := 0
        sub %o0,8,%o0
        sub %o1,8,%o1
        sllx %o3,3,%o3          // %o3 = 8*count
        sub %o0,%o3,%o0         // %o0 = &sourceptr1[-count-1]
        sub %o1,%o3,%o1         // %o1 = &sourceptr2[-count-1]
        sub %o2,%o3,%o2         // %o2 = &destptr[-count]
1:        ldx [%o0+%o3],%o4     // source1-digit
          ldx [%o1+%o3],%o5     // source2-digit
          addcc %o5,%g1,%o5
          movcc %xcc,0,%g1      // %g1|%o5 := %o5 + alter Carry %g1
          subcc %o4,%o5,%o4
          movcs %xcc,1,%g1      // %o4-2^64*%g1 := %o4 - %o5 - alter Carry %g1
          subcc %o3,8,%o3
          bne,pt %xcc,1b
         _ stx %o4,[%o2+%o3]    // Digit ablegen
2:      retl
       _ mov %g1,%o0
#endif

// extern uintD subx_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);
        DECLARE_FUNCTION(subx_loop_down)
C(subx_loop_down:) // Input in %o0,%o1,%o2,%o3,%o4, verändert %g1, Output in %o0
#if STANDARD_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %o4,%g1            // Carry (0 oder -1)
        sub %o1,8,%o1
1:        ldx [%o1],%o5         // source2-digit
          sub %o0,8,%o0
          ldx [%o0],%o4         // source1-digit
          subcc %o5,%g1,%o5
          movcc %xcc,0,%g1      // %o5-2^64*%g1 := %o5 - alter Carry %g1
          subcc %o4,%o5,%o4
          movcs %xcc,-1,%g1     // %o4+2^64*%g1 := %o4 - %o5 + alter Carry %g1
          sub %o2,8,%o2
          stx %o4,[%o2]         // Digit ablegen
          subcc %o3,1,%o3
          bne,pt %xcc,1b
         _ sub %o1,8,%o1
2:      retl
       _ mov %g1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %o4,%g1            // Carry (0 oder -1)
        sub %o0,8,%o0
        sub %o1,8,%o1
        sllx %o3,3,%o3          // %o3 = 8*count
        sub %o0,%o3,%o0         // %o0 = &sourceptr1[-count-1]
        sub %o1,%o3,%o1         // %o1 = &sourceptr2[-count-1]
        sub %o2,%o3,%o2         // %o2 = &destptr[-count]
1:        ldx [%o1+%o3],%o5     // source2-digit
          ldx [%o0+%o3],%o4     // source1-digit
          subcc %o5,%g1,%o5
          movcc %xcc,0,%g1      // %o5-2^64*%g1 := %o5 - alter Carry %g1
          subcc %o4,%o5,%o4
          movcs %xcc,-1,%g1     // %o4+2^64*%g1 := %o4 - %o5 + alter Carry %g1
          subcc %o3,8,%o3
          bne,pt %xcc,1b
         _ stx %o4,[%o2+%o3]    // Digit ablegen
2:      retl
       _ mov %g1,%o0
#endif

// extern uintD subfrom_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
        DECLARE_FUNCTION(subfrom_loop_down)
C(subfrom_loop_down:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ mov %g0,%o5            // Carry := 0
        sub %o0,8,%o0
1:        ldx [%o0],%o3         // source-digit
          sub %o1,8,%o1
          ldx [%o1],%o4         // dest-digit
          addcc %o3,%o5,%o3
          movcc %xcc,0,%o5      // %o5|%o3 := %o3 + alter Carry %o5
          subcc %o4,%o3,%o4
          movcs %xcc,1,%o5      // %o4-2^64*%o5 := %o4 - %o3 - alter Carry %o5
          stx %o4,[%o1]         // Digit ablegen
          subcc %o2,1,%o2
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov %o5,%o0
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ mov %g0,%o5            // Carry := 0
        sub %o0,8,%o0
        sub %o1,8,%o1
        sllx %o2,3,%o2          // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &sourceptr[-count-1]
        sub %o1,%o2,%o1         // %o1 = &destptr[-count-1]
          ldx [%o0+%o2],%o3     // source-digit
1:        ldx [%o1+%o2],%o4     // dest-digit
          addcc %o3,%o5,%o3
          movcc %xcc,0,%o5      // %o5|%o3 := %o3 + alter Carry %o5
          subcc %o4,%o3,%o4
          movcs %xcc,1,%o5      // %o4-2^64*%o5 := %o4 - %o3 - alter Carry %o5
          stx %o4,[%o1+%o2]     // Digit ablegen
          subcc %o2,8,%o2
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o2],%o3   // source-digit
2:      retl
       _ mov %o5,%o0
#endif

// extern uintD dec_loop_down (uintD* ptr, uintC count);
        DECLARE_FUNCTION(dec_loop_down)
C(dec_loop_down:) // Input in %o0,%o1, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %o0,8,%o0
1:        ldx [%o0],%o2
          subcc %o2,1,%o2
          bcc,pn %xcc,3f
         _ stx %o2,[%o0]
          subcc %o1,1,%o1
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov -1,%o0
3:      retl
       _ mov 0,%o0
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %o0,8,%o0
        sllx %o1,3,%o1          // %o1 = 8*count
        sub %o0,%o1,%o0         // %o0 = &ptr[-count-1]
          ldx [%o0+%o1],%o2     // digit holen
1:        subcc %o2,1,%o2       // decrementieren
          bcc,pn %xcc,3f
         _ stx %o2,[%o0+%o1]    // ablegen
          subcc %o1,8,%o1       // Zähler erniedrigen, Pointer erniedrigen
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o1],%o2
2:      retl
       _ mov -1,%o0
3:      retl
       _ mov 0,%o0
#endif

// extern uintD neg_loop_down (uintD* ptr, uintC count);
        DECLARE_FUNCTION(neg_loop_down)
C(neg_loop_down:) // Input in %o0,%o1, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        // erstes Digit /=0 suchen:
        brz,pn %o1,2f
       _ sub %o0,8,%o0
1:        ldx [%o0],%o2
          subcc %g0,%o2,%o2
          bne,pn %xcc,3f
         _ subcc %o1,1,%o1
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov 0,%o0
3:      // erstes Digit /=0 gefunden, ab jetzt gibt's Carrys
        stx %o2,[%o0]           // 1 Digit negieren
        // alle anderen Digits invertieren:
        be,pn %xcc,5f
       _ sub %o0,8,%o0
4:        ldx [%o0],%o2
          subcc %o1,1,%o1
          xnor %g0,%o2,%o2
          stx %o2,[%o0]
          bne,pt %xcc,4b
         _ sub %o0,8,%o0
5:      retl
       _ mov -1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        // erstes Digit /=0 suchen:
        brz,pn %o1,2f
       _ sub %o0,8,%o0
        sllx %o1,3,%o1          // %o1 = 8*count
        sub %o0,%o1,%o0         // %o0 = &ptr[-count-1]
          ldx [%o0+%o1],%o2     // digit holen
1:        subcc %g0,%o2,%o2     // negieren, testen
          bne,pn %xcc,3f
         _ subcc %o1,8,%o1      // Zähler erniedrigen, Pointer erniedrigen
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o1],%o2
2:      retl
       _ mov 0,%o0
3:      // erstes Digit /=0 gefunden, ab jetzt gibt's Carrys
        // alle anderen Digits invertieren:
        add %o1,8,%o1
        stx %o2,[%o0+%o1]       // ablegen
        subcc %o1,8,%o1
        be,pn %xcc,5f
       _ nop
          ldx [%o0+%o1],%o2
4:        xnor %g0,%o2,%o2
          stx %o2,[%o0+%o1]
          subcc %o1,8,%o1
          bne,a,pt %xcc,4b
         __ ldx [%o0+%o1],%o2
5:      retl
       _ mov -1,%o0
#endif

// extern uintD shift1left_loop_down (uintD* ptr, uintC count);
        DECLARE_FUNCTION(shift1left_loop_down)
C(shift1left_loop_down:) // Input in %o0,%o1, Output in %o0
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ mov 0,%o3              // Carry := 0
        sub %o0,8,%o0
1:        ldx [%o0],%o2         // Digit
          addcc %o2,%o2,%o4     // shiften
          add %o4,%o3,%o4       // und carry
          srlx %o2,63,%o3       // neues Carry
          stx %o4,[%o0]         // Digit ablegen
          subcc %o1,1,%o1
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov %o3,%o0

// extern uintD shiftleft_loop_down (uintD* ptr, uintC count, uintC i, uintD carry);
        DECLARE_FUNCTION(shiftleft_loop_down)
C(shiftleft_loop_down:) // Input in %o0,%o1,%o2,%o3, verändert %g1, Output in %o0
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %g0,%o2,%g1        // 64-i (mod 64)
        sub %o0,8,%o0
1:        ldx [%o0],%o4         // Digit
          subcc %o1,1,%o1
          sllx %o4,%o2,%o5      // dessen niedere (64-i) Bits
          or %o3,%o5,%o5        // mit dem alten Carry kombinieren
          stx %o5,[%o0]         // Digit ablegen
          srlx %o4,%g1,%o3      // dessen höchste i Bits liefern den neuen Carry
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov %o3,%o0

// extern uintD shiftleftcopy_loop_down (uintD* sourceptr, uintD* destptr, uintC count, uintC i);
        DECLARE_FUNCTION(shiftleftcopy_loop_down)
C(shiftleftcopy_loop_down:) // Input in %o0,%o1,%o2,%o3, verändert %g1,%g2, Output in %o0
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ mov 0,%o4              // Carry := 0
        sub %g0,%o3,%g1         // 64-i (mod 64)
        sub %o0,8,%o0
1:        ldx [%o0],%o5         // Digit
          subcc %o2,1,%o2
          sllx %o5,%o3,%g2      // dessen niedere (64-i) Bits
          or %o4,%g2,%g2        // mit dem alten Carry kombinieren
          sub %o1,8,%o1
          stx %g2,[%o1]         // Digit ablegen
          srlx %o5,%g1,%o4      // dessen höchste i Bits liefern den neuen Carry
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov %o4,%o0

// extern uintD shift1right_loop_up (uintD* ptr, uintC count, uintD carry);
        DECLARE_FUNCTION(shift1right_loop_up)
C(shift1right_loop_up:) // Input in %o0,%o1,%o2, Output in %o0
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sllx %o2,63,%o2        // Carry
1:        ldx [%o0],%o3         // Digit
          subcc %o1,1,%o1
          srlx %o3,1,%o4        // shiften
          or %o2,%o4,%o4        // und mit altem Carry kombinieren
          stx %o4,[%o0]         // und ablegen
          sllx %o3,63,%o2       // neuer Carry
          bne,pt %xcc,1b
         _ add %o0,8,%o0
2:      retl
       _ mov %o2,%o0

// extern uintD shiftright_loop_up (uintD* ptr, uintC count, uintC i);
        DECLARE_FUNCTION(shiftright_loop_up)
C(shiftright_loop_up:) // Input in %o0,%o1,%o2, verändert %g1, Output in %o0
//      srl %o1,0,%o1           // zero-extend %o1 = count
        sub %g0,%o2,%g1         // 64-i (mod 64)
        brz,pn %o1,2f
       _ or %g0,%g0,%o3         // Carry := 0
1:        ldx [%o0],%o4         // Digit
          subcc %o1,1,%o1
          srlx %o4,%o2,%o5      // shiften
          or %o3,%o5,%o5        // und mit altem Carry kombinieren
          stx %o5,[%o0]         // und ablegen
          sllx %o4,%g1,%o3      // neuer Carry
          bne,pt %xcc,1b
         _ add %o0,8,%o0
2:      retl
       _ mov %o3,%o0

// extern uintD shiftrightsigned_loop_up (uintD* ptr, uintC count, uintC i);
        DECLARE_FUNCTION(shiftrightsigned_loop_up)
C(shiftrightsigned_loop_up:) // Input in %o0,%o1,%o2, verändert %g1, Output in %o0
//      srl %o1,0,%o1           // zero-extend %o1 = count
        ldx [%o0],%o4           // erstes Digit
        sub %g0,%o2,%g1         // 64-i (mod 64)
        srax %o4,%o2,%o5        // shiften
        stx %o5,[%o0]           // und ablegen
        sllx %o4,%g1,%o3        // neuer Carry
        subcc %o1,1,%o1
        be,pn %xcc,2f
       _ add %o0,8,%o0
1:        ldx [%o0],%o4         // Digit
          subcc %o1,1,%o1
          srlx %o4,%o2,%o5      // shiften
          or %o3,%o5,%o5        // und mit altem Carry kombinieren
          stx %o5,[%o0]         // und ablegen
          sllx %o4,%g1,%o3      // neuer Carry
          bne,pt %xcc,1b
         _ add %o0,8,%o0
2:      retl
       _ mov %o3,%o0

// extern uintD shiftrightcopy_loop_up (uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry);
        DECLARE_FUNCTION(shiftrightcopy_loop_up)
C(shiftrightcopy_loop_up:) // Input in %o0,%o1,%o2,%o3,%o4, verändert %g1,%g2, Output in %o0
//      srl %o2,0,%o2           // zero-extend %o2 = count
        sub %g0,%o3,%g1         // 64-i (mod 64)
        brz,pn %o2,2f
       _ sllx %o4,%g1,%g2       // erster Carry
1:        ldx [%o0],%o4         // Digit
          add %o0,8,%o0
          srlx %o4,%o3,%o5      // shiften
          or %g2,%o5,%o5        // und mit altem Carry kombinieren
          stx %o5,[%o1]         // und ablegen
          sllx %o4,%g1,%g2      // neuer Carry
          subcc %o2,1,%o2
          bne,pt %xcc,1b
         _ add %o1,8,%o1
2:      retl
       _ mov %g2,%o0

// extern uintD mulusmall_loop_down (uintD digit, uintD* ptr, uintC len, uintD newdigit);
        DECLARE_FUNCTION(mulusmall_loop_down)
C(mulusmall_loop_down:) // Input in %o0,%o1,%o2,%o3, Output in %o0, verändert %g1
//      srl %o2,0,%o2           // zero-extend %o2 = len
        brz,pn %o2,2f
       _ sub %o1,8,%o1
1:        // nächstes Digit [%o1] mit der 6-Bit-Zahl %o0 multiplizieren
          // und kleinen Carry %o3 dazu:
          ldx [%o1],%o4
          sub %o2,1,%o2
          srlx %o4,32,%o5       // high32(x)
          srl %o4,0,%o4         // low32(x)
          mulx %o4,%o0,%o4      // low32(x)*digit
          mulx %o5,%o0,%o5      // high32(x)*digit
          sllx %o5,32,%g1       // low32(high32(x)*digit)*2^32
          add %g1,%o3,%g1       // plus carry
          addcc %o4,%g1,%o4     // plus low32(x)*digit
          srlx %o5,32,%o3       // high32(high32(x)*digit)
          add %o3,1,%g1
          movcs %xcc,%g1,%o3    // neuer Carry
          stx %o4,[%o1]         // neues Digit ablegen
          brnz,pt %o2,1b
         _ sub %o1,8,%o1
2:      retl
       _ mov %o3,%o0

// extern void mulu_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
        DECLARE_FUNCTION(mulu_loop_down)
C(mulu_loop_down:) // Input in %i0,%i1,%i2,%i3
        save %sp,-192,%sp
        mov 0,%l0               // Carry
        srlx %i0,32,%l1         // %l1 = high32(digit)
        srl %i0,0,%l2           // %l2 = low32(digit)
        mov 1,%l3
        sllx %l3,32,%l3         // %l3 = 2^32
        sub %i1,%i2,%i1         // %i1 = sourceptr - destptr
1:        sub %i2,8,%i2
          ldx [%i1+%i2],%o0     // nächstes Digit
          subcc %i3,1,%i3
          // mit digit multiplizieren: (%l1*2^32+%l2) * %o0 + %l0 -> %l0|%o0
          srlx %o0,32,%o1
          srl %o0,0,%o2
          mulx %l1,%o1,%o3      // high part
          mulx %l1,%o2,%o4      // first mid part
          mulx %l2,%o1,%o1      // second mid part
          mulx %l2,%o2,%o2      // low part
          srlx %o2,32,%o5       // low part's upper half
          add %o4,%o5,%o4       // add to one of the mid parts, no carry
          addcc %o4,%o1,%o4     // add other mid part
          add %o3,%l3,%o5
          movcs %xcc,%o5,%o3    // if carry, add 2^32 to the high part
          srlx %o4,32,%o5
          sllx %o4,32,%o4
          srl %o2,0,%o2
          add %o2,%o4,%o0       // combine low32(midparts) and low32(lowpart)
          addcc %o0,%l0,%o0     // alten Carry addieren
          add %o3,%o5,%l0       // add high32(midparts) to high part
          add %l0,1,%o5
          movcs %xcc,%o5,%l0    // neuer Carry
          // Multiplikation fertig
          brnz,pt %i3,1b
         _ stx %o0,[%i2]        // Low-Digit ablegen
        stx %l0,[%i2-8]         // letzten Carry ablegen
        ret
       _ restore

// extern uintD muluadd_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
        DECLARE_FUNCTION(muluadd_loop_down)
C(muluadd_loop_down:) // Input in %i0,%i1,%i2,%i3, Output in %i0
        save %sp,-192,%sp
        mov 0,%l0               // Carry
        srlx %i0,32,%l1         // %l1 = high32(digit)
        srl %i0,0,%l2           // %l2 = low32(digit)
        mov 1,%l3
        sllx %l3,32,%l3         // %l3 = 2^32
        sub %i1,%i2,%i1         // %i1 = sourceptr - destptr
1:        sub %i2,8,%i2
          ldx [%i1+%i2],%o0     // nächstes Digit
          ldx [%i2],%i4         // *destptr
          subcc %i3,1,%i3
          // mit digit multiplizieren: (%l1*2^32+%l2) * %o0 + %l0 -> %l0|%o0
          srlx %o0,32,%o1
          srl %o0,0,%o2
          mulx %l1,%o1,%o3      // high part
          mulx %l1,%o2,%o4      // first mid part
          mulx %l2,%o1,%o1      // second mid part
          mulx %l2,%o2,%o2      // low part
          srlx %o2,32,%o5       // low part's upper half
          add %o4,%o5,%o4       // add to one of the mid parts, no carry
          addcc %o4,%o1,%o4     // add other mid part
          add %o3,%l3,%o5
          movcs %xcc,%o5,%o3    // if carry, add 2^32 to the high part
          srlx %o4,32,%o5
          sllx %o4,32,%o4
          srl %o2,0,%o2
          add %o2,%o4,%o0       // combine low32(midparts) and low32(lowpart)
          addcc %o0,%l0,%o0     // alten Carry addieren
          add %o3,%o5,%l0       // add high32(midparts) to high part
          add %l0,1,%o5
          movcs %xcc,%o5,%l0    // neuer Carry
          // Multiplikation fertig
          addcc %i4,%o0,%o0     // alten *destptr addieren
          add %l0,1,%o2
          movcs %xcc,%o2,%l0    // neuer Carry
          brnz,pt %i3,1b
         _ stx %o0,[%i2]        // Low-Digit ablegen
        mov %l0,%i0             // letzter Carry
        ret
       _ restore

// extern uintD mulusub_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
        DECLARE_FUNCTION(mulusub_loop_down)
C(mulusub_loop_down:) // Input in %i0,%i1,%i2,%i3, Output in %i0
        save %sp,-192,%sp
        mov 0,%l0               // Carry
        srlx %i0,32,%l1         // %l1 = high32(digit)
        srl %i0,0,%l2           // %l2 = low32(digit)
        mov 1,%l3
        sllx %l3,32,%l3         // %l3 = 2^32
        sub %i1,%i2,%i1         // %i1 = sourceptr - destptr
1:        sub %i2,8,%i2
          ldx [%i1+%i2],%o0     // nächstes Digit
          ldx [%i2],%i4         // *destptr
          subcc %i3,1,%i3
          // mit digit multiplizieren: (%l1*2^32+%l2) * %o0 + %l0 -> %l0|%o0
          srlx %o0,32,%o1
          srl %o0,0,%o2
          mulx %l1,%o1,%o3      // high part
          mulx %l1,%o2,%o4      // first mid part
          mulx %l2,%o1,%o1      // second mid part
          mulx %l2,%o2,%o2      // low part
          srlx %o2,32,%o5       // low part's upper half
          add %o4,%o5,%o4       // add to one of the mid parts, no carry
          addcc %o4,%o1,%o4     // add other mid part
          add %o3,%l3,%o5
          movcs %xcc,%o5,%o3    // if carry, add 2^32 to the high part
          srlx %o4,32,%o5
          sllx %o4,32,%o4
          srl %o2,0,%o2
          add %o2,%o4,%o0       // combine low32(midparts) and low32(lowpart)
          addcc %o0,%l0,%o0     // alten Carry addieren
          add %o3,%o5,%l0       // add high32(midparts) to high part
          add %l0,1,%o5
          movcs %xcc,%o5,%l0    // neuer Carry
          // Multiplikation fertig
          subcc %i4,%o0,%o0     // vom alten *destptr subtrahieren
          add %l0,1,%o2
          movcs %xcc,%o2,%l0    // neuer Carry
          brnz,pt %i3,1b
         _ stx %o0,[%i2]        // Low-Digit ablegen
        mov %l0,%i0             // letzter Carry
        ret
       _ restore

#endif

#if !CL_DS_BIG_ENDIAN_P

// extern void or_loop_down (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(or_loop_down)
C(or_loop_down:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
        sub %o0,8,%o0
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          or %o3,%o4,%o3        // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ sub %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sllx %o2,3,%o2         // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[-count]
        sub %o1,%o2,%o1         // %o1 = &yptr[-count]
1:        subcc %o2,8,%o2       // Zähler erniedrigen, Pointer erniedrigen
          ldx [%o1+%o2],%o3     // nächstes Digit holen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          or %o4,%o3,%o3        // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void xor_loop_down (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(xor_loop_down)
C(xor_loop_down:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
        sub %o0,8,%o0
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          xor %o3,%o4,%o3       // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ sub %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sllx %o2,3,%o2         // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[-count]
        sub %o1,%o2,%o1         // %o1 = &yptr[-count]
1:        subcc %o2,8,%o2       // Zähler erniedrigen, Pointer erniedrigen
          ldx [%o1+%o2],%o3     // nächstes Digit holen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          xor %o4,%o3,%o3       // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void and_loop_down (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(and_loop_down)
C(and_loop_down:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
        sub %o0,8,%o0
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          and %o3,%o4,%o3       // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ sub %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sllx %o2,3,%o2         // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[-count]
        sub %o1,%o2,%o1         // %o1 = &yptr[-count]
1:        subcc %o2,8,%o2       // Zähler erniedrigen, Pointer erniedrigen
          ldx [%o1+%o2],%o3     // nächstes Digit holen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          and %o4,%o3,%o3       // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void eqv_loop_down (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(eqv_loop_down)
C(eqv_loop_down:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
        sub %o0,8,%o0
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          xnor %o3,%o4,%o3      // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ sub %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sllx %o2,3,%o2         // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[-count]
        sub %o1,%o2,%o1         // %o1 = &yptr[-count]
1:        subcc %o2,8,%o2       // Zähler erniedrigen, Pointer erniedrigen
          ldx [%o1+%o2],%o3     // nächstes Digit holen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          xnor %o4,%o3,%o3      // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void nand_loop_down (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(nand_loop_down)
C(nand_loop_down:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
        sub %o0,8,%o0
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          and %o3,%o4,%o3       // verknüpfen
          xnor %g0,%o3,%o3
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ sub %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sllx %o2,3,%o2         // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[-count]
        sub %o1,%o2,%o1         // %o1 = &yptr[-count]
1:        subcc %o2,8,%o2       // Zähler erniedrigen, Pointer erniedrigen
          ldx [%o1+%o2],%o3     // nächstes Digit holen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          and %o4,%o3,%o3       // beide verknüpfen
          xnor %g0,%o3,%o3
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void nor_loop_down (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(nor_loop_down)
C(nor_loop_down:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
        sub %o0,8,%o0
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          or %o3,%o4,%o3        // verknüpfen
          xnor %g0,%o3,%o3
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ sub %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sllx %o2,3,%o2         // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[-count]
        sub %o1,%o2,%o1         // %o1 = &yptr[-count]
1:        subcc %o2,8,%o2       // Zähler erniedrigen, Pointer erniedrigen
          ldx [%o1+%o2],%o3     // nächstes Digit holen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          or %o4,%o3,%o3        // beide verknüpfen
          xnor %g0,%o3,%o3
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void andc2_loop_down (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(andc2_loop_down)
C(andc2_loop_down:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
        sub %o0,8,%o0
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          andn %o3,%o4,%o3      // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ sub %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sllx %o2,3,%o2         // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[-count]
        sub %o1,%o2,%o1         // %o1 = &yptr[-count]
1:        subcc %o2,8,%o2       // Zähler erniedrigen, Pointer erniedrigen
          ldx [%o1+%o2],%o3     // nächstes Digit holen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          andn %o4,%o3,%o3      // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void orc2_loop_down (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(orc2_loop_down)
C(orc2_loop_down:) // Input in %o0,%o1,%o2
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %o1,%o0,%o1        // %o1 = yptr-xptr
        sub %o0,8,%o0
1:        ldx [%o0],%o3         // *xptr
          ldx [%o0+%o1],%o4     // *yptr
          subcc %o2,1,%o2
          orn %o3,%o4,%o3       // verknüpfen
          stx %o3,[%o0]         // =: *xptr
          bne,pt %xcc,1b
         _ sub %o0,8,%o0        // xptr++, yptr++
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sllx %o2,3,%o2         // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[-count]
        sub %o1,%o2,%o1         // %o1 = &yptr[-count]
1:        subcc %o2,8,%o2       // Zähler erniedrigen, Pointer erniedrigen
          ldx [%o1+%o2],%o3     // nächstes Digit holen
          ldx [%o0+%o2],%o4     // noch ein Digit holen
          orn %o4,%o3,%o3       // beide verknüpfen
          bne,pt %xcc,1b
         _ stx %o3,[%o1+%o2]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern void not_loop_down (uintD* xptr, uintC count);
        DECLARE_FUNCTION(not_loop_down)
C(not_loop_down:) // Input in %o0,%o1
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %o0,8,%o0
1:        ldx [%o0],%o2
          subcc %o1,1,%o1
          xnor %g0,%o2,%o2
          stx %o2,[%o0]
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ nop
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sllx %o1,3,%o1         // %o1 = 8*count
        sub %o0,%o1,%o0         // %o0 = &destptr[-count]
1:        subcc %o1,8,%o1       // Zähler erniedrigen, Pointer erniedrigen
          ldx [%o0+%o1],%o2     // nächstes Digit holen
          xnor %g0,%o2,%o2
          bne,pt %xcc,1b
         _ stx %o2,[%o0+%o1]    // Digit ablegen
2:      retl
       _ nop
#endif

// extern boolean and_test_loop_down (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(and_test_loop_down)
C(and_test_loop_down:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,4f
       _ sub %o0,8,%o0
1:        ldx [%o0],%o3
          sub %o1,8,%o1
          ldx [%o1],%o4
          subcc %o2,1,%o2
          be,pn %xcc,3f
         _ andcc %o3,%o4,%g0
          be,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov 1,%o0
3:      bne 2b
       _ nop
4:      retl
       _ mov 0,%o0
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        sllx %o2,3,%o2          // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[-count]
        sub %o1,%o2,%o1         // %o1 = &yptr[-count]
        subcc %o2,8,%o2
        bcs,pn %xcc,2f
       _ nop
          ldx [%o0+%o2],%o3     // nächstes Digit holen
1:        ldx [%o1+%o2],%o4     // noch ein Digit holen
          andcc %o3,%o4,%g0     // beide verknüpfen
          bne,pn %xcc,3f
         _ subcc %o2,8,%o2      // Zähler erniedrigen, Pointer erniedrigen
          bcc,a,pt %xcc,1b
         __ ldx [%o0+%o2],%o3   // nächstes Digit holen
2:      retl
       _ mov 0,%o0
3:      retl
       _ mov 1,%o0
#endif

// extern cl_signean compare_loop_down (uintD* xptr, uintD* yptr, uintC count);
        DECLARE_FUNCTION(compare_loop_down)
C(compare_loop_down:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ nop
1:        ldx [%o0-8],%o3
          ldx [%o1-8],%o4
          subcc %o3,%o4,%g0
          bne,pn %xcc,3f
         _ sub %o0,8,%o0
          subcc %o2,1,%o2
          bne,pn %xcc,1b
         _ sub %o1,8,%o1
2:      retl
       _ mov 0,%o0
3:      mov 1,%o0
        movlu %xcc,-1,%o0
        retl
       _ sra %o0,0,%o0          // sign-extend %o0
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        sllx %o2,3,%o2          // %o2 = 8*count
        sub %o0,%o2,%o0         // %o0 = &xptr[-count]
        sub %o1,%o2,%o1         // %o1 = &yptr[-count]
        subcc %o2,8,%o2
        bcs,pn %xcc,4f
       _ nop
          ldx [%o0+%o2],%o3     // nächstes Digit holen
1:        ldx [%o1+%o2],%o4     // noch ein Digit holen
          subcc %o2,8,%o2       // Zähler erniedrigen, Pointer erniedrigen
          bcs,pn %xcc,3f
         _ subcc %o3,%o4,%g0    // vergleichen
          be,a,pt %xcc,1b
         __ ldx [%o0+%o2],%o3   // nächstes Digit holen
2:      mov 1,%o0
        movlu %xcc,-1,%o0
        retl
       _ sra %o0,0,%o0          // sign-extend %o0
3:      bne 2b
       _ nop
4:      retl
       _ mov 0,%o0
#endif

// extern uintD add_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
        DECLARE_FUNCTION(add_loop_up)
C(add_loop_up:) // Input in %o0,%o1,%o2,%o3, verändert %g1, Output in %o0
#if STANDARD_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %g0,%g1            // Carry := 0
1:        ldx [%o0],%o4         // source1-digit
          add %o0,8,%o0
          ldx [%o1],%o5         // source2-digit
          add %o1,8,%o1
          addcc %o4,%g1,%o4
          movcc %xcc,0,%g1      // %g1|%o4 := %o4 + alter Carry %g1
          addcc %o4,%o5,%o4
          movcs %xcc,1,%g1      // %g1|%o4 := %o4 + alter Carry %g1 + %o5
          stx %o4,[%o2]         // Digit ablegen
          subcc %o3,1,%o3
          bne,pt %xcc,1b
         _ add %o2,8,%o2
2:      retl
       _ mov %g1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %g0,%g1            // Carry := 0
        sub %g0,%o3,%o3         // %o3 = -count
        sllx %o3,3,%o3          // %o3 = -8*count
        sub %o2,8,%o2
        sub %o0,%o3,%o0         // %o0 = &sourceptr1[count]
        sub %o1,%o3,%o1         // %o1 = &sourceptr2[count]
        sub %o2,%o3,%o2         // %o2 = &destptr[count-1]
1:        ldx [%o0+%o3],%o4     // source1-digit
          ldx [%o1+%o3],%o5     // source2-digit
          addcc %o4,%g1,%o4
          movcc %xcc,0,%g1      // %g1|%o4 := %o4 + alter Carry %g1
          addcc %o4,%o5,%o4
          movcs %xcc,1,%g1      // %g1|%o4 := %o4 + alter Carry %g1 + %o5
          addcc %o3,8,%o3       // Zähler erniedrigen, Pointer erhöhen
          bne,pt %xcc,1b
         _ stx %o4,[%o2+%o3]    // Digit ablegen
2:      retl
       _ mov %g1,%o0
#endif

// extern uintD addto_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
        DECLARE_FUNCTION(addto_loop_up)
C(addto_loop_up:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ mov %g0,%o5            // Carry := 0
1:        ldx [%o0],%o3         // source-digit
          add %o0,8,%o0
          ldx [%o1],%o4         // dest-digit
          addcc %o3,%o5,%o3
          movcc %xcc,0,%o5      // %o5|%o3 := %o3 + alter Carry %o5
          addcc %o3,%o4,%o4
          movcs %xcc,1,%o5      // %o5|%o4 := %o3 + alter Carry %o5 + %o4
          stx %o4,[%o1]         // Digit ablegen
          subcc %o2,1,%o2
          bne,pt %xcc,1b
         _ add %o1,8,%o1
2:      retl
       _ mov %o5,%o0
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ mov %g0,%o5            // Carry := 0
        sub %g0,%o2,%o2         // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &sourceptr[count]
        sub %o1,%o2,%o1         // %o1 = &destptr[count]
          ldx [%o0+%o2],%o3     // source-digit
1:        ldx [%o1+%o2],%o4     // dest-digit
          addcc %o3,%o5,%o3
          movcc %xcc,0,%o5      // %o5|%o3 := %o3 + alter Carry %o5
          addcc %o3,%o4,%o4
          movcs %xcc,1,%o5      // %o5|%o4 := %o3 + alter Carry %o5 + %o4
          stx %o4,[%o1+%o2]     // Digit ablegen
          addcc %o2,8,%o2       // Zähler erniedrigen, Pointer erhöhen
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o2],%o3   // source-digit
2:      retl
       _ mov %o5,%o0
#endif

// extern uintD inc_loop_up (uintD* ptr, uintC count);
        DECLARE_FUNCTION(inc_loop_up)
C(inc_loop_up:) // Input in %o0,%o1, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ nop
          ldx [%o0],%o2
1:        add %o0,8,%o0
          addcc %o2,1,%o2
          bne,pn %xcc,3f
         _ stx %o2,[%o0-8]
          subcc %o1,1,%o1
          bne,a,pt %xcc,1b
         __ ldx [%o0],%o2
2:      retl
       _ mov 1,%o0
3:      retl
       _ mov 0,%o0
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %g0,%o1,%o1        // %o1 = -count
        sllx %o1,3,%o1          // %o1 = -8*count
        sub %o0,%o1,%o0         // %o0 = &ptr[count]
          ldx [%o0+%o1],%o2     // digit holen
1:        addcc %o2,1,%o2       // incrementieren
          bne,pn %xcc,3f
         _ stx %o2,[%o0+%o1]    // ablegen
          addcc %o1,8,%o1       // Zähler erniedrigen, Pointer erhöhen
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o1],%o2
2:      retl
       _ mov 1,%o0
3:      retl
       _ mov 0,%o0
#endif

// extern uintD sub_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
        DECLARE_FUNCTION(sub_loop_up)
C(sub_loop_up:) // Input in %o0,%o1,%o2,%o3, verändert %g1, Output in %o0
#if STANDARD_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %g0,%g1            // Carry := 0
1:        ldx [%o0],%o4         // source1-digit
          add %o0,8,%o0
          ldx [%o1],%o5         // source2-digit
          add %o1,8,%o1
          addcc %o5,%g1,%o5
          movcc %xcc,0,%g1      // %g1|%o5 := %o5 + alter Carry %g1
          subcc %o4,%o5,%o4
          movcs %xcc,1,%g1      // %o4-2^64*%g1 := %o4 - %o5 - alter Carry %g1
          stx %o4,[%o2]         // Digit ablegen
          subcc %o3,1,%o3
          bne,pt %xcc,1b
         _ add %o2,8,%o2
2:      retl
       _ mov %g1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %g0,%g1            // Carry := 0
        sub %g0,%o3,%o3         // %o3 = -count
        sllx %o3,3,%o3          // %o3 = -8*count
        sub %o2,8,%o2
        sub %o0,%o3,%o0         // %o0 = &sourceptr1[count]
        sub %o1,%o3,%o1         // %o1 = &sourceptr2[count]
        sub %o2,%o3,%o2         // %o2 = &destptr[count-1]
1:        ldx [%o1+%o3],%o5     // source2-digit
          ldx [%o0+%o3],%o4     // source1-digit
          addcc %o5,%g1,%o5
          movcc %xcc,0,%g1      // %g1|%o5 := %o5 + alter Carry %g1
          subcc %o4,%o5,%o4
          movcs %xcc,1,%g1      // %o4-2^64*%g1 := %o4 - %o5 - alter Carry %g1
          addcc %o3,8,%o3
          bne,pt %xcc,1b
         _ stx %o4,[%o2+%o3]    // Digit ablegen
2:      retl
       _ mov %g1,%o0
#endif

// extern uintD subx_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);
        DECLARE_FUNCTION(subx_loop_up)
C(subx_loop_up:) // Input in %o0,%o1,%o2,%o3,%o4, verändert %g1, Output in %o0
#if STANDARD_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %o4,%g1            // Carry (0 oder -1)
1:        ldx [%o0],%o4         // source1-digit
          add %o0,8,%o0
          ldx [%o1],%o5         // source2-digit
          add %o1,8,%o1
          subcc %o5,%g1,%o5
          movcc %xcc,0,%g1      // %o5-2^64*%g1 := %o5 - alter Carry %g1
          subcc %o4,%o5,%o4
          movcs %xcc,-1,%g1     // %o4+2^64*%g1 := %o4 - %o5 + alter Carry %g1
          stx %o4,[%o2]         // Digit ablegen
          subcc %o3,1,%o3
          bne,pt %xcc,1b
         _ add %o2,8,%o2
2:      retl
       _ mov %g1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o3,0,%o3           // zero-extend %o3 = count
        brz,pn %o3,2f
       _ mov %o4,%g1            // Carry (0 oder -1)
        sub %g0,%o3,%o3         // %o3 = -count
        sllx %o3,3,%o3          // %o3 = -8*count
        sub %o2,8,%o2
        sub %o0,%o3,%o0         // %o0 = &sourceptr1[count]
        sub %o1,%o3,%o1         // %o1 = &sourceptr2[count]
        sub %o2,%o3,%o2         // %o2 = &destptr[count-1]
1:        ldx [%o1+%o3],%o5     // source2-digit
          ldx [%o0+%o3],%o4     // source1-digit
          subcc %o5,%g1,%o5
          movcc %xcc,0,%g1      // %o5-2^64*%g1 := %o5 - alter Carry %g1
          subcc %o4,%o5,%o4
          movcs %xcc,-1,%g1     // %o4+2^64*%g1 := %o4 - %o5 + alter Carry %g1
          addcc %o3,8,%o3
          bne,pt %xcc,1b
         _ stx %o4,[%o2+%o3]    // Digit ablegen
2:      retl
       _ mov %g1,%o0
#endif

// extern uintD subfrom_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
        DECLARE_FUNCTION(subfrom_loop_up)
C(subfrom_loop_up:) // Input in %o0,%o1,%o2, Output in %o0
#if STANDARD_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ mov %g0,%o5            // Carry := 0
1:        ldx [%o0],%o3         // source-digit
          add %o0,8,%o0
          ldx [%o1],%o4         // dest-digit
          addcc %o3,%o5,%o3
          movcc %xcc,0,%o5      // %o5|%o3 := %o3 + alter Carry %o5
          subcc %o4,%o3,%o4
          movcs %xcc,1,%o5      // %o4-2^64*%o5 := %o4 - %o3 - alter Carry %o5
          stx %o4,[%o1]         // Digit ablegen
          subcc %o2,1,%o2
          bne,pt %xcc,1b
         _ add %o1,8,%o1
2:      retl
       _ mov %o5,%o0
#endif
#if COUNTER_LOOPS
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ mov %g0,%o5            // Carry := 0
        sub %g0,%o2,%o2         // %o2 = -count
        sllx %o2,3,%o2          // %o2 = -8*count
        sub %o0,%o2,%o0         // %o0 = &sourceptr[count]
        sub %o1,%o2,%o1         // %o1 = &destptr[count]
          ldx [%o0+%o2],%o3     // source-digit
1:        ldx [%o1+%o2],%o4     // dest-digit
          addcc %o3,%o5,%o3
          movcc %xcc,0,%o5      // %o5|%o3 := %o3 + alter Carry %o5
          subcc %o4,%o3,%o4
          movcs %xcc,1,%o5      // %o4-2^64*%o5 := %o4 - %o3 - alter Carry %o5
          stx %o4,[%o1+%o2]     // Digit ablegen
          addcc %o2,8,%o2
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o2],%o3   // source-digit
2:      retl
       _ mov %o5,%o0
#endif

// extern uintD dec_loop_up (uintD* ptr, uintC count);
        DECLARE_FUNCTION(dec_loop_up)
C(dec_loop_up:) // Input in %o0,%o1, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ nop
          ldx [%o0],%o2
1:        add %o0,8,%o0
          subcc %o2,1,%o2
          bcc,pn %xcc,3f
         _ stx %o2,[%o0-8]
          subcc %o1,1,%o1
          bne,a,pt %xcc,1b
         __ ldx [%o0],%o2
2:      retl
       _ mov -1,%o0
3:      retl
       _ mov 0,%o0
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %g0,%o1,%o1        // %o1 = -count
        sllx %o1,3,%o1          // %o1 = -8*count
        sub %o0,%o1,%o0         // %o0 = &ptr[count]
          ldx [%o0+%o1],%o2     // digit holen
1:        subcc %o2,1,%o2       // decrementieren
          bcc,pn %xcc,3f
         _ stx %o2,[%o0+%o1]    // ablegen
          addcc %o1,8,%o1       // Zähler erniedrigen, Pointer erhöhen
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o1],%o2
2:      retl
       _ mov -1,%o0
3:      retl
       _ mov 0,%o0
#endif

// extern uintD neg_loop_up (uintD* ptr, uintC count);
        DECLARE_FUNCTION(neg_loop_up)
C(neg_loop_up:) // Input in %o0,%o1, Output in %o0
#if STANDARD_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        // erstes Digit /=0 suchen:
        brz,pn %o1,2f
       _ add %o0,8,%o0
1:        ldx [%o0-8],%o2
          subcc %g0,%o2,%o2
          bne,pn %xcc,3f
         _ subcc %o1,1,%o1
          bne,pt %xcc,1b
         _ add %o0,8,%o0
2:      retl
       _ mov 0,%o0
3:      // erstes Digit /=0 gefunden, ab jetzt gibt's Carrys
        // 1 Digit negieren, alle anderen Digits invertieren:
        be,pn %xcc,5f
       _ stx %o2,[%o0-8]
4:        ldx [%o0],%o2
          subcc %o1,1,%o1
          xnor %g0,%o2,%o2
          stx %o2,[%o0]
          bne,pt %xcc,4b
         _ add %o0,8,%o0
5:      retl
       _ mov -1,%o0
#endif
#if COUNTER_LOOPS
//      srl %o1,0,%o1           // zero-extend %o1 = count
        // erstes Digit /=0 suchen:
        brz,pn %o1,2f
       _ sub %g0,%o1,%o1        // %o1 = -count
        sllx %o1,3,%o1          // %o1 = -8*count
        sub %o0,%o1,%o0         // %o0 = &ptr[count]
          ldx [%o0+%o1],%o2     // digit holen
1:        subcc %g0,%o2,%o2     // negieren, testen
          bne,pn %xcc,3f
         _ addcc %o1,8,%o1      // Zähler erniedrigen, Pointer erhöhen
          bne,a,pt %xcc,1b
         __ ldx [%o0+%o1],%o2
2:      retl
       _ mov 0,%o0
3:      // erstes Digit /=0 gefunden, ab jetzt gibt's Carrys
        // alle anderen Digits invertieren:
        sub %o1,8,%o1
        stx %o2,[%o0+%o1]       // ablegen
        addcc %o1,8,%o1
        be,pn %xcc,5f
       _ nop
          ldx [%o0+%o1],%o2
4:        xnor %g0,%o2,%o2
          stx %o2,[%o0+%o1]
          addcc %o1,8,%o1
          bne,a,pt %xcc,4b
         __ ldx [%o0+%o1],%o2
5:      retl
       _ mov -1,%o0
#endif

// extern uintD shift1left_loop_up (uintD* ptr, uintC count);
        DECLARE_FUNCTION(shift1left_loop_up)
C(shift1left_loop_up:) // Input in %o0,%o1, Output in %o0
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ mov 0,%o3              // Carry := 0
1:        ldx [%o0],%o2         // Digit
          addcc %o2,%o2,%o4     // shiften
          add %o4,%o3,%o4       // und carry
          srlx %o2,63,%o3       // neues Carry
          stx %o4,[%o0]         // Digit ablegen
          subcc %o1,1,%o1
          bne,pt %xcc,1b
         _ add %o0,8,%o0
2:      retl
       _ mov %o3,%o0

// extern uintD shiftleft_loop_up (uintD* ptr, uintC count, uintC i, uintD carry);
        DECLARE_FUNCTION(shiftleft_loop_up)
C(shiftleft_loop_up:) // Input in %o0,%o1,%o2,%o3, verändert %g1, Output in %o0
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sub %g0,%o2,%g1        // 64-i (mod 64)
1:        ldx [%o0],%o4         // Digit
          subcc %o1,1,%o1
          sllx %o4,%o2,%o5      // dessen niedere (64-i) Bits
          or %o3,%o5,%o5        // mit dem alten Carry kombinieren
          stx %o5,[%o0]         // Digit ablegen
          srlx %o4,%g1,%o3      // dessen höchste i Bits liefern den neuen Carry
          bne,pt %xcc,1b
         _ add %o0,8,%o0
2:      retl
       _ mov %o3,%o0

#endif

// extern uintD shiftleftcopy_loop_up (uintD* sourceptr, uintD* destptr, uintC count, uintC i);
        DECLARE_FUNCTION(shiftleftcopy_loop_up)
C(shiftleftcopy_loop_up:) // Input in %o0,%o1,%o2,%o3, verändert %g1,%g2, Output in %o0
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ mov 0,%o4              // Carry := 0
        sub %g0,%o3,%g1         // 64-i (mod 64)
1:        ldx [%o0],%o5         // Digit
          subcc %o2,1,%o2
          sllx %o5,%o3,%g2      // dessen niedere (64-i) Bits
          or %o4,%g2,%g2        // mit dem alten Carry kombinieren
          stx %g2,[%o1]         // Digit ablegen
          add %o1,8,%o1
          srlx %o5,%g1,%o4      // dessen höchste i Bits liefern den neuen Carry
          bne,pt %xcc,1b
         _ add %o0,8,%o0
2:      retl
       _ mov %o4,%o0

#if !CL_DS_BIG_ENDIAN_P

// extern uintD shift1right_loop_down (uintD* ptr, uintC count, uintD carry);
        DECLARE_FUNCTION(shift1right_loop_down)
C(shift1right_loop_down:) // Input in %o0,%o1,%o2, Output in %o0
//      srl %o1,0,%o1           // zero-extend %o1 = count
        brz,pn %o1,2f
       _ sllx %o2,63,%o2        // Carry
        sub %o0,8,%o0
1:        ldx [%o0],%o3         // Digit
          subcc %o1,1,%o1
          srlx %o3,1,%o4        // shiften
          or %o2,%o4,%o4        // und mit altem Carry kombinieren
          stx %o4,[%o0]         // und ablegen
          sllx %o3,63,%o2       // neuer Carry
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov %o2,%o0

// extern uintD shiftright_loop_down (uintD* ptr, uintC count, uintC i);
        DECLARE_FUNCTION(shiftright_loop_down)
C(shiftright_loop_down:) // Input in %o0,%o1,%o2, verändert %g1, Output in %o0
//      srl %o1,0,%o1           // zero-extend %o1 = count
        sub %g0,%o2,%g1         // 64-i (mod 64)
        brz,pn %o1,2f
       _ or %g0,%g0,%o3         // Carry := 0
        sub %o0,8,%o0
1:        ldx [%o0],%o4         // Digit
          subcc %o1,1,%o1
          srlx %o4,%o2,%o5      // shiften
          or %o3,%o5,%o5        // und mit altem Carry kombinieren
          stx %o5,[%o0]         // und ablegen
          sllx %o4,%g1,%o3      // neuer Carry
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov %o3,%o0

// extern uintD shiftrightsigned_loop_down (uintD* ptr, uintC count, uintC i);
        DECLARE_FUNCTION(shiftrightsigned_loop_down)
C(shiftrightsigned_loop_down:) // Input in %o0,%o1,%o2, verändert %g1, Output in %o0
//      srl %o1,0,%o1           // zero-extend %o1 = count
        ldx [%o0-8],%o4         // erstes Digit
        sub %g0,%o2,%g1         // 64-i (mod 64)
        srax %o4,%o2,%o5        // shiften
        stx %o5,[%o0-8]         // und ablegen
        sllx %o4,%g1,%o3        // neuer Carry
        subcc %o1,1,%o1
        be,pn %xcc,2f
       _ sub %o0,16,%o0
1:        ldx [%o0],%o4         // Digit
          subcc %o1,1,%o1
          srlx %o4,%o2,%o5      // shiften
          or %o3,%o5,%o5        // und mit altem Carry kombinieren
          stx %o5,[%o0]         // und ablegen
          sllx %o4,%g1,%o3      // neuer Carry
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov %o3,%o0

// extern uintD shiftrightcopy_loop_down (uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry);
        DECLARE_FUNCTION(shiftrightcopy_loop_down)
C(shiftrightcopy_loop_down:) // Input in %o0,%o1,%o2,%o3,%o4, verändert %g1,%g2, Output in %o0
//      srl %o2,0,%o2           // zero-extend %o2 = count
        sub %g0,%o3,%g1         // 64-i (mod 64)
        brz,pn %o2,2f
       _ sllx %o4,%g1,%g2       // erster Carry
          sub %o0,8,%o0
1:        ldx [%o0],%o4         // Digit
          sub %o1,8,%o1
          srlx %o4,%o3,%o5      // shiften
          or %g2,%o5,%o5        // und mit altem Carry kombinieren
          stx %o5,[%o1]         // und ablegen
          sllx %o4,%g1,%g2      // neuer Carry
          subcc %o2,1,%o2
          bne,pt %xcc,1b
         _ sub %o0,8,%o0
2:      retl
       _ mov %g2,%o0

// extern uintD mulusmall_loop_up (uintD digit, uintD* ptr, uintC len, uintD newdigit);
        DECLARE_FUNCTION(mulusmall_loop_up)
C(mulusmall_loop_up:) // Input in %o0,%o1,%o2,%o3, Output in %o0, verändert %g1
//      srl %o2,0,%o2           // zero-extend %o2 = len
        brz,pn %o2,2f
       _ nop
1:        // nächstes Digit [%o1] mit der 6-Bit-Zahl %o0 multiplizieren
          // und kleinen Carry %o3 dazu:
          ldx [%o1],%o4
          sub %o2,1,%o2
          srlx %o4,32,%o5       // high32(x)
          srl %o4,0,%o4         // low32(x)
          mulx %o4,%o0,%o4      // low32(x)*digit
          mulx %o5,%o0,%o5      // high32(x)*digit
          sllx %o5,32,%g1       // low32(high32(x)*digit)*2^32
          add %g1,%o3,%g1       // plus carry
          addcc %o4,%g1,%o4     // plus low32(x)*digit
          srlx %o5,32,%o3       // high32(high32(x)*digit)
          add %o3,1,%g1
          movcs %xcc,%g1,%o3    // neuer Carry
          stx %o4,[%o1]         // neues Digit ablegen
          brnz,pt %o2,1b
         _ add %o1,8,%o1
2:      retl
       _ mov %o3,%o0

// extern void mulu_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
        DECLARE_FUNCTION(mulu_loop_up)
C(mulu_loop_up:) // Input in %i0,%i1,%i2,%i3
        save %sp,-192,%sp
        mov 0,%l0               // Carry
        srlx %i0,32,%l1         // %l1 = high32(digit)
        srl %i0,0,%l2           // %l2 = low32(digit)
        mov 1,%l3
        sllx %l3,32,%l3         // %l3 = 2^32
        sub %i1,%i2,%i1         // %i1 = sourceptr - destptr
1:        ldx [%i1+%i2],%o0     // nächstes Digit
          subcc %i3,1,%i3
          // mit digit multiplizieren: (%l1*2^32+%l2) * %o0 + %l0 -> %l0|%o0
          srlx %o0,32,%o1
          srl %o0,0,%o2
          mulx %l1,%o1,%o3      // high part
          mulx %l1,%o2,%o4      // first mid part
          mulx %l2,%o1,%o1      // second mid part
          mulx %l2,%o2,%o2      // low part
          srlx %o2,32,%o5       // low part's upper half
          add %o4,%o5,%o4       // add to one of the mid parts, no carry
          addcc %o4,%o1,%o4     // add other mid part
          add %o3,%l3,%o5
          movcs %xcc,%o5,%o3    // if carry, add 2^32 to the high part
          srlx %o4,32,%o5
          sllx %o4,32,%o4
          srl %o2,0,%o2
          add %o2,%o4,%o0       // combine low32(midparts) and low32(lowpart)
          addcc %o0,%l0,%o0     // alten Carry addieren
          add %o3,%o5,%l0       // add high32(midparts) to high part
          add %l0,1,%o5
          movcs %xcc,%o5,%l0    // neuer Carry
          // Multiplikation fertig
          stx %o0,[%i2]         // Low-Digit ablegen
          brnz,pt %i3,1b
         _ add %i2,8,%i2
        stx %l0,[%i2]           // letzten Carry ablegen
        ret
       _ restore

// extern uintD muluadd_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
        DECLARE_FUNCTION(muluadd_loop_up)
C(muluadd_loop_up:) // Input in %i0,%i1,%i2,%i3, Output in %i0
        save %sp,-192,%sp
        mov 0,%l0               // Carry
        srlx %i0,32,%l1         // %l1 = high32(digit)
        srl %i0,0,%l2           // %l2 = low32(digit)
        mov 1,%l3
        sllx %l3,32,%l3         // %l3 = 2^32
        sub %i1,%i2,%i1         // %i1 = sourceptr - destptr
1:        ldx [%i1+%i2],%o0     // nächstes Digit
          ldx [%i2],%i4         // *destptr
          subcc %i3,1,%i3
          // mit digit multiplizieren: (%l1*2^32+%l2) * %o0 + %l0 -> %l0|%o0
          srlx %o0,32,%o1
          srl %o0,0,%o2
          mulx %l1,%o1,%o3      // high part
          mulx %l1,%o2,%o4      // first mid part
          mulx %l2,%o1,%o1      // second mid part
          mulx %l2,%o2,%o2      // low part
          srlx %o2,32,%o5       // low part's upper half
          add %o4,%o5,%o4       // add to one of the mid parts, no carry
          addcc %o4,%o1,%o4     // add other mid part
          add %o3,%l3,%o5
          movcs %xcc,%o5,%o3    // if carry, add 2^32 to the high part
          srlx %o4,32,%o5
          sllx %o4,32,%o4
          srl %o2,0,%o2
          add %o2,%o4,%o0       // combine low32(midparts) and low32(lowpart)
          addcc %o0,%l0,%o0     // alten Carry addieren
          add %o3,%o5,%l0       // add high32(midparts) to high part
          add %l0,1,%o5
          movcs %xcc,%o5,%l0    // neuer Carry
          // Multiplikation fertig
          addcc %i4,%o0,%o0     // alten *destptr addieren
          add %l0,1,%o2
          movcs %xcc,%o2,%l0    // neuer Carry
          stx %o0,[%i2]         // Low-Digit ablegen
          brnz,pt %i3,1b
         _ add %i2,8,%i2
        mov %l0,%i0             // letzter Carry
        ret
       _ restore

// extern uintD mulusub_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
        DECLARE_FUNCTION(mulusub_loop_up)
C(mulusub_loop_up:) // Input in %i0,%i1,%i2,%i3, Output in %i0
        save %sp,-192,%sp
        mov 0,%l0               // Carry
        srlx %i0,32,%l1         // %l1 = high32(digit)
        srl %i0,0,%l2           // %l2 = low32(digit)
        mov 1,%l3
        sllx %l3,32,%l3         // %l3 = 2^32
        sub %i1,%i2,%i1         // %i1 = sourceptr - destptr
1:        ldx [%i1+%i2],%o0     // nächstes Digit
          ldx [%i2],%i4         // *destptr
          subcc %i3,1,%i3
          // mit digit multiplizieren: (%l1*2^32+%l2) * %o0 + %l0 -> %l0|%o0
          srlx %o0,32,%o1
          srl %o0,0,%o2
          mulx %l1,%o1,%o3      // high part
          mulx %l1,%o2,%o4      // first mid part
          mulx %l2,%o1,%o1      // second mid part
          mulx %l2,%o2,%o2      // low part
          srlx %o2,32,%o5       // low part's upper half
          add %o4,%o5,%o4       // add to one of the mid parts, no carry
          addcc %o4,%o1,%o4     // add other mid part
          add %o3,%l3,%o5
          movcs %xcc,%o5,%o3    // if carry, add 2^32 to the high part
          srlx %o4,32,%o5
          sllx %o4,32,%o4
          srl %o2,0,%o2
          add %o2,%o4,%o0       // combine low32(midparts) and low32(lowpart)
          addcc %o0,%l0,%o0     // alten Carry addieren
          add %o3,%o5,%l0       // add high32(midparts) to high part
          add %l0,1,%o5
          movcs %xcc,%o5,%l0    // neuer Carry
          // Multiplikation fertig
          subcc %i4,%o0,%o0     // vom alten *destptr subtrahieren
          add %l0,1,%o2
          movcs %xcc,%o2,%l0    // neuer Carry
          stx %o0,[%i2]         // Low-Digit ablegen
          brnz,pt %i3,1b
         _ add %i2,8,%i2
        mov %l0,%i0             // letzter Carry
        ret
       _ restore

#endif

// extern void shiftxor_loop_up (uintD* xptr, const uintD* yptr, uintC count, uintC i);
        DECLARE_FUNCTION(shiftxor_loop_up)
C(shiftxor_loop_up:) // Input in %o0,%o1,%o2,%o3, verändert %g1,%g2
//      srl %o2,0,%o2           // zero-extend %o2 = count
        brz,pn %o2,2f
       _ sub %g0,%o3,%g1        // 64-i (mod 64)
        sub %o1,%o0,%o1
        ldx [%o0],%o4           // *xptr holen
1:        ldx [%o0+%o1],%o5     // *yptr holen
          subcc %o2,1,%o2
          sllx %o5,%o3,%g2      // dessen niedere (64-i) Bits
          xor %o4,%g2,%o4       // mit dem modifizierten *xptr kombinieren
          stx %o4,[%o0]         // und ablegen
          add %o0,8,%o0
          srlx %o5,%g1,%g2      // höchste i Bits von *yptr
          ldx [%o0],%o4         // schon mal mit dem nächsten *xptr
          bne,pt %xcc,1b
         _ xor %o4,%g2,%o4      // verknüpfen
        stx %o4,[%o0]           // und ablegen
2:      retl
       _ nop

