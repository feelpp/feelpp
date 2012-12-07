// Externe Routinen zu ARILEV1.D
// Prozessor: 80386 im native mode
// Assembler-Syntax: GNU oder SUN, Moves von links nach rechts
// Compiler: GNU-C oder SUN-C
// Parameter-Übergabe: auf dem Stack 4(%esp),8(%esp),...
// Register: %eax,%edx,%ecx dürfen stets verändert werden, alles andere retten.
// Ergebnis-Übergabe: in %eax
// Einstellungen: intCsize=32, intDsize=32.

// Bruno Haible 14.8.1992
// Zum Teil abgeschrieben von Bernhard Degels "v-i386.s"

  #ifdef ASM_UNDERSCORE
    #if defined(__STDC__) || defined (__cplusplus)
      #define C(entrypoint) _##entrypoint
    #else
      #define C(entrypoint) _/**/entrypoint
    #endif
  #else
    #define C(entrypoint) entrypoint
  #endif
  #ifdef ASM_UNDERSCORE
    #if defined(__STDC__) || defined (__cplusplus)
      #define L(label) L##label
    #else
      #define L(label) L/**/label
    #endif
  #else
    #if defined(__STDC__) || defined (__cplusplus)
      #define L(label) .L##label
    #else
      #define L(label) .L/**/label
    #endif
  #endif
  #if defined(ASM_UNDERSCORE) || defined(COHERENT) /* defined(__EMX__) || defined(__GO32__) || defined(linux) || defined(__386BSD__) || defined(__NetBSD__) || defined(__FreeBSD__) || defined(COHERENT) || ... */
    // GNU-Assembler oder MWC-Assembler
    #define repz     repe
    #define shcl     %cl,
  #else /* defined(sun) || ... */
    // SUN-Assembler oder Consensys-Assembler
    #define jecxz    orl %ecx,%ecx ; jz
    #define shcl
  #endif
  #if defined(__EMX__)
    // Direction-Flag ist defaultmäßig gelöscht
    #define dir0start
    #define dir0end
    #define dir1start  std
    #define dir1end    cld
  #elif 1
    // Wir gehen auf Nummer sicher.
    #define dir0start  cld
    #define dir0end
    #define dir1start  std
    #define dir1end    cld
  #else
    // Direction-Flag darf nach Belieben modifiziert werden
    #define dir0start  cld
    #define dir0end
    #define dir1start  std
    #define dir1end
  #endif
  // Alignment. Note that some assemblers need ".align 3,0x90" whereas other
  // assemblers don't like this syntax. So we put in the "nop"s by hand.
  #if defined(ASM_UNDERSCORE) && !(defined(__CYGWIN32__) || defined(__MINGW32__))
    // BSD syntax assembler
    #define ALIGN  .align 3
  #else
    // ELF syntax assembler
    #define ALIGN  .align 8
  #endif
  // When this file is compiled into a shared library, ELF linkers need to
  // know which symbols are functions.
  #if defined(__svr4__) || defined(__ELF__) || defined(__NetBSD__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__ROSE__) || defined(_SEQUENT_) || defined(DGUX) || defined(_SCO_COFF) || defined(_SCO_ELF)
    #define DECLARE_FUNCTION(name) .type C(name),@function
  #else
    #define DECLARE_FUNCTION(name)
  #endif

            .text

            .globl C(copy_loop_up)
            .globl C(copy_loop_down)
            .globl C(fill_loop_up)
            .globl C(fill_loop_down)
            .globl C(clear_loop_up)
            .globl C(clear_loop_down)
            .globl C(test_loop_up)
            .globl C(test_loop_down)
            .globl C(xor_loop_up)
            .globl C(compare_loop_up)
            .globl C(shiftleftcopy_loop_up)
            .globl C(shiftxor_loop_up)
#if CL_DS_BIG_ENDIAN_P
            .globl C(or_loop_up)
            .globl C(and_loop_up)
            .globl C(eqv_loop_up)
            .globl C(nand_loop_up)
            .globl C(nor_loop_up)
            .globl C(andc2_loop_up)
            .globl C(orc2_loop_up)
            .globl C(not_loop_up)
            .globl C(and_test_loop_up)
            .globl C(add_loop_down)
            .globl C(addto_loop_down)
            .globl C(inc_loop_down)
            .globl C(sub_loop_down)
            .globl C(subx_loop_down)
            .globl C(subfrom_loop_down)
            .globl C(dec_loop_down)
            .globl C(neg_loop_down)
            .globl C(shift1left_loop_down)
            .globl C(shiftleft_loop_down)
            .globl C(shiftleftcopy_loop_down)
            .globl C(shift1right_loop_up)
            .globl C(shiftright_loop_up)
            .globl C(shiftrightsigned_loop_up)
            .globl C(shiftrightcopy_loop_up)
            .globl C(mulusmall_loop_down)
            .globl C(mulu_loop_down)
            .globl C(muluadd_loop_down)
            .globl C(mulusub_loop_down)
            .globl C(divu_loop_up)
            .globl C(divucopy_loop_up)
#else
            .globl C(or_loop_down)
            .globl C(xor_loop_down)
            .globl C(and_loop_down)
            .globl C(eqv_loop_down)
            .globl C(nand_loop_down)
            .globl C(nor_loop_down)
            .globl C(andc2_loop_down)
            .globl C(orc2_loop_down)
            .globl C(not_loop_down)
            .globl C(and_test_loop_down)
            .globl C(compare_loop_down)
            .globl C(add_loop_up)
            .globl C(addto_loop_up)
            .globl C(inc_loop_up)
            .globl C(sub_loop_up)
            .globl C(subx_loop_up)
            .globl C(subfrom_loop_up)
            .globl C(dec_loop_up)
            .globl C(neg_loop_up)
            .globl C(shift1left_loop_up)
            .globl C(shiftleft_loop_up)
            .globl C(shift1right_loop_down)
            .globl C(shiftright_loop_down)
            .globl C(shiftrightsigned_loop_down)
            .globl C(shiftrightcopy_loop_down)
            .globl C(mulusmall_loop_up)
            .globl C(mulu_loop_up)
            .globl C(muluadd_loop_up)
            .globl C(mulusub_loop_up)
            .globl C(divu_loop_down)
            .globl C(divucopy_loop_down)
#endif

#ifndef __GNUC__ /* mit GNU-C machen wir mulu32() als Macro, der inline multipliziert */

// extern struct { uint32 lo; uint32 hi; } mulu32_ (uint32 arg1, uint32 arg2);
// 2^32*hi+lo := arg1*arg2.
            .globl C(mulu32_)
            ALIGN
            DECLARE_FUNCTION(mulu32_)
C(mulu32_:)
            movl    4(%esp),%eax    // arg1
            mull    8(%esp)         // %edx|%eax := arg1 * arg2
            movl    %edx,C(mulu32_high) // %edx = hi abspeichern
            ret                     // %eax = lo als Ergebnis

#endif

#ifndef __GNUC__ /* mit GNU-C machen wir divu_6432_3232() als Macro, der inline dividiert */

// extern struct { uint32 q; uint32 r; } divu_6432_3232_ (uint32 xhi, uint32 xlo, uint32 y);
// x = 2^32*xhi+xlo = q*y+r schreiben. Sei bekannt, daß 0 <= x < 2^32*y .
            .globl C(divu_6432_3232_)
            ALIGN
            DECLARE_FUNCTION(divu_6432_3232_)
C(divu_6432_3232_:)
            movl    4(%esp),%edx
            movl    8(%esp),%eax
            divl    12(%esp)       // x = %edx|%eax durch dividieren
            movl    %edx,C(divu_32_rest) // Rest %edx = r abspeichern
            ret                    // Quotient %eax = q als Ergebnis

#endif

// extern uintD* copy_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(copy_loop_up)
C(copy_loop_up:)
            movl    %edi,%edx       // %edi retten
            movl    %esi,%eax       // %esi retten
            movl    4(%esp),%esi    // %esi = sourceptr
            movl    8(%esp),%edi    // %edi = destptr
            movl    12(%esp),%ecx   // %ecx = count
            dir0start
            rep
              movsl                 // %ecx mal aufwärts (%edi) := (%esi)
            dir0end
            movl    %eax,%esi       // %esi zurück
            movl    %edi,%eax       // %edi als Ergebnis
            movl    %edx,%edi       // %edi zurück
            ret

// extern uintD* copy_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(copy_loop_down)
C(copy_loop_down:)
            movl    %edi,%edx       // %edi retten
            movl    %esi,%eax       // %esi retten
            movl    4(%esp),%esi    // %esi = sourceptr
            movl    8(%esp),%edi    // %edi = destptr
            movl    12(%esp),%ecx   // %ecx = count
            leal    -4(%esi),%esi
            leal    -4(%edi),%edi
            dir1start
            rep
              movsl                 // %ecx mal abwärts (%edi) := (%esi)
            dir1end
            movl    %eax,%esi       // %esi zurück
            leal    4(%edi),%eax    // %edi als Ergebnis
            movl    %edx,%edi       // %edi zurück
            ret

// extern uintD* fill_loop_up (uintD* destptr, uintC count, uintD filler);
            ALIGN
            DECLARE_FUNCTION(fill_loop_up)
C(fill_loop_up:)
            movl    %edi,%edx       // %edi retten
            movl    4(%esp),%edi    // %edi = destptr
            movl    8(%esp),%ecx    // %ecx = count
            movl    12(%esp),%eax   // %eax = filler
            dir0start
            rep
              stosl                 // %ecx mal aufwärts (%edi) := %eax
            dir0end
            movl    %edi,%eax       // %edi als Ergebnis
            movl    %edx,%edi       // %edi zurück
            ret

// extern uintD* fill_loop_down (uintD* destptr, uintC count, uintD filler);
            ALIGN
            DECLARE_FUNCTION(fill_loop_down)
C(fill_loop_down:)
            movl    %edi,%edx       // %edi retten
            movl    4(%esp),%edi    // %edi = destptr
            movl    8(%esp),%ecx    // %ecx = count
            movl    12(%esp),%eax   // %eax = filler
            leal    -4(%edi),%edi
            dir1start
            rep
              stosl                 // %ecx mal abwärts (%edi) := %eax
            dir1end
            leal    4(%edi),%eax    // %edi als Ergebnis
            movl    %edx,%edi       // %edi zurück
            ret

// extern uintD* clear_loop_up (uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(clear_loop_up)
C(clear_loop_up:)
            movl    %edi,%edx       // %edi retten
            movl    4(%esp),%edi    // %edi = destptr
            movl    8(%esp),%ecx    // %ecx = count
            xorl    %eax,%eax       // %eax = 0
            dir0start
            rep
              stosl                 // %ecx mal aufwärts (%edi) := %eax
            dir0end
            movl    %edi,%eax       // %edi als Ergebnis
            movl    %edx,%edi       // %edi zurück
            ret

// extern uintD* clear_loop_down (uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(clear_loop_down)
C(clear_loop_down:)
            movl    %edi,%edx       // %edi retten
            movl    4(%esp),%edi    // %edi = destptr
            movl    8(%esp),%ecx    // %ecx = count
            leal    -4(%edi),%edi
            xorl    %eax,%eax       // %eax = 0
            dir1start
            rep
              stosl                 // %ecx mal abwärts (%edi) := %eax
            dir1end
            leal    4(%edi),%eax    // %edi als Ergebnis
            movl    %edx,%edi       // %edi zurück
            ret

// extern boolean test_loop_up (uintD* ptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(test_loop_up)
C(test_loop_up:)
            movl    %edi,%edx       // %edi retten
            movl    4(%esp),%edi    // %edi = ptr
            movl    8(%esp),%ecx    // %ecx = count
            xorl    %eax,%eax       // %eax = 0
            dir0start
            repz                    // Falls %ecx > 0:
              scasl                 // %ecx mal aufwärts (%edi) testen
                                    // und weiterschleifen, falls Z, d.h. (%edi)=0.
            dir0end
            // Noch ist %eax = 0.
            jz      L(tlu1)         // alles =0 -> Ergebnis 0
            incl    %eax            // Ergebnis 1
L(tlu1:)    movl    %edx,%edi       // %edi zurück
            ret

// extern boolean test_loop_down (uintD* ptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(test_loop_down)
C(test_loop_down:)
            movl    %edi,%edx       // %edi retten
            movl    4(%esp),%edi    // %edi = ptr
            movl    8(%esp),%ecx    // %ecx = count
            xorl    %eax,%eax       // %eax = 0
            leal    -4(%edi),%edi
            dir1start
            repz                    // Falls %ecx > 0:
              scasl                 // %ecx mal aufwärts (%edi) testen
                                    // und weiterschleifen, falls Z, d.h. (%edi)=0.
            dir1end
            // Noch ist %eax = 0.
            jz      L(tld1)         // alles =0 -> Ergebnis 0
            incl    %eax            // Ergebnis 1
L(tld1:)    movl    %edx,%edi       // %edi zurück
            ret

#if CL_DS_BIG_ENDIAN_P

// extern void or_loop_up (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(or_loop_up)
C(or_loop_up:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(olu2)         // %ecx = 0 ?
L(olu1:)      movl    (%edx,%esi),%eax // *yptr
              orl     %eax,(%edx)      // *xptr |= ...
              leal    4(%edx),%edx     // xptr++, yptr++
              decl    %ecx
              jnz     L(olu1)
L(olu2:)    popl    %esi            // %esi zurück
            ret

#endif

// extern void xor_loop_up (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(xor_loop_up)
C(xor_loop_up:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(xlu2)         // %ecx = 0 ?
L(xlu1:)      movl    (%edx,%esi),%eax // *yptr
              xorl    %eax,(%edx)      // *xptr ^= ...
              leal    4(%edx),%edx     // xptr++, yptr++
              decl    %ecx
              jnz     L(xlu1)
L(xlu2:)    popl    %esi            // %esi zurück
            ret

#if CL_DS_BIG_ENDIAN_P

// extern void and_loop_up (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(and_loop_up)
C(and_loop_up:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(alu2)         // %ecx = 0 ?
L(alu1:)      movl    (%edx,%esi),%eax // *yptr
              andl    %eax,(%edx)      // *xptr &= ...
              leal    4(%edx),%edx     // xptr++, yptr++
              decl    %ecx
              jnz     L(alu1)
L(alu2:)    popl    %esi            // %esi zurück
            ret

// extern void eqv_loop_up (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(eqv_loop_up)
C(eqv_loop_up:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(elu2)         // %ecx = 0 ?
L(elu1:)      movl    (%edx),%eax      // *xptr
              xorl    (%edx,%esi),%eax // ^ *yptr
              notl    %eax             // ~(...)
              movl    %eax,(%edx)      // =: *xptr
              leal    4(%edx),%edx     // xptr++, yptr++
              decl    %ecx
              jnz     L(elu1)
L(elu2:)    popl    %esi            // %esi zurück
            ret

// extern void nand_loop_up (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(nand_loop_up)
C(nand_loop_up:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(nalu2)        // %ecx = 0 ?
L(nalu1:)     movl    (%edx),%eax      // *xptr
              andl    (%edx,%esi),%eax // & *yptr
              notl    %eax             // ~(...)
              movl    %eax,(%edx)      // =: *xptr
              leal    4(%edx),%edx     // xptr++, yptr++
              decl    %ecx
              jnz     L(nalu1)
L(nalu2:)   popl    %esi            // %esi zurück
            ret

// extern void nor_loop_up (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(nor_loop_up)
C(nor_loop_up:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(nolu2)        // %ecx = 0 ?
L(nolu1:)     movl    (%edx),%eax      // *xptr
              orl     (%edx,%esi),%eax // | *yptr
              notl    %eax             // ~(...)
              movl    %eax,(%edx)      // =: *xptr
              leal    4(%edx),%edx     // xptr++, yptr++
              decl    %ecx
              jnz     L(nolu1)
L(nolu2:)   popl    %esi            // %esi zurück
            ret

// extern void andc2_loop_up (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(andc2_loop_up)
C(andc2_loop_up:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(aclu2)        // %ecx = 0 ?
L(aclu1:)     movl    (%edx,%esi),%eax // *yptr
              notl    %eax             // ~ *yptr
              andl    %eax,(%edx)      // *xptr &= ...
              leal    4(%edx),%edx     // xptr++, yptr++
              decl    %ecx
              jnz     L(aclu1)
L(aclu2:)   popl    %esi            // %esi zurück
            ret

// extern void orc2_loop_up (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(orc2_loop_up)
C(orc2_loop_up:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(oclu2)        // %ecx = 0 ?
L(oclu1:)     movl    (%edx,%esi),%eax // *yptr
              notl    %eax             // ~ *yptr
              orl     %eax,(%edx)      // *xptr |= ...
              leal    4(%edx),%edx     // xptr++, yptr++
              decl    %ecx
              jnz     L(oclu1)
L(oclu2:)   popl    %esi            // %esi zurück
            ret

// extern void not_loop_up (uintD* xptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(not_loop_up)
C(not_loop_up:)
            movl    4(%esp),%edx    // %edx = xptr
            movl    8(%esp),%ecx    // %ecx = count
            jecxz   L(nlu2)         // %ecx = 0 ?
            nop ; nop ; nop ; nop ; nop ; nop
L(nlu1:)      notl    (%edx)           // ~= *xptr
              leal    4(%edx),%edx     // xptr++
              decl    %ecx
              jnz     L(nlu1)
L(nlu2:)    ret

// extern boolean and_test_loop_up (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(and_test_loop_up)
C(and_test_loop_up:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            jecxz   L(atlu2)        // %ecx = 0 ?
            subl    %edx,%esi
L(atlu1:)     movl    (%edx,%esi),%eax // *yptr
              andl    (%edx),%eax      // *xptr & ...
              jnz     L(atlu3)
              leal    4(%edx),%edx     // xptr++, yptr++
              decl    %ecx
              jnz     L(atlu1)
L(atlu2:)   xorl    %eax,%eax       // Ergebnis 0
            popl    %esi            // %esi zurück
            ret
L(atlu3:)   movl    $1,%eax         // Ergebnis 1 (nicht irgendwas /=0 !)
            popl    %esi            // %esi zurück
            ret

#endif

// extern cl_signean compare_loop_up (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(compare_loop_up)
C(compare_loop_up:)
            movl    %esi,%edx       // %esi retten
            movl    %edi,%eax       // %edi retten
            movl    4(%esp),%esi    // %esi = xptr
            movl    8(%esp),%edi    // %edi = yptr
            movl    12(%esp),%ecx   // %ecx = count
            cmpl    %ecx,%ecx       // initialize flags for the case %ecx is 0
            dir0start
            repz                    // Falls %ecx > 0:
              cmpsl                 // %ecx mal aufwärts (%edi) und (%esi) vergleichen
                                    // und weiterschleifen, falls Z, d.h. (%edi)=(%esi).
            dir0end
            // Flags -> Ergebnis:
            // Z,NC -> bis zum Schluß (%esi)-(%edi) = 0 -> x=y -> Ergebnis 0
            // NZ,C -> schließlich (%esi)-(%edi) < 0 -> x<y -> Ergebnis -1
            // NZ,NC -> schließlich (%esi)-(%edi) > 0 -> x>y -> Ergebnis +1
            movl    %eax,%edi       // %edi zurück
            movl    %edx,%esi       // %esi zurück
            jbe     L(cmlu1)        // "be" = Z oder C
            movl    $1,%eax         // Ergebnis +1
            ret
L(cmlu1:)   sbbl    %eax,%eax       // Ergebnis -1 (falls C) oder 0 (falls NC)
            ret

#if CL_DS_BIG_ENDIAN_P

// extern uintD add_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(add_loop_down)
C(add_loop_down:)
            pushl   %esi            // %esi retten
            pushl   %edi            // %edi retten
            movl    12(%esp),%edx   // %edx = sourceptr1
            movl    16(%esp),%esi   // %esi = sourceptr2
            movl    20(%esp),%edi   // %edi = destptr
            movl    24(%esp),%ecx   // %ecx = count
            subl    %edi,%edx
            subl    %edi,%esi
            orl     %ecx,%ecx       // %ecx = 0 ?, Carry löschen
            jz      L(ald2)
L(ald1:)      leal    -4(%edi),%edi   // sourceptr1--, sourceptr2--, destptr--
              movl    (%edx,%edi),%eax // *sourceptr1
              adcl    (%esi,%edi),%eax // + *sourceptr2 + carry
              movl    %eax,(%edi)     // =: *destptr, neuen Carry behalten
              decl    %ecx
              jnz     L(ald1)
L(ald2:)    sbbl    %eax,%eax      // Ergebnis := - Carry
            popl    %edi           // %edi zurück
            popl    %esi           // %esi zurück
            ret

// extern uintD addto_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(addto_loop_down)
C(addto_loop_down:)
            pushl   %edi            // %edi retten
            movl    8(%esp),%edx    // %edx = sourceptr
            movl    12(%esp),%edi   // %edi = destptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edi,%edx
            orl     %ecx,%ecx       // %ecx = 0 ?, Carry löschen
            jz      L(atld2)
L(atld1:)     leal    -4(%edi),%edi   // sourceptr--, destptr--
              movl    (%edx,%edi),%eax // *sourceptr
              adcl    %eax,(%edi)     // + *destptr + carry =: *destptr, neuer Carry
              decl    %ecx
              jnz     L(atld1)
L(atld2:)   sbbl    %eax,%eax       // Ergebnis := - Carry
            popl    %edi            // %edi zurück
            ret

// extern uintD inc_loop_down (uintD* ptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(inc_loop_down)
C(inc_loop_down:)
            movl    4(%esp),%edx    // %edx = ptr
            movl    8(%esp),%ecx    // %ecx = count
            jecxz   L(ild2)         // %ecx = 0 ?
L(ild1:)      leal    -4(%edx),%edx
              addl    $1,(%edx)       // (*ptr)++
              jnc     L(ild3)         // kein Carry -> fertig
              decl    %ecx
              jnz     L(ild1)
L(ild2:)    movl    $1,%eax         // Ergebnis := 1
            ret
L(ild3:)    xorl    %eax,%eax       // Ergebnis := 0
            ret

// extern uintD sub_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(sub_loop_down)
C(sub_loop_down:)
            pushl   %esi            // %esi retten
            pushl   %edi            // %edi retten
            movl    12(%esp),%edx   // %edx = sourceptr1
            movl    16(%esp),%esi   // %esi = sourceptr2
            movl    20(%esp),%edi   // %edi = destptr
            movl    24(%esp),%ecx   // %ecx = count
            subl    %edi,%edx
            subl    %edi,%esi
            orl     %ecx,%ecx       // %ecx = 0 ?, Carry löschen
            jz      L(sld2)
L(sld1:)      leal    -4(%edi),%edi   // sourceptr1--, sourceptr2--, destptr--
              movl    (%edx,%edi),%eax // *sourceptr1
              sbbl    (%esi,%edi),%eax // - *sourceptr2 - carry
              movl    %eax,(%edi)     // =: *destptr, neuen Carry behalten
              decl    %ecx
              jnz     L(sld1)
L(sld2:)    sbbl    %eax,%eax      // Ergebnis := - Carry
            popl    %edi           // %edi zurück
            popl    %esi           // %esi zurück
            ret

// extern uintD subx_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);
            ALIGN
            DECLARE_FUNCTION(subx_loop_down)
C(subx_loop_down:)
            pushl   %esi            // %esi retten
            pushl   %edi            // %edi retten
            movl    12(%esp),%edx   // %edx = sourceptr1
            movl    16(%esp),%esi   // %esi = sourceptr2
            movl    20(%esp),%edi   // %edi = destptr
            movl    24(%esp),%ecx   // %ecx = count
            jecxz   L(sxld2)        // %ecx = 0 ?
            subl    %edi,%edx
            subl    %edi,%esi
            movl    28(%esp),%eax   // carry, 0 oder -1
            addl    %eax,%eax       // Bit 31 davon in den Carry
            nop ; nop
L(sxld1:)     leal    -4(%edi),%edi   // sourceptr1--, sourceptr2--, destptr--
              movl    (%edx,%edi),%eax // *sourceptr1
              sbbl    (%esi,%edi),%eax // - *sourceptr2 - carry
              movl    %eax,(%edi)     // =: *destptr, neuen Carry behalten
              decl    %ecx
              jnz     L(sxld1)
            sbbl    %eax,%eax      // Ergebnis := - Carry
            popl    %edi           // %edi zurück
            popl    %esi           // %esi zurück
            ret
L(sxld2:)   movl    28(%esp),%eax  // Ergebnis := carry
            popl    %edi           // %edi zurück
            popl    %esi           // %esi zurück
            ret

// extern uintD subfrom_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(subfrom_loop_down)
C(subfrom_loop_down:)
            pushl   %edi            // %edi retten
            movl    8(%esp),%edx    // %edx = sourceptr
            movl    12(%esp),%edi   // %edi = destptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edi,%edx
            orl     %ecx,%ecx       // %ecx = 0 ?, Carry löschen
            jz      L(sfld2)
L(sfld1:)     leal    -4(%edi),%edi   // sourceptr--, destptr--
              movl    (%edx,%edi),%eax // *sourceptr
              sbbl    %eax,(%edi)     // *destptr - *sourceptr - carry =: *destptr, neuer Carry
              decl    %ecx
              jnz     L(sfld1)
L(sfld2:)   sbbl    %eax,%eax       // Ergebnis := - Carry
            popl    %edi            // %edi zurück
            ret

// extern uintD dec_loop_down (uintD* ptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(dec_loop_down)
C(dec_loop_down:)
            movl    4(%esp),%edx    // %edx = ptr
            movl    8(%esp),%ecx    // %ecx = count
            jecxz   L(dld2)         // %ecx = 0 ?
L(dld1:)      leal    -4(%edx),%edx
              subl    $1,(%edx)       // (*ptr)--
              jnc     L(dld3)         // kein Carry -> fertig
              decl    %ecx
              jnz     L(dld1)
L(dld2:)    movl    $-1,%eax        // Ergebnis := -1
            ret
L(dld3:)    xorl    %eax,%eax       // Ergebnis := 0
            ret

// extern uintD neg_loop_down (uintD* ptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(neg_loop_down)
C(neg_loop_down:)
            movl    4(%esp),%edx    // %edx = ptr
            movl    8(%esp),%ecx    // %ecx = count
            // erstes Digit /=0 suchen:
            jecxz   L(nld2)         // %ecx = 0 ?
L(nld1:)      leal    -4(%edx),%edx
              negl    (%edx)
              jnz     L(nld3)
              decl    %ecx
              jnz     L(nld1)
L(nld2:)    xorl    %eax,%eax       // Ergebnis := 0
            ret
            nop ; nop ; nop ; nop ; nop ; nop
L(nld3:)    // erstes Digit /=0 gefunden, ab jetzt gibt's Carrys
            // alle anderen Digits invertieren:
            decl    %ecx
            jz      L(nld5)
L(nld4:)      leal    -4(%edx),%edx
              notl    (%edx)
              decl    %ecx
              jnz     L(nld4)
L(nld5:)    movl    $-1,%eax        // Ergebnis := -1
            ret

// extern uintD shift1left_loop_down (uintD* ptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(shift1left_loop_down)
C(shift1left_loop_down:)
            movl    4(%esp),%edx    // %edx = ptr
            movl    8(%esp),%ecx    // %ecx = count
            orl     %ecx,%ecx       // %ecx = 0 ?, Carry löschen
            jz      L(s1lld2)
            nop ; nop ; nop ; nop
L(s1lld1:)    leal    -4(%edx),%edx   // ptr--
              rcll    $1,(%edx)       // *ptr und Carry um 1 Bit links rotieren
              decl    %ecx
              jnz     L(s1lld1)
L(s1lld2:)  sbbl    %eax,%eax       // Ergebnis := - Carry
            ret

// extern uintD shiftleft_loop_down (uintD* ptr, uintC count, uintC i, uintD carry);
            ALIGN
            DECLARE_FUNCTION(shiftleft_loop_down)
C(shiftleft_loop_down:)
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    12(%esp),%edi   // %edi = ptr
            movl    16(%esp),%edx   // %edx = count
            movb    20(%esp),%cl    // %cl = i
            orl     %edx,%edx       // count = 0 ?
            jz      L(slld4)
            // erstes Digit shiften:
            leal    -4(%edi),%edi
            movl    (%edi),%eax     // Digit in %eax halten
            movl    %eax,%ebx       // und in %ebx rechnen:
            shll    %cl,%ebx        // um i Bits links shiften
            orl     24(%esp),%ebx   // und die unteren i Bits eintragen
            movl    %ebx,(%edi)     // und wieder ablegen
            // Letztes Digit in %eax.
            decl    %edx
            jz      L(slld2)
            nop ; nop ; nop ; nop
L(slld1:)     // weiteres Digit shiften:
              leal    -4(%edi),%edi
              movl    (%edi),%ebx
              shldl   shcl %eax,(%edi) // (%edi) um %cl=i Bits links shiften, %eax von rechts reinshiften
              // Letztes Digit in %ebx.
              decl    %edx
              jz      L(slld3)
              // weiteres Digit shiften:
              leal    -4(%edi),%edi
              movl    (%edi),%eax
              shldl   shcl %ebx,(%edi) // (%edi) um %cl=i Bits links shiften, %ebx von rechts reinshiften
              // Letztes Digit in %eax.
              decl    %edx
              jnz     L(slld1)
L(slld2:)   movl    %eax,%ebx
L(slld3:)   xorl    %eax,%eax       // %eax := 0
            shldl   shcl %ebx,%eax  // %eax := höchste %cl=i Bits von %ebx
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret
L(slld4:)   movl    24(%esp),%eax   // %eax := carry
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret

// extern uintD shiftleftcopy_loop_down (uintD* sourceptr, uintD* destptr, uintC count, uintC i);
            ALIGN
            DECLARE_FUNCTION(shiftleftcopy_loop_down)
C(shiftleftcopy_loop_down:)
            pushl   %esi            // %esi retten
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    16(%esp),%esi   // %esi = sourceptr
            movl    20(%esp),%edi   // %edi = destptr
            movl    24(%esp),%edx   // count
            movb    28(%esp),%cl    // i
            orl     %edx,%edx       // count = 0 ?
            jz      L(slcld4)
            subl    %edi,%esi
            // erstes Digit shiften:
            leal    -4(%edi),%edi   // sourceptr--, destptr--
            movl    (%edi,%esi),%ebx // *sourceptr in %ebx halten
            movl    %ebx,%eax       // und in %eax rechnen:
            shll    %cl,%eax        // um i Bits links shiften, rechts Nullen rein
            movl    %eax,(%edi)     // und als *destptr ablegen
            // Letztes Digit in %ebx.
            negb    %cl             // 32-i
            decl    %edx
            jz      L(slcld2)
L(slcld1:)    // weiteres Digit shiften:
              leal    -4(%edi),%edi   // sourceptr--, destptr--
              movl    (%edi,%esi),%eax // nächstes Digit nach %eax
              shrdl   shcl %eax,%ebx  // %ebx um %cl=32-i Bits rechts shiften, %eax von links reinshiften
              movl    %ebx,(%edi)     // %ebx als *destptr ablegen
              // Letztes Digit in %eax.
              decl    %edx
              jz      L(slcld3)
              // weiteres Digit shiften:
              leal    -4(%edi),%edi   // sourceptr--, destptr--
              movl    (%edi,%esi),%ebx // nächstes Digit nach %ebx
              shrdl   shcl %ebx,%eax  // %eax um %cl=32-i Bits rechts shiften, %ebx von links reinshiften
              movl    %eax,(%edi)     // %eax als *destptr ablegen
              // Letztes Digit in %ebx.
              decl    %edx
              jnz     L(slcld1)
L(slcld2:)  movl    %ebx,%eax
L(slcld3:)  shrl    %cl,%eax        // %eax um 32-i Bits nach rechts shiften
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            popl    %esi            // %esi zurück
            ret
L(slcld4:)  xorl    %eax,%eax       // %eax := 0
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            popl    %esi            // %esi zurück
            ret

// extern uintD shift1right_loop_up (uintD* ptr, uintC count, uintD carry);
            ALIGN
            DECLARE_FUNCTION(shift1right_loop_up)
C(shift1right_loop_up:)
            movl    4(%esp),%edx    // %edx = ptr
            movl    8(%esp),%ecx    // %ecx = count
            movl    12(%esp),%eax   // %eax = carry (0 oder -1)
            jecxz   L(s1rld3)       // %ecx = 0 ?
            addl    %eax,%eax       // Carry := Bit 31 von carry
L(s1rld1:)    rcrl    $1,(%edx)       // *ptr und Carry um 1 Bit rechts rotieren
              leal    4(%edx),%edx    // ptr++
              decl    %ecx
              jnz     L(s1rld1)
L(s1rld2:)  sbbl    %eax,%eax       // Ergebnis := - Carry
L(s1rld3:)  ret

// extern uintD shiftright_loop_up (uintD* ptr, uintC count, uintC i);
            ALIGN
            DECLARE_FUNCTION(shiftright_loop_up)
C(shiftright_loop_up:)
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    12(%esp),%edi   // %edi = ptr
            movl    16(%esp),%edx   // %edx = count
            movb    20(%esp),%cl    // %cl = i
            orl     %edx,%edx       // count = 0 ?
            jz      L(srlu4)
            // erstes Digit shiften:
            movl    (%edi),%eax     // Digit in %eax halten
            movl    %eax,%ebx       // und in %ebx rechnen:
            shrl    %cl,%ebx        // um i Bits rechts shiften
            movl    %ebx,(%edi)     // und wieder ablegen
            // Letztes Digit in %eax.
            decl    %edx
            jz      L(srlu2)
            nop ; nop ; nop
L(srlu1:)     // weiteres Digit shiften:
              leal    4(%edi),%edi
              movl    (%edi),%ebx
              shrdl   shcl %eax,(%edi) // (%edi) um %cl=i Bits rechts shiften, %eax von links reinshiften
              // Letztes Digit in %ebx.
              decl    %edx
              jz      L(srlu3)
              // weiteres Digit shiften:
              leal    4(%edi),%edi
              movl    (%edi),%eax
              shrdl   shcl %ebx,(%edi) // (%edi) um %cl=i Bits rechts shiften, %ebx von links reinshiften
              // Letztes Digit in %eax.
              decl    %edx
              jnz     L(srlu1)
L(srlu2:)   movl    %eax,%ebx
L(srlu3:)   xorl    %eax,%eax       // %eax := 0
            shrdl   shcl %ebx,%eax  // %eax := niedrigste %cl=i Bits von %ebx, als Bits 31..32-i
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret
L(srlu4:)   xorl    %eax,%eax       // %eax := 0
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret

// extern uintD shiftrightsigned_loop_up (uintD* ptr, uintC count, uintC i);
            ALIGN
            DECLARE_FUNCTION(shiftrightsigned_loop_up)
C(shiftrightsigned_loop_up:)
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    12(%esp),%edi   // %edi = ptr
            movl    16(%esp),%edx   // %edx = count
            movb    20(%esp),%cl    // %cl = i
            // erstes Digit shiften:
            movl    (%edi),%eax     // Digit in %eax halten
            movl    %eax,%ebx       // und in %ebx rechnen:
            sarl    %cl,%ebx        // um i Bits rechts shiften, Vorzeichen vervielfachen
            movl    %ebx,(%edi)     // und wieder ablegen
            // Letztes Digit in %eax.
            decl    %edx
            jz      L(srslu2)
L(srslu1:)    // weiteres Digit shiften:
              leal    4(%edi),%edi
              movl    (%edi),%ebx
              shrdl   shcl %eax,(%edi) // (%edi) um %cl=i Bits rechts shiften, %eax von links reinshiften
              // Letztes Digit in %ebx.
              decl    %edx
              jz      L(srslu3)
              // weiteres Digit shiften:
              leal    4(%edi),%edi
              movl    (%edi),%eax
              shrdl   shcl %ebx,(%edi) // (%edi) um %cl=i Bits rechts shiften, %ebx von links reinshiften
              // Letztes Digit in %eax.
              decl    %edx
              jnz     L(srslu1)
L(srslu2:)  movl    %eax,%ebx
L(srslu3:)  xorl    %eax,%eax       // %eax := 0
            shrdl   shcl %ebx,%eax  // %eax := niedrigste %cl=i Bits von %ebx, als Bits 31..32-i
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret

// extern uintD shiftrightcopy_loop_up (uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry);
            ALIGN
            DECLARE_FUNCTION(shiftrightcopy_loop_up)
C(shiftrightcopy_loop_up:)
            pushl   %esi            // %esi retten
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    16(%esp),%esi   // %esi = sourceptr
            movl    20(%esp),%edi   // %edi = destptr
            movl    24(%esp),%edx   // count
            movb    28(%esp),%cl    // i
            negb    %cl             // 32-i
            movl    32(%esp),%eax   // %eax = carry
            orl     %edx,%edx       // count = 0 ?
            jz      L(srcld3)
            subl    %edi,%esi
            // erstes Digit shiften:
            movl    (%edi,%esi),%ebx // *sourceptr in %ebx halten
            shldl   shcl %ebx,%eax  // carry um %cl=32-i Bits links shiften, dabei *sourceptr rein
            movl    %eax,(%edi)     // und als *destptr ablegen
            // Letztes Digit in %ebx.
            decl    %edx
            jz      L(srcld2)
L(srcld1:)    // weiteres Digit shiften:
              leal    4(%edi),%edi    // sourceptr++, destptr++
              movl    (%edi,%esi),%eax // nächstes Digit nach %eax
              shldl   shcl %eax,%ebx  // %ebx um %cl=32-i Bits links shiften, %eax von rechts reinshiften
              movl    %ebx,(%edi)     // %ebx als *destptr ablegen
              // Letztes Digit in %eax.
              decl    %edx
              jz      L(srcld3)
              // weiteres Digit shiften:
              leal    4(%edi),%edi    // sourceptr++, destptr++
              movl    (%edi,%esi),%ebx // nächstes Digit nach %ebx
              shldl   shcl %ebx,%eax  // %eax um %cl=32-i Bits links shiften, %ebx von rechts reinshiften
              movl    %eax,(%edi)     // %eax als *destptr ablegen
              // Letztes Digit in %ebx.
              decl    %edx
              jnz     L(srcld1)
L(srcld2:)  movl    %ebx,%eax
L(srcld3:)  shll    %cl,%eax        // %eax um 32-i Bits nach links shiften
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            popl    %esi            // %esi zurück
            ret

// extern uintD mulusmall_loop_down (uintD digit, uintD* ptr, uintC len, uintD newdigit);
            ALIGN
            DECLARE_FUNCTION(mulusmall_loop_down)
C(mulusmall_loop_down:)
            pushl   %ebp            // %ebp retten
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    16(%esp),%ebx   // %ebx = digit
            movl    20(%esp),%edi   // %edi = ptr
            movl    24(%esp),%ecx   // %ecx = len
            movl    28(%esp),%ebp   // %ebp = carry := newdigit
            movl    %ecx,%eax
            negl    %eax            // %eax = -len
            jz      L(msld2)
            leal    -4(%edi,%eax,4),%edi // %edi = &ptr[-1-len]
            nop ; nop ; nop
L(msld1:)     movl    (%edi,%ecx,4),%eax // *ptr
              mull    %ebx               // %edx|%eax := digit * *ptr
              addl    %ebp,%eax          // carry und Low-Teil des Produktes addieren
              movl    $0,%ebp
              adcl    %edx,%ebp          // Übertrag zum High-Teil %edx dazu, gibt neuen carry
              movl    %eax,(%edi,%ecx,4) // Low-Teil als *ptr ablegen
              decl    %ecx               // count--, ptr--
              jnz     L(msld1)
L(msld2:)   movl    %ebp,%eax       // Ergebnis := letzter Übertrag
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            popl    %ebp            // %ebp zurück
            ret

// extern void mulu_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
            ALIGN
            DECLARE_FUNCTION(mulu_loop_down)
C(mulu_loop_down:)
            pushl   %ebp            // %ebp retten
            pushl   %edi            // %edi retten
            pushl   %esi            // %esi retten
            pushl   %ebx            // %ebx retten
            movl    20(%esp),%ebx   // %ebx = digit
            movl    24(%esp),%esi   // %esi = sourceptr
            movl    28(%esp),%edi   // %edi = destptr
            movl    32(%esp),%ecx   // %ecx = len
            movl    %ecx,%eax
            notl    %eax            // %eax = -1-len
            leal    (%esi,%eax,4),%esi // %esi = &sourceptr[-1-len]
            leal    (%edi,%eax,4),%edi // %edi = &destptr[-1-len]
            xorl    %ebp,%ebp       // %epb = carry := 0
L(muld1:)     movl    (%esi,%ecx,4),%eax // *sourceptr
              mull    %ebx               // %edx|%eax := digit * *sourceptr
              addl    %ebp,%eax          // carry und Low-Teil des Produktes addieren
              movl    $0,%ebp
              adcl    %edx,%ebp          // Übertrag zum High-Teil %edx dazu, gibt neuen carry
              movl    %eax,(%edi,%ecx,4) // Low-Teil als *destptr ablegen
              decl    %ecx               // count--, sourceptr--, destptr--
              jnz     L(muld1)
            movl    %ebp,(%edi)     // letzten Übertrag ablegen
            popl    %ebx            // %ebx zurück
            popl    %esi            // %esi zurück
            popl    %edi            // %edi zurück
            popl    %ebp            // %ebp zurück
            ret

// extern uintD muluadd_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
            ALIGN
            DECLARE_FUNCTION(muluadd_loop_down)
C(muluadd_loop_down:)
            pushl   %ebp            // %ebp retten
            pushl   %edi            // %edi retten
            pushl   %esi            // %esi retten
            pushl   %ebx            // %ebx retten
            movl    20(%esp),%ebx   // %ebx = digit
            movl    24(%esp),%esi   // %esi = sourceptr
            movl    28(%esp),%edi   // %edi = destptr
            movl    32(%esp),%ecx   // %ecx = len
            movl    %ecx,%eax
            notl    %eax            // %eax = -1-len
            leal    (%esi,%eax,4),%esi // %esi = &sourceptr[-1-len]
            leal    (%edi,%eax,4),%edi // %edi = &destptr[-1-len]
            xorl    %ebp,%ebp       // %epb = carry := 0
L(muald1:)    movl    (%esi,%ecx,4),%eax // *sourceptr
              mull    %ebx               // %edx|%eax := digit * *sourceptr
              addl    %ebp,%eax          // carry und Low-Teil des Produktes addieren
              movl    $0,%ebp
              adcl    %ebp,%edx          // Übertrag zum High-Teil %edx dazu
              addl    %eax,(%edi,%ecx,4) // Low-Teil zu *destptr addieren
              adcl    %edx,%ebp          // zweiten Übertrag zu %edx addieren, gibt neuen carry
              decl    %ecx               // count--, sourceptr--, destptr--
              jnz     L(muald1)
            movl    %ebp,%eax       // Ergebnis := letzter Übertrag
            popl    %ebx            // %ebx zurück
            popl    %esi            // %esi zurück
            popl    %edi            // %edi zurück
            popl    %ebp            // %ebp zurück
            ret

// extern uintD mulusub_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
            ALIGN
            DECLARE_FUNCTION(mulusub_loop_down)
C(mulusub_loop_down:)
            pushl   %ebp            // %ebp retten
            pushl   %edi            // %edi retten
            pushl   %esi            // %esi retten
            pushl   %ebx            // %ebx retten
            movl    20(%esp),%ebx   // %ebx = digit
            movl    24(%esp),%esi   // %esi = sourceptr
            movl    28(%esp),%edi   // %edi = destptr
            movl    32(%esp),%ecx   // %ecx = len
            movl    %ecx,%eax
            notl    %eax            // %eax = -1-len
            leal    (%esi,%eax,4),%esi // %esi = &sourceptr[-1-len]
            leal    (%edi,%eax,4),%edi // %edi = &destptr[-1-len]
            xorl    %ebp,%ebp       // %epb = carry := 0
L(musld1:)    movl    (%esi,%ecx,4),%eax // *sourceptr
              mull    %ebx               // %edx|%eax := digit * *sourceptr
              addl    %ebp,%eax          // carry und Low-Teil des Produktes addieren
              movl    $0,%ebp
              adcl    %ebp,%edx          // Übertrag zum High-Teil %edx dazu
              subl    %eax,(%edi,%ecx,4) // Low-Teil von *destptr subtrahieren
              adcl    %edx,%ebp          // zweiten Übertrag zu %edx addieren, gibt neuen carry
              decl    %ecx               // count--, sourceptr--, destptr--
              jnz     L(musld1)
            movl    %ebp,%eax       // Ergebnis := letzter Übertrag
            popl    %ebx            // %ebx zurück
            popl    %esi            // %esi zurück
            popl    %edi            // %edi zurück
            popl    %ebp            // %ebp zurück
            ret

// extern uintD divu_loop_up (uintD digit, uintD* ptr, uintC len);
            ALIGN
            DECLARE_FUNCTION(divu_loop_up)
C(divu_loop_up:)
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    12(%esp),%ebx   // %ebx = digit
            movl    16(%esp),%edi   // %edi = ptr
            movl    20(%esp),%ecx   // %ecx = len
            xorl    %edx,%edx       // %edx = Rest := 0
            jecxz   L(dlu2)         // %ecx = 0 ?
L(dlu1:)      movl    (%edi),%eax     // nächstes Digit *ptr
              divl    %ebx            // Division von %edx|%eax durch %ebx
              movl    %eax,(%edi)     // Quotient %eax ablegen, Rest in %edx behalten
              leal    4(%edi),%edi    // ptr++
              decl    %ecx
              jnz     L(dlu1)
L(dlu2:)    movl    %edx,%eax       // Ergebnis := letzter Rest
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret

// extern uintD divucopy_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
            ALIGN
            DECLARE_FUNCTION(divucopy_loop_up)
C(divucopy_loop_up:)
            pushl   %edi            // %edi retten
            pushl   %esi            // %esi retten
            pushl   %ebx            // %ebx retten
            movl    16(%esp),%ebx   // %ebx = digit
            movl    20(%esp),%esi   // %esi = sourceptr
            movl    24(%esp),%edi   // %edi = destptr
            movl    28(%esp),%ecx   // %ecx = len
            xorl    %edx,%edx       // %edx = Rest := 0
            jecxz   L(dclu2)        // %ecx = 0 ?
            subl    %edi,%esi
L(dclu1:)     movl    (%esi,%edi),%eax // nächstes Digit *ptr
              divl    %ebx            // Division von %edx|%eax durch %ebx
              movl    %eax,(%edi)     // Quotient %eax ablegen, Rest in %edx behalten
              leal    4(%edi),%edi    // sourceptr++, destptr++
              decl    %ecx
              jnz     L(dclu1)
L(dclu2:)   movl    %edx,%eax       // Ergebnis := letzter Rest
            popl    %ebx            // %ebx zurück
            popl    %esi            // %esi zurück
            popl    %edi            // %edi zurück
            ret

#endif

#if !CL_DS_BIG_ENDIAN_P

// extern void or_loop_down (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(or_loop_down)
C(or_loop_down:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(old2)         // %ecx = 0 ?
L(old1:)      leal    -4(%edx),%edx    // xptr--, yptr--
              movl    (%edx,%esi),%eax // *yptr
              orl     %eax,(%edx)      // *xptr |= ...
              decl    %ecx
              jnz     L(old1)
L(old2:)    popl    %esi            // %esi zurück
            ret

// extern void xor_loop_down (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(xor_loop_down)
C(xor_loop_down:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(xld2)         // %ecx = 0 ?
L(xld1:)      leal    -4(%edx),%edx    // xptr--, yptr--
              movl    (%edx,%esi),%eax // *yptr
              xorl    %eax,(%edx)      // *xptr ^= ...
              decl    %ecx
              jnz     L(xld1)
L(xld2:)    popl    %esi            // %esi zurück
            ret

// extern void and_loop_down (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(and_loop_down)
C(and_loop_down:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(ald2)         // %ecx = 0 ?
L(ald1:)      leal    -4(%edx),%edx    // xptr--, yptr--
              movl    (%edx,%esi),%eax // *yptr
              andl    %eax,(%edx)      // *xptr &= ...
              decl    %ecx
              jnz     L(ald1)
L(ald2:)    popl    %esi            // %esi zurück
            ret

// extern void eqv_loop_down (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(eqv_loop_down)
C(eqv_loop_down:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(eld2)         // %ecx = 0 ?
L(eld1:)      leal    -4(%edx),%edx    // xptr--, yptr--
              movl    (%edx),%eax      // *xptr
              xorl    (%edx,%esi),%eax // ^ *yptr
              notl    %eax             // ~(...)
              movl    %eax,(%edx)      // =: *xptr
              decl    %ecx
              jnz     L(eld1)
L(eld2:)    popl    %esi            // %esi zurück
            ret

// extern void nand_loop_down (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(nand_loop_down)
C(nand_loop_down:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(nald2)        // %ecx = 0 ?
L(nald1:)     leal    -4(%edx),%edx    // xptr--, yptr--
              movl    (%edx),%eax      // *xptr
              andl    (%edx,%esi),%eax // & *yptr
              notl    %eax             // ~(...)
              movl    %eax,(%edx)      // =: *xptr
              decl    %ecx
              jnz     L(nald1)
L(nald2:)   popl    %esi            // %esi zurück
            ret

// extern void nor_loop_down (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(nor_loop_down)
C(nor_loop_down:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(nold2)        // %ecx = 0 ?
L(nold1:)     leal    -4(%edx),%edx    // xptr--, yptr--
              movl    (%edx),%eax      // *xptr
              orl     (%edx,%esi),%eax // | *yptr
              notl    %eax             // ~(...)
              movl    %eax,(%edx)      // =: *xptr
              decl    %ecx
              jnz     L(nold1)
L(nold2:)   popl    %esi            // %esi zurück
            ret

// extern void andc2_loop_down (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(andc2_loop_down)
C(andc2_loop_down:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(acld2)        // %ecx = 0 ?
L(acld1:)     leal    -4(%edx),%edx    // xptr--, yptr--
              movl    (%edx,%esi),%eax // *yptr
              notl    %eax             // ~ *yptr
              andl    %eax,(%edx)      // *xptr &= ...
              decl    %ecx
              jnz     L(acld1)
L(acld2:)   popl    %esi            // %esi zurück
            ret

// extern void orc2_loop_down (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(orc2_loop_down)
C(orc2_loop_down:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edx,%esi
            jecxz   L(ocld2)        // %ecx = 0 ?
L(ocld1:)     leal    -4(%edx),%edx    // xptr--, yptr--
              movl    (%edx,%esi),%eax // *yptr
              notl    %eax             // ~ *yptr
              orl     %eax,(%edx)      // *xptr |= ...
              decl    %ecx
              jnz     L(ocld1)
L(ocld2:)   popl    %esi            // %esi zurück
            ret

// extern void not_loop_down (uintD* xptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(not_loop_down)
C(not_loop_down:)
            movl    4(%esp),%edx    // %edx = xptr
            movl    8(%esp),%ecx    // %ecx = count
            jecxz   L(nld2)         // %ecx = 0 ?
            nop ; nop ; nop ; nop ; nop ; nop
L(nld1:)      leal    -4(%edx),%edx    // xptr--
              notl    (%edx)           // ~= *xptr
              decl    %ecx
              jnz     L(nld1)
L(nld2:)    ret

// extern boolean and_test_loop_down (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(and_test_loop_down)
C(and_test_loop_down:)
            pushl   %esi            // %esi retten
            movl    8(%esp),%edx    // %edx = xptr
            movl    12(%esp),%esi   // %esi = yptr
            movl    16(%esp),%ecx   // %ecx = count
            jecxz   L(atld2)        // %ecx = 0 ?
            subl    %edx,%esi
L(atld1:)     leal    -4(%edx),%edx    // xptr--, yptr--
              movl    (%edx,%esi),%eax // *yptr
              andl    (%edx),%eax      // *xptr & ...
              jnz     L(atld3)
              decl    %ecx
              jnz     L(atld1)
L(atld2:)   xorl    %eax,%eax       // Ergebnis 0
            popl    %esi            // %esi zurück
            ret
L(atld3:)   movl    $1,%eax         // Ergebnis 1 (nicht irgendwas /=0 !)
            popl    %esi            // %esi zurück
            ret

// extern cl_signean compare_loop_down (uintD* xptr, uintD* yptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(compare_loop_down)
C(compare_loop_down:)
            movl    %esi,%edx       // %esi retten
            movl    %edi,%eax       // %edi retten
            movl    4(%esp),%esi    // %esi = xptr
            movl    8(%esp),%edi    // %edi = yptr
            movl    12(%esp),%ecx   // %ecx = count
            leal    -4(%esi),%esi
            leal    -4(%edi),%edi
            cmpl    %ecx,%ecx       // initialize flags for the case %ecx is 0
            dir1start
            repz                    // Falls %ecx > 0:
              cmpsl                 // %ecx mal aufwärts (%edi) und (%esi) vergleichen
                                    // und weiterschleifen, falls Z, d.h. (%edi)=(%esi).
            dir1end
            // Flags -> Ergebnis:
            // Z,NC -> bis zum Schluß (%esi)-(%edi) = 0 -> x=y -> Ergebnis 0
            // NZ,C -> schließlich (%esi)-(%edi) < 0 -> x<y -> Ergebnis -1
            // NZ,NC -> schließlich (%esi)-(%edi) > 0 -> x>y -> Ergebnis +1
            movl    %eax,%edi       // %edi zurück
            movl    %edx,%esi       // %esi zurück
            jbe     L(cmld1)        // "be" = Z oder C
            movl    $1,%eax         // Ergebnis +1
            ret
L(cmld1:)   sbbl    %eax,%eax       // Ergebnis -1 (falls C) oder 0 (falls NC)
            ret

// extern uintD add_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(add_loop_up)
C(add_loop_up:)
            pushl   %esi            // %esi retten
            pushl   %edi            // %edi retten
            movl    12(%esp),%edx   // %edx = sourceptr1
            movl    16(%esp),%esi   // %esi = sourceptr2
            movl    20(%esp),%edi   // %edi = destptr
            movl    24(%esp),%ecx   // %ecx = count
            subl    %edi,%edx
            subl    %edi,%esi
            orl     %ecx,%ecx       // %ecx = 0 ?, Carry löschen
            jz      L(alu2)
L(alu1:)      movl    (%edx,%edi),%eax // *sourceptr1
              adcl    (%esi,%edi),%eax // + *sourceptr2 + carry
              movl    %eax,(%edi)     // =: *destptr, neuen Carry behalten
              leal    4(%edi),%edi    // sourceptr1++, sourceptr2++, destptr++
              decl    %ecx
              jnz     L(alu1)
L(alu2:)    sbbl    %eax,%eax      // Ergebnis := - Carry
            popl    %edi           // %edi zurück
            popl    %esi           // %esi zurück
            ret

// extern uintD addto_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(addto_loop_up)
C(addto_loop_up:)
            pushl   %edi            // %edi retten
            movl    8(%esp),%edx    // %edx = sourceptr
            movl    12(%esp),%edi   // %edi = destptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edi,%edx
            orl     %ecx,%ecx       // %ecx = 0 ?, Carry löschen
            jz      L(atlu2)
L(atlu1:)     movl    (%edx,%edi),%eax // *sourceptr
              adcl    %eax,(%edi)     // + *destptr + carry =: *destptr, neuer Carry
              leal    4(%edi),%edi    // sourceptr++, destptr++
              decl    %ecx
              jnz     L(atlu1)
L(atlu2:)   sbbl    %eax,%eax       // Ergebnis := - Carry
            popl    %edi            // %edi zurück
            ret

// extern uintD inc_loop_up (uintD* ptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(inc_loop_up)
C(inc_loop_up:)
            movl    4(%esp),%edx    // %edx = ptr
            movl    8(%esp),%ecx    // %ecx = count
            jecxz   L(ilu2)         // %ecx = 0 ?
L(ilu1:)      addl    $1,(%edx)       // (*ptr)++
              jnc     L(ilu3)         // kein Carry -> fertig
              leal    4(%edx),%edx
              decl    %ecx
              jnz     L(ilu1)
L(ilu2:)    movl    $1,%eax         // Ergebnis := 1
            ret
L(ilu3:)    xorl    %eax,%eax       // Ergebnis := 0
            ret

// extern uintD sub_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(sub_loop_up)
C(sub_loop_up:)
            pushl   %esi            // %esi retten
            pushl   %edi            // %edi retten
            movl    12(%esp),%edx   // %edx = sourceptr1
            movl    16(%esp),%esi   // %esi = sourceptr2
            movl    20(%esp),%edi   // %edi = destptr
            movl    24(%esp),%ecx   // %ecx = count
            subl    %edi,%edx
            subl    %edi,%esi
            orl     %ecx,%ecx       // %ecx = 0 ?, Carry löschen
            jz      L(slu2)
L(slu1:)      movl    (%edx,%edi),%eax // *sourceptr1
              sbbl    (%esi,%edi),%eax // - *sourceptr2 - carry
              movl    %eax,(%edi)     // =: *destptr, neuen Carry behalten
              leal    4(%edi),%edi    // sourceptr1++, sourceptr2++, destptr++
              decl    %ecx
              jnz     L(slu1)
L(slu2:)    sbbl    %eax,%eax      // Ergebnis := - Carry
            popl    %edi           // %edi zurück
            popl    %esi           // %esi zurück
            ret

// extern uintD subx_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);
            ALIGN
            DECLARE_FUNCTION(subx_loop_up)
C(subx_loop_up:)
            pushl   %esi            // %esi retten
            pushl   %edi            // %edi retten
            movl    12(%esp),%edx   // %edx = sourceptr1
            movl    16(%esp),%esi   // %esi = sourceptr2
            movl    20(%esp),%edi   // %edi = destptr
            movl    24(%esp),%ecx   // %ecx = count
            jecxz   L(sxlu2)        // %ecx = 0 ?
            subl    %edi,%edx
            subl    %edi,%esi
            movl    28(%esp),%eax   // carry, 0 oder -1
            addl    %eax,%eax       // Bit 31 davon in den Carry
            nop ; nop
L(sxlu1:)     movl    (%edx,%edi),%eax // *sourceptr1
              sbbl    (%esi,%edi),%eax // - *sourceptr2 - carry
              movl    %eax,(%edi)     // =: *destptr, neuen Carry behalten
              leal    4(%edi),%edi    // sourceptr1++, sourceptr2++, destptr++
              decl    %ecx
              jnz     L(sxlu1)
            sbbl    %eax,%eax      // Ergebnis := - Carry
            popl    %edi           // %edi zurück
            popl    %esi           // %esi zurück
            ret
L(sxlu2:)   movl    28(%esp),%eax  // Ergebnis := carry
            popl    %edi           // %edi zurück
            popl    %esi           // %esi zurück
            ret

// extern uintD subfrom_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(subfrom_loop_up)
C(subfrom_loop_up:)
            pushl   %edi            // %edi retten
            movl    8(%esp),%edx    // %edx = sourceptr
            movl    12(%esp),%edi   // %edi = destptr
            movl    16(%esp),%ecx   // %ecx = count
            subl    %edi,%edx
            orl     %ecx,%ecx       // %ecx = 0 ?, Carry löschen
            jz      L(sflu2)
L(sflu1:)     movl    (%edx,%edi),%eax // *sourceptr
              sbbl    %eax,(%edi)     // *destptr - *sourceptr - carry =: *destptr, neuer Carry
              leal    4(%edi),%edi    // sourceptr++, destptr++
              decl    %ecx
              jnz     L(sflu1)
L(sflu2:)   sbbl    %eax,%eax       // Ergebnis := - Carry
            popl    %edi            // %edi zurück
            ret

// extern uintD dec_loop_up (uintD* ptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(dec_loop_up)
C(dec_loop_up:)
            movl    4(%esp),%edx    // %edx = ptr
            movl    8(%esp),%ecx    // %ecx = count
            jecxz   L(dlu2)         // %ecx = 0 ?
L(dlu1:)      subl    $1,(%edx)       // (*ptr)--
              jnc     L(dlu3)         // kein Carry -> fertig
              leal    4(%edx),%edx
              decl    %ecx
              jnz     L(dlu1)
L(dlu2:)    movl    $-1,%eax        // Ergebnis := -1
            ret
L(dlu3:)    xorl    %eax,%eax       // Ergebnis := 0
            ret

// extern uintD neg_loop_up (uintD* ptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(neg_loop_up)
C(neg_loop_up:)
            movl    4(%esp),%edx    // %edx = ptr
            movl    8(%esp),%ecx    // %ecx = count
            // erstes Digit /=0 suchen:
            jecxz   L(nlu2)         // %ecx = 0 ?
L(nlu1:)      negl    (%edx)
              jnz     L(nlu3)
              leal    4(%edx),%edx
              decl    %ecx
              jnz     L(nlu1)
L(nlu2:)    xorl    %eax,%eax       // Ergebnis := 0
            ret
            nop ; nop ; nop ; nop ; nop ; nop
L(nlu3:)    // erstes Digit /=0 gefunden, ab jetzt gibt's Carrys
            // alle anderen Digits invertieren:
            decl    %ecx
            jz      L(nlu5)
L(nlu4:)      leal    4(%edx),%edx
              notl    (%edx)
              decl    %ecx
              jnz     L(nlu4)
L(nlu5:)    movl    $-1,%eax        // Ergebnis := -1
            ret

// extern uintD shift1left_loop_up (uintD* ptr, uintC count);
            ALIGN
            DECLARE_FUNCTION(shift1left_loop_up)
C(shift1left_loop_up:)
            movl    4(%esp),%edx    // %edx = ptr
            movl    8(%esp),%ecx    // %ecx = count
            orl     %ecx,%ecx       // %ecx = 0 ?, Carry löschen
            jz      L(s1llu2)
            nop ; nop ; nop ; nop
L(s1llu1:)    rcll    $1,(%edx)       // *ptr und Carry um 1 Bit links rotieren
              leal    4(%edx),%edx    // ptr++
              decl    %ecx
              jnz     L(s1llu1)
L(s1llu2:)  sbbl    %eax,%eax       // Ergebnis := - Carry
            ret

// extern uintD shiftleft_loop_up (uintD* ptr, uintC count, uintC i, uintD carry);
            ALIGN
            DECLARE_FUNCTION(shiftleft_loop_up)
C(shiftleft_loop_up:)
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    12(%esp),%edi   // %edi = ptr
            movl    16(%esp),%edx   // %edx = count
            movb    20(%esp),%cl    // %cl = i
            orl     %edx,%edx       // count = 0 ?
            jz      L(sllu4)
            // erstes Digit shiften:
            movl    (%edi),%eax     // Digit in %eax halten
            movl    %eax,%ebx       // und in %ebx rechnen:
            shll    %cl,%ebx        // um i Bits links shiften
            orl     24(%esp),%ebx   // und die unteren i Bits eintragen
            movl    %ebx,(%edi)     // und wieder ablegen
            leal    4(%edi),%edi
            // Letztes Digit in %eax.
            decl    %edx
            jz      L(sllu2)
            nop ; nop ; nop ; nop
L(sllu1:)     // weiteres Digit shiften:
              movl    (%edi),%ebx
              shldl   shcl %eax,(%edi) // (%edi) um %cl=i Bits links shiften, %eax von rechts reinshiften
              leal    4(%edi),%edi
              // Letztes Digit in %ebx.
              decl    %edx
              jz      L(sllu3)
              // weiteres Digit shiften:
              movl    (%edi),%eax
              shldl   shcl %ebx,(%edi) // (%edi) um %cl=i Bits links shiften, %ebx von rechts reinshiften
              leal    4(%edi),%edi
              // Letztes Digit in %eax.
              decl    %edx
              jnz     L(sllu1)
L(sllu2:)   movl    %eax,%ebx
L(sllu3:)   xorl    %eax,%eax       // %eax := 0
            shldl   shcl %ebx,%eax  // %eax := höchste %cl=i Bits von %ebx
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret
L(sllu4:)   movl    24(%esp),%eax   // %eax := carry
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret

#endif

// extern uintD shiftleftcopy_loop_up (uintD* sourceptr, uintD* destptr, uintC count, uintC i);
            ALIGN
            DECLARE_FUNCTION(shiftleftcopy_loop_up)
C(shiftleftcopy_loop_up:)
            pushl   %esi            // %esi retten
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    16(%esp),%esi   // %esi = sourceptr
            movl    20(%esp),%edi   // %edi = destptr
            movl    24(%esp),%edx   // count
            movb    28(%esp),%cl    // i
            orl     %edx,%edx       // count = 0 ?
            jz      L(slclu4)
            subl    %edi,%esi
            // erstes Digit shiften:
            movl    (%edi,%esi),%ebx // *sourceptr in %ebx halten
            movl    %ebx,%eax       // und in %eax rechnen:
            shll    %cl,%eax        // um i Bits links shiften, rechts Nullen rein
            movl    %eax,(%edi)     // und als *destptr ablegen
            leal    4(%edi),%edi    // sourceptr++, destptr++
            // Letztes Digit in %ebx.
            negb    %cl             // 32-i
            decl    %edx
            jz      L(slclu2)
L(slclu1:)    // weiteres Digit shiften:
              movl    (%edi,%esi),%eax // nächstes Digit nach %eax
              shrdl   shcl %eax,%ebx  // %ebx um %cl=32-i Bits rechts shiften, %eax von links reinshiften
              movl    %ebx,(%edi)     // %ebx als *destptr ablegen
              leal    4(%edi),%edi    // sourceptr++, destptr++
              // Letztes Digit in %eax.
              decl    %edx
              jz      L(slclu3)
              // weiteres Digit shiften:
              movl    (%edi,%esi),%ebx // nächstes Digit nach %ebx
              shrdl   shcl %ebx,%eax  // %eax um %cl=32-i Bits rechts shiften, %ebx von links reinshiften
              movl    %eax,(%edi)     // %eax als *destptr ablegen
              leal    4(%edi),%edi    // sourceptr++, destptr++
              // Letztes Digit in %ebx.
              decl    %edx
              jnz     L(slclu1)
L(slclu2:)  movl    %ebx,%eax
L(slclu3:)  shrl    %cl,%eax        // %eax um 32-i Bits nach rechts shiften
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            popl    %esi            // %esi zurück
            ret
L(slclu4:)  xorl    %eax,%eax       // %eax := 0
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            popl    %esi            // %esi zurück
            ret

#if !CL_DS_BIG_ENDIAN_P

// extern uintD shift1right_loop_down (uintD* ptr, uintC count, uintD carry);
            ALIGN
            DECLARE_FUNCTION(shift1right_loop_down)
C(shift1right_loop_down:)
            movl    4(%esp),%edx    // %edx = ptr
            movl    8(%esp),%ecx    // %ecx = count
            movl    12(%esp),%eax   // %eax = carry (0 oder -1)
            jecxz   L(s1rlu3)       // %ecx = 0 ?
            addl    %eax,%eax       // Carry := Bit 31 von carry
L(s1rlu1:)    leal    -4(%edx),%edx   // ptr--
              rcrl    $1,(%edx)       // *ptr und Carry um 1 Bit rechts rotieren
              decl    %ecx
              jnz     L(s1rlu1)
L(s1rlu2:)  sbbl    %eax,%eax       // Ergebnis := - Carry
L(s1rlu3:)  ret

// extern uintD shiftright_loop_down (uintD* ptr, uintC count, uintC i);
            ALIGN
            DECLARE_FUNCTION(shiftright_loop_down)
C(shiftright_loop_down:)
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    12(%esp),%edi   // %edi = ptr
            movl    16(%esp),%edx   // %edx = count
            movb    20(%esp),%cl    // %cl = i
            orl     %edx,%edx       // count = 0 ?
            jz      L(srld4)
            // erstes Digit shiften:
            leal    -4(%edi),%edi
            movl    (%edi),%eax     // Digit in %eax halten
            movl    %eax,%ebx       // und in %ebx rechnen:
            shrl    %cl,%ebx        // um i Bits rechts shiften
            movl    %ebx,(%edi)     // und wieder ablegen
            // Letztes Digit in %eax.
            decl    %edx
            jz      L(srld2)
L(srld1:)     // weiteres Digit shiften:
              leal    -4(%edi),%edi
              movl    (%edi),%ebx
              shrdl   shcl %eax,(%edi) // (%edi) um %cl=i Bits rechts shiften, %eax von links reinshiften
              // Letztes Digit in %ebx.
              decl    %edx
              jz      L(srld3)
              // weiteres Digit shiften:
              leal    -4(%edi),%edi
              movl    (%edi),%eax
              shrdl   shcl %ebx,(%edi) // (%edi) um %cl=i Bits rechts shiften, %ebx von links reinshiften
              // Letztes Digit in %eax.
              decl    %edx
              jnz     L(srld1)
L(srld2:)   movl    %eax,%ebx
L(srld3:)   xorl    %eax,%eax       // %eax := 0
            shrdl   shcl %ebx,%eax  // %eax := niedrigste %cl=i Bits von %ebx, als Bits 31..32-i
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret
L(srld4:)   xorl    %eax,%eax       // %eax := 0
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret

// extern uintD shiftrightsigned_loop_down (uintD* ptr, uintC count, uintC i);
            ALIGN
            DECLARE_FUNCTION(shiftrightsigned_loop_down)
C(shiftrightsigned_loop_down:)
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    12(%esp),%edi   // %edi = ptr
            movl    16(%esp),%edx   // %edx = count
            movb    20(%esp),%cl    // %cl = i
            // erstes Digit shiften:
            leal    -4(%edi),%edi
            movl    (%edi),%eax     // Digit in %eax halten
            movl    %eax,%ebx       // und in %ebx rechnen:
            sarl    %cl,%ebx        // um i Bits rechts shiften, Vorzeichen vervielfachen
            movl    %ebx,(%edi)     // und wieder ablegen
            // Letztes Digit in %eax.
            decl    %edx
            jz      L(srsld2)
            nop ; nop ; nop ; nop
L(srsld1:)    // weiteres Digit shiften:
              leal    -4(%edi),%edi
              movl    (%edi),%ebx
              shrdl   shcl %eax,(%edi) // (%edi) um %cl=i Bits rechts shiften, %eax von links reinshiften
              // Letztes Digit in %ebx.
              decl    %edx
              jz      L(srsld3)
              // weiteres Digit shiften:
              leal    -4(%edi),%edi
              movl    (%edi),%eax
              shrdl   shcl %ebx,(%edi) // (%edi) um %cl=i Bits rechts shiften, %ebx von links reinshiften
              // Letztes Digit in %eax.
              decl    %edx
              jnz     L(srsld1)
L(srsld2:)  movl    %eax,%ebx
L(srsld3:)  xorl    %eax,%eax       // %eax := 0
            shrdl   shcl %ebx,%eax  // %eax := niedrigste %cl=i Bits von %ebx, als Bits 31..32-i
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret

// extern uintD shiftrightcopy_loop_down (uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry);
            ALIGN
            DECLARE_FUNCTION(shiftrightcopy_loop_down)
C(shiftrightcopy_loop_down:)
            pushl   %esi            // %esi retten
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    16(%esp),%esi   // %esi = sourceptr
            movl    20(%esp),%edi   // %edi = destptr
            movl    24(%esp),%edx   // count
            movb    28(%esp),%cl    // i
            negb    %cl             // 32-i
            movl    32(%esp),%eax   // %eax = carry
            orl     %edx,%edx       // count = 0 ?
            jz      L(srclu3)
            subl    %edi,%esi
            // erstes Digit shiften:
            leal    -4(%edi),%edi   // sourceptr--, destptr--
            movl    (%edi,%esi),%ebx // *sourceptr in %ebx halten
            shldl   shcl %ebx,%eax  // carry um %cl=32-i Bits links shiften, dabei *sourceptr rein
            movl    %eax,(%edi)     // und als *destptr ablegen
            // Letztes Digit in %ebx.
            decl    %edx
            jz      L(srclu2)
            nop ; nop ; nop
L(srclu1:)    // weiteres Digit shiften:
              leal    -4(%edi),%edi   // sourceptr--, destptr--
              movl    (%edi,%esi),%eax // nächstes Digit nach %eax
              shldl   shcl %eax,%ebx  // %ebx um %cl=32-i Bits links shiften, %eax von rechts reinshiften
              movl    %ebx,(%edi)     // %ebx als *destptr ablegen
              // Letztes Digit in %eax.
              decl    %edx
              jz      L(srclu3)
              // weiteres Digit shiften:
              leal    -4(%edi),%edi   // sourceptr--, destptr--
              movl    (%edi,%esi),%ebx // nächstes Digit nach %ebx
              shldl   shcl %ebx,%eax  // %eax um %cl=32-i Bits links shiften, %ebx von rechts reinshiften
              movl    %eax,(%edi)     // %eax als *destptr ablegen
              // Letztes Digit in %ebx.
              decl    %edx
              jnz     L(srclu1)
L(srclu2:)  movl    %ebx,%eax
L(srclu3:)  shll    %cl,%eax        // %eax um 32-i Bits nach links shiften
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            popl    %esi            // %esi zurück
            ret

// extern uintD mulusmall_loop_up (uintD digit, uintD* ptr, uintC len, uintD newdigit);
            ALIGN
            DECLARE_FUNCTION(mulusmall_loop_up)
C(mulusmall_loop_up:)
            pushl   %ebp            // %ebp retten
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    16(%esp),%ebx   // %ebx = digit
            movl    20(%esp),%edi   // %edi = ptr
            movl    24(%esp),%ecx   // %ecx = len
            movl    28(%esp),%ebp   // %ebp = carry := newdigit
            leal    (%edi,%ecx,4),%edi // %edi = &ptr[len]
            negl    %ecx            // %ecx = -count
            jz      L(mslu2)
L(mslu1:)     movl    (%edi,%ecx,4),%eax // *ptr
              mull    %ebx               // %edx|%eax := digit * *ptr
              addl    %ebp,%eax          // carry und Low-Teil des Produktes addieren
              movl    $0,%ebp
              adcl    %edx,%ebp          // Übertrag zum High-Teil %edx dazu, gibt neuen carry
              movl    %eax,(%edi,%ecx,4) // Low-Teil als *ptr ablegen
              incl    %ecx               // count--, ptr++
              jnz     L(mslu1)
L(mslu2:)   movl    %ebp,%eax       // Ergebnis := letzter Übertrag
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            popl    %ebp            // %ebp zurück
            ret

// extern void mulu_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
            ALIGN
            DECLARE_FUNCTION(mulu_loop_up)
C(mulu_loop_up:)
            pushl   %ebp            // %ebp retten
            pushl   %edi            // %edi retten
            pushl   %esi            // %esi retten
            pushl   %ebx            // %ebx retten
            movl    20(%esp),%ebx   // %ebx = digit
            movl    24(%esp),%esi   // %esi = sourceptr
            movl    28(%esp),%edi   // %edi = destptr
            movl    32(%esp),%ecx   // %ecx = len
            leal    (%esi,%ecx,4),%esi // %esi = &sourceptr[len]
            leal    (%edi,%ecx,4),%edi // %edi = &destptr[len]
            negl    %ecx            // %ecx = -count
            xorl    %ebp,%ebp       // %epb = carry := 0
            nop ; nop
L(mulu1:)     movl    (%esi,%ecx,4),%eax // *sourceptr
              mull    %ebx               // %edx|%eax := digit * *sourceptr
              addl    %ebp,%eax          // carry und Low-Teil des Produktes addieren
              movl    $0,%ebp
              adcl    %edx,%ebp          // Übertrag zum High-Teil %edx dazu, gibt neuen carry
              movl    %eax,(%edi,%ecx,4) // Low-Teil als *destptr ablegen
              incl    %ecx               // count--, sourceptr++, destptr++
              jnz     L(mulu1)
            movl    %ebp,(%edi)     // letzten Übertrag ablegen
            popl    %ebx            // %ebx zurück
            popl    %esi            // %esi zurück
            popl    %edi            // %edi zurück
            popl    %ebp            // %ebp zurück
            ret

// extern uintD muluadd_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
            ALIGN
            DECLARE_FUNCTION(muluadd_loop_up)
C(muluadd_loop_up:)
            pushl   %ebp            // %ebp retten
            pushl   %edi            // %edi retten
            pushl   %esi            // %esi retten
            pushl   %ebx            // %ebx retten
            movl    20(%esp),%ebx   // %ebx = digit
            movl    24(%esp),%esi   // %esi = sourceptr
            movl    28(%esp),%edi   // %edi = destptr
            movl    32(%esp),%ecx   // %ecx = len
            leal    (%esi,%ecx,4),%esi // %esi = &sourceptr[len]
            leal    (%edi,%ecx,4),%edi // %edi = &destptr[len]
            negl    %ecx            // %ecx = -count
            xorl    %ebp,%ebp       // %epb = carry := 0
            nop ; nop
L(mualu1:)    movl    (%esi,%ecx,4),%eax // *sourceptr
              mull    %ebx               // %edx|%eax := digit * *sourceptr
              addl    %ebp,%eax          // carry und Low-Teil des Produktes addieren
              movl    $0,%ebp
              adcl    %ebp,%edx          // Übertrag zum High-Teil %edx dazu
              addl    %eax,(%edi,%ecx,4) // Low-Teil zu *destptr addieren
              adcl    %edx,%ebp          // zweiten Übertrag zu %edx addieren, gibt neuen carry
              incl    %ecx               // count--, sourceptr++, destptr++
              jnz     L(mualu1)
            movl    %ebp,%eax       // Ergebnis := letzter Übertrag
            popl    %ebx            // %ebx zurück
            popl    %esi            // %esi zurück
            popl    %edi            // %edi zurück
            popl    %ebp            // %ebp zurück
            ret

// extern uintD mulusub_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
            ALIGN
            DECLARE_FUNCTION(mulusub_loop_up)
C(mulusub_loop_up:)
            pushl   %ebp            // %ebp retten
            pushl   %edi            // %edi retten
            pushl   %esi            // %esi retten
            pushl   %ebx            // %ebx retten
            movl    20(%esp),%ebx   // %ebx = digit
            movl    24(%esp),%esi   // %esi = sourceptr
            movl    28(%esp),%edi   // %edi = destptr
            movl    32(%esp),%ecx   // %ecx = len
            leal    (%esi,%ecx,4),%esi // %esi = &sourceptr[len]
            leal    (%edi,%ecx,4),%edi // %edi = &destptr[len]
            negl    %ecx            // %ecx = -count
            xorl    %ebp,%ebp       // %epb = carry := 0
            nop ; nop
L(muslu1:)    movl    (%esi,%ecx,4),%eax // *sourceptr
              mull    %ebx               // %edx|%eax := digit * *sourceptr
              addl    %ebp,%eax          // carry und Low-Teil des Produktes addieren
              movl    $0,%ebp
              adcl    %ebp,%edx          // Übertrag zum High-Teil %edx dazu
              subl    %eax,(%edi,%ecx,4) // Low-Teil von *destptr subtrahieren
              adcl    %edx,%ebp          // zweiten Übertrag zu %edx addieren, gibt neuen carry
              incl    %ecx               // count--, sourceptr++, destptr++
              jnz     L(muslu1)
            movl    %ebp,%eax       // Ergebnis := letzter Übertrag
            popl    %ebx            // %ebx zurück
            popl    %esi            // %esi zurück
            popl    %edi            // %edi zurück
            popl    %ebp            // %ebp zurück
            ret

// extern uintD divu_loop_down (uintD digit, uintD* ptr, uintC len);
            ALIGN
            DECLARE_FUNCTION(divu_loop_down)
C(divu_loop_down:)
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    12(%esp),%ebx   // %ebx = digit
            movl    16(%esp),%edi   // %edi = ptr
            movl    20(%esp),%ecx   // %ecx = len
            xorl    %edx,%edx       // %edx = Rest := 0
            jecxz   L(dld2)         // %ecx = 0 ?
L(dld1:)      leal    -4(%edi),%edi   // ptr--
              movl    (%edi),%eax     // nächstes Digit *ptr
              divl    %ebx            // Division von %edx|%eax durch %ebx
              movl    %eax,(%edi)     // Quotient %eax ablegen, Rest in %edx behalten
              decl    %ecx
              jnz     L(dld1)
L(dld2:)    movl    %edx,%eax       // Ergebnis := letzter Rest
            popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            ret

// extern uintD divucopy_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
            ALIGN
            DECLARE_FUNCTION(divucopy_loop_down)
C(divucopy_loop_down:)
            pushl   %edi            // %edi retten
            pushl   %esi            // %esi retten
            pushl   %ebx            // %ebx retten
            movl    16(%esp),%ebx   // %ebx = digit
            movl    20(%esp),%esi   // %esi = sourceptr
            movl    24(%esp),%edi   // %edi = destptr
            movl    28(%esp),%ecx   // %ecx = len
            xorl    %edx,%edx       // %edx = Rest := 0
            jecxz   L(dcld2)        // %ecx = 0 ?
            subl    %edi,%esi
L(dcld1:)     leal    -4(%edi),%edi   // sourceptr--, destptr--
              movl    (%esi,%edi),%eax // nächstes Digit *ptr
              divl    %ebx            // Division von %edx|%eax durch %ebx
              movl    %eax,(%edi)     // Quotient %eax ablegen, Rest in %edx behalten
              decl    %ecx
              jnz     L(dcld1)
L(dcld2:)   movl    %edx,%eax       // Ergebnis := letzter Rest
            popl    %ebx            // %ebx zurück
            popl    %esi            // %esi zurück
            popl    %edi            // %edi zurück
            ret

#endif

// extern void shiftxor_loop_up (uintD* xptr, const uintD* yptr, uintC count, uintC i);
            ALIGN
            DECLARE_FUNCTION(shiftxor_loop_up)
C(shiftxor_loop_up:)
            pushl   %esi            // %esi retten
            pushl   %edi            // %edi retten
            pushl   %ebx            // %ebx retten
            movl    16(%esp),%esi   // %esi = xptr
            movl    20(%esp),%edi   // %edi = yptr
            movl    24(%esp),%edx   // count
            movb    28(%esp),%cl    // i
            orl     %edx,%edx       // count = 0 ?
            jz      L(shxlu4)
            subl    %esi,%edi
            // erstes Digit shiften:
            movl    (%esi,%edi),%ebx // *yptr in %ebx halten
            movl    %ebx,%eax       // und in %eax rechnen:
            shll    %cl,%eax        // um i Bits links shiften, rechts Nullen rein
            xorl    %eax,(%esi)     // und mit *xptr verknüpfen und ablegen
            leal    4(%esi),%esi    // sourceptr++, destptr++
            // Letztes Digit in %ebx.
            negb    %cl             // 32-i
            decl    %edx
            jz      L(shxlu2)
L(shxlu1:)    // weiteres Digit shiften:
              movl    (%esi,%edi),%eax // nächstes Digit nach %eax
              shrdl   shcl %eax,%ebx  // %ebx um %cl=32-i Bits rechts shiften, %eax von links reinshiften
              xorl    %ebx,(%esi)     // %ebx mit *xptr verknüpfen und ablegen
              leal    4(%esi),%esi    // xptr++, yptr++
              // Letztes Digit in %eax.
              decl    %edx
              jz      L(shxlu3)
              // weiteres Digit shiften:
              movl    (%esi,%edi),%ebx // nächstes Digit nach %ebx
              shrdl   shcl %ebx,%eax  // %eax um %cl=32-i Bits rechts shiften, %ebx von links reinshiften
              xorl    %eax,(%esi)     // %eax mit *xptr verknüpfen und ablegen
              leal    4(%esi),%esi    // xptr++, yptr++
              // Letztes Digit in %ebx.
              decl    %edx
              jnz     L(shxlu1)
L(shxlu2:)  movl    %ebx,%eax
L(shxlu3:)  shrl    %cl,%eax        // %eax um 32-i Bits nach rechts shiften
            xorl    %eax,(%esi)     // und mit *xptr verknüpfen und ablegen
L(shxlu4:)  popl    %ebx            // %ebx zurück
            popl    %edi            // %edi zurück
            popl    %esi            // %esi zurück
            ret

