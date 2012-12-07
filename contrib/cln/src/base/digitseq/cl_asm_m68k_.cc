// Externe Routinen zu ARILEV1.D
// Prozessor: 680x0 mit x>=2
// Assembler-Syntax: meist "$" streichen, auf A/UX "$" durch "%" ersetzen
// Compiler: CC oder GNU-C auf SUN3 oder AMIGA oder A/UX
// Parameter-Übergabe:
//   auf dem Stack: sp@(4), sp@(8), ... (.W-Größen belegen 4 Byte!),
//   Rückgabewert in d0.
// Register a0-a1,d0-d1 frei verwendbar,
// Register a2-a4,d2-d7 müssen gerettet werden.
// Einstellungen: intCsize=16, intDsize=32.

  #ifdef ASM_UNDERSCORE
    #if defined(__STDC__) || defined (__cplusplus)
      #define C(entrypoint) _##entrypoint
    #else
      #define C(entrypoint) _/**/entrypoint
    #endif
  #else
    #define C(entrypoint) entrypoint
  #endif

  // Befehl, der das X- und das C-Bit löscht (und evtl. d1 modifiziert):
  #if defined(sun)
    // SUN-Assembler
    #define clrx   subw $d1,$d1
  #else
    // GNU-Assembler
    #if defined(__STDC__) || defined (__cplusplus)
      // Some preprocessors keep the backslash in place, some don't.
      // In any case, we will filter it out later.
      #define clrx   andb \#0x0e,$ccr
    #else
      #define clrx   andb #0x0e,$ccr
    #endif
  #endif

  // When this file is compiled into a shared library, ELF linkers need to
  // know which symbols are functions.
  #if defined(__NetBSD__) || defined(__OpenBSD__) || defined(__ELF__) || defined(__svr4__)
    #define DECLARE_FUNCTION(name) .type C(name),@function
  #else
    #define DECLARE_FUNCTION(name)
  #endif

           .text

           .globl C(copy_loop_up),C(copy_loop_down),C(fill_loop_up),C(fill_loop_down)
           .globl C(clear_loop_up),C(clear_loop_down)
           .globl C(test_loop_up),C(test_loop_down)
           .globl C(xor_loop_up),C(compare_loop_up),C(shiftleftcopy_loop_up),C(shiftxor_loop_up)
#if CL_DS_BIG_ENDIAN_P
           .globl C(or_loop_up),C(and_loop_up),C(eqv_loop_up)
           .globl C(nand_loop_up),C(nor_loop_up),C(andc2_loop_up),C(orc2_loop_up)
           .globl C(not_loop_up)
           .globl C(and_test_loop_up)
           .globl C(add_loop_down),C(addto_loop_down),C(inc_loop_down)
           .globl C(sub_loop_down),C(subx_loop_down),C(subfrom_loop_down),C(dec_loop_down)
           .globl C(neg_loop_down)
           .globl C(shift1left_loop_down),C(shiftleft_loop_down),C(shiftleftcopy_loop_down)
           .globl C(shift1right_loop_up),C(shiftright_loop_up),C(shiftrightsigned_loop_up),C(shiftrightcopy_loop_up)
           .globl C(mulusmall_loop_down),C(mulu_loop_down),C(muluadd_loop_down),C(mulusub_loop_down)
           .globl C(divu_loop_up),C(divucopy_loop_up)
#else
           .globl C(or_loop_down),C(xor_loop_down),C(and_loop_down),C(eqv_loop_down)
           .globl C(nand_loop_down),C(nor_loop_down),C(andc2_loop_down),C(orc2_loop_down)
           .globl C(not_loop_down)
           .globl C(and_test_loop_down),C(compare_loop_down)
           .globl C(add_loop_up),C(addto_loop_up),C(inc_loop_up)
           .globl C(sub_loop_up),C(subx_loop_up),C(subfrom_loop_up),C(dec_loop_up)
           .globl C(neg_loop_up)
           .globl C(shift1left_loop_up),C(shiftleft_loop_up)
           .globl C(shift1right_loop_down),C(shiftright_loop_down),C(shiftrightsigned_loop_down),C(shiftrightcopy_loop_down)
           .globl C(mulusmall_loop_up),C(mulu_loop_up),C(muluadd_loop_up),C(mulusub_loop_up)
           .globl C(divu_loop_down),C(divucopy_loop_down)
#endif

#ifndef __GNUC__ /* mit GNU-C machen wir mulu32() als Macro, der inline multipliziert */
           .globl C(mulu32_)
| extern struct { uint32 lo; uint32 hi; } mulu32_ (uint32 arg1, uint32 arg2);
| 2^32*hi+lo := arg1*arg2.
           DECLARE_FUNCTION(mulu32_)
C(mulu32_:) | Input in d0,d1, Output in d0,mulu32_high
           movel $sp@(4),$d0
           movel $sp@(8),$d1
           mulul $d1,$d1:$d0
           movel $d1,(C(mulu32_high)) | Adressierung?? Deklaration??
           rts
#endif

#ifndef __GNUC__ /* mit GNU-C machen wir divu_6432_3232() als Macro, der inline dividiert */
           .globl C(divu_6432_3232_)
| extern struct { uint32 q; uint32 r; } divu_6432_3232_ (uint32 xhi, uint32 xlo, uint32 y);
| x = 2^32*xhi+xlo = q*y+r schreiben. Sei bekannt, daß 0 <= x < 2^32*y .
           DECLARE_FUNCTION(divu_6432_3232_)
C(divu_6432_3232_:) | Input in d1,d0,d2, Output in d0,divu_32_rest
           movel $sp@(4),$d1
           movel $sp@(8),$d0
           divul $sp@(12),$d1:$d0 | x = d1|d0 durch y dividieren
           movel $d1,(C(divu_32_rest)) | Rest ablegen | Adressierung?? Deklaration??
           rts
#endif

| extern uintD* copy_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
           DECLARE_FUNCTION(copy_loop_up)
C(copy_loop_up:) | Input in a0,a1,d0.W, Output in d0
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a0@+,$a1@+
    2:       dbra $d0,1b
           movel $a1,$d0
           rts

| extern uintD* copy_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
           DECLARE_FUNCTION(copy_loop_down)
C(copy_loop_down:) | Input in a0,a1,d0.W, Output in d0
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a0@-,$a1@-
    2:       dbra $d0,1b
           movel $a1,$d0
           rts

| extern uintD* fill_loop_up (uintD* destptr, uintC count, uintD filler);
           DECLARE_FUNCTION(fill_loop_up)
C(fill_loop_up:) | Input in a0,d0.W,d1, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           movel $sp@(12),$d1
           bras 2f
    1:       movel $d1,$a0@+
    2:       dbra $d0,1b
           movel $a0,$d0
           rts

| extern uintD* fill_loop_down (uintD* destptr, uintC count, uintD filler);
           DECLARE_FUNCTION(fill_loop_down)
C(fill_loop_down:) | Input in a0,d0.W,d1, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           movel $sp@(12),$d1
           bras 2f
    1:       movel $d1,$a0@-
    2:       dbra $d0,1b
           movel $a0,$d0
           rts

| extern uintD* clear_loop_up (uintD* destptr, uintC count);
           DECLARE_FUNCTION(clear_loop_up)
C(clear_loop_up:) | Input in a0,d0.W, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           bras 2f
    1:       clrl $a0@+
    2:       dbra $d0,1b
           movel $a0,$d0
           rts

| extern uintD* clear_loop_down (uintD* destptr, uintC count);
           DECLARE_FUNCTION(clear_loop_down)
C(clear_loop_down:) | Input in a0,d0.W, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           bras 2f
    1:       clrl $a0@-
    2:       dbra $d0,1b
           movel $a0,$d0
           rts

| extern boolean test_loop_up (uintD* ptr, uintC count);
           DECLARE_FUNCTION(test_loop_up)
C(test_loop_up:) | Input in a0,d0.W, Output in d0.W=d0.L
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           bras 2f
    1:       tstl $a0@+
             bnes 3f
    2:       dbra $d0,1b
           clrl $d0
           rts
    3:     moveq #1,$d0
           rts

| extern boolean test_loop_down (uintD* ptr, uintC count);
           DECLARE_FUNCTION(test_loop_down)
C(test_loop_down:) | Input in a0,d0.W, Output in d0.W=d0.L
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           bras 2f
    1:       tstl $a0@-
             bnes 3f
    2:       dbra $d0,1b
           clrl $d0
           rts
    3:     moveq #1,$d0
           rts

#if CL_DS_BIG_ENDIAN_P

| extern void or_loop_up (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(or_loop_up)
C(or_loop_up:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@+,$d1
             orl $d1,$a0@+
    2:       dbra $d0,1b
           rts

#endif

| extern void xor_loop_up (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(xor_loop_up)
C(xor_loop_up:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@+,$d1
             eorl $d1,$a0@+
    2:       dbra $d0,1b
           rts

#if CL_DS_BIG_ENDIAN_P

| extern void and_loop_up (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(and_loop_up)
C(and_loop_up:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@+,$d1
             andl $d1,$a0@+
    2:       dbra $d0,1b
           rts

| extern void eqv_loop_up (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(eqv_loop_up)
C(eqv_loop_up:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@+,$d1
             eorl $d1,$a0@
             notl $a0@+
    2:       dbra $d0,1b
           rts

| extern void nand_loop_up (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(nand_loop_up)
C(nand_loop_up:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@+,$d1
             andl $d1,$a0@
             notl $a0@+
    2:       dbra $d0,1b
           rts

| extern void nor_loop_up (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(nor_loop_up)
C(nor_loop_up:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@+,$d1
             orl $d1,$a0@
             notl $a0@+
    2:       dbra $d0,1b
           rts

| extern void andc2_loop_up (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(andc2_loop_up)
C(andc2_loop_up:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@+,$d1
             notl $d1
             andl $d1,$a0@+
    2:       dbra $d0,1b
           rts

| extern void orc2_loop_up (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(orc2_loop_up)
C(orc2_loop_up:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@+,$d1
             notl $d1
             orl $d1,$a0@+
    2:       dbra $d0,1b
           rts

| extern void not_loop_up (uintD* xptr, uintC count);
           DECLARE_FUNCTION(not_loop_up)
C(not_loop_up:) | Input in a0,d0.W
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           bras 2f
    1:       notl $a0@+
    2:       dbra $d0,1b
           rts

| extern boolean and_test_loop_up (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(and_test_loop_up)
C(and_test_loop_up:) | Input in a0,a1,d0.W, verändert d1, Output in d0.W=d0.L
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a0@+,$d1
             andl $a1@+,$d1
             bnes 3f
    2:       dbra $d0,1b
           clrl $d0
           rts
    3:     moveq #1,$d0
           rts

#endif

| extern cl_signean compare_loop_up (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(compare_loop_up)
C(compare_loop_up:) | Input in a0,a1,d0.W, Output in d0.W=d0.L
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       cmpml $a1@+,$a0@+
             bnes 3f
    2:       dbra $d0,1b
           clrl $d0
           rts
    3:     bcss 4f
           moveq #1,$d0
           rts
    4:     moveq #-1,$d0
           rts

#if CL_DS_BIG_ENDIAN_P

| extern uintD add_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
           DECLARE_FUNCTION(add_loop_down)
C(add_loop_down:) | Input in a0,a1,a2,d0.W, verändert d1,d2, Output in d0
           moveml $a2/$d2,$sp@-
           movel $sp@(8+4),$a0
           movel $sp@(8+8),$a1
           movel $sp@(8+12),$a2
           movew $sp@(8+16+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       movel $a0@-,$d1
             movel $a1@-,$d2
             addxl $d2,$d1
             movel $d1,$a2@-
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           moveml $sp@+,$a2/$d2
           rts

| extern uintD addto_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
           DECLARE_FUNCTION(addto_loop_down)
C(addto_loop_down:) | Input in a0,a1,d0.W, Output in d0
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       addxl $a0@-,$a1@-
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           rts

| extern uintD inc_loop_down (uintD* ptr, uintC count);
           DECLARE_FUNCTION(inc_loop_down)
C(inc_loop_down:) | Input in a0,d0.W, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           dbra $d0,1f          | simuliere gesetzten Carry
           moveq #-1,$d0     | d0.L=-1 für Übertrag
           rts
    1:       addql #1,$a0@-
             dbcc $d0,1b
           subxl $d0,$d0       | kein Carry -> d0.L=0, sonst d0.L=-1 für Übertrag
           rts

| extern uintD sub_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
           DECLARE_FUNCTION(sub_loop_down)
C(sub_loop_down:) | Input in a0,a1,a2,d0.W, verändert d1,d2, Output in d0
           moveml $a2/$d2,$sp@-
           movel $sp@(8+4),$a0
           movel $sp@(8+8),$a1
           movel $sp@(8+12),$a2
           movew $sp@(8+16+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       movel $a0@-,$d1
             movel $a1@-,$d2
             subxl $d2,$d1
             movel $d1,$a2@-
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           moveml $sp@+,$a2/$d2
           rts

| extern uintD subx_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);
           DECLARE_FUNCTION(subx_loop_down)
C(subx_loop_down:) | Input in a0,a1,a2,d0.W,d1, verändert d2, Output in d0
           moveml $a2/$d2,$sp@-
           movel $sp@(8+4),$a0
           movel $sp@(8+8),$a1
           movel $sp@(8+12),$a2
           movew $sp@(8+16+2),$d0
           movel $sp@(8+20),$d1
           roxrl #1,$d1      | X-Bit initialisieren
           bras 2f
    1:       movel $a0@-,$d1
             movel $a1@-,$d2
             subxl $d2,$d1
             movel $d1,$a2@-
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           moveml $sp@+,$a2/$d2
           rts

| extern uintD subfrom_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
           DECLARE_FUNCTION(subfrom_loop_down)
C(subfrom_loop_down:) | Input in a0,a1,d0.W, Output in d0
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       subxl $a0@-,$a1@-
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           rts

| extern uintD dec_loop_down (uintD* ptr, uintC count);
           DECLARE_FUNCTION(dec_loop_down)
C(dec_loop_down:) | Input in a0,d0.W, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           dbra $d0,1f          | simuliere gesetzten Carry
           moveq #-1,$d0     | d0.L=-1 als Übertrag
           rts
    1:       subql #1,$a0@-
             dbcc $d0,1b       | kein Carry -> Schleife abbrechen
           subxl $d0,$d0       | kein Carry -> d0.L=0, sonst d0.L=-1 als Übertrag
           rts

| extern uintD neg_loop_down (uintD* ptr, uintC count);
           DECLARE_FUNCTION(neg_loop_down)
C(neg_loop_down:) | Input in a0,d0.W, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       negxl $a0@-
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           rts

| extern uintD shift1left_loop_down (uintD* ptr, uintC count);
           DECLARE_FUNCTION(shift1left_loop_down)
C(shift1left_loop_down:) | Input in a0,d0.W, Output in d0.L
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       roxlw $a0@-     | Digit a0@- um 1 Bit links schieben, X-Bit als Buffer
             roxlw $a0@-
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           rts

| extern uintD shiftleft_loop_down (uintD* ptr, uintC count, uintC i, uintD carry);
           DECLARE_FUNCTION(shiftleft_loop_down)
C(shiftleft_loop_down:) | Input in a0,d0.W,d1.W,d2, Output in d0
#if 1
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           movel $sp@(16+16),$d2
           moveq #32,$d5
           subw $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.W = 32-i,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       movel $a0@-,$d3  | d3.L = neues Digit
             movel $d3,$d4
             lsll $d1,$d4      | um i Bits nach links schieben
             orl $d2,$d4       | mit vorigem Übertrag kombinieren
             movel $d4,$a0@   | 32 Bits ablegen
             movel $d3,$d2
             lsrl $d5,$d2      | neuen Übertrag bilden
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           movel $sp@(16+16),$d2
           moveq #1,$d5
           lsll $d1,$d5
           subql #1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.L = 2^i-1
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       movel $a0@-,$d3  | d3.L = neues Digit
             roll $d1,$d3      | um i Bits links rotieren
             movel $d3,$d4
             andl $d5,$d3      | untere i Bits in d3
             eorl $d3,$d4      | obere 32-i Bits in d4
             orl $d2,$d4       | mit vorigem übertrag kombinieren
             movel $d4,$a0@   | 32 Bits ablegen
             movel $d3,$d2     | neuer Übertrag
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#endif

| extern uintD shiftleftcopy_loop_down (uintD* sourceptr, uintD* destptr, uintC count, uintC i);
           DECLARE_FUNCTION(shiftleftcopy_loop_down)
C(shiftleftcopy_loop_down:) | Input in a0,a1,d0.W,d1.W, Output in d0
#if 1
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movel $sp@(16+8),$a1
           movew $sp@(16+12+2),$d0
           movew $sp@(16+16+2),$d1
           moveq #32,$d5
           subw $d1,$d5
           clrl $d2
           | a0 = sourceptr, a1 = destptr, d0.W = count, d1.W = i, d5.W = 32-i,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       movel $a0@-,$d3  | d3.L = neues Digit
             movel $d3,$d4
             lsll $d1,$d4      | um i Bits nach links schieben
             orl $d2,$d4       | mit vorigem Übertrag kombinieren
             movel $d4,$a1@-  | 32 Bits ablegen
             movel $d3,$d2
             lsrl $d5,$d2      | neuen Übertrag bilden
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movel $sp@(16+8),$a1
           movew $sp@(16+12+2),$d0
           movew $sp@(16+16+2),$d1
           moveq #1,$d5
           lsll $d1,$d5
           subql #1,$d5
           | a0 = sourceptr, a1 = destptr, d0.W = count, d1.W = i, d5.L = 2^i-1
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       movel $a0@-,$d3  | d3.L = neues Digit
             roll $d1,$d3      | um i Bits links rotieren
             movel $d3,$d4
             andl $d5,$d3      | untere i Bits in d3
             eorl $d3,$d4      | obere 32-i Bits in d4
             orl $d2,$d4       | mit vorigem übertrag kombinieren
             movel $d4,$a1@-  | 32 Bits ablegen
             movel $d3,$d2     | neuer Übertrag
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#endif

| extern uintD shift1right_loop_up (uintD* ptr, uintC count, uintD carry);
           DECLARE_FUNCTION(shift1right_loop_up)
C(shift1right_loop_up:) | Input in a0,d0.W,d1, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           movel $sp@(12),$d1
           roxrl #1,$d1       | X-Bit löschen oder setzen, je nach d1
           bras 2f
    1:       roxrw $a0@+     | Digit a0@+ um 1 Bit rechts schieben, X-Bit als Buffer
             roxrw $a0@+
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           rts

| extern uintD shiftright_loop_up (uintD* ptr, uintC count, uintC i);
           DECLARE_FUNCTION(shiftright_loop_up)
C(shiftright_loop_up:) | Input in a0,d0.W,d1.W, Output in d0
#if 1
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           moveq #32,$d5
           subw $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.W = 32-i,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           clrl $d2
           bras 2f
    1:       | a0 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.W = 32-i,
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@,$d3   | neue Daten
             movel $d3,$d4
             lsrl $d1,$d3      | um i Bits rechts schieben
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
             movel $d3,$a0@+  | ablegen
             lsll $d5,$d4      | um (32-i) Bits links geschoben
             movel $d4,$d2     | liefert neuen Übertrag
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           moveq #-1,$d5
           lsrl $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.L = 2^(32-i)-1,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           clrl $d2
           bras 2f
    1:       | a0 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.L = 2^(32-i)-1,
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@,$d3   | neue Daten
             rorl $d1,$d3      | um i Bits rechts rotieren
             movel $d3,$d4
             andl $d5,$d3      | untere 32-i Bits
             eorl $d3,$d4      | obere i Bits
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
             movel $d4,$d2     | neuer Übertrag
             movel $d3,$a0@+  | ablegen
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#endif

| extern uintD shiftrightsigned_loop_up (uintD* ptr, uintC count, uintC i);
           DECLARE_FUNCTION(shiftrightsigned_loop_up)
C(shiftrightsigned_loop_up:) | Input in a0,d0.W,d1.W, Output in d0
#if 1
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           moveq #32,$d5
           subw $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.W = 32-i,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           subqw #1,$d0
           movel $a0@,$d3     | erstes Digit
           movel $d3,$d4
           asrl $d1,$d3        | um i Bits rechts schieben
           bras 2f
    1:       | a0 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.W = 32-i,
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@,$d3   | neue Daten
             movel $d3,$d4
             lsrl $d1,$d3      | um i Bits rechts schieben
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
    2:       movel $d3,$a0@+  | ablegen
             lsll $d5,$d4      | um (32-i) Bits links geschoben
             movel $d4,$d2     | liefert neuen Übertrag
             dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           moveq #-1,$d5
           lsrl $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.L = 2^(32-i)-1,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           subqw #1,$d0
           movel $a0@,$d3     | erstes Digit
           movel $d3,$d4
           rorl $d1,$d4        | um i Bits rechts rotieren
           movel $d5,$d2
           notl $d2
           andl $d4,$d2        | obere 32-i Bits
           asrl $d1,$d3        | erstes Digit um i Bits rechts shiften
           bras 2f
    1:       | a0 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.L = 2^(32-i)-1,
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@,$d3   | neue Daten
             rorl $d1,$d3      | um i Bits rechts rotieren
             movel $d3,$d4
             andl $d5,$d3      | untere 32-i Bits
             eorl $d3,$d4      | obere i Bits
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
             movel $d4,$d2     | neuer Übertrag
    2:       movel $d3,$a0@+  | ablegen
             dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#endif

| extern uintD shiftrightcopy_loop_up (uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry);
           DECLARE_FUNCTION(shiftrightcopy_loop_up)
C(shiftrightcopy_loop_up:) | Input in a0,a1,d0.W,d1.W,d2, Output in d0
#if 1
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movel $sp@(16+8),$a1
           movew $sp@(16+12+2),$d0
           movew $sp@(16+16+2),$d1
           movel $sp@(16+20),$d2
           moveq #32,$d5
           subw $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.W = 32-i
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       | a0,a1 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.W = 32-i
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@+,$d3  | neue Daten
             movel $d3,$d4
             lsrl $d1,$d3      | um i Bits rechts schieben
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
             movel $d3,$a1@+  | ablegen
             movel $d4,$d2
    2:       lsll $d5,$d2      | um (32-i) Bits links geschoben, gibt neuen Übertrag
             dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movel $sp@(16+8),$a1
           movew $sp@(16+12+2),$d0
           movew $sp@(16+16+2),$d1
           movel $sp@(16+20),$d2
           moveq #-1,$d5
           lsrl $d1,$d5
           rorl $d1,$d2
           | a0 = ptr, d0.W = count, d1.W = i, d5.L = 2^(32-i)-1
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       | a0,a1 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.L = 2^(32-i)-1
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@+,$d3  | neue Daten
             rorl $d1,$d3      | um i Bits rechts rotieren
             movel $d3,$d4
             andl $d5,$d3      | untere 32-i Bits
             eorl $d3,$d4      | obere i Bits
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
             movel $d4,$d2     | neuer Übertrag
             movel $d3,$a1@+  | ablegen
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#endif

| extern uintD mulusmall_loop_down (uintD digit, uintD* ptr, uintC len, uintD newdigit);
           DECLARE_FUNCTION(mulusmall_loop_down)
C(mulusmall_loop_down:) | Input in d0,a0,d1.W,d2, Output in d0
           moveml $d2-$d4,$sp@-
           movel $sp@(12+4),$d0
           movel $sp@(12+8),$a0
           movew $sp@(12+12+2),$d1
           movel $sp@(12+16),$d2
           addw #0,$d1        | X-Bit löschen
           bras 2f
    1:       movel $a0@-,$d3  | nächstes Digit
             mulul $d0,$d4:$d3  | mit digit multiplizieren
             addxl $d2,$d3     | und bisherigen Carry und X-Bit addieren
             movel $d3,$a0@   | Low-Digit ablegen
             movel $d4,$d2     | High-Digit gibt neuen Carry
    2:       dbra $d1,1b
           clrl $d0
           addxl $d2,$d0       | letzter Carry (incl. X-Bit)
           moveml $sp@+,$d2-$d4
           rts

| extern void mulu_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
           DECLARE_FUNCTION(mulu_loop_down)
C(mulu_loop_down:) | Input in d0,a0,a1,d1.W
#if 1
           moveml $d2-$d4,$sp@-
           movel $sp@(12+4),$d0
           movel $sp@(12+8),$a0
           movel $sp@(12+12),$a1
           movew $sp@(12+16+2),$d1
           subl $d2,$d2        | carry := 0, X-Bit löschen
           bras 2f
    1:       movel $a0@-,$d3  | nächstes Digit
             mulul $d0,$d4:$d3  | mit digit multiplizieren
             addxl $d2,$d3     | und bisherigen Carry und X-Bit addieren
             movel $d3,$a1@-  | Low-Digit ablegen
             movel $d4,$d2     | High-Digit gibt neuen Carry
    2:       dbra $d1,1b
           clrl $d3
           addxl $d3,$d2       | letztes X-Bit verarbeiten
           movel $d2,$a1@-    | letzten Carry ablegen
           moveml $sp@+,$d2-$d4
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$d0
           movel $sp@(16+8),$a0
           movel $sp@(16+12),$a1
           movew $sp@(16+16+2),$d1
           clrl $d5           | 0
           clrl $d2           | carry
           bras 2f
    1:       movel $a0@-,$d3  | nächstes Digit
             mulul $d0,$d4:$d3  | mit digit multiplizieren
             addl $d2,$d3      | und bisherigen Carry addieren
             addxl $d5,$d4
             movel $d3,$a1@-  | Low-Digit ablegen
             movel $d4,$d2     | High-Digit gibt neuen Carry
    2:       dbra $d1,1b
           movel $d2,$a1@-    | letzten Carry ablegen
           moveml $sp@+,$d2-$d5
           rts
#endif

| extern uintD muluadd_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
           DECLARE_FUNCTION(muluadd_loop_down)
C(muluadd_loop_down:) | Input in d0,a0,a1,d1.W, Output in d0
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$d0
           movel $sp@(16+8),$a0
           movel $sp@(16+12),$a1
           movew $sp@(16+16+2),$d1
           clrl $d5           | 0
           subl $d2,$d2        | carry := 0, X-Bit löschen
           bras 2f
    1:       movel $a0@-,$d3  | nächstes Digit
             mulul $d0,$d4:$d3  | mit digit multiplizieren
             addxl $d2,$d3     | und bisherigen Carry und X-Bit addieren
             addxl $d5,$d4
             addl $d3,$a1@-   | Low-Digit zum dest-Digit addieren, X als Übertrag
             movel $d4,$d2     | High-Digit gibt neuen Carry
    2:       dbra $d1,1b
           addxl $d5,$d2       | letztes X-Bit addieren
           movel $d2,$d0       | letzten Carry als Ergebnis
           moveml $sp@+,$d2-$d5
           rts

| extern uintD mulusub_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
           DECLARE_FUNCTION(mulusub_loop_down)
C(mulusub_loop_down:) | Input in d0,a0,a1,d1.W, Output in d0
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$d0
           movel $sp@(16+8),$a0
           movel $sp@(16+12),$a1
           movew $sp@(16+16+2),$d1
           clrl $d5           | 0
           subl $d2,$d2        | carry := 0, X-Bit löschen
           bras 2f
    1:       movel $a0@-,$d3  | nächstes Digit
             mulul $d0,$d4:$d3  | mit digit multiplizieren
             addxl $d2,$d3     | und bisherigen Carry und X-Bit addieren
             addxl $d5,$d4
             subl $d3,$a1@-   | Low-Digit vom dest-Digit subtrahieren, X als Übertrag
             movel $d4,$d2     | High-Digit gibt neuen Carry
    2:       dbra $d1,1b
           clrl $d0
           addxl $d2,$d0       | letzter Carry und letztes X-Bit
           moveml $sp@+,$d2-$d5
           rts

| extern uintD divu_loop_up (uintD digit, uintD* ptr, uintC len);
           DECLARE_FUNCTION(divu_loop_up)
C(divu_loop_up:) | Input in d0,a0,d1.W, Output in d0
           moveml $d2-$d3,$sp@-
           movel $sp@(8+4),$d0
           movel $sp@(8+8),$a0
           movew $sp@(8+12+2),$d1
           clrl $d2           | Rest := 0
           bras 2f
    1:       movel $a0@,$d3   | nächst-niedriges Digit
             divul $d0,$d2:$d3  | mit Rest kombinieren und durch digit dividieren
             movel $d3,$a0@+  | Quotient ablegen, Rest in d2
    2:       dbra $d1,1b
           movel $d2,$d0       | Rest
           moveml $sp@+,$d2-$d3
           rts

| extern uintD divucopy_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
           DECLARE_FUNCTION(divucopy_loop_up)
C(divucopy_loop_up:) | Input in d0,a0,a1,d1.W, Output in d0
           moveml $d2-$d3,$sp@-
           movel $sp@(8+4),$d0
           movel $sp@(8+8),$a0
           movel $sp@(8+12),$a1
           movew $sp@(8+16+2),$d1
           clrl $d2           | Rest := 0
           bras 2f
    1:       movel $a0@+,$d3  | nächst-niedriges Digit
             divul $d0,$d2:$d3  | mit Rest kombinieren und durch digit dividieren
             movel $d3,$a1@+  | Quotient ablegen, Rest in d2
    2:       dbra $d1,1b
           movel $d2,$d0       | Rest
           moveml $sp@+,$d2-$d3
           rts

#endif

#if !CL_DS_BIG_ENDIAN_P

| extern void or_loop_down (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(or_loop_down)
C(or_loop_down:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@-,$d1
             orl $d1,$a0@-
    2:       dbra $d0,1b
           rts

| extern void xor_loop_down (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(xor_loop_down)
C(xor_loop_down:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@-,$d1
             eorl $d1,$a0@-
    2:       dbra $d0,1b
           rts

| extern void and_loop_down (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(and_loop_down)
C(and_loop_down:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@-,$d1
             andl $d1,$a0@-
    2:       dbra $d0,1b
           rts

| extern void eqv_loop_down (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(eqv_loop_down)
C(eqv_loop_down:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@-,$d1
             eorl $d1,$a0@-
             notl $a0@
    2:       dbra $d0,1b
           rts

| extern void nand_loop_down (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(nand_loop_down)
C(nand_loop_down:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@-,$d1
             andl $d1,$a0@-
             notl $a0@
    2:       dbra $d0,1b
           rts

| extern void nor_loop_down (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(nor_loop_down)
C(nor_loop_down:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@-,$d1
             orl $d1,$a0@-
             notl $a0@
    2:       dbra $d0,1b
           rts

| extern void andc2_loop_down (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(andc2_loop_down)
C(andc2_loop_down:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@-,$d1
             notl $d1
             andl $d1,$a0@-
    2:       dbra $d0,1b
           rts

| extern void orc2_loop_down (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(orc2_loop_down)
C(orc2_loop_down:) | Input in a0,a1,d0.W, verändert d1
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a1@-,$d1
             notl $d1
             orl $d1,$a0@-
    2:       dbra $d0,1b
           rts

| extern void not_loop_down (uintD* xptr, uintC count);
           DECLARE_FUNCTION(not_loop_down)
C(not_loop_down:) | Input in a0,d0.W
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           bras 2f
    1:       notl $a0@-
    2:       dbra $d0,1b
           rts

| extern boolean and_test_loop_down (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(and_test_loop_down)
C(and_test_loop_down:) | Input in a0,a1,d0.W, verändert d1, Output in d0.W=d0.L
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           bras 2f
    1:       movel $a0@-,$d1
             andl $a1@-,$d1
             bnes 3f
    2:       dbra $d0,1b
           clrl $d0
           rts
    3:     moveq #1,$d0
           rts

| extern cl_signean compare_loop_down (uintD* xptr, uintD* yptr, uintC count);
           DECLARE_FUNCTION(compare_loop_down)
C(compare_loop_down:) | Input in a0,a1,d0.W, verändert d1,d2, Output in d0.W=d0.L
           movel $d2,$sp@-
           movel $sp@(4+4),$a0
           movel $sp@(4+8),$a1
           movew $sp@(4+12+2),$d0
           bras 2f
    1:       movel $a0@-,$d1
             movel $a1@-,$d2
             cmpl $d2,$d1
             bnes 3f
    2:       dbra $d0,1b
           clrl $d0
           movel $sp@+,$d2
           rts
    3:     bcss 4f
           moveq #1,$d0
           movel $sp@+,$d2
           rts
    4:     moveq #-1,$d0
           movel $sp@+,$d2
           rts

| extern uintD add_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
           DECLARE_FUNCTION(add_loop_up)
C(add_loop_up:) | Input in a0,a1,a2,d0.W, verändert d1,d2, Output in d0
           moveml $a2/$d2,$sp@-
           movel $sp@(8+4),$a0
           movel $sp@(8+8),$a1
           movel $sp@(8+12),$a2
           movew $sp@(8+16+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       movel $a0@+,$d1
             movel $a1@+,$d2
             addxl $d2,$d1
             movel $d1,$a2@+
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           moveml $sp@+,$a2/$d2
           rts

| extern uintD addto_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
           DECLARE_FUNCTION(addto_loop_up)
C(addto_loop_up:) | Input in a0,a1,d0.W, verändert d1,d2, Output in d0
           movel $d2,$sp@-
           movel $sp@(4+4),$a0
           movel $sp@(4+8),$a1
           movew $sp@(4+12+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       movel $a0@+,$d1
             movel $a1@,$d2
             addxl $d1,$d2
             movel $d2,$a1@+
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           movel $sp@+,$d2
           rts

| extern uintD inc_loop_up (uintD* ptr, uintC count);
           DECLARE_FUNCTION(inc_loop_up)
C(inc_loop_up:) | Input in a0,d0.W, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           dbra $d0,1f          | simuliere gesetzten Carry
           moveq #-1,$d0     | d0.L=-1 für Übertrag
           rts
    1:       addql #1,$a0@+
             dbcc $d0,1b
           subxl $d0,$d0       | kein Carry -> d0.L=0, sonst d0.L=-1 für Übertrag
           rts

| extern uintD sub_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
           DECLARE_FUNCTION(sub_loop_up)
C(sub_loop_up:) | Input in a0,a1,a2,d0.W, verändert d1,d2, Output in d0
           moveml $a2/$d2,$sp@-
           movel $sp@(8+4),$a0
           movel $sp@(8+8),$a1
           movel $sp@(8+12),$a2
           movew $sp@(8+16+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       movel $a0@+,$d1
             movel $a1@+,$d2
             subxl $d2,$d1
             movel $d1,$a2@+
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           moveml $sp@+,$a2/$d2
           rts

| extern uintD subx_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);
           DECLARE_FUNCTION(subx_loop_up)
C(subx_loop_up:) | Input in a0,a1,a2,d0.W,d1, verändert d2, Output in d0
           moveml $a2/$d2,$sp@-
           movel $sp@(8+4),$a0
           movel $sp@(8+8),$a1
           movel $sp@(8+12),$a2
           movew $sp@(8+16+2),$d0
           movel $sp@(8+20),$d1
           roxrl #1,$d1      | X-Bit initialisieren
           bras 2f
    1:       movel $a0@+,$d1
             movel $a1@+,$d2
             subxl $d2,$d1
             movel $d1,$a2@+
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           moveml $sp@+,$a2/$d2
           rts

| extern uintD subfrom_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
           DECLARE_FUNCTION(subfrom_loop_up)
C(subfrom_loop_up:) | Input in a0,a1,d0.W, Output in d0
           movel $sp@(4),$a0
           movel $sp@(8),$a1
           movew $sp@(12+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       movel $a0@+,$d1
             movel $a1@,$d2
             subxl $d1,$d2
             movel $d2,$a1@+
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           rts

| extern uintD dec_loop_up (uintD* ptr, uintC count);
           DECLARE_FUNCTION(dec_loop_up)
C(dec_loop_up:) | Input in a0,d0.W, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           dbra $d0,1f          | simuliere gesetzten Carry
           moveq #-1,$d0     | d0.L=-1 als Übertrag
           rts
    1:       subql #1,$a0@+
             dbcc $d0,1b       | kein Carry -> Schleife abbrechen
           subxl $d0,$d0       | kein Carry -> d0.L=0, sonst d0.L=-1 als Übertrag
           rts

| extern uintD neg_loop_up (uintD* ptr, uintC count);
           DECLARE_FUNCTION(neg_loop_up)
C(neg_loop_up:) | Input in a0,d0.W, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       negxl $a0@+
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           rts

| extern uintD shift1left_loop_up (uintD* ptr, uintC count);
           DECLARE_FUNCTION(shift1left_loop_up)
C(shift1left_loop_up:) | Input in a0,d0.W, verändert d1, Output in d0.L
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           clrx   | X-Bit löschen
           bras 2f
    1:       movel $a0@,$d1
             roxll #1,$d1
             movel $d1,$a0@+
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           rts

| extern uintD shiftleft_loop_up (uintD* ptr, uintC count, uintC i, uintD carry);
           DECLARE_FUNCTION(shiftleft_loop_up)
C(shiftleft_loop_up:) | Input in a0,d0.W,d1.W,d2, Output in d0
#if 1
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           movel $sp@(16+16),$d2
           moveq #32,$d5
           subw $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.W = 32-i,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       movel $a0@,$d3  | d3.L = neues Digit
             movel $d3,$d4
             lsll $d1,$d4      | um i Bits nach links schieben
             orl $d2,$d4       | mit vorigem Übertrag kombinieren
             movel $d4,$a0@+   | 32 Bits ablegen
             movel $d3,$d2
             lsrl $d5,$d2      | neuen Übertrag bilden
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           movel $sp@(16+16),$d2
           moveq #1,$d5
           lsll $d1,$d5
           subql #1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.L = 2^i-1
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       movel $a0@,$d3  | d3.L = neues Digit
             roll $d1,$d3      | um i Bits links rotieren
             movel $d3,$d4
             andl $d5,$d3      | untere i Bits in d3
             eorl $d3,$d4      | obere 32-i Bits in d4
             orl $d2,$d4       | mit vorigem übertrag kombinieren
             movel $d4,$a0@+   | 32 Bits ablegen
             movel $d3,$d2     | neuer Übertrag
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#endif

#endif

| extern uintD shiftleftcopy_loop_up (uintD* sourceptr, uintD* destptr, uintC count, uintC i);
           DECLARE_FUNCTION(shiftleftcopy_loop_up)
C(shiftleftcopy_loop_up:) | Input in a0,a1,d0.W,d1.W, Output in d0
#if 1
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movel $sp@(16+8),$a1
           movew $sp@(16+12+2),$d0
           movew $sp@(16+16+2),$d1
           moveq #32,$d5
           subw $d1,$d5
           clrl $d2
           | a0 = sourceptr, a1 = destptr, d0.W = count, d1.W = i, d5.W = 32-i,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       movel $a0@+,$d3  | d3.L = neues Digit
             movel $d3,$d4
             lsll $d1,$d4      | um i Bits nach links schieben
             orl $d2,$d4       | mit vorigem Übertrag kombinieren
             movel $d4,$a1@+  | 32 Bits ablegen
             movel $d3,$d2
             lsrl $d5,$d2      | neuen Übertrag bilden
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movel $sp@(16+8),$a1
           movew $sp@(16+12+2),$d0
           movew $sp@(16+16+2),$d1
           moveq #1,$d5
           lsll $d1,$d5
           subql #1,$d5
           | a0 = sourceptr, a1 = destptr, d0.W = count, d1.W = i, d5.L = 2^i-1
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       movel $a0@+,$d3  | d3.L = neues Digit
             roll $d1,$d3      | um i Bits links rotieren
             movel $d3,$d4
             andl $d5,$d3      | untere i Bits in d3
             eorl $d3,$d4      | obere 32-i Bits in d4
             orl $d2,$d4       | mit vorigem übertrag kombinieren
             movel $d4,$a1@+  | 32 Bits ablegen
             movel $d3,$d2     | neuer Übertrag
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#endif

#if !CL_DS_BIG_ENDIAN_P

| extern uintD shift1right_loop_down (uintD* ptr, uintC count, uintD carry);
           DECLARE_FUNCTION(shift1right_loop_down)
C(shift1right_loop_down:) | Input in a0,d0.W,d1, Output in d0
           movel $sp@(4),$a0
           movew $sp@(8+2),$d0
           movel $sp@(12),$d1
           roxrl #1,$d1       | X-Bit löschen oder setzen, je nach d1
           bras 2f
    1:       movel $a0@-,$d1
             roxrl #1,$d1
             movel $d1,$a0@
    2:       dbra $d0,1b
           subxl $d0,$d0       | -1 falls X gesetzt, 0 falls X gelöscht
           rts

| extern uintD shiftright_loop_down (uintD* ptr, uintC count, uintC i);
           DECLARE_FUNCTION(shiftright_loop_down)
C(shiftright_loop_down:) | Input in a0,d0.W,d1.W, Output in d0
#if 1
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           moveq #32,$d5
           subw $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.W = 32-i,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           clrl $d2
           bras 2f
    1:       | a0 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.W = 32-i,
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@-,$d3   | neue Daten
             movel $d3,$d4
             lsrl $d1,$d3      | um i Bits rechts schieben
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
             movel $d3,$a0@  | ablegen
             lsll $d5,$d4      | um (32-i) Bits links geschoben
             movel $d4,$d2     | liefert neuen Übertrag
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           moveq #-1,$d5
           lsrl $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.L = 2^(32-i)-1,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           clrl $d2
           bras 2f
    1:       | a0 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.L = 2^(32-i)-1,
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@-,$d3   | neue Daten
             rorl $d1,$d3      | um i Bits rechts rotieren
             movel $d3,$d4
             andl $d5,$d3      | untere 32-i Bits
             eorl $d3,$d4      | obere i Bits
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
             movel $d4,$d2     | neuer Übertrag
             movel $d3,$a0@  | ablegen
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#endif

| extern uintD shiftrightsigned_loop_down (uintD* ptr, uintC count, uintC i);
           DECLARE_FUNCTION(shiftrightsigned_loop_down)
C(shiftrightsigned_loop_down:) | Input in a0,d0.W,d1.W, Output in d0
#if 1
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           moveq #32,$d5
           subw $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.W = 32-i,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           subqw #1,$d0
           movel $a0@-,$d3     | erstes Digit
           movel $d3,$d4
           asrl $d1,$d3        | um i Bits rechts schieben
           bras 2f
    1:       | a0 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.W = 32-i,
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@-,$d3   | neue Daten
             movel $d3,$d4
             lsrl $d1,$d3      | um i Bits rechts schieben
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
    2:       movel $d3,$a0@  | ablegen
             lsll $d5,$d4      | um (32-i) Bits links geschoben
             movel $d4,$d2     | liefert neuen Übertrag
             dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movew $sp@(16+8+2),$d0
           movew $sp@(16+12+2),$d1
           moveq #-1,$d5
           lsrl $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.L = 2^(32-i)-1,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           subqw #1,$d0
           movel $a0@-,$d3     | erstes Digit
           movel $d3,$d4
           rorl $d1,$d4        | um i Bits rechts rotieren
           movel $d5,$d2
           notl $d2
           andl $d4,$d2        | obere 32-i Bits
           asrl $d1,$d3        | erstes Digit um i Bits rechts shiften
           bras 2f
    1:       | a0 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.L = 2^(32-i)-1,
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@-,$d3   | neue Daten
             rorl $d1,$d3      | um i Bits rechts rotieren
             movel $d3,$d4
             andl $d5,$d3      | untere 32-i Bits
             eorl $d3,$d4      | obere i Bits
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
             movel $d4,$d2     | neuer Übertrag
    2:       movel $d3,$a0@  | ablegen
             dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#endif

| extern uintD shiftrightcopy_loop_down (uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry);
           DECLARE_FUNCTION(shiftrightcopy_loop_down)
C(shiftrightcopy_loop_down:) | Input in a0,a1,d0.W,d1.W,d2, Output in d0
#if 1
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movel $sp@(16+8),$a1
           movew $sp@(16+12+2),$d0
           movew $sp@(16+16+2),$d1
           movel $sp@(16+20),$d2
           moveq #32,$d5
           subw $d1,$d5
           | a0 = ptr, d0.W = count, d1.W = i, d5.W = 32-i
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       | a0,a1 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.W = 32-i
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@-,$d3  | neue Daten
             movel $d3,$d4
             lsrl $d1,$d3      | um i Bits rechts schieben
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
             movel $d3,$a1@-  | ablegen
             movel $d4,$d2
    2:       lsll $d5,$d2      | um (32-i) Bits links geschoben, gibt neuen Übertrag
             dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movel $sp@(16+8),$a1
           movew $sp@(16+12+2),$d0
           movew $sp@(16+16+2),$d1
           movel $sp@(16+20),$d2
           moveq #-1,$d5
           lsrl $d1,$d5
           rorl $d1,$d2
           | a0 = ptr, d0.W = count, d1.W = i, d5.L = 2^(32-i)-1
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       | a0,a1 = Aufwärtszähler Adresse, d0.W = Herabzähler, d1.W = i, d5.L = 2^(32-i)-1
             | d2.L = Schiebe-Übertrag (obere i Bits, restliche 32-i Bits sind 0)
             | d3.L = Schiebe-Akku
             movel $a0@-,$d3  | neue Daten
             rorl $d1,$d3      | um i Bits rechts rotieren
             movel $d3,$d4
             andl $d5,$d3      | untere 32-i Bits
             eorl $d3,$d4      | obere i Bits
             orl $d2,$d3       | und mit vorigem Übertrag kombinieren
             movel $d4,$d2     | neuer Übertrag
             movel $d3,$a1@-  | ablegen
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           movel $d2,$d0
           moveml $sp@+,$d2-$d5
           rts
#endif

| extern uintD mulusmall_loop_up (uintD digit, uintD* ptr, uintC len, uintD newdigit);
           DECLARE_FUNCTION(mulusmall_loop_up)
C(mulusmall_loop_up:) | Input in d0,a0,d1.W,d2, Output in d0
           moveml $d2-$d4,$sp@-
           movel $sp@(12+4),$d0
           movel $sp@(12+8),$a0
           movew $sp@(12+12+2),$d1
           movel $sp@(12+16),$d2
           addw #0,$d1        | X-Bit löschen
           bras 2f
    1:       movel $a0@,$d3  | nächstes Digit
             mulul $d0,$d4:$d3  | mit digit multiplizieren
             addxl $d2,$d3     | und bisherigen Carry und X-Bit addieren
             movel $d3,$a0@+   | Low-Digit ablegen
             movel $d4,$d2     | High-Digit gibt neuen Carry
    2:       dbra $d1,1b
           clrl $d0
           addxl $d2,$d0       | letzter Carry (incl. X-Bit)
           moveml $sp@+,$d2-$d4
           rts

| extern void mulu_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
           DECLARE_FUNCTION(mulu_loop_up)
C(mulu_loop_up:) | Input in d0,a0,a1,d1.W
#if 1
           moveml $d2-$d4,$sp@-
           movel $sp@(12+4),$d0
           movel $sp@(12+8),$a0
           movel $sp@(12+12),$a1
           movew $sp@(12+16+2),$d1
           subl $d2,$d2        | carry := 0, X-Bit löschen
           bras 2f
    1:       movel $a0@+,$d3  | nächstes Digit
             mulul $d0,$d4:$d3  | mit digit multiplizieren
             addxl $d2,$d3     | und bisherigen Carry und X-Bit addieren
             movel $d3,$a1@+  | Low-Digit ablegen
             movel $d4,$d2     | High-Digit gibt neuen Carry
    2:       dbra $d1,1b
           clrl $d3
           addxl $d3,$d2       | letztes X-Bit verarbeiten
           movel $d2,$a1@+    | letzten Carry ablegen
           moveml $sp@+,$d2-$d4
           rts
#else
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$d0
           movel $sp@(16+8),$a0
           movel $sp@(16+12),$a1
           movew $sp@(16+16+2),$d1
           clrl $d5           | 0
           clrl $d2           | carry
           bras 2f
    1:       movel $a0@+,$d3  | nächstes Digit
             mulul $d0,$d4:$d3  | mit digit multiplizieren
             addl $d2,$d3      | und bisherigen Carry addieren
             addxl $d5,$d4
             movel $d3,$a1@+  | Low-Digit ablegen
             movel $d4,$d2     | High-Digit gibt neuen Carry
    2:       dbra $d1,1b
           movel $d2,$a1@+    | letzten Carry ablegen
           moveml $sp@+,$d2-$d5
           rts
#endif

| extern uintD muluadd_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
           DECLARE_FUNCTION(muluadd_loop_up)
C(muluadd_loop_up:) | Input in d0,a0,a1,d1.W, Output in d0
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$d0
           movel $sp@(16+8),$a0
           movel $sp@(16+12),$a1
           movew $sp@(16+16+2),$d1
           clrl $d5           | 0
           subl $d2,$d2        | carry := 0, X-Bit löschen
           bras 2f
    1:       movel $a0@+,$d3  | nächstes Digit
             mulul $d0,$d4:$d3  | mit digit multiplizieren
             addxl $d2,$d3     | und bisherigen Carry und X-Bit addieren
             addxl $d5,$d4
             addl $d3,$a1@+   | Low-Digit zum dest-Digit addieren, X als Übertrag
             movel $d4,$d2     | High-Digit gibt neuen Carry
    2:       dbra $d1,1b
           addxl $d5,$d2       | letztes X-Bit addieren
           movel $d2,$d0       | letzten Carry als Ergebnis
           moveml $sp@+,$d2-$d5
           rts

| extern uintD mulusub_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
           DECLARE_FUNCTION(mulusub_loop_up)
C(mulusub_loop_up:) | Input in d0,a0,a1,d1.W, Output in d0
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$d0
           movel $sp@(16+8),$a0
           movel $sp@(16+12),$a1
           movew $sp@(16+16+2),$d1
           clrl $d5           | 0
           subl $d2,$d2        | carry := 0, X-Bit löschen
           bras 2f
    1:       movel $a0@+,$d3  | nächstes Digit
             mulul $d0,$d4:$d3  | mit digit multiplizieren
             addxl $d2,$d3     | und bisherigen Carry und X-Bit addieren
             addxl $d5,$d4
             subl $d3,$a1@+   | Low-Digit vom dest-Digit subtrahieren, X als Übertrag
             movel $d4,$d2     | High-Digit gibt neuen Carry
    2:       dbra $d1,1b
           clrl $d0
           addxl $d2,$d0       | letzter Carry und letztes X-Bit
           moveml $sp@+,$d2-$d5
           rts

| extern uintD divu_loop_down (uintD digit, uintD* ptr, uintC len);
           DECLARE_FUNCTION(divu_loop_down)
C(divu_loop_down:) | Input in d0,a0,d1.W, Output in d0
           moveml $d2-$d3,$sp@-
           movel $sp@(8+4),$d0
           movel $sp@(8+8),$a0
           movew $sp@(8+12+2),$d1
           clrl $d2           | Rest := 0
           bras 2f
    1:       movel $a0@-,$d3   | nächst-niedriges Digit
             divul $d0,$d2:$d3  | mit Rest kombinieren und durch digit dividieren
             movel $d3,$a0@  | Quotient ablegen, Rest in d2
    2:       dbra $d1,1b
           movel $d2,$d0       | Rest
           moveml $sp@+,$d2-$d3
           rts

| extern uintD divucopy_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
           DECLARE_FUNCTION(divucopy_loop_down)
C(divucopy_loop_down:) | Input in d0,a0,a1,d1.W, Output in d0
           moveml $d2-$d3,$sp@-
           movel $sp@(8+4),$d0
           movel $sp@(8+8),$a0
           movel $sp@(8+12),$a1
           movew $sp@(8+16+2),$d1
           clrl $d2           | Rest := 0
           bras 2f
    1:       movel $a0@-,$d3  | nächst-niedriges Digit
             divul $d0,$d2:$d3  | mit Rest kombinieren und durch digit dividieren
             movel $d3,$a1@-  | Quotient ablegen, Rest in d2
    2:       dbra $d1,1b
           movel $d2,$d0       | Rest
           moveml $sp@+,$d2-$d3
           rts

#endif

| extern void shiftxor_loop_up (uintD* xptr, const uintD* yptr, uintC count, uintC i);
           DECLARE_FUNCTION(shiftxor_loop_up)
C(shiftxor_loop_up:) | Input in a0,a1,d0.W,d1.W
           moveml $d2-$d5,$sp@-
           movel $sp@(16+4),$a0
           movel $sp@(16+8),$a1
           movew $sp@(16+12+2),$d0
           movew $sp@(16+16+2),$d1
           moveq #32,$d5
           subw $d1,$d5
           clrl $d2
           | a0 = xptr, a1 = yptr, d0.W = count, d1.W = i, d5.W = 32-i,
           | d2.L = Schiebe-Übertrag (i Bits), d3.L = Schiebe-Akku
           bras 2f
    1:       movel $a1@+,$d3  | d3.L = neues Digit
             movel $d3,$d4
             lsll $d1,$d4      | um i Bits nach links schieben
             orl $d2,$d4       | mit vorigem Übertrag kombinieren
             eorl $d4,$a0@+  | 32 Bits ablegen
             movel $d3,$d2
             lsrl $d5,$d2      | neuen Übertrag bilden
    2:       dbra $d0,1b        | Schleife d0.W mal durchlaufen
           eorl $d2,$a0@
           moveml $sp@+,$d2-$d5
           rts
