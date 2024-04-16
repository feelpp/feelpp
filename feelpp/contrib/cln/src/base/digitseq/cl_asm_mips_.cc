// Externe Routinen zu ARILEV1.D
// Prozessor: MIPS
// Endianness: irrelevant
// Compiler: GNU-C oder ...
// Parameter-Übergabe:
//   o32: in Registern $4,$5,$6,$7, und auf dem Stack 16($sp),...
//   n32: in Registern $4,$5,$6,$7,$8,$9,$10,$11, und auf dem Stack 4($sp),...
// Rückgabewert: in Register $2
// Einstellungen: intCsize=32, intDsize=32.
// Besonderheiten: Nach jedem Ladebefehl ein Wartetakt nötig, bevor der
//   geholte Wert benutzt werden darf.

// Strictly speaking, the MIPS ABI (-32 or -n32) is independent from the CPU
// identification (-mips[12] or -mips[34]). But -n32 is commonly used together
// with -mips3, and it's easier to test the CPU identification.
#if __mips >= 3
  #define ABI_N32 1
#else
  #define ABI_O32 1
#endif

// When this file is compiled into a shared library, ELF linkers need to
// know which symbols are functions.
#if defined(__GNU__) || defined(__NetBSD__)
  #define DECLARE_FUNCTION(name) .type name,@function
#else
  #define DECLARE_FUNCTION(name)
#endif

        .text

        .globl copy_loop_up
        .globl copy_loop_down
        .globl fill_loop_up
        .globl fill_loop_down
        .globl clear_loop_up
        .globl clear_loop_down
        .globl test_loop_up
        .globl test_loop_down
        .globl xor_loop_up
        .globl compare_loop_up
#if CL_DS_BIG_ENDIAN_P
        .globl or_loop_up
        .globl and_loop_up
        .globl eqv_loop_up
        .globl nand_loop_up
        .globl nor_loop_up
        .globl andc2_loop_up
        .globl orc2_loop_up
        .globl not_loop_up
        .globl and_test_loop_up
        .globl add_loop_down
        .globl addto_loop_down
        .globl inc_loop_down
        .globl sub_loop_down
        .globl subx_loop_down
        .globl subfrom_loop_down
        .globl dec_loop_down
        .globl neg_loop_down
#else
        .globl or_loop_down
        .globl xor_loop_down
        .globl and_loop_down
        .globl eqv_loop_down
        .globl nand_loop_down
        .globl nor_loop_down
        .globl andc2_loop_down
        .globl orc2_loop_down
        .globl not_loop_down
        .globl and_test_loop_down
        .globl compare_loop_down
        .globl add_loop_up
        .globl addto_loop_up
        .globl inc_loop_up
        .globl sub_loop_up
        .globl subx_loop_up
        .globl subfrom_loop_up
        .globl dec_loop_up
        .globl neg_loop_up
#endif

#ifndef __GNUC__ /* mit GNU-C machen wir mulu32() als Macro, der inline multipliziert */

// extern struct { uint32 lo; uint32 hi; } mulu32_ (uint32 arg1, uint32 arg2);
// 2^32*hi+lo := arg1*arg2.
        .globl mulu32_
        .align 2
        DECLARE_FUNCTION(mulu32_)
        .ent mulu32_ // Input in $4,$5, Output in $2,mulu32_high
mulu32_:
#if __mips_isa_rev >= 6
        mulu $2,$5,$4           // arg1 * arg2, lo
        muhu $6,$5,$4           // arg1 * arg2, hi
#else
        multu $5,$4             // arg1 * arg2
        mfhi $6                 // hi
        mflo $2                 // lo
#endif
        sw $6,mulu32_high       // hi abspeichern // Adressierung?? Deklaration??
        j $31                   // return
        .end mulu32_

#endif

// extern uintD* copy_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(copy_loop_up)
        .ent copy_loop_up // Input in $4,$5,$6, Output in $2
colu1:    lw $12,($4)           // d = *sourceptr
          addu $4,4             // sourceptr++
          sw $12,($5)           // *destptr = d
          addu $5,4             // destptr++
          subu $6,1             // count--
copy_loop_up:
          bnez $6,colu1         // until (count==0)
        move $2,$5              // destptr
        j $31                   // return
        .end copy_loop_up

// extern uintD* copy_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(copy_loop_down)
        .ent copy_loop_down // Input in $4,$5,$6, Output in $2
cold1:    subu $4,4             // sourceptr--
          lw $12,($4)           // d = *sourceptr
          subu $5,4             // destptr--
          sw $12,($5)           // *destptr = d
          subu $6,1             // count--
copy_loop_down:
          bnez $6,cold1         // until (count==0)
        move $2,$5              // destptr
        j $31                   // return
        .end copy_loop_down

// extern uintD* fill_loop_up (uintD* destptr, uintC count, uintD filler);
        .align 2
        DECLARE_FUNCTION(fill_loop_up)
        .ent fill_loop_up // Input in $4,$5,$6, Output in $2
flu1:     sw $6,($4)            // *destptr = filler
          addu $4,4             // destptr++
          subu $5,1             // count--
fill_loop_up:
          bnez $5,flu1          // until (count==0)
        move $2,$4              // destptr
        j $31                   // return
        .end fill_loop_up

// extern uintD* fill_loop_down (uintD* destptr, uintC count, uintD filler);
        .align 2
        DECLARE_FUNCTION(fill_loop_down)
        .ent fill_loop_down // Input in $4,$5,$6, Output in $2
fld1:     subu $4,4             // destptr--
          sw $6,($4)            // *destptr = filler
          subu $5,1             // count--
fill_loop_down:
          bnez $5,fld1          // until (count==0)
        move $2,$4              // destptr
        j $31                   // return
        .end fill_loop_down

// extern uintD* clear_loop_up (uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(clear_loop_up)
        .ent clear_loop_up // Input in $4,$5, Output in $2
cllu1:    sw $0,($4)            // *destptr = 0
          addu $4,4             // destptr++
          subu $5,1             // count--
clear_loop_up:
          bnez $5,cllu1         // until (count==0)
        move $2,$4              // destptr
        j $31                   // return
        .end clear_loop_up

// extern uintD* clear_loop_down (uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(clear_loop_down)
        .ent clear_loop_down // Input in $4,$5, Output in $2
clld1:    subu $4,4             // destptr--
          sw $0,($4)            // *destptr = 0
          subu $5,1             // count--
clear_loop_down:
          bnez $5,clld1         // until (count==0)
        move $2,$4              // destptr
        j $31                   // return
        .end clear_loop_down

// extern boolean test_loop_up (uintD* ptr, uintC count);
        .align 2
        DECLARE_FUNCTION(test_loop_up)
        .ent test_loop_up // Input in $4,$5
tlu1:     lw $12,($4)           // x = *ptr
          addu $4,4             // ptr++
          bnez $12,tlu3
          subu $5,1             // count--
test_loop_up:
          bnez $5,tlu1          // until (count==0)
        move $2,$0              // 0
        j $31                   // return
tlu3:   li $2,1                 // 1
        j $31                   // return
        .end test_loop_up

// extern boolean test_loop_down (uintD* ptr, uintC count);
        .align 2
        DECLARE_FUNCTION(test_loop_down)
        .ent test_loop_down // Input in $4,$5
tld1:     subu $4,4             // ptr--
          lw $12,($4)           // x = *ptr
          subu $5,1             // count--
          bnez $12,tld3
test_loop_down:
          bnez $5,tld1          // until (count==0)
        move $2,$0              // 0
        j $31                   // return
tld3:   li $2,1                 // 1
        j $31                   // return
        .end test_loop_down

#if CL_DS_BIG_ENDIAN_P

// extern void or_loop_up (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(or_loop_up)
        .ent or_loop_up // Input in $4,$5,$6
olu1:     lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          addu $5,4             // yptr++
          or $12,$13            // x |= y
          sw $12,($4)           // *xptr = x
          addu $4,4             // xptr++
          subu $6,1             // count--
or_loop_up:
          bnez $6,olu1          // until (count==0)
        j $31                   // return
        .end or_loop_up

#endif

// extern void xor_loop_up (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(xor_loop_up)
        .ent xor_loop_up // Input in $4,$5,$6
xlu1:     lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          addu $5,4             // yptr++
          xor $12,$13           // x ^= y
          sw $12,($4)           // *xptr = x
          addu $4,4             // xptr++
          subu $6,1             // count--
xor_loop_up:
          bnez $6,xlu1          // until (count==0)
        j $31                   // return
        .end xor_loop_up

#if CL_DS_BIG_ENDIAN_P

// extern void and_loop_up (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(and_loop_up)
        .ent and_loop_up // Input in $4,$5,$6
alu1:     lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          addu $5,4             // yptr++
          and $12,$13           // x &= y
          sw $12,($4)           // *xptr = x
          addu $4,4             // xptr++
          subu $6,1             // count--
and_loop_up:
          bnez $6,alu1          // until (count==0)
        j $31                   // return
        .end and_loop_up

// extern void eqv_loop_up (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(eqv_loop_up)
        .ent eqv_loop_up // Input in $4,$5,$6
nxlu1:    lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          addu $5,4             // yptr++
          xor $12,$13           // x ^= y
          nor $12,$0            // x = ~x
          sw $12,($4)           // *xptr = x
          addu $4,4             // xptr++
          subu $6,1             // count--
eqv_loop_up:
          bnez $6,nxlu1         // until (count==0)
        j $31                   // return
        .end eqv_loop_up

// extern void nand_loop_up (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(nand_loop_up)
        .ent nand_loop_up // Input in $4,$5,$6
nalu1:    lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          addu $5,4             // yptr++
          and $12,$13           // x &= y        // Gibt es 'nand $12,$13' ??
          nor $12,$0            // x = ~x
          sw $12,($4)           // *xptr = x
          addu $4,4             // xptr++
          subu $6,1             // count--
nand_loop_up:
          bnez $6,nalu1         // until (count==0)
        j $31                   // return
        .end nand_loop_up

// extern void nor_loop_up (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(nor_loop_up)
        .ent nor_loop_up // Input in $4,$5,$6
nolu1:    lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          addu $5,4             // yptr++
          nor $12,$13           // x = ~(x|y)
          sw $12,($4)           // *xptr = x
          addu $4,4             // xptr++
          subu $6,1             // count--
nor_loop_up:
          bnez $6,nolu1         // until (count==0)
        j $31                   // return
        .end nor_loop_up

// extern void andc2_loop_up (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(andc2_loop_up)
        .ent andc2_loop_up // Input in $4,$5,$6
aclu1:    lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          addu $5,4             // yptr++
          nor $13,$0            // y = ~y
          and $12,$13           // x &= y
          sw $12,($4)           // *xptr = x
          addu $4,4             // xptr++
          subu $6,1             // count--
andc2_loop_up:
          bnez $6,aclu1         // until (count==0)
        j $31                   // return
        .end andc2_loop_up

// extern void orc2_loop_up (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(orc2_loop_up)
        .ent orc2_loop_up // Input in $4,$5,$6
oclu1:    lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          addu $5,4             // yptr++
          nor $13,$0            // y = ~y
          or $12,$13            // x |= y
          sw $12,($4)           // *xptr = x
          addu $4,4             // xptr++
          subu $6,1             // count--
orc2_loop_up:
          bnez $6,oclu1         // until (count==0)
        j $31                   // return
        .end orc2_loop_up

// extern void not_loop_up (uintD* xptr, uintC count);
        .align 2
        DECLARE_FUNCTION(not_loop_up)
        .ent not_loop_up // Input in $4,$5
nlu1:     lw $12,($4)           // x = *xptr
          subu $5,1             // count--
          nor $12,$0            // x = ~x
          sw $12,($4)           // *xptr = x
          addu $4,4             // xptr++
not_loop_up:
          bnez $5,nlu1          // until (count==0)
        j $31                   // return
        .end not_loop_up

// extern boolean and_test_loop_up (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(and_test_loop_up)
        .ent and_test_loop_up // Input in $4,$5,$6
atlu1:    lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          addu $5,4             // yptr++
          and $12,$13           // x &= y
          bnez $12,atlu3        // if (x) ...
          addu $4,4             // xptr++
          subu $6,1             // count--
and_test_loop_up:
          bnez $6,atlu1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
atlu3:  li $2,1                 // 1
        j $31                   // return
        .end and_test_loop_up

#endif

// extern cl_signean compare_loop_up (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(compare_loop_up)
        .ent compare_loop_up // Input in $4,$5,$6
cmlu1:    lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          addu $5,4             // yptr++
          bne $12,$13,cmlu3     // if (!(x==y)) ...
          addu $4,4             // xptr++
          subu $6,1             // count--
compare_loop_up:
          bnez $6,cmlu1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
cmlu3:  bltu $12,$13,cmlu4      // if (x<y) ...
        li $2,1                 // 1
        j $31                   // return
cmlu4:  li $2,-1                // -1
        j $31                   // return
        .end compare_loop_up

#if CL_DS_BIG_ENDIAN_P

// extern uintD add_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(add_loop_down)
        .ent add_loop_down // Input in $4,$5,$6,$7, Output in $2
ald1:     // kein Carry
          subu $4,4             // sourceptr1--
          subu $5,4             // sourceptr2--
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          subu $6,4             // destptr--
          addu $12,$13          // dest = source1 + source2
          sw $12,($6)           // *destptr = dest
          bltu $12,$13,ald4     // if (dest < source2) [also Carry] ...
ald2:
          subu $7,1             // count--
add_loop_down:
          bnez $7,ald1          // until (count==0)
        move $2,$0              // 0
        j $31                   // return
ald3:   // Hier Carry
          subu $4,4             // sourceptr1--
          subu $5,4             // sourceptr2--
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          subu $6,4             // destptr--
          addu $12,$13          // dest = source1 + source2
          addu $12,1            //        + 1
          sw $12,($6)           // *destptr = dest
          bgtu $12,$13,ald2     // if (dest > source2) [also kein Carry] ...
ald4:     subu $7,1             // count--
          bnez $7,ald3          // until (count==0)
        li $2,1                 // 1
        j $31                   // return
        .end add_loop_down

// extern uintD addto_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(addto_loop_down)
        .ent addto_loop_down // Input in $4,$5,$6, Output in $2
atld1:    // kein Carry
          subu $4,4             // sourceptr--
          subu $5,4             // destptr--
          lw $12,($4)           // source1 = *sourceptr
          lw $13,($5)           // source2 = *destptr
          subu $6,1             // count--
          addu $12,$13          // dest = source1 + source2
          sw $12,($5)           // *destptr = dest
          bltu $12,$13,atld4    // if (dest < source2) [also Carry] ...
addto_loop_down:
atld2:    bnez $6,atld1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
atld3:  // Hier Carry
          subu $4,4             // sourceptr--
          subu $5,4             // destptr--
          lw $12,($4)           // source1 = *sourceptr
          lw $13,($5)           // source2 = *destptr
          subu $6,1             // count--
          addu $12,$13          // dest = source1 + source2
          addu $12,1            //        + 1
          sw $12,($5)           // *destptr = dest
          bgtu $12,$13,atld2    // if (dest > source2) [also kein Carry] ...
atld4:    bnez $6,atld3         // until (count==0)
        li $2,1                 // 1
        j $31                   // return
        .end addto_loop_down

// extern uintD inc_loop_down (uintD* ptr, uintC count);
        .align 2
        DECLARE_FUNCTION(inc_loop_down)
        .ent inc_loop_down // Input in $4,$5, Output in $2
ild1:     subu $4,4             // ptr--
          lw $12,($4)           // x = *ptr
          subu $5,1             // count--
          addu $12,1            // x++;
          sw $12,($4)           // *ptr = x
          bnez $12,ild3         // if (!(x==0)) ...
inc_loop_down:
          bnez $5,ild1          // until (count==0)
        li $2,1                 // 1
        j $31                   // return
ild3:   move $2,$0              // 0
        j $31                   // return
        .end inc_loop_down

// extern uintD sub_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(sub_loop_down)
        .ent sub_loop_down // Input in $4,$5,$6,$7, Output in $2
sld1:     // kein Carry
          subu $4,4             // sourceptr1--
          subu $5,4             // sourceptr2--
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          subu $6,4             // destptr--
          bltu $12,$13,sld2     // if (source1 < source2) [also Carry] ...
          subu $12,$13          // dest = source1 - source2
          sw $12,($6)           // *destptr = dest
          subu $7,1             // count--
sub_loop_down:
          bnez $7,sld1          // until (count==0)
        move $2,$0              // 0
        j $31                   // return
sld2:     subu $12,$13          // dest = source1 - source2
          sw $12,($6)           // *destptr = dest
          subu $7,1             // count--
          bnez $7,sld3          // until (count==0)
        li $2,-1                // -1
        j $31                   // return
sld3:   // Hier Carry
          subu $4,4             // sourceptr1--
          subu $5,4             // sourceptr2--
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          subu $6,4             // destptr--
          bgtu $12,$13,sld4     // if (source1 > source2) [also kein Carry] ...
          subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($6)           // *destptr = dest
          subu $7,1             // count--
          bnez $7,sld3          // until (count==0)
        li $2,-1                // -1
        j $31                   // return
sld4:     subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($6)           // *destptr = dest
          subu $7,1             // count--
          bnez $7,sld1          // until (count==0)
        move $2,$0              // 0
        j $31                   // return
        .end sub_loop_down

// extern uintD subx_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);
        .align 2
        DECLARE_FUNCTION(subx_loop_down)
        .ent subx_loop_down // Input in $4,$5,$6,$7,$8 Output in $2
subx_loop_down:
#if ABI_N32
        move $12,$8             // carry
#else
        lw $12,16($sp)          // carry
#endif
        bnez $12,sxld5          // !(carry==0) ?
        b sxld2
sxld1:    // kein Carry
          subu $4,4             // sourceptr1--
          subu $5,4             // sourceptr2--
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          subu $6,4             // destptr--
          bltu $12,$13,sxld3    // if (source1 < source2) [also Carry] ...
          subu $12,$13          // dest = source1 - source2
          sw $12,($6)           // *destptr = dest
          subu $7,1             // count--
sxld2:    bnez $7,sxld1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
sxld3:    subu $12,$13          // dest = source1 - source2
          sw $12,($6)           // *destptr = dest
          subu $7,1             // count--
          bnez $7,sxld4         // until (count==0)
        li $2,-1                // -1
        j $31                   // return
sxld4:  // Hier Carry
          subu $4,4             // sourceptr1--
          subu $5,4             // sourceptr2--
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          subu $6,4             // destptr--
          bgtu $12,$13,sxld6    // if (source1 > source2) [also kein Carry] ...
          subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($6)           // *destptr = dest
          subu $7,1             // count--
sxld5:    bnez $7,sxld4         // until (count==0)
        li $2,-1                // -1
        j $31                   // return
sxld6:    subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($6)           // *destptr = dest
          subu $7,1             // count--
          bnez $7,sxld1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
        .end subx_loop_down

// extern uintD subfrom_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(subfrom_loop_down)
        .ent subfrom_loop_down // Input in $4,$5,$6,$7, Output in $2
sfld1:    // kein Carry
          subu $4,4             // sourceptr--
          subu $5,4             // destptr--
          lw $12,($5)           // source1 = *destptr
          lw $13,($4)           // source2 = *sourceptr
          subu $6,1             // count--
          bltu $12,$13,sfld2    // if (source1 < source2) [also Carry] ...
          subu $12,$13          // dest = source1 - source2
          sw $12,($5)           // *destptr = dest
subfrom_loop_down:
          bnez $6,sfld1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
sfld2:    subu $12,$13          // dest = source1 - source2
          sw $12,($5)           // *destptr = dest
          bnez $6,sfld3         // until (count==0)
        li $2,-1                // -1
        j $31                   // return
sfld3:  // Hier Carry
          subu $4,4             // sourceptr--
          subu $5,4             // destptr--
          lw $12,($5)           // source1 = *destptr
          lw $13,($4)           // source2 = *sourceptr
          subu $6,1             // count--
          bgtu $12,$13,sfld4    // if (source1 > source2) [also kein Carry] ...
          subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($5)           // *destptr = dest
          bnez $6,sfld3         // until (count==0)
        li $2,-1                // -1
        j $31                   // return
sfld4:    subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($5)           // *destptr = dest
          bnez $6,sfld1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
        .end subfrom_loop_down

// extern uintD dec_loop_down (uintD* ptr, uintC count);
        .align 2
        DECLARE_FUNCTION(dec_loop_down)
        .ent dec_loop_down // Input in $4,$5, Output in $2
dld1:     subu $4,4             // ptr--
          lw $12,($4)           // x = *ptr
          subu $5,1             // count--
          bnez $12,dld3         // if (!(x==0)) ...
          subu $12,1            // x--;
          sw $12,($4)           // *ptr = x
dec_loop_down:
          bnez $5,dld1          // until (count==0)
        li $2,-1                // -1
        j $31                   // return
dld3:   subu $12,1              // x--;
        sw $12,($4)             // *ptr = x
        move $2,$0              // 0
        j $31                   // return
        .end dec_loop_down

// extern uintD neg_loop_down (uintD* ptr, uintC count);
        .align 2
        DECLARE_FUNCTION(neg_loop_down)
        .ent neg_loop_down // Input in $4,$5, Output in $2
        // erstes Digit /=0 suchen:
nld1:     subu $4,4             // ptr--
          lw $12,($4)           // x = *ptr
          subu $5,1             // count--
          bnez $12,nld3         // if (!(x==0)) ...
neg_loop_down:
          bnez $5,nld1          // until (count==0)
        move $2,$0              // 0
        j $31                   // return
nld3:   // erstes Digit /=0 gefunden, ab jetzt gibt's Carrys
        // 1 Digit negieren:
        subu $12,$0,$12         // x = -x
        sw $12,($4)             // *ptr = x
        // alle anderen Digits invertieren:
        b nld5
nld4:     subu $4,4             // xptr--
          lw $12,($4)           // x = *xptr
          subu $5,1             // count--
          nor $12,$0            // x = ~x
          sw $12,($4)           // *xptr = x
nld5:     bnez $5,nld4          // until (count==0)
        li $2,-1                // -1
        j $31                   // return
        .end neg_loop_down

#endif

#if !CL_DS_BIG_ENDIAN_P

// extern void or_loop_down (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(or_loop_down)
        .ent or_loop_down // Input in $4,$5,$6
old1:     subu $4,4             // xptr--
          subu $5,4             // yptr--
          lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          subu $6,1             // count--
          or $12,$13            // x |= y
          sw $12,($4)           // *xptr = x
or_loop_down:
          bnez $6,old1          // until (count==0)
        j $31                   // return
        .end or_loop_down

// extern void xor_loop_down (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(xor_loop_down)
        .ent xor_loop_down // Input in $4,$5,$6
xld1:     subu $4,4             // xptr--
          subu $5,4             // yptr--
          lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          subu $6,1             // count--
          xor $12,$13           // x ^= y
          sw $12,($4)           // *xptr = x
xor_loop_down:
          bnez $6,xld1          // until (count==0)
        j $31                   // return
        .end xor_loop_down

// extern void and_loop_down (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(and_loop_down)
        .ent and_loop_down // Input in $4,$5,$6
ald1:     subu $4,4             // xptr--
          subu $5,4             // yptr--
          lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          subu $6,1             // count--
          and $12,$13           // x &= y
          sw $12,($4)           // *xptr = x
and_loop_down:
          bnez $6,ald1          // until (count==0)
        j $31                   // return
        .end and_loop_down

// extern void eqv_loop_down (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(eqv_loop_down)
        .ent eqv_loop_down // Input in $4,$5,$6
nxld1:    subu $4,4             // xptr--
          subu $5,4             // yptr--
          lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          subu $6,1             // count--
          xor $12,$13           // x ^= y
          nor $12,$0            // x = ~x
          sw $12,($4)           // *xptr = x
eqv_loop_down:
          bnez $6,nxld1         // until (count==0)
        j $31                   // return
        .end eqv_loop_down

// extern void nand_loop_down (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(nand_loop_down)
        .ent nand_loop_down // Input in $4,$5,$6
nald1:    subu $4,4             // xptr--
          subu $5,4             // yptr--
          lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          subu $6,1             // count--
          and $12,$13           // x &= y        // Gibt es 'nand $12,$13' ??
          nor $12,$0            // x = ~x
          sw $12,($4)           // *xptr = x
nand_loop_down:
          bnez $6,nald1         // until (count==0)
        j $31                   // return
        .end nand_loop_down

// extern void nor_loop_down (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(nor_loop_down)
        .ent nor_loop_down // Input in $4,$5,$6
nold1:    subu $4,4             // xptr--
          subu $5,4             // yptr--
          lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          subu $6,1             // count--
          nor $12,$13           // x = ~(x|y)
          sw $12,($4)           // *xptr = x
nor_loop_down:
          bnez $6,nold1         // until (count==0)
        j $31                   // return
        .end nor_loop_down

// extern void andc2_loop_down (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(andc2_loop_down)
        .ent andc2_loop_down // Input in $4,$5,$6
acld1:    subu $4,4             // xptr--
          subu $5,4             // yptr--
          lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          subu $6,1             // count--
          nor $13,$0            // y = ~y
          and $12,$13           // x &= y
          sw $12,($4)           // *xptr = x
andc2_loop_down:
          bnez $6,acld1         // until (count==0)
        j $31                   // return
        .end andc2_loop_down

// extern void orc2_loop_down (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(orc2_loop_down)
        .ent orc2_loop_down // Input in $4,$5,$6
ocld1:    subu $4,4             // xptr--
          subu $5,4             // yptr--
          lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          subu $6,1             // count--
          nor $13,$0            // y = ~y
          or $12,$13            // x |= y
          sw $12,($4)           // *xptr = x
orc2_loop_down:
          bnez $6,ocld1         // until (count==0)
        j $31                   // return
        .end orc2_loop_down

// extern void not_loop_down (uintD* xptr, uintC count);
        .align 2
        DECLARE_FUNCTION(not_loop_down)
        .ent not_loop_down // Input in $4,$5
nld1:     subu $4,4             // xptr--
          lw $12,($4)           // x = *xptr
          subu $5,1             // count--
          nor $12,$0            // x = ~x
          sw $12,($4)           // *xptr = x
not_loop_down:
          bnez $5,nld1          // until (count==0)
        j $31                   // return
        .end not_loop_down

// extern boolean and_test_loop_down (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(and_test_loop_down)
        .ent and_test_loop_down // Input in $4,$5,$6
atld1:    subu $4,4             // xptr--
          subu $5,4             // yptr--
          lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          and $12,$13           // x &= y
          bnez $12,atld3        // if (x) ...
          subu $6,1             // count--
and_test_loop_down:
          bnez $6,atld1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
atld3:  li $2,1                 // 1
        j $31                   // return
        .end and_test_loop_down

// extern cl_signean compare_loop_down (uintD* xptr, uintD* yptr, uintC count);
        .align 2
        DECLARE_FUNCTION(compare_loop_down)
        .ent compare_loop_down // Input in $4,$5,$6
cmld1:    subu $4,4             // xptr--
          subu $5,4             // yptr--
          lw $12,($4)           // x = *xptr
          lw $13,($5)           // y = *yptr
          subu $6,1             // count--
          bne $12,$13,cmld3     // if (!(x==y)) ...
compare_loop_down:
          bnez $6,cmld1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
cmld3:  bltu $12,$13,cmld4      // if (x<y) ...
        li $2,1                 // 1
        j $31                   // return
cmld4:  li $2,-1                // -1
        j $31                   // return
        .end compare_loop_down

// extern uintD add_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(add_loop_up)
        .ent add_loop_up // Input in $4,$5,$6,$7, Output in $2
alu1:     // kein Carry
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          addu $4,4             // sourceptr1++
          addu $5,4             // sourceptr2++
          addu $12,$13          // dest = source1 + source2
          sw $12,($6)           // *destptr = dest
          addu $6,4             // destptr++
          bltu $12,$13,alu4     // if (dest < source2) [also Carry] ...
alu2:
          subu $7,1             // count--
add_loop_up:
          bnez $7,alu1          // until (count==0)
        move $2,$0              // 0
        j $31                   // return
alu3:   // Hier Carry
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          addu $4,4             // sourceptr1++
          addu $5,4             // sourceptr2++
          addu $12,$13          // dest = source1 + source2
          addu $12,1            //        + 1
          sw $12,($6)           // *destptr = dest
          addu $6,4             // destptr++
          bgtu $12,$13,alu2     // if (dest > source2) [also kein Carry] ...
alu4:     subu $7,1             // count--
          bnez $7,alu3          // until (count==0)
        li $2,1                 // 1
        j $31                   // return
        .end add_loop_up

// extern uintD addto_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(addto_loop_up)
        .ent addto_loop_up // Input in $4,$5,$6, Output in $2
atlu1:    // kein Carry
          lw $12,($4)           // source1 = *sourceptr
          lw $13,($5)           // source2 = *destptr
          addu $4,4             // sourceptr++
          subu $6,1             // count--
          addu $12,$13          // dest = source1 + source2
          sw $12,($5)           // *destptr = dest
          addu $5,4             // destptr++
          bltu $12,$13,atlu4    // if (dest < source2) [also Carry] ...
addto_loop_up:
atlu2:    bnez $6,atlu1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
atlu3:  // Hier Carry
          lw $12,($4)           // source1 = *sourceptr
          lw $13,($5)           // source2 = *destptr
          addu $4,4             // sourceptr++
          subu $6,1             // count--
          addu $12,$13          // dest = source1 + source2
          addu $12,1            //        + 1
          sw $12,($5)           // *destptr = dest
          addu $5,4             // destptr++
          bgtu $12,$13,atlu2    // if (dest > source2) [also kein Carry] ...
atlu4:    bnez $6,atlu3         // until (count==0)
        li $2,1                 // 1
        j $31                   // return
        .end addto_loop_up

// extern uintD inc_loop_up (uintD* ptr, uintC count);
        .align 2
        DECLARE_FUNCTION(inc_loop_up)
        .ent inc_loop_up // Input in $4,$5, Output in $2
ilu1:     lw $12,($4)           // x = *ptr
          subu $5,1             // count--
          addu $12,1            // x++;
          sw $12,($4)           // *ptr = x
          addu $4,4             // ptr++
          bnez $12,ilu3         // if (!(x==0)) ...
inc_loop_up:
          bnez $5,ilu1          // until (count==0)
        li $2,1                 // 1
        j $31                   // return
ilu3:   move $2,$0              // 0
        j $31                   // return
        .end inc_loop_up

// extern uintD sub_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(sub_loop_up)
        .ent sub_loop_up // Input in $4,$5,$6,$7, Output in $2
slu1:     // kein Carry
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          addu $4,4             // sourceptr1++
          addu $5,4             // sourceptr2++
          subu $7,1             // count--
          bltu $12,$13,slu2     // if (source1 < source2) [also Carry] ...
          subu $12,$13          // dest = source1 - source2
          sw $12,($6)           // *destptr = dest
          addu $6,4             // destptr++
sub_loop_up:
          bnez $7,slu1          // until (count==0)
        move $2,$0              // 0
        j $31                   // return
slu2:     subu $12,$13          // dest = source1 - source2
          sw $12,($6)           // *destptr = dest
          addu $6,4             // destptr++
          bnez $7,slu3          // until (count==0)
        li $2,-1                // -1
        j $31                   // return
slu3:   // Hier Carry
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          addu $4,4             // sourceptr1++
          addu $5,4             // sourceptr2++
          subu $7,1             // count--
          bgtu $12,$13,slu4     // if (source1 > source2) [also kein Carry] ...
          subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($6)           // *destptr = dest
          addu $6,4             // destptr++
          bnez $7,slu3          // until (count==0)
        li $2,-1                // -1
        j $31                   // return
slu4:     subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($6)           // *destptr = dest
          addu $6,4             // destptr++
          bnez $7,slu1          // until (count==0)
        move $2,$0              // 0
        j $31                   // return
        .end sub_loop_up

// extern uintD subx_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);
        .align 2
        DECLARE_FUNCTION(subx_loop_up)
        .ent subx_loop_up // Input in $4,$5,$6,$7,$8, Output in $2
subx_loop_up:
#if ABI_N32
        move $12,$8             // carry
#else
        lw $12,16($sp)          // carry
#endif
        bnez $12,sxlu5          // !(carry==0) ?
        b sxlu2
sxlu1:    // kein Carry
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          addu $4,4             // sourceptr1++
          addu $5,4             // sourceptr2++
          subu $7,1             // count--
          bltu $12,$13,sxlu3    // if (source1 < source2) [also Carry] ...
          subu $12,$13          // dest = source1 - source2
          sw $12,($6)           // *destptr = dest
          addu $6,4             // destptr++
sxlu2:    bnez $7,sxlu1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
sxlu3:    subu $12,$13          // dest = source1 - source2
          sw $12,($6)           // *destptr = dest
          addu $6,4             // destptr++
          bnez $7,sxlu4         // until (count==0)
        li $2,-1                // -1
        j $31                   // return
sxlu4:  // Hier Carry
          lw $12,($4)           // source1 = *sourceptr1
          lw $13,($5)           // source2 = *sourceptr2
          addu $4,4             // sourceptr1++
          addu $5,4             // sourceptr2++
          subu $7,1             // count--
          bgtu $12,$13,sxlu6    // if (source1 > source2) [also kein Carry] ...
          subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($6)           // *destptr = dest
          addu $6,4             // destptr++
sxlu5:    bnez $7,sxlu4         // until (count==0)
        li $2,-1                // -1
        j $31                   // return
sxlu6:    subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($6)           // *destptr = dest
          addu $6,4             // destptr++
          bnez $7,sxlu1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
        .end subx_loop_up

// extern uintD subfrom_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
        .align 2
        DECLARE_FUNCTION(subfrom_loop_up)
        .ent subfrom_loop_up // Input in $4,$5,$6,$7, Output in $2
sflu1:    // kein Carry
          lw $12,($5)           // source1 = *destptr
          lw $13,($4)           // source2 = *sourceptr
          addu $4,4             // sourceptr++
          subu $6,1             // count--
          bltu $12,$13,sflu2    // if (source1 < source2) [also Carry] ...
          subu $12,$13          // dest = source1 - source2
          sw $12,($5)           // *destptr = dest
          addu $5,4             // destptr++
subfrom_loop_up:
          bnez $6,sflu1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
sflu2:    subu $12,$13          // dest = source1 - source2
          sw $12,($5)           // *destptr = dest
          addu $5,4             // destptr++
          bnez $6,sflu3         // until (count==0)
        li $2,-1                // -1
        j $31                   // return
sflu3:  // Hier Carry
          lw $12,($5)           // source1 = *destptr
          lw $13,($4)           // source2 = *sourceptr
          addu $4,4             // sourceptr++
          subu $6,1             // count--
          bgtu $12,$13,sflu4    // if (source1 > source2) [also kein Carry] ...
          subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($5)           // *destptr = dest
          addu $5,4             // destptr++
          bnez $6,sflu3         // until (count==0)
        li $2,-1                // -1
        j $31                   // return
sflu4:    subu $12,$13          // dest = source1 - source2
          subu $12,1            //        - 1
          sw $12,($5)           // *destptr = dest
          addu $5,4             // destptr++
          bnez $6,sflu1         // until (count==0)
        move $2,$0              // 0
        j $31                   // return
        .end subfrom_loop_up

// extern uintD dec_loop_up (uintD* ptr, uintC count);
        .align 2
        DECLARE_FUNCTION(dec_loop_up)
        .ent dec_loop_up // Input in $4,$5, Output in $2
dlu1:     lw $12,($4)           // x = *ptr
          subu $5,1             // count--
          bnez $12,dlu3         // if (!(x==0)) ...
          subu $12,1            // x--;
          sw $12,($4)           // *ptr = x
          addu $4,4             // ptr++
dec_loop_up:
          bnez $5,dlu1          // until (count==0)
        li $2,-1                // -1
        j $31                   // return
dlu3:   subu $12,1              // x--;
        sw $12,($4)             // *ptr = x
        move $2,$0              // 0
        j $31                   // return
        .end dec_loop_up

// extern uintD neg_loop_up (uintD* ptr, uintC count);
        .align 2
        DECLARE_FUNCTION(neg_loop_up)
        .ent neg_loop_up // Input in $4,$5, Output in $2
        // erstes Digit /=0 suchen:
nlu1:     lw $12,($4)           // x = *ptr
          subu $5,1             // count--
          bnez $12,nlu3         // if (!(x==0)) ...
          addu $4,4             // ptr++
neg_loop_up:
          bnez $5,nlu1          // until (count==0)
        move $2,$0              // 0
        j $31                   // return
nlu3:   // erstes Digit /=0 gefunden, ab jetzt gibt's Carrys
        // 1 Digit negieren:
        subu $12,$0,$12         // x = -x
        sw $12,($4)             // *ptr = x
        // alle anderen Digits invertieren:
        b nlu5
nlu4:     lw $12,($4)           // x = *xptr
          subu $5,1             // count--
          nor $12,$0            // x = ~x
          sw $12,($4)           // *xptr = x
nlu5:     addu $4,4             // xptr++
          bnez $5,nlu4          // until (count==0)
        li $2,-1                // -1
        j $31                   // return
        .end neg_loop_up

#endif

