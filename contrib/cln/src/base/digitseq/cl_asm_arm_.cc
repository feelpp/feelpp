// ariarm.d (c) Copyright 1994, 1997 P.J.Burwood
// little-endian modifications (c) Copyright 1996 B. Haible
// external routines for arilev1.d
// Processor: ARM in APCS mode
// Assembler-Syntax: ObjAsm under RISC OS, GAS otherwise
// Assumptions: intCsize=32, intDsize=32.
// Parameter passing conventions: APCS means that registers a1-a4 and ip
//   do not have to be preserved across function calls.
// Note: A sequence of up to 4 conditional instructions is used in preference
//   to a branch.


#ifdef __riscos

// ObjAsm syntax

a1      RN      0
a2      RN      1
a3      RN      2
a4      RN      3
v1      RN      4
v2      RN      5
v3      RN      6
v4      RN      7
v5      RN      8
v6      RN      9
sl      RN      10
fp      RN      11
ip      RN      12
sp      RN      13
lr      RN      14
pc      RN      15

f0      FN      0
f1      FN      1
f2      FN      2
f3      FN      3
f4      FN      4
f5      FN      5
f6      FN      6
f7      FN      7

#define C(x) _##x
#define EXPORT(x) EXPORT x
#define DECLARE_FUNCTION(x)
#define GLABEL(x) _##x
#define LABEL(x) _##x

        AREA    |C$$code|,CODE,READONLY

#else

// GAS syntax

a1      .req    r0
a2      .req    r1
a3      .req    r2
a4      .req    r3
v1      .req    r4
v2      .req    r5
v3      .req    r6
v4      .req    r7
v5      .req    r8
v6      .req    r9
rfp     .req    r9
sl      .req    r10
fp      .req    r11
ip      .req    r12
sp      .req    r13
lr      .req    r14
pc      .req    r15

#define C(x) _##x
#define EXPORT(x) .global _##x
#if defined(__NetBSD__)
#define DECLARE_FUNCTION(x) .type _##x,%function
#else
#define DECLARE_FUNCTION(x)
#endif
#define GLABEL(x) _##x##:
#define LABEL(x) x##:
#define RRX rrx
#define END

#endif


#if defined(__arm7m__) || defined(__arm8__) || defined(__arm9__) || defined(__strongarm__)
  // ARM7M and later have 32x32 -> 64 multiplies which execute in 2-4 clocks.
  #define HAVE_umull
#endif


#if defined(__GNUC__) && 0
  // With GNU C, we would like to pass the second return value in a2, don't
  // need a global variable. Unfortunately, the current Acorn gcc crashes if
  // we declare an appropriate local register variable with __asm__.
  // It would be possible to declare the functions as returning a 64-bit
  // result, but given the quality of gcc code dealing with 64-bit entities
  // and the subtleties of 64-bit returns values (passed in register or in
  // memory?) we now let it be.
#else
  // Use three global variables.
  #define MULU32_HIGH
  #define DIVU_16_REST
  #define DIVU_32_REST
#endif

#ifdef __riscos

#ifdef MULU32_HIGH
ptr_mulu32_high
        IMPORT  mulu32_high
        DCD     mulu32_high
#endif
#ifdef DIVU_16_REST
ptr_divu_16_rest
        IMPORT  divu_16_rest
        DCD     divu_16_rest
#endif
#ifdef DIVU_32_REST
ptr_divu_32_rest
        IMPORT  divu_32_rest
        DCD     divu_32_rest
#endif

#else

#ifdef MULU32_HIGH
ptr_mulu32_high:
        .word   _mulu32_high
        .align  0
#endif
#ifdef DIVU_16_REST
ptr_divu_16_rest:
        .word   _divu_16_rest
        .align  0
#endif
#ifdef DIVU_32_REST
ptr_divu_32_rest:
        .word   _divu_32_rest
        .align  0
#endif

#endif


// extern uint32 mulu32_ (uint32 x, uint32 y);
//       entry
//               a1 = x
//               a2 = y
//       exit
//               a1 = low32(x*y)
//               a2 = high32(x*y)
//               mulu32_high = high32(x*y)
//               a3,a4,ip destroyed
        EXPORT(mulu32_)
        DECLARE_FUNCTION(mulu32_)
GLABEL(mulu32_)
#ifdef HAVE_umull
        MOV     a3,a2
        UMULL   a1,a2,a3,a1
#else
        MOV     ip,a1,LSR #16           // temp := top half of x
        MOV     a3,a2,LSR #16           // hi := top half of y
        BIC     a1,a1,ip,LSL #16        // x  := bottom half of x
        BIC     a2,a2,a3,LSL #16        // y  := bottom half of y
        MUL     a4,a1,a2                // low section of result
        MUL     a2,ip,a2                // ) middle sections
        MUL     a1,a3,a1                // )   of result
        MUL     a3,ip,a3                // high section of result
        ADDS    a2,a2,a1                // add middle sections
                                        // (can't use mla as we need carry)
        ADDCS   a3,a3,#0x10000          // carry from above add
        ADDS    a1,a4,a2,LSL #16        // x is now bottom 32 bits of result
        ADC     a2,a3,a2,LSR #16        // hi is top 32 bits
#endif
#ifdef MULU32_HIGH
        LDR     a3,[pc,#ptr_mulu32_high-.-8]
        STR     a2,[a3,#0]
#endif
        MOVS    pc,lr

// extern uint16 divu_3216_1616_ (uint32 x, uint16 y);
//       entry
//               a1 = x
//               a2 = y
//       exit
//               a1 = q = floor(x/y)
//               a2 = r = x-q*y
//               divu_16_rest = r = x-q*y
//               a3 destroyed
        EXPORT(divu_3216_1616_)
        DECLARE_FUNCTION(divu_3216_1616_)
GLABEL(divu_3216_1616_)
        // see cl_low_div.cc for algorithm
        // in that notation: a1 = r, a2 = -s.
        MOV     a2,a2,LSL#15            // multiply divisor by 2^15
        RSB     a2,a2,#0                // negate divisor
        ADDS    a1,a2,a1                // dividend = dividend + -divisor/2
        SUBCC   a1,a1,a2                // dividend = dividend - -divisor/2
        ADCS    a1,a2,a1,LSL#1          // dividend = dividend*2 + -divisor
                                        // and shift quotient
        SUBCC   a1,a1,a2                // do this another 14 times
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2
        ADCS    a1,a2,a1,LSL#1
        SUBCC   a1,a1,a2                // do the last conditional subtraction
        MOV     a2,a1,LSR#15            // move remainder into a2 and shift
        ADC     a1,a1,a1                // move last bit of quotient in
        MOV     a1,a1,LSL#16            // AND out top 16 bits by shifting up
        MOV     a1,a1,LSR#16            // and back down again
#ifdef DIVU_16_REST
        LDR     a3,[pc,#ptr_divu_16_rest-.-8]   // save rest so can be picked up later
        STR     a2,[a3,#0]              // the result is 16 bits
#endif
        MOVS    pc, lr

// extern uint32 divu_6432_3232_ (uint32 xhi, uint32 xlo, uint32 y); // -> Quotient q
// extern uint32 divu_32_rest;                                       // -> Rest r
//       see cl_low_div.cc for algorithm
//       entry
//               a1 = xhi (dividend)
//               a2 = xlo (dividend)
//               a3 = y (divisor)
//       exit
//               a1 = 32 bit quotient
//               a2 = 32 bit remainder
//               a3, a4 destroyed
        EXPORT(divu_6432_3232_)
        DECLARE_FUNCTION(divu_6432_3232_)
GLABEL(divu_6432_3232_)
        STMFD   sp!, {v1,v2,v3,v4,v5,v6,lr}
        MOV     v2, a2                  // = xlo
        MOV     v1, a3                  // = y
        CMP     a3,#0x10000             // y <= (uint32)(bit(16)-1)
        BCS     divu_6432_3232_l1
        MOV     a2, v2, LSR #16
        ORR     a1, a2, a1, ASL #16     // = highlow32(low16(xhi),high16(xlo))
        MOV     a2, v1
        BL      C(divu_3216_1616_)
        MOV     v3, a1                  // = q1
        MOV     a1, v2, ASL #16
        MOV     a1, a1, LSR #16
        ORR     a1, a1, a2, ASL #16     // = highlow32(r1,low16(xlo))
        MOV     a2, v1
        BL      C(divu_3216_1616_)
        ORR     a1, a1, v3, ASL #16     // = highlow32(q1,q0)
#ifdef DIVU_32_REST
        LDR     a4,[pc,#ptr_divu_32_rest-.-8]
        STR     a2,[a4,#0]              // divu_32_rest = remainder
#endif
        LDMFD   sp!, {v1,v2,v3,v4,v5,v6,pc}^

LABEL(divu_6432_3232_l1)
        MOV     v3, #0                  // s = 0
        MOVS    a4, v1, LSR #16         // while ((sint32)y >= 0)
        ADDEQ   v3, v3, #16             //   { y = y<<1; s++; }
        MOVEQ   v1, v1, ASL #16
        MOVS    a4, v1, LSR #24
        ADDEQ   v3, v3, #8
        MOVEQ   v1, v1, ASL #8
        MOVS    a4, v1, LSR #28
        ADDEQ   v3, v3, #4
        MOVEQ   v1, v1, ASL #4
        MOVS    a4, v1, LSR #30
        ADDEQ   v3, v3, #2
        MOVEQ   v1, v1, ASL #2
        MOVS    a4, v1, LSR #31
        ADDEQ   v3, v3, #1
        MOVEQ   v1, v1, ASL #1

        CMPS    v3, #0
        MOVNE   a2, a1, ASL v3          // if (!(s==0))
        RSBNE   a1, v3, #32             //   { xhi = (xhi << s)
        ORRNE   a1, a2, v2, LSR a1      //         | (xlo >> (32-s));
        MOVNE   v2, v2, ASL v3          //     xlo = xlo << s; }
        ADD     a2, v1, #0x10000        // y1_1 = high16(y)+1
        MOVS    v5, a2, LSR #16         // if (y1_1 = 0)
        MOVEQ   v4, a1, ASL #16         // r16 = low16(xhi) * 2^16
        MOVEQ   a1, a1, LSR #16         // q1 = high16(xhi)
        MOVNE   a2, v5
        BLNE    C(divu_3216_1616_)      // divu_3216_1616(xhi,y1_1, q1=,r16=)
        MOVNE   v4, a2, ASL #16         // r16 = r16 * 2^16
        ORR     v4, v4, v2, LSR #16     // r = highlow32(r16,high16(xlo))
        MOV     a4, v1, ASL #16         // tmp = mulu16(low16(y),q1)
        MOV     a4, a4, LSR #16
        MUL     a3, a4, a1
        RSB     a3, a3, a1, ASL #16     // r2 = highlow32_0(q1) - tmp
        MOV     v6, a1                  // = q1
        ADDS    a1, v4, a3              // r += r2
        ADDCS   v6, v6, #1              // if ( r < r2 ) { q1 += 1
        SUBCS   a1, a1, v1              //                 r -= y }
        CMP     a1, v1                  // if (r >= y)
        ADDCS   v6, v6, #1              //     { q1 += 1
        SUBCS   a1, a1, v1              //       r -= y }
        CMP     v5, #0                  // if (y1_1 = 0)
        MOVEQ   v4, a1, ASL #16         //    { r16 = low16(r) * 2^16
        MOVEQ   a1, a1, LSR #16         //      q0  = high16(r) }
        MOVNE   a2, v5
        BLNE    C(divu_3216_1616_)      // divu_3216_1616(r,y1_1, q0=,r16=)
        MOVNE   v4, a2, ASL #16         // r16 = r16 * 2^16
        MOV     v2, v2, ASL #16
        ORR     v4, v4, v2, LSR #16     // r = highlow32(r16,low16(xlo))
        MOV     a4, v1, ASL #16         // tmp = mulu16(low16(y),q0)
        MOV     a4, a4, LSR #16
        MUL     a3, a4, a1
        RSB     a3, a3, a1, ASL #16     // r2 = highlow32_0(q0) - tmp
        ADDS    v4, v4, a3              // r += r2
        ADDCS   a1, a1, #1              // if ( r < r2 ) { q0 += 1
        SUBCS   v4, v4, v1              //                 r -= y }
        CMP     v4, v1                  // if (r >= y)
        ADDCS   a1, a1, #1              //     { q0 += 1
        SUBCS   v4, v4, v1              //       r -= y }
        MOV     a2, v4, LSR v3          // remainder = r >> s
        ORR     a1, a1, v6, ASL #16     // return highlow32(q1,q0)
#ifdef DIVU_32_REST
        LDR     a3,[pc,#ptr_divu_32_rest-.-8]
        STR     a2,[a3,#0]              // divu_32_rest = remainder
#endif
        LDMFD   sp!, {v1,v2,v3,v4,v5,v6,pc}^

// extern uintD* copy_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
//       entry
//               a1 = source pointer
//               a2 = destination pointer
//               a3 = count of words to store
//       exit
//               a1 = address of last word stored + 1
//               a2 - a4, ip destroyed
        EXPORT(copy_loop_up)            // word aligned copy loop up
        DECLARE_FUNCTION(copy_loop_up)
GLABEL(copy_loop_up)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     copy_loop_up_l1         // yup, so branch
        CMP     a4,#2                   // copy the first 1-3 words
        LDR     a4,[a1],#4              // to align the total to a multiple
        STR     a4,[a2],#4              // of 4 words
        LDRGE   a4,[a1],#4
        STRGE   a4,[a2],#4
        LDRGT   a4,[a1],#4
        STRGT   a4,[a2],#4
LABEL(copy_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,a2                   // return addr of last word stored
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1,lr}             // save work regs
LABEL(copy_loop_up_l2)
        LDMIA   a1!,{a3,v1,ip,lr}       // copy 4 words in one go
        STMIA   a2!,{a3,v1,ip,lr}
        SUBS    a4,a4,#8                // decrement counter by 8
        LDMGEIA a1!,{a3,v1,ip,lr}       // if count still positive then copy
        STMGEIA a2!,{a3,v1,ip,lr}       // 4 more words
        BGT     copy_loop_up_l2         // and loop
        MOV     a1,a2                   // return addr of last word stored
        LDMFD   sp!,{v1,pc}^            // restore work regs and return

// extern uintD* copy_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
//       entry
//               a1 = source pointer
//               a2 = destination pointer
//               a3 = count of words to store
//       exit
//               a1 = address of last word stored
//               a2 - a4, ip destroyed
        EXPORT(copy_loop_down)          // word aligned copy loop down
        DECLARE_FUNCTION(copy_loop_down)
GLABEL(copy_loop_down)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     copy_loop_down_l1       // yup, so branch
        CMP     a4,#2                   // copy the first 1-3 words
        LDR     a4,[a1,#-4]!            // to align the total to a multiple
        STR     a4,[a2,#-4]!            // of 4 words
        LDRGE   a4,[a1,#-4]!
        STRGE   a4,[a2,#-4]!
        LDRGT   a4,[a1,#-4]!
        STRGT   a4,[a2,#-4]!
LABEL(copy_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,a2                   // return addr of last word stored
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1,lr}             // save work regs
LABEL(copy_loop_down_l2)
        LDMDB   a1!,{a3,v1,ip,lr}       // copy 4 words in one go
        STMDB   a2!,{a3,v1,ip,lr}
        SUBS    a4,a4,#8                // decrement counter by 8
        LDMGEDB a1!,{a3,v1,ip,lr}       // if count still positive then copy
        STMGEDB a2!,{a3,v1,ip,lr}       // 4 more words
        BGT     copy_loop_down_l2       // and loop
        MOV     a1,a2                   // return addr of last word stored
        LDMFD   sp!,{v1,pc}^            // restore work regs and return

// extern uintD* clear_loop_up (uintD* destptr, uintC count);
//       entry
//               a1 = destination pointer
//               a2 = count of words to store
//       exit
//               a1 = address of last word stored + 1
//               a2 - a4, ip destroyed
        EXPORT(clear_loop_up)           // word aligned clear loop up
        DECLARE_FUNCTION(clear_loop_up)
GLABEL(clear_loop_up)
        MOV     a3,#0                   // set filler to 0
                                        // and drop into fill_loop_up

// extern uintD* fill_loop_up (uintD* destptr, uintC count, uintD filler);
//       entry
//               a1 = destination pointer
//               a2 = count of words to store
//               a3 = word to store
//       exit
//               a1 = address of last word stored + 1
//               a2 - a4, ip destroyed
        EXPORT(fill_loop_up)            // word aligned fill loop up
        DECLARE_FUNCTION(fill_loop_up)
GLABEL(fill_loop_up)
        ANDS    a4,a2,#3                // multiple of 4 words ?
        BEQ     fill_loop_up_l1         // yup, so branch
        CMP     a4,#2                   // store the first 1-3 words
        STR     a3,[a1],#4              // to align the total to a multiple
        STRGE   a3,[a1],#4              // of 4 words
        STRGT   a3,[a1],#4
LABEL(fill_loop_up_l1)
        BICS    a4,a2,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1,lr}             // save work regs
        MOV     v1,a3                   // copy filler to three other
        MOV     ip,a3                   // registers
        MOV     lr,a3
LABEL(fill_loop_up_l2)
        STMIA   a1!,{a3,v1,ip,lr}       // store 4 fillers in one go
        SUBS    a4,a4,#8                // decrement counter by 8
        STMGEIA a1!,{a3,v1,ip,lr}       // if count still positive then store 4
        BGT     fill_loop_up_l2         // more and loop
        LDMFD   sp!,{v1,pc}^            // restore work regs and return


// extern uintD* clear_loop_down (uintD* destptr, uintC count);
//       entry
//               a1 = destination pointer
//               a2 = count of words to store
//       exit
//               a1 = address of last word stored + 1
//               a2 - a4, ip destroyed
        EXPORT(clear_loop_down)         // word aligned clear loop down
        DECLARE_FUNCTION(clear_loop_down)
GLABEL(clear_loop_down)
        MOV     a3,#0                   // set filler to 0
                                        // and drop into fill_loop_down

// extern uintD* fill_loop_down (uintD* destptr, uintC count, uintD filler);
//       entry
//               a1 = destination pointer
//               a2 = count of words to store
//               a3 = word to store
//       exit
//               a1 = address of last word stored
//               a2 - a4, ip destroyed
        EXPORT(fill_loop_down)          // word aligned fill loop down
        DECLARE_FUNCTION(fill_loop_down)
GLABEL(fill_loop_down)
        ANDS    a4,a2,#3                // multiple of 4 words ?
        BEQ     fill_loop_down_l1       // yup, so branch
        CMP     a4,#2                   // store the first 1-3 words
        STR     a3,[a1,#-4]!            // to align the total to a multiple
        STRGE   a3,[a1,#-4]!            // of 4 words
        STRGT   a3,[a1,#-4]!
LABEL(fill_loop_down_l1)
        BICS    a4,a2,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1,lr}             // save work regs
        MOV     v1,a3                   // copy filler to three other
        MOV     ip,a3                   // registers
        MOV     lr,a3
LABEL(fill_loop_down_l2)
        STMDB   a1!,{a3,v1,ip,lr}       // store 4 fillers in one go
        SUBS    a4,a4,#8                // decrement counter by 8
        STMGEDB a1!,{a3,v1,ip,lr}       // if count still positive then store 4
        BGT     fill_loop_down_l2       // more and loop
        LDMFD   sp!,{v1,pc}^            // restore work regs and return

// extern void test_loop_up (uintD* xptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = count of words to be TESTed
//       exit
//               a1 = TRUE if any words are non-zero else FALSE
//               a2 - a4, ip destroyed
        EXPORT(test_loop_up)            // word aligned test loop up
        DECLARE_FUNCTION(test_loop_up)
GLABEL(test_loop_up)
        MOV     ip,a1                   // move xptr to ip
        MOV     a1,#1                   // set result to TRUE
        ANDS    a3,a2,#3                // multiple of 4 words ?
        BEQ     test_loop_up_l1         // yup, so branch
        LDR     a4,[ip],#4              // TEST the first 1-3 words
        TEQ     a4,#0                   // align the total to a multiple of 4
        MOVNES  pc,lr                   // return TRUE if AND_TEST ok
        CMP     a3,#2
        BLT     test_loop_up_l1         // need to branch 'cos PSR set
        LDRGE   a4,[ip],#4              // when checking against zero
        TEQGE   a4,#0
        MOVNES  pc,lr
        CMP     a3,#2
        BLE     test_loop_up_l1         // need to branch 'cos PSR set
        LDRGT   a4,[ip],#4              // when checking against zero
        TEQGT   a4,#0
        MOVNES  pc,lr
LABEL(test_loop_up_l1)
        BICS    a4,a2,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // return FALSE
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1,lr}             // save work regs
LABEL(test_loop_up_l2)
        LDMIA   ip!,{a2,a3,v1,lr}       // load 4 words in one go
        TEQ     a2,#0                   // TEST the four words
        TEQEQ   a3,#0
        TEQEQ   v1,#0
        TEQEQ   lr,#0
        LDMNEFD sp!,{v1,pc}^
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     test_loop_up_l2         // if count still positive then loop
        MOV     a1,#0
        LDMFD   sp!,{v1,pc}^            // restore work regs and return

// extern void test_loop_down (uintD* xptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = count of words to be TESTed
//       exit
//               a1 = TRUE if any words are non-zero else FALSE
//               a2 - a4, ip destroyed
        EXPORT(test_loop_down)          // word aligned test loop down
        DECLARE_FUNCTION(test_loop_down)
GLABEL(test_loop_down)
        MOV     ip,a1                   // move xptr to ip
        MOV     a1,#1                   // set result to TRUE
        ANDS    a3,a2,#3                // multiple of 4 words ?
        BEQ     test_loop_down_l1       // yup, so branch
        LDR     a4,[ip,#-4]!            // TEST the first 1-3 words
        TEQ     a4,#0                   // align the total to a multiple of 4
        MOVNES  pc,lr                   // return TRUE if AND_TEST ok
        CMP     a3,#2
        BLT     test_loop_down_l1       // need to branch 'cos PSR set
        LDRGE   a4,[ip,#-4]!            // when checking against zero
        TEQGE   a4,#0
        MOVNES  pc,lr
        CMP     a3,#2
        BLE     test_loop_down_l1       // need to branch 'cos PSR set
        LDRGT   a4,[ip,#-4]!            // when checking against zero
        TEQGT   a4,#0
        MOVNES  pc,lr
LABEL(test_loop_down_l1)
        BICS    a4,a2,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // return FALSE
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1,lr}             // save work regs
LABEL(test_loop_down_l2)
        LDMDB   ip!,{a2,a3,v1,lr}       // load 4 words in one go
        TEQ     a2,#0                   // TEST the four words
        TEQEQ   a3,#0
        TEQEQ   v1,#0
        TEQEQ   lr,#0
        LDMNEFD sp!,{v1,pc}^
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     test_loop_down_l2       // if count still positive then loop
        MOV     a1,#0
        LDMFD   sp!,{v1,pc}^            // restore work regs and return

#if CL_DS_BIG_ENDIAN_P

// extern void or_loop_up (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be ORed
//       exit
//               xptr |= yptr for count words
//               a1 - a4, ip destroyed
        EXPORT(or_loop_up)              // word aligned or loop up
        DECLARE_FUNCTION(or_loop_up)
GLABEL(or_loop_up)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     or_loop_up_l1           // yup, so branch
        CMP     a4,#2                   // OR the first 1-3 words
        LDR     a4,[a2],#4              // to align the total to a multiple
        LDR     ip,[a1]                 // of 4 words
        ORR     ip,ip,a4
        STR     ip,[a1],#4
        BLT     or_loop_up_l1           // better to branch than skip instrs.
        LDRGE   a4,[a2],#4
        LDRGE   ip,[a1]
        ORRGE   ip,ip,a4
        STRGE   ip,[a1],#4
        LDRGT   a4,[a2],#4
        LDRGT   ip,[a1]
        ORRGT   ip,ip,a4
        STRGT   ip,[a1],#4
LABEL(or_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(or_loop_up_l2)
        LDMIA   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMIA   a1,{v3,v4,v5,lr}        // load target words
        ORR     v3,v3,a3                // OR the four words
        ORR     v4,v4,v1
        ORR     v5,v5,v2
        ORR     lr,lr,ip
        STMIA   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     or_loop_up_l2           // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

#endif

// extern void xor_loop_up (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be XORed
//       exit
//               xptr ^= yptr for count words
//               a1 - a4, ip destroyed
        EXPORT(xor_loop_up)             // word aligned xor loop up
        DECLARE_FUNCTION(xor_loop_up)
GLABEL(xor_loop_up)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     xor_loop_up_l1          // yup, so branch
        CMP     a4,#2                   // XOR the first 1-3 words
        LDR     a4,[a2],#4              // to align the total to a multiple
        LDR     ip,[a1]                 // of 4 words
        EOR     ip,ip,a4
        STR     ip,[a1],#4
        BLT     xor_loop_up_l1          // better to branch than skip instrs.
        LDRGE   a4,[a2],#4
        LDRGE   ip,[a1]
        EORGE   ip,ip,a4
        STRGE   ip,[a1],#4
        LDRGT   a4,[a2],#4
        LDRGT   ip,[a1]
        EORGT   ip,ip,a4
        STRGT   ip,[a1],#4
LABEL(xor_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(xor_loop_up_l2)
        LDMIA   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMIA   a1,{v3,v4,v5,lr}        // load target words
        EOR     v3,v3,a3                // XOR the four words
        EOR     v4,v4,v1
        EOR     v5,v5,v2
        EOR     lr,lr,ip
        STMIA   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     xor_loop_up_l2          // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

#if CL_DS_BIG_ENDIAN_P

// extern void and_loop_up (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be ANDed
//       exit
//               xptr &= yptr for count words
//               a1 - a4, ip destroyed
        EXPORT(and_loop_up)             // word aligned and loop up
        DECLARE_FUNCTION(and_loop_up)
GLABEL(and_loop_up)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     and_loop_up_l1          // yup, so branch
        CMP     a4,#2                   // AND the first 1-3 words
        LDR     a4,[a2],#4              // to align the total to a multiple
        LDR     ip,[a1]                 // of 4 words
        AND     ip,ip,a4
        STR     ip,[a1],#4
        BLT     and_loop_up_l1          // better to branch than skip instrs.
        LDRGE   a4,[a2],#4
        LDRGE   ip,[a1]
        ANDGE   ip,ip,a4
        STRGE   ip,[a1],#4
        LDRGT   a4,[a2],#4
        LDRGT   ip,[a1]
        ANDGT   ip,ip,a4
        STRGT   ip,[a1],#4
LABEL(and_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(and_loop_up_l2)
        LDMIA   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMIA   a1,{v3,v4,v5,lr}        // load target words
        AND     v3,v3,a3                // AND the four words
        AND     v4,v4,v1
        AND     v5,v5,v2
        AND     lr,lr,ip
        STMIA   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     and_loop_up_l2          // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void eqv_loop_up (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be XORed
//       exit
//               xptr = ~(xptr ^ yptr) for count words
//               a1 - a4, ip destroyed
        EXPORT(eqv_loop_up)             // word aligned eqv loop up
        DECLARE_FUNCTION(eqv_loop_up)
GLABEL(eqv_loop_up)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     eqv_loop_up_l1          // yup, so branch
        CMP     a4,#2                   // EQV the first 1-3 words
        LDR     a4,[a2],#4              // to align the total to a multiple
        LDR     ip,[a1]                 // of 4 words
        EOR     ip,ip,a4
        MVN     ip,ip
        STR     ip,[a1],#4
        BLT     eqv_loop_up_l1          // better to branch than skip instrs.
        LDRGE   a4,[a2],#4
        LDRGE   ip,[a1]
        EORGE   ip,ip,a4
        MVNGE   ip,ip
        STRGE   ip,[a1],#4
        BLE     eqv_loop_up_l1          // better to branch than skip instrs.
        LDRGT   a4,[a2],#4
        LDRGT   ip,[a1]
        EORGT   ip,ip,a4
        MVNGT   ip,ip
        STRGT   ip,[a1],#4
LABEL(eqv_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(eqv_loop_up_l2)
        LDMIA   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMIA   a1,{v3,v4,v5,lr}        // load target words
        EOR     v3,v3,a3                // EVQ the four words
        MVN     v3,v3
        EOR     v4,v4,v1
        MVN     v4,v4
        EOR     v5,v5,v2
        MVN     v5,v5
        EOR     lr,lr,ip
        MVN     lr,lr
        STMIA   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     eqv_loop_up_l2          // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void nand_loop_up (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be NANDed
//       exit
//               xptr = ~(xptr & yptr) for count words
//               a1 - a4, ip destroyed
        EXPORT(nand_loop_up)            // word aligned nand loop up
        DECLARE_FUNCTION(nand_loop_up)
GLABEL(nand_loop_up)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     nand_loop_up_l1         // yup, so branch
        CMP     a4,#2                   // NAND the first 1-3 words
        LDR     a4,[a2],#4              // to align the total to a multiple
        LDR     ip,[a1]                 // of 4 words
        AND     ip,ip,a4
        MVN     ip,ip
        STR     ip,[a1],#4
        BLT     nand_loop_up_l1         // better to branch than skip instrs.
        LDRGE   a4,[a2],#4
        LDRGE   ip,[a1]
        ANDGE   ip,ip,a4
        MVNGE   ip,ip
        STRGE   ip,[a1],#4
        BLE     nand_loop_up_l1         // better to branch than skip instrs.
        LDRGT   a4,[a2],#4
        LDRGT   ip,[a1]
        ANDGT   ip,ip,a4
        MVNGT   ip,ip
        STRGT   ip,[a1],#4
LABEL(nand_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(nand_loop_up_l2)
        LDMIA   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMIA   a1,{v3,v4,v5,lr}        // load target words
        AND     v3,v3,a3                // NAND the four words
        MVN     v3,v3
        AND     v4,v4,v1
        MVN     v4,v4
        AND     v5,v5,v2
        MVN     v5,v5
        AND     lr,lr,ip
        MVN     lr,lr
        STMIA   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     nand_loop_up_l2         // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void nor_loop_up (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be NORed
//       exit
//               xptr = ~(xptr | yptr) for count words
//               a1 - a4, ip destroyed
        EXPORT(nor_loop_up)             // word aligned nor loop up
        DECLARE_FUNCTION(nor_loop_up)
GLABEL(nor_loop_up)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     nor_loop_up_l1          // yup, so branch
        CMP     a4,#2                   // NOR the first 1-3 words
        LDR     a4,[a2],#4              // to align the total to a multiple
        LDR     ip,[a1]                 // of 4 words
        ORR     ip,ip,a4
        MVN     ip,ip
        STR     ip,[a1],#4
        BLT     nor_loop_up_l1          // better to branch than skip instrs.
        LDRGE   a4,[a2],#4
        LDRGE   ip,[a1]
        ORRGE   ip,ip,a4
        MVNGE   ip,ip
        STRGE   ip,[a1],#4
        BLE     nor_loop_up_l1          // better to branch than skip instrs.
        LDRGT   a4,[a2],#4
        LDRGT   ip,[a1]
        ORRGT   ip,ip,a4
        MVNGT   ip,ip
        STRGT   ip,[a1],#4
LABEL(nor_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(nor_loop_up_l2)
        LDMIA   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMIA   a1,{v3,v4,v5,lr}        // load target words
        ORR     v3,v3,a3                // NOR the four words
        MVN     v3,v3
        ORR     v4,v4,v1
        MVN     v4,v4
        ORR     v5,v5,v2
        MVN     v5,v5
        ORR     lr,lr,ip
        MVN     lr,lr
        STMIA   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     nor_loop_up_l2          // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void andc2_loop_up (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be ANDC2ed
//       exit
//               xptr = xptr & ~yptr for count words
//               a1 - a4, ip destroyed
        EXPORT(andc2_loop_up)           // word aligned andc2 loop up
        DECLARE_FUNCTION(andc2_loop_up)
GLABEL(andc2_loop_up)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     andc2_loop_up_l1        // yup, so branch
        CMP     a4,#2                   // ANDC2 the first 1-3 words
        LDR     a4,[a2],#4              // to align the total to a multiple
        LDR     ip,[a1]                 // of 4 words
        BIC     ip,ip,a4
        STR     ip,[a1],#4
        BLT     andc2_loop_up_l1        // better to branch than skip instrs.
        LDRGE   a4,[a2],#4
        LDRGE   ip,[a1]
        BICGE   ip,ip,a4
        STRGE   ip,[a1],#4
        LDRGT   a4,[a2],#4
        LDRGT   ip,[a1]
        BICGT   ip,ip,a4
        STRGT   ip,[a1],#4
LABEL(andc2_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(andc2_loop_up_l2)
        LDMIA   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMIA   a1,{v3,v4,v5,lr}        // load target words
        BIC     v3,v3,a3                // ANDC2 the four words
        BIC     v4,v4,v1
        BIC     v5,v5,v2
        BIC     lr,lr,ip
        STMIA   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     andc2_loop_up_l2        // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void orc2_loop_up (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be XORed
//       exit
//               xptr = xptr | ~yptr for count words
//               a1 - a4, ip destroyed
        EXPORT(orc2_loop_up)            // word aligned orc2 loop up
        DECLARE_FUNCTION(orc2_loop_up)
GLABEL(orc2_loop_up)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     orc2_loop_up_l1         // yup, so branch
        CMP     a4,#2                   // ORC2 the first 1-3 words
        LDR     a4,[a2],#4              // to align the total to a multiple
        LDR     ip,[a1]                 // of 4 words
        MVN     a4,a4
        ORR     ip,ip,a4
        STR     ip,[a1],#4
        BLT     orc2_loop_up_l1         // better to branch than skip instrs.
        LDRGE   a4,[a2],#4
        LDRGE   ip,[a1]
        MVNGE   a4,a4
        ORRGE   ip,ip,a4
        STRGE   ip,[a1],#4
        BLE     orc2_loop_up_l1         // better to branch than skip instrs.
        LDRGT   a4,[a2],#4
        LDRGT   ip,[a1]
        MVNGT   a4,a4
        ORRGT   ip,ip,a4
        STRGT   ip,[a1],#4
LABEL(orc2_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(orc2_loop_up_l2)
        LDMIA   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMIA   a1,{v3,v4,v5,lr}        // load target words
        MVN     a3,a3                   // ORC2 the four words
        ORR     v3,v3,a3
        MVN     v1,v1
        ORR     v4,v4,v1
        MVN     v2,v2
        ORR     v5,v5,v2
        MVN     ip,ip
        ORR     lr,lr,ip
        STMIA   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     orc2_loop_up_l2         // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void not_loop_up (uintD* xptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = count of words to be NOTed
//       exit
//               xptr = ~xptr for count words
//               a1 - a4, ip destroyed
        EXPORT(not_loop_up)             // word aligned not loop up
        DECLARE_FUNCTION(not_loop_up)
GLABEL(not_loop_up)
        ANDS    a3,a2,#3                // multiple of 4 words ?
        BEQ     not_loop_up_l1          // yup, so branch
        CMP     a3,#2                   // NOT the first 1-3 words
        LDR     a3,[a1]                 // to align the total to a multiple
        MVN     a3,a3                   // of 4 words
        STR     a3,[a1],#4
        BLT     not_loop_up_l1          // better to branch than skip instrs.
        LDRGE   a3,[a1]
        MVNGE   a3,a3
        STRGE   a3,[a1],#4
        LDRGT   a3,[a1]
        MVNGT   a3,a3
        STRGT   a3,[a1],#4
LABEL(not_loop_up_l1)
        BICS    a4,a2,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{lr}                // save work regs
LABEL(not_loop_up_l2)
        LDMIA   a1,{a2,a3,ip,lr}        // load 4 words in one go,NO writeback
        MVN     a2,a2                   // NOT the four words
        MVN     a3,a3
        MVN     ip,ip
        MVN     lr,lr
        STMIA   a1!,{a2,a3,ip,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     not_loop_up_l2          // if count still positive then loop
        LDMFD   sp!,{pc}^               // restore work regs and return

// extern void and_test_loop_up (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be AND_TESTed
//       exit
//               a1 = TRUE if any words ANDed together are non-zero else FALSE
//               a2 - a4, ip destroyed
        EXPORT(and_test_loop_up)        // word aligned and_test loop up
        DECLARE_FUNCTION(and_test_loop_up)
GLABEL(and_test_loop_up)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     and_test_loop_up_l1     // yup, so branch
        CMP     a4,#2
        LDR     a4,[a2],#4              // AND_TEST the first 1-3 words
        LDR     ip,[a1],#4              // to align the total to a multiple
        TST     ip,a4                   // of 4 words
        MOVNE   a1,#1                   // return TRUE if AND_TEST ok
        MOVNES  pc,lr
        BCC     and_test_loop_up_l1     // better to branch than skip instrs.
        LDRGE   a4,[a2],#4
        LDRGE   ip,[a1],#4
        TSTGE   ip,a4
        MOVNE   a1,#1
        MOVNES  pc,lr
        ANDS    a4,a3,#3
        CMP     a4,#2
        BLE     and_test_loop_up_l1     // better to branch than skip instrs.
        LDRGT   a4,[a2],#4
        LDRGT   ip,[a1],#4
        TSTGT   ip,a4
        MOVNE   a1,#1
        MOVNES  pc,lr
LABEL(and_test_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // return FALSE
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v6,lr}          // save work regs
        MOV     v6,a1                   // move xptr to v6
        MOV     a1,#1                   // set result to TRUE
LABEL(and_test_loop_up_l2)
        LDMIA   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMIA   v6!,{v3,v4,v5,lr}       // load target words
        TST     v3,a3                   // AND_TEST the four words
        TSTEQ   v4,v1
        TSTEQ   v5,v2
        TSTEQ   lr,ip
        LDMNEFD sp!,{v1-v6,pc}^
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     and_test_loop_up_l2     // if count still positive then loop
        MOV     a1,#0
        LDMFD   sp!,{v1-v6,pc}^         // restore work regs and return

#endif

// extern void compare_loop_up (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be COMPAREd
//       exit
//               a1 = +1 if first non-equal word in xptr[] and yptr[]
//                       xptr[i] > yptr[i]
//                    -1 if xptr[i] < yptr[i]
//                     0 otherwise
//               a2 - a4, ip destroyed
        EXPORT(compare_loop_up)         // word aligned compare loop up
        DECLARE_FUNCTION(compare_loop_up)
GLABEL(compare_loop_up)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     compare_loop_up_l1      // yup, so branch
        LDR     a4,[a2],#4              // COMPARE the first 1-3 words
        LDR     ip,[a1],#4              // to align the total to a multiple
        CMP     ip,a4                   // of 4 words
        MVNLO   a1,#0                   // x < y -> -1
        MOVHI   a1,#1                   // x > y -> +1
        MOVNES  pc,lr                   // and return result if not equal
        ANDS    a4,a3,#3
        CMP     a4,#2
        BLT     compare_loop_up_l1      // need to branch 'cos PSR used
        LDR     a4,[a2],#4
        LDR     ip,[a1],#4
        CMP     ip,a4
        MVNLO   a1,#0
        MOVHI   a1,#1
        MOVNES  pc,lr
        ANDS    a4,a3,#3
        CMP     a4,#2
        BLE     compare_loop_up_l1      // need to branch 'cos PSR used
        LDR     a4,[a2],#4
        LDR     ip,[a1],#4
        CMP     ip,a4
        MVNLO   a1,#0
        MOVHI   a1,#1
        MOVNES  pc,lr
LABEL(compare_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // xptr[] == yptr[] -> 0
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v6,lr}          // save work regs
        MOV     v6,a1                   // move xptr to v6
        MOV     a1,#1                   // set result to +1
LABEL(compare_loop_up_l2)
        LDMIA   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMIA   v6!,{v3,v4,v5,lr}       // load test words
        CMP     v3,a3                   // COMPARE the four words
        CMPEQ   v4,v1
        CMPEQ   v5,v2
        CMPEQ   lr,ip
        MVNLO   a1,#0                   // x < y -> -1 (a1 already holds +1)
        LDMNEFD sp!,{v1-v6,pc}^
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     compare_loop_up_l2      // if count still positive then loop
        MOV     a1,#0
        LDMFD   sp!,{v1-v6,pc}^         // restore work regs and return

#if CL_DS_BIG_ENDIAN_P

// extern uintD addto_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
//       entry
//               a1 = sourceptr
//               a2 = destptr
//               a3 = count of words to be added
//       exit
//               destptr[] = sourceptr[] + destptr[]
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(addto_loop_down)         // word aligned addto loop down
        DECLARE_FUNCTION(addto_loop_down)
GLABEL(addto_loop_down)
                MOV     a4,a3           // set regs for a call
                MOV     a3,a2           // to add_loop_down
                                        // and drop into add_loop_down

// extern uintD add_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
//       entry
//               a1 = sourceptr1
//               a2 = sourceptr2
//               a3 = destptr
//               a4 = count of words to be added
//       exit
//               destptr[] = sourceptr1[] + sourceptr2[]
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(add_loop_down)           // word aligned add loop down
        DECLARE_FUNCTION(add_loop_down)
GLABEL(add_loop_down)
        ANDS    ip,a4,#3                // multiple of 4 words ?
        BEQ     add_loop_down_l1        // yup, so branch
        STMFD   sp!,{v6,lr}
        LDR     v6,[a2,#-4]!            // add the first 1-3 words
        LDR     lr,[a1,#-4]!            // to align the total to a multiple
        ADDS    lr,lr,v6                // of 4 words
        STR     lr,[a3,#-4]!
        TEQ     ip,#1
        BEQ     add_loop_down_l0        // need to branch 'cos PSR used
        LDR     v6,[a2,#-4]!
        LDR     lr,[a1,#-4]!
        ADCS    lr,lr,v6
        STR     lr,[a3,#-4]!
        TEQ     ip,#2
        BEQ     add_loop_down_l0        // need to branch 'cos PSR used
        LDR     v6,[a2,#-4]!
        LDR     lr,[a1,#-4]!
        ADCS    lr,lr,v6
        STR     lr,[a3,#-4]!
LABEL(add_loop_down_l0)                 // at least one add has happened
        BICS    a4,a4,#3                // set counter to multiple of 4
        BNE     add_loop_down_l3        // branch if more adds to do
        ADCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        LDMEQFD sp!,{v6,pc}^            // and return
LABEL(add_loop_down_l1)
        BICS    a4,a4,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // no adds, so C = 0
        MOVEQS  pc,lr                   // if zero then we're done
        CMN     a4,#0                   // clear carry bit
        STMFD   sp!,{v6,lr}
LABEL(add_loop_down_l3)
        STMFD   sp!,{v1-v5}             // save work regs
LABEL(add_loop_down_l2)
        LDMDB   a2!,{v1,v2,v3,ip}       // load 4 words in one go
        LDMDB   a1!,{v4,v5,v6,lr}       // and from source2
        ADCS    lr,lr,ip                // add the four words with carry
        ADCS    v6,v6,v3
        ADCS    v5,v5,v2
        ADCS    v4,v4,v1
        STMDB   a3!,{v4,v5,v6,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4, preserve C
        TEQ     a4,#0                   // are we done ?
        BNE     add_loop_down_l2        // if count non-zero then loop
        ADC     a1,a4,a4                // set result to Carry (a4 is 0)
        LDMFD   sp!,{v1-v6,pc}^         // restore work regs and return

// extern uintD inc_loop_down (uintD* ptr, uintC count);
//       entry
//               a1 = ptr
//               a2 = count of words to be INCed
//       exit
//               a1 = 0 if any words are non-zero after increment else 1
//                      stop incrementing when first word becomes non-zero
//               a2 - a4, ip destroyed
        EXPORT(inc_loop_down)           // word aligned inc loop down
        DECLARE_FUNCTION(inc_loop_down)
GLABEL(inc_loop_down)
        ANDS    a3,a2,#1                // multiple of 2 words ?
        BEQ     inc_loop_down_l1        // yup, so branch
        LDR     a4,[a1,#-4]!            // INC the first word
        ADDS    a4,a4,#1                // align the total to a multiple of 2
        STR     a4,[a1]
        MOVNE   a1,#0                   // set result to 0
        MOVNES  pc,lr                   // return 0 if non-zero result
LABEL(inc_loop_down_l1)
        BICS    a4,a2,#1                // set counter to multiple of 2
        MOVEQ   a1,#1                   // return 1
        MOVEQS  pc,lr                   // if zero then we're done
        MOV     ip,a1                   // move ptr to ip
        MOV     a1,#0                   // set result to 0
        ANDS    a3,a4,#3
        BEQ     inc_loop_down_l3
        LDMDB   ip,{a2,a3}              // load 2 words in one go
        ADDS    a3,a3,#1                // INC the two words
        ADDEQS  a2,a2,#1                // stopping when first word non-zero
        STMDB   ip!,{a2,a3}             // store 2 results
        MOVNES  pc,lr                   // return 0 if any result non-zero
        SUBS    a4,a4,#2                // decrement counter by 2
        MOVEQ   a1,#1                   // if finished loop then
        MOVEQS  pc,lr                   // return 1
LABEL(inc_loop_down_l3)                 // now a multiple of 4 words
        STMFD   sp!,{v1,lr}             // save work regs
LABEL(inc_loop_down_l2)
        LDMDB   ip,{a2,a3,v1,lr}        // load 4 words in one go
        ADDS    lr,lr,#1                // INC the four words
        ADDEQS  v1,v1,#1                // stopping when first word non-zero
        ADDEQS  a3,a3,#1
        ADDEQS  a2,a2,#1
        STMDB   ip!,{a2,a3,v1,lr}       // store 4 results
        LDMNEFD sp!,{v1,pc}^            // return 0 if any result non-zero
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     inc_loop_down_l2        // if count still positive then loop
        MOV     a1,#1
        LDMFD   sp!,{v1,pc}^            // restore work regs and return 1

// extern uintD sub_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
//       entry
//               a1 = sourceptr1
//               a2 = sourceptr2
//               a3 = destptr
//               a4 = count of words to be subtracted
//       exit
//               destptr[] = sourceptr1[] -  sourceptr2[]
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(sub_loop_down)           // word aligned sub loop down
LABEL(sub_loop_down)
        ANDS    ip,a4,#3                // multiple of 4 words ?
        BEQ     sub_loop_down_l1        // yup, so branch
        STMFD   sp!,{v6,lr}
        LDR     v6,[a2,#-4]!            // subtract the first 1-3 words
        LDR     lr,[a1,#-4]!            // to align the total to a multiple
        SUBS    lr,lr,v6                // of 4 words
        STR     lr,[a3,#-4]!
        TEQ     ip,#1
        BNE     sub_loop_down_l0        // branch if more than one subtract
LABEL(sub_loop_down_l4)                 // drop through for better instr. timings
        BICS    a4,a4,#3                // set counter to multiple of 4
        SBCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        LDMEQFD sp!,{v6,pc}^            // and return
        STMFD   sp!,{v1-v5}             // save work regs
        B       sub_loop_down_l2        // branch if more subtracts to do
LABEL(sub_loop_down_l0)
        LDR     v6,[a2,#-4]!
        LDR     lr,[a1,#-4]!
        SBCS    lr,lr,v6
        STR     lr,[a3,#-4]!
        TEQ     ip,#2
        BEQ     sub_loop_down_l4        // need to branch 'cos PSR used
        LDR     v6,[a2,#-4]!
        LDR     lr,[a1,#-4]!
        SBCS    lr,lr,v6
        STR     lr,[a3,#-4]!
        B       sub_loop_down_l4
LABEL(sub_loop_down_l1)
        BICS    a4,a4,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // no subtracts, so C = 0
        MOVEQS  pc,lr                   // if zero then we're done
        CMP     a4,#0                   // set carry bit, since a4 > 0
        STMFD   sp!,{v1-v6,lr}          // save work regs
LABEL(sub_loop_down_l2)
        LDMDB   a2!,{v1,v2,v3,ip}       // load 4 words in one go
        LDMDB   a1!,{v4,v5,v6,lr}       // and from source2
        SBCS    lr,lr,ip                // subtract the four words with carry
        SBCS    v6,v6,v3
        SBCS    v5,v5,v2
        SBCS    v4,v4,v1
        STMDB   a3!,{v4,v5,v6,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4, preserve C
        TEQ     a4,#0                   // are we done ?
        BNE     sub_loop_down_l2        // if count non-zero then loop
        SBC     a1,a4,a4                // set result to Carry (a4 is 0)
        LDMFD   sp!,{v1-v6,pc}^         // restore work regs and return

// extern uintD subx_loop_down (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);
//       entry
//               a1 = sourceptr1
//               a2 = sourceptr2
//               a3 = destptr
//               a4 = count of words to be subtracted
//               [sp] = carry
//       exit
//               destptr[] = sourceptr1[] -  sourceptr2[]
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(subx_loop_down)          // word aligned xsub loop down
LABEL(subx_loop_down)
        LDR     ip,[sp]                 // get starting value of carry
LABEL(subx_loop_down_lsub)
        RSBS    ip,ip,#0                // set carry in PSR
        ANDS    ip,a4,#3                // multiple of 4 words ?
        BEQ     subx_loop_down_l1       // yup, so branch
        STMFD   sp!,{v6,lr}
        LDR     v6,[a2,#-4]!            // subtract the first 1-3 words
        LDR     lr,[a1,#-4]!            // to align the total to a multiple
        SBCS    lr,lr,v6                // of 4 words
        STR     lr,[a3,#-4]!
        TEQ     ip,#1
        BNE     subx_loop_down_l0       // branch if more than one subtract
LABEL(subx_loop_down_l4)                // drop through for better instr. timings
        BICS    a4,a4,#3                // set counter to multiple of 4
        SBCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        LDMEQFD sp!,{v6,pc}^            // and return
        STMFD   sp!,{v1-v5}             // save work regs
        B       subx_loop_down_l2       // branch if more subtracts to do
LABEL(subx_loop_down_l0)
        LDR     v6,[a2,#-4]!
        LDR     lr,[a1,#-4]!
        SBCS    lr,lr,v6
        STR     lr,[a3,#-4]!
        TEQ     ip,#2
        BEQ     subx_loop_down_l4       // need to branch 'cos PSR used
        LDR     v6,[a2,#-4]!
        LDR     lr,[a1,#-4]!
        SBCS    lr,lr,v6
        STR     lr,[a3,#-4]!
        B       subx_loop_down_l4
LABEL(subx_loop_down_l1)
        BICS    a4,a4,#3                // set counter to multiple of 4
        SBCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v6,lr}          // save work regs
LABEL(subx_loop_down_l2)
        LDMDB   a2!,{v1,v2,v3,ip}       // load 4 words in one go
        LDMDB   a1!,{v4,v5,v6,lr}       // and from source2
        SBCS    lr,lr,ip                // subtract the four words with carry
        SBCS    v6,v6,v3
        SBCS    v5,v5,v2
        SBCS    v4,v4,v1
        STMDB   a3!,{v4,v5,v6,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4, preserve C
        TEQ     a4,#0                   // are we done ?
        BNE     subx_loop_down_l2       // if count non-zero then loop
        SBC     a1,a4,a4                // set result to Carry (a4 is 0)
        LDMFD   sp!,{v1-v6,pc}^         // restore work regs and return

// extern uintD subfrom_loop_down (uintD* sourceptr, uintD* destptr, uintC count);
//       entry
//               a1 = sourceptr
//               a2 = destptr
//               a3 = count of words to be subtracted
//       exit
//               destptr[] = destptr[] - sourceptr[]
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(subfrom_loop_down)       // word aligned subfrom loop down
        DECLARE_FUNCTION(subfrom_loop_down)
GLABEL(subfrom_loop_down)
        ANDS    ip,a3,#3                // multiple of 4 words ?
        BEQ     subfrom_loop_down_l1    // yup, so branch
        STMFD   sp!,{lr}
        LDR     a4,[a1,#-4]!            // subtract the first 1-3 words
        LDR     lr,[a2,#-4]!            // to align the total to a multiple
        SUBS    lr,lr,a4                // of 4 words
        STR     lr,[a2]
        TEQ     ip,#1
        BNE     subfrom_loop_down_l0    // branch if more than one subtract
LABEL(subfrom_loop_down_l4)             // drop through for better instr. timings
        BICS    a4,a3,#3                // set counter to multiple of 4
        SBCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        LDMEQFD sp!,{pc}^               // and return
        STMFD   sp!,{v1-v5}             // save work regs
        B       subfrom_loop_down_l2    // branch if more subtracts to do
LABEL(subfrom_loop_down_l0)
        LDR     a4,[a1,#-4]!
        LDR     lr,[a2,#-4]!
        SBCS    lr,lr,a4
        STR     lr,[a2]
        TEQ     ip,#2
        BEQ     subfrom_loop_down_l4    // need to branch 'cos PSR used
        LDR     a4,[a1,#-4]!
        LDR     lr,[a2,#-4]!
        SBCS    lr,lr,a4
        STR     lr,[a2]
        B       subfrom_loop_down_l4
LABEL(subfrom_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // no subtracts, so C = 0
        MOVEQS  pc,lr                   // if zero then we're done
        CMP     a4,#0                   // set carry bit, since a4 > 0
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(subfrom_loop_down_l2)
        LDMDB   a1!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMDB   a2,{v3,v4,v5,lr}        // and from destptr
        SBCS    lr,lr,ip                // subtract the four words with carry
        SBCS    v5,v5,v2
        SBCS    v4,v4,v1
        SBCS    v3,v3,a3
        STMDB   a2!,{v3,v4,v5,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4, preserve C
        TEQ     a4,#0                   // are we done ?
        BNE     subfrom_loop_down_l2    // if count non-zero then loop
        SBC     a1,a4,a4                // set result to Carry (a4 is 0)
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern uintD dec_loop_down (uintD* ptr, uintC count);
//       entry
//               a1 = ptr
//               a2 = count of words to be DECed
//       exit
//               a1 = 0 if any words are non-zero before decrement else -1
//                      stop decrementing when first word is non-zero
//               a2 - a4, ip destroyed
        EXPORT(dec_loop_down)           // word aligned dec loop down
        DECLARE_FUNCTION(dec_loop_down)
GLABEL(dec_loop_down)
        ANDS    a3,a2,#1                // multiple of 2 words ?
        BEQ     dec_loop_down_l1        // yup, so branch
        LDR     a4,[a1,#-4]!            // DEC the first word
        SUBS    a4,a4,#1                // align the total to a multiple of 2
        STR     a4,[a1]
        MOVCS   a1,#0                   // set result to 0
        MOVCSS  pc,lr                   // return 0 if non-zero result
LABEL(dec_loop_down_l1)
        BICS    a4,a2,#1                // set counter to multiple of 2
        MVNEQ   a1,#0                   // return -1
        MOVEQS  pc,lr                   // if zero then we're done
        MOV     ip,a1                   // move ptr to ip
        MOV     a1,#0                   // set result to 0
        ANDS    a3,a4,#3
        BEQ     dec_loop_down_l3
        LDMDB   ip,{a2,a3}              // load 2 words in one go
        SUBS    a3,a3,#1                // DEC the two words
        SUBCCS  a2,a2,#1                // stopping when first word non-zero
        STMDB   ip!,{a2,a3}             // store 2 results
        MOVCSS  pc,lr                   // return 0 if any result non-zero
        SUBS    a4,a4,#2                // decrement counter by 2
        MVNEQ   a1,#0                   // if finished loop then
        MOVEQS  pc,lr                   // return -1
LABEL(dec_loop_down_l3)                 // now a multiple of 4 words
        STMFD   sp!,{v1,lr}             // save work regs
LABEL(dec_loop_down_l2)
        LDMDB   ip,{a2,a3,v1,lr}        // load 4 words in one go
        SUBS    lr,lr,#1                // DEC the four words
        SUBCCS  v1,v1,#1                // stopping when first word non-zero
        SUBCCS  a3,a3,#1
        SUBCCS  a2,a2,#1
        STMDB   ip!,{a2,a3,v1,lr}       // store 4 results
        LDMCSFD sp!,{v1,pc}^            // return 0 if any carry
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     dec_loop_down_l2        // if count still positive then loop
        MVN     a1,#0
        LDMFD   sp!,{v1,pc}^            // restore work regs and return -1

// extern void neg_loop_down (uintD* ptr, uintC count);
//       entry
//               a1 = ptr
//               a2 = count of words. The long integer is to be NEGated
//       exit
//               ptr[] = -ptr[] for count words
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(neg_loop_down)           // word aligned neg loop down
        DECLARE_FUNCTION(neg_loop_down)
GLABEL(neg_loop_down)
        CMPS    a2,#0                   // count = 0 ?
        MOVEQ   a1,#0                   // yup, so return 0
        MOVEQS  pc,lr
LABEL(neg_loop_down_l1)                 // skip all the zero words first
        LDR     a3,[a1,#-4]!            // compare words against zero
        CMPS    a3,#0                   // downwards in memory
        BNE     neg_loop_down_l2        // non-zero, so negate rest of words
        SUBS    a2,a2,#1                // reduce count of words
        BNE     neg_loop_down_l1        // more ?, so loop
        MOV     a1,#0                   // return 0
        MOVS    pc,lr
LABEL(neg_loop_down_l2)
        RSB     a3,a3,#0                // first non-zero word = -word
        STR     a3,[a1]
        SUBS    a2,a2,#1
        MVNEQ   a1,#0                   // done ? -> return -1
        MOVEQS  pc,lr
                                        // now NOT rest of the words
        ANDS    a3,a2,#3                // multiple of 4 words ?
        BEQ     neg_loop_down_l3        // yup, so branch
        CMP     a3,#2                   // NOT the first 1-3 words
        LDR     a3,[a1,#-4]!            // to align the total to a multiple
        MVN     a3,a3                   // of 4 words
        STR     a3,[a1]
        BLT     neg_loop_down_l3        // better to branch than skip instrs.
        LDRGE   a3,[a1,#-4]!
        MVNGE   a3,a3
        STRGE   a3,[a1]
        LDRGT   a3,[a1,#-4]!
        MVNGT   a3,a3
        STRGT   a3,[a1]
LABEL(neg_loop_down_l3)
        BICS    a4,a2,#3                // set counter to multiple of 4
        MVNEQ   a1,#0                   // set result to -1
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{lr}                // save work regs
LABEL(neg_loop_down_l4)
        LDMDB   a1,{a2,a3,ip,lr}        // load 4 words in one go,NO writeback
        MVN     a2,a2                   // NOT the four words
        MVN     a3,a3
        MVN     ip,ip
        MVN     lr,lr
        STMDB   a1!,{a2,a3,ip,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     neg_loop_down_l4        // if count still positive then loop
        MVN     a1,#0                   // set result to -1
        LDMFD   sp!,{pc}^               // restore work regs and return -1

// extern uintD shift1left_loop_down (uintD* ptr, uintC count);
//       entry
//               a1 = ptr
//               a2 = count of words to be shifted left
//       exit
//               a1 = carry out from last shift left
//               a2 - a4, ip destroyed
        EXPORT(shift1left_loop_down)    // word aligned shift1left loop down
        DECLARE_FUNCTION(shift1left_loop_down)
GLABEL(shift1left_loop_down)
        CMN     a1,#0                   // clear carry bit, since a1 > 0
        ANDS    a3,a2,#1                // multiple of 2 words ?
        BEQ     shift1left_loop_down_l1 // yup, so branch
        LDR     a4,[a1,#-4]!            // shift left the first word
        ADDS    a4,a4,a4
        STR     a4,[a1]
LABEL(shift1left_loop_down_l1)
        BICS    a4,a2,#1                // set counter to multiple of 2
        ADCEQ   a1,a4,a4                // if zero set result to C (a4 is 0)
        MOVEQS  pc,lr                   // and return
        ANDS    a3,a4,#3                // multiple of 4 words ?
        BEQ     shift1left_loop_down_l3 // yup, so branch
        LDMDB   a1,{a2,a3}              // load 2 words in one go
        ADCS    a3,a3,a3                // shift left the two words
        ADCS    a2,a2,a2
        STMDB   a1!,{a2,a3}             // store 2 results
        BICS    a4,a4,#2                // decrement counter by 2
        ADCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        MOVEQS  pc,lr                   // and return
LABEL(shift1left_loop_down_l3)          // now a multiple of 4 words
        STMFD   sp!,{lr}                // save work regs
LABEL(shift1left_loop_down_l2)
        LDMDB   a1,{a2,a3,ip,lr}        // load 4 words in one go
        ADCS    lr,lr,lr                // shift left the four words
        ADCS    ip,ip,ip
        ADCS    a3,a3,a3
        ADCS    a2,a2,a2
        STMDB   a1!,{a2,a3,ip,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4
        TEQ     a4,#0                   // are we done ?
        BNE     shift1left_loop_down_l2 // if count non-zero then loop
        ADC     a1,a4,a4                // set result to Carry (a4 is 0)
        LDMFD   sp!,{pc}^               // restore work regs and return 1

// extern uintD shiftleft_loop_down (uintD* ptr, uintC count, uintC i, uintD carry);
//       entry
//               a1 = ptr
//               a2 = count of words to be shifted left
//               a3 = size of left shift
//               a4 = value to ORR in for first shift
//       exit
//               a1 = shift out from last shift left
//               a2 - a4, ip destroyed
        EXPORT(shiftleft_loop_down)     // word aligned shiftleft loop down
        DECLARE_FUNCTION(shiftleft_loop_down)
GLABEL(shiftleft_loop_down)
        STMFD   sp!,{v6,lr}
        RSB     v6,a3,#32               // size of complementary right shift
        ANDS    ip,a2,#3                // multiple of 4 words ?
        BEQ     shiftleft_loop_down_l1  // yup, so branch
        LDR     lr,[a1,#-4]!            // shiftleft the first 1-3 words
        ORR     a4,a4,lr,ASL a3         // to align the total to a multiple
        STR     a4,[a1,#0]              // of 4 words
        MOV     a4,lr,LSR v6
        CMP     ip,#2
        BLT     shiftleft_loop_down_l1  // better to branch than skip instrs.
        LDRGE   lr,[a1,#-4]!
        ORRGE   a4,a4,lr,ASL a3
        STRGE   a4,[a1,#0]
        MOVGE   a4,lr,LSR v6
        LDRGT   lr,[a1,#-4]!
        ORRGT   a4,a4,lr,ASL a3
        STRGT   a4,[a1,#0]
        MOVGT   a4,lr,LSR v6
LABEL(shiftleft_loop_down_l1)
        BICS    ip,a2,#3                // set counter to multiple of 4
        MOVEQ   a1,a4                   // if zero then we're done
        LDMEQFD sp!,{v6,pc}^            // so return last shift out
        STMFD   sp!,{v1-v3}             // save work regs
LABEL(shiftleft_loop_down_l2)
        LDMDB   a1,{a2,v1,v2,v3}        // load 4 words in one go
        ORR     lr,a4,v3,ASL a3         // shiftleft the four words
        MOV     a4,v3,LSR v6            // keep carry in a4
        ORR     v3,a4,v2,ASL a3         // and store results up a register
        MOV     a4,v2,LSR v6            // to regs v1-v3,lr
        ORR     v2,a4,v1,ASL a3
        MOV     a4,v1,LSR v6
        ORR     v1,a4,a2,ASL a3
        MOV     a4,a2,LSR v6
        STMDB   a1!,{v1,v2,v3,lr}       // store 4 results
        SUBS    ip,ip,#4                // decrement counter by 4
        BGT     shiftleft_loop_down_l2  // if count still positive then loop
        MOV     a1,a4                   // result = last shift out
        LDMFD   sp!,{v1-v3,v6,pc}^      // restore work regs and return

// extern uintD shiftleftcopy_loop_down (uintD* sourceptr, uintD* destptr, uintC count, uintC i);
//       entry
//               a1 = sourceptr
//               a2 = destptr
//               a3 = count of words to be shifted left
//               a4 = size of left shift
//       exit
//               a1 = shift out from last shift left
//               a2 - a4, ip destroyed
        EXPORT(shiftleftcopy_loop_down) // word aligned shiftleftcopy loop down
        DECLARE_FUNCTION(shiftleftcopy_loop_down)
GLABEL(shiftleftcopy_loop_down)
        STMFD   sp!,{v5,v6,lr}
        MOV     v5,#0                   // initial shift carry
        RSB     v6,a4,#32               // size of complementary right shift
        ANDS    ip,a3,#3                // multiple of 4 words ?
        BEQ     shiftleftcopy_loop_down_l1      // yup, so branch
        LDR     lr,[a1,#-4]!            // shiftleft the first 1-3 words
        ORR     v5,v5,lr,ASL a4         // to align the total to a multiple
        STR     v5,[a2,#-4]!            // of 4 words
        MOV     v5,lr,LSR v6
        CMP     ip,#2
        BLT     shiftleftcopy_loop_down_l1      // better to branch than skip instrs.
        LDRGE   lr,[a1,#-4]!
        ORRGE   v5,v5,lr,ASL a4
        STRGE   v5,[a2,#-4]!
        MOVGE   v5,lr,LSR v6
        LDRGT   lr,[a1,#-4]!
        ORRGT   v5,v5,lr,ASL a4
        STRGT   v5,[a2,#-4]!
        MOVGT   v5,lr,LSR v6
LABEL(shiftleftcopy_loop_down_l1)
        BICS    ip,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,v5                   // if zero then we're done
        LDMEQFD sp!,{v5,v6,pc}^         // so return last shift out
        STMFD   sp!,{v1-v3}             // save work regs
LABEL(shiftleftcopy_loop_down_l2)
        LDMDB   a1!,{a3,v1,v2,v3}       // load 4 words in one go
        ORR     lr,v5,v3,ASL a4         // shiftleft the four words
        MOV     v5,v3,LSR v6            // keep carry in v5
        ORR     v3,v5,v2,ASL a4         // and store results up a register
        MOV     v5,v2,LSR v6            // to regs v1-v3,lr
        ORR     v2,v5,v1,ASL a4
        MOV     v5,v1,LSR v6
        ORR     v1,v5,a3,ASL a4
        MOV     v5,a3,LSR v6
        STMDB   a2!,{v1,v2,v3,lr}       // store 4 results
        SUBS    ip,ip,#4                // decrement counter by 4
        BGT     shiftleftcopy_loop_down_l2      // if count still positive then loop
        MOV     a1,v5                   // result = last shift out
        LDMFD   sp!,{v1-v3,v5,v6,pc}^   // restore work regs and return

// extern uintD shift1right_loop_up (uintD* ptr, uintC count, uintD carry);
//       entry
//               a1 = ptr
//               a2 = count of words to be shifted right
//               a3 = carry
//       exit
//               a1 = carry out from last shift right
//               a2 - a4, ip destroyed
        EXPORT(shift1right_loop_up)     // word aligned shift1right loop up
        DECLARE_FUNCTION(shift1right_loop_up)
GLABEL(shift1right_loop_up)
        MOVS    a3,a3,LSR #1            // set carry
        ANDS    a3,a2,#1                // multiple of 2 words ?
        BEQ     shift1right_loop_up_l1  // yup, so branch
        LDR     a4,[a1]                 // shift right the first word
        MOVS    a4,a4,RRX
        STR     a4,[a1],#4
LABEL(shift1right_loop_up_l1)
        BICS    a4,a2,#1                // set counter to multiple of 2
        MOVEQ   a1,a4,RRX               // if zero set result to C (a4 is 0)
        MOVEQS  pc,lr                   // and return
        ANDS    a3,a4,#3                // multiple of 4 words ?
        BEQ     shift1right_loop_up_l3  // yup, so branch
        LDMIA   a1,{a2,a3}              // load 2 words in one go
        MOVS    a2,a2,RRX               // shift right the two words
        MOVS    a3,a3,RRX
        STMIA   a1!,{a2,a3}             // store 2 results
        BICS    a4,a4,#2                // decrement counter by 2
        ADCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        MOVEQS  pc,lr                   // and return
LABEL(shift1right_loop_up_l3)           // now a multiple of 4 words
        STMFD   sp!,{lr}                // save work regs
LABEL(shift1right_loop_up_l2)
        LDMIA   a1,{a2,a3,ip,lr}        // load 4 words in one go
        MOVS    a2,a2,RRX               // shift right the four words
        MOVS    a3,a3,RRX
        MOVS    ip,ip,RRX
        MOVS    lr,lr,RRX
        STMIA   a1!,{a2,a3,ip,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4
        TEQ     a4,#0                   // are we done ?
        BNE     shift1right_loop_up_l2  // if count non-zero then loop
        MOV     a1,a4,RRX               // set result to Carry (a4 is 0)
        LDMFD   sp!,{pc}^               // restore work regs and return 1

// extern uintD shiftright_loop_up (uintD* ptr, uintC count, uintC i);
//       entry
//               a1 = ptr
//               a2 = count of words to be shifted right
//               a3 = size of right shift
//       exit
//               a1 = shift out from last shift right
//               a2 - a4, ip destroyed
        EXPORT(shiftright_loop_up)      // word aligned shiftright loop up
        DECLARE_FUNCTION(shiftright_loop_up)
GLABEL(shiftright_loop_up)
        STMFD   sp!,{v6,lr}
        MOV     a4,#0                   // initial shift carry
        RSB     v6,a3,#32               // size of complementary left shift
LABEL(shiftright_loop_up_l0)
        ANDS    ip,a2,#3                // multiple of 4 words ?
        BEQ     shiftright_loop_up_l1   // yup, so branch
        LDR     lr,[a1]                 // shiftright the first 1-3 words
        ORR     a4,a4,lr,LSR a3         // to align the total to a multiple
        STR     a4,[a1],#4              // of 4 words
        MOV     a4,lr,ASL v6
        CMP     ip,#2
        BLT     shiftright_loop_up_l1   // better to branch than skip instrs.
        LDRGE   lr,[a1]
        ORRGE   a4,a4,lr,LSR a3
        STRGE   a4,[a1],#4
        MOVGE   a4,lr,ASL v6
        LDRGT   lr,[a1]
        ORRGT   a4,a4,lr,LSR a3
        STRGT   a4,[a1],#4
        MOVGT   a4,lr,ASL v6
LABEL(shiftright_loop_up_l1)
        BICS    ip,a2,#3                // set counter to multiple of 4
        MOVEQ   a1,a4                   // if zero then we're done
        LDMEQFD sp!,{v6,pc}^            // so return last shift out
        STMFD   sp!,{v1-v3}             // save work regs
LABEL(shiftright_loop_up_l2)
        LDMIA   a1,{v1,v2,v3,lr}        // load 4 words in one go
        ORR     a2,a4,v1,LSR a3         // shiftright the four words
        MOV     a4,v1,ASL v6            // keep carry in a4
        ORR     v1,a4,v2,LSR a3         // and store results down a register
        MOV     a4,v2,ASL v6            // to regs a2,v1-v3
        ORR     v2,a4,v3,LSR a3
        MOV     a4,v3,ASL v6
        ORR     v3,a4,lr,LSR a3
        MOV     a4,lr,ASL v6
        STMIA   a1!,{a2,v1,v2,v3}       // store 4 results
        SUBS    ip,ip,#4                // decrement counter by 4
        BGT     shiftright_loop_up_l2   // if count still positive then loop
        MOV     a1,a4                   // result = last shift out
        LDMFD   sp!,{v1-v3,v6,pc}^      // restore work regs and return

// extern uintD shiftrightsigned_loop_up (uintD* ptr, uintC count, uintC i);
//       entry
//               a1 = ptr
//               a2 = count of words to be shifted right signed
//               a3 = size of right shift
//       exit
//               a1 = shift out from last shift right
//               a2 - a4, ip destroyed
        EXPORT(shiftrightsigned_loop_up)// word aligned shiftrightsigned loop up
        DECLARE_FUNCTION(shiftrightsigned_loop_up)
GLABEL(shiftrightsigned_loop_up)
        STMFD   sp!,{v6,lr}
        RSB     v6,a3,#32               // size of complementary left shift
        LDR     lr,[a1]                 // setup carry for first shift.
        MOV     a4,lr,ASR #31           // this is the sign extended bits
        AND     a4,a4,a4,LSL v6         // 31->(32-i) of the first word
        B       shiftright_loop_up_l0   // use right shift code now

// extern uintD shiftrightcopy_loop_up (uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry);
//       entry
//               a1 = sourceptr
//               a2 = destptr
//               a3 = count of words to be shifted right
//               a4 = size of right shift
//               [sp] = carry for first shift
//       exit
//               a1 = shift out from last shift right
//               a2 - a4, ip destroyed
        EXPORT(shiftrightcopy_loop_up)  // word aligned shiftrightcopy loop up
        DECLARE_FUNCTION(shiftrightcopy_loop_up)
GLABEL(shiftrightcopy_loop_up)
        STMFD   sp!,{v5,v6,lr}
        LDR     v5,[sp,#12]             // initial shift carry
        RSB     v6,a4,#32               // size of complementary left shift
        MOV     v5,v5,ASL v6
LABEL(shiftrightcopy_loop_up_l0)
        ANDS    ip,a3,#3                // multiple of 4 words ?
        BEQ     shiftrightcopy_loop_up_l1       // yup, so branch
        LDR     lr,[a1],#4              // shiftright the first 1-3 words
        ORR     v5,v5,lr,LSR a4         // to align the total to a multiple
        STR     v5,[a2],#4              // of 4 words
        MOV     v5,lr,ASL v6
        CMP     ip,#2
        BLT     shiftrightcopy_loop_up_l1       // better to branch than skip instrs.
        LDRGE   lr,[a1],#4
        ORRGE   v5,v5,lr,LSR a4
        STRGE   v5,[a2],#4
        MOVGE   v5,lr,ASL v6
        LDRGT   lr,[a1],#4
        ORRGT   v5,v5,lr,LSR a4
        STRGT   v5,[a2],#4
        MOVGT   v5,lr,ASL v6
LABEL(shiftrightcopy_loop_up_l1)
        BICS    ip,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,v5                   // if zero then we're done
        LDMEQFD sp!,{v5,v6,pc}^         // so return last shift out
        STMFD   sp!,{v1-v3}             // save work regs
LABEL(shiftrightcopy_loop_up_l2)
        LDMIA   a1!,{v1,v2,v3,lr}       // load 4 words in one go
        ORR     a3,v5,v1,LSR a4         // shiftright the four words
        MOV     v5,v1,ASL v6            // keep carry in v5
        ORR     v1,v5,v2,LSR a4         // and store results down a register
        MOV     v5,v2,ASL v6            // to regs a2,v1-v3
        ORR     v2,v5,v3,LSR a4
        MOV     v5,v3,ASL v6
        ORR     v3,v5,lr,LSR a4
        MOV     v5,lr,ASL v6
        STMIA   a2!,{a3,v1,v2,v3}       // store 4 results
        SUBS    ip,ip,#4                // decrement counter by 4
        BGT     shiftrightcopy_loop_up_l2       // if count still positive then loop
        MOV     a1,v5                   // result = last shift out
        LDMFD   sp!,{v1-v3,v5,v6,pc}^   // restore work regs and return

#ifndef HAVE_umull
// mulu32_64_vregs
//       entry
//               a1 = x
//               ip = y
//       exit
//               v1 = low32(x*y)
//               ip = high32(x*y)
//               v2,v3,v4 destroyed
LABEL(mulu32_64_vregs)
        MOV     v1,a1,LSR #16           // temp := top half of x
        MOV     v2,ip,LSR #16           // hi := top half of y
        BIC     v3,a1,v1,LSL #16        // x  := bottom half of x
        BIC     ip,ip,v2,LSL #16        // y  := bottom half of y
        MUL     v4,v3,ip                // low section of result
        MUL     ip,v1,ip                // ) middle sections
        MUL     v3,v2,v3                // )   of result
        MUL     v2,v1,v2                // high section of result
        ADDS    ip,ip,v3                // add middle sections
                                        // (can't use mla as we need carry)
        ADDCS   v2,v2,#0x10000          // carry from above add
        ADDS    v1,v4,ip,LSL #16        // x is now bottom 32 bits of result
        ADC     ip,v2,ip,LSR #16        // hi is top 32 bits
        MOVS    pc,lr
#endif

// extern uintD mulusmall_loop_down (uintD digit, uintD* ptr, uintC len, uintD newdigit);
//       entry
//               a1 = digit
//               a2 = ptr
//               a3 = count of words to be multiplied down
//               a4 = new digit = carry
//       exit
//               a1 = final carry of multiply
//               a2 - a4, ip destroyed
        EXPORT(mulusmall_loop_down)
        DECLARE_FUNCTION(mulusmall_loop_down)
GLABEL(mulusmall_loop_down)
        CMP     a3,#0
        MOVEQ   a1,a4
        MOVEQS  pc,lr
#ifdef HAVE_umull
        STMFD   sp!,{v1,lr}
LABEL(mulusmall_loop_down_l1)
        LDR     ip,[a2,#-4]!
        UMULL   v1,ip,a1,ip             // muluD(digit,*--ptr,hi=,lo=)
        ADDS    v1,v1,a4                // lo += carry
        ADC     a4,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a2,#0]              // *ptr = lo
        SUBS    a3,a3,#1                // len--
        BNE     mulusmall_loop_down_l1  // until len==0
        MOV     a1,a4                   // return carry
        LDMFD   sp!,{v1,pc}^
#else
        STMFD   sp!,{v1-v2,lr}
LABEL(mulusmall_loop_down_l1)
        LDR     ip,[a2,#-4]!

//      BL      mulu32_64_vregs         // muluD(digit,*--ptr,hi=,lo=)
// replaced by multiplication of a small x = a1 and a big y = ip :
        MOV     v1,ip,LSR #16           // top half of y
        BIC     ip,ip,v1,LSL #16        // bottom half of y
        MUL     v2,a1,v1                // middle section of result
        MUL     v1,a1,ip                // low section of result
        MOV     ip,#0                   // high section of result
        ADDS    v1,v1,v2,LSL #16        // bottom 32 bits of result
        ADC     ip,ip,v2,LSR #16        // top 32 bits of result

        ADDS    v1,v1,a4                // lo += carry
        ADC     a4,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a2,#0]              // *ptr = lo
        SUBS    a3,a3,#1                // len--
        BNE     mulusmall_loop_down_l1  // until len==0
        MOV     a1,a4                   // return carry
        LDMFD   sp!,{v1-v2,pc}^
#endif

// extern void mulu_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
//       entry
//               a1 = digit
//               a2 = sourceptr
//               a3 = destptr
//               a4 = count of words to be multiplied down
//       exit
//               a1 - a4, ip destroyed
        EXPORT(mulu_loop_down)
        DECLARE_FUNCTION(mulu_loop_down)
GLABEL(mulu_loop_down)
#ifdef HAVE_umull
        STMFD   sp!,{v1,v5,lr}
        MOV     v5,#0
LABEL(mulu_loop_down_l1)
        LDR     ip,[a2,#-4]!
        UMULL   v1,ip,a1,ip             // muluD(digit,*--sourceptr,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a3,#-4]!            // *--destptr = lo
        SUBS    a4,a4,#1                // len--
        BNE     mulu_loop_down_l1       // until len==0
        STR     v5,[a3,#-4]!            // *--destptr = carry
        LDMFD   sp!,{v1,v5,pc}^
#else
        STMFD   sp!,{v1-v5,lr}
        MOV     v5,#0
LABEL(mulu_loop_down_l1)
        LDR     ip,[a2,#-4]!
        BL      mulu32_64_vregs         // muluD(digit,*--sourceptr,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a3,#-4]!            // *--destptr = lo
        SUBS    a4,a4,#1                // len--
        BNE     mulu_loop_down_l1       // until len==0
        STR     v5,[a3,#-4]!            // *--destptr = carry
        LDMFD   sp!,{v1-v5,pc}^
#endif

// extern void muluadd_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
//       entry
//               a1 = digit
//               a2 = sourceptr
//               a3 = destptr
//               a4 = count of words to be multiplied added down
//       exit
//               a1 - a4, ip destroyed
        EXPORT(muluadd_loop_down)
        DECLARE_FUNCTION(muluadd_loop_down)
GLABEL(muluadd_loop_down)
#ifdef HAVE_umull
        STMFD   sp!,{v1,v5,lr}
        MOV     v5,#0
LABEL(muluadd_loop_down_l1)
        LDR     ip,[a2,#-4]!
        UMULL   v1,ip,a1,ip             // muluD(digit,*--sourceptr,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADCCS   ip,ip,#0                // if (lo<carry) { hi += 1 };
        LDR     v5,[a3,#-4]!            // carry = *--destptr
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a3,#0]              // *destptr = lo
        SUBS    a4,a4,#1                // len--
        BNE     muluadd_loop_down_l1    // until len==0
        MOV     a1,v5                   // return carry
        LDMFD   sp!,{v1,v5,pc}^
#else
        STMFD   sp!,{v1-v5,lr}
        MOV     v5,#0
LABEL(muluadd_loop_down_l1)
        LDR     ip,[a2,#-4]!
        BL      mulu32_64_vregs         // muluD(digit,*--sourceptr,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADCCS   ip,ip,#0                // if (lo<carry) { hi += 1 };
        LDR     v5,[a3,#-4]!            // carry = *--destptr
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a3,#0]              // *destptr = lo
        SUBS    a4,a4,#1                // len--
        BNE     muluadd_loop_down_l1    // until len==0
        MOV     a1,v5                   // return carry
        LDMFD   sp!,{v1-v5,pc}^
#endif

// extern void mulusub_loop_down (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
//       entry
//               a1 = digit
//               a2 = sourceptr
//               a3 = destptr
//               a4 = count of words to be multiplied subtracted down
//       exit
//               a1 - a4, ip destroyed
        EXPORT(mulusub_loop_down)
        DECLARE_FUNCTION(mulusub_loop_down)
GLABEL(mulusub_loop_down)
#ifdef HAVE_umull
        STMFD   sp!,{v1,v5,lr}
        MOV     v5,#0
LABEL(mulusub_loop_down_l1)
        LDR     ip,[a2,#-4]!
        UMULL   v1,ip,a1,ip             // muluD(digit,*--sourceptr,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 };
        LDR     ip,[a3,#-4]!            // carry = *--destptr
        SUBS    ip,ip,v1
        STR     ip,[a3,#0]              // *destptr = carry - lo
        ADDCC   v5,v5,#1                // if (carry<lo) { hi += 1 }; carry=hi
        SUBS    a4,a4,#1                // len--
        BNE     mulusub_loop_down_l1    // until len==0
        MOV     a1,v5                   // return carry
        LDMFD   sp!,{v1,v5,pc}^
#else
        STMFD   sp!,{v1-v5,lr}
        MOV     v5,#0
LABEL(mulusub_loop_down_l1)
        LDR     ip,[a2,#-4]!
        BL      mulu32_64_vregs         // muluD(digit,*--sourceptr,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 };
        LDR     ip,[a3,#-4]!            // carry = *--destptr
        SUBS    ip,ip,v1
        STR     ip,[a3,#0]              // *destptr = carry - lo
        ADDCC   v5,v5,#1                // if (carry<lo) { hi += 1 }; carry=hi
        SUBS    a4,a4,#1                // len--
        BNE     mulusub_loop_down_l1    // until len==0
        MOV     a1,v5                   // return carry
        LDMFD   sp!,{v1-v5,pc}^
#endif

#endif

#if !CL_DS_BIG_ENDIAN_P

// extern void or_loop_down (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be ORed
//       exit
//               xptr |= yptr for count words
//               a1 - a4, ip destroyed
        EXPORT(or_loop_down)            // word aligned or loop down
        DECLARE_FUNCTION(or_loop_down)
GLABEL(or_loop_down)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     or_loop_down_l1         // yup, so branch
        CMP     a4,#2                   // OR the first 1-3 words
        LDR     a4,[a2,#-4]!            // to align the total to a multiple
        LDR     ip,[a1,#-4]!            // of 4 words
        ORR     ip,ip,a4
        STR     ip,[a1]
        BLT     or_loop_down_l1         // better to branch than skip instrs.
        LDRGE   a4,[a2,#-4]!
        LDRGE   ip,[a1,#-4]!
        ORRGE   ip,ip,a4
        STRGE   ip,[a1]
        LDRGT   a4,[a2,#-4]!
        LDRGT   ip,[a1,#-4]!
        ORRGT   ip,ip,a4
        STRGT   ip,[a1]
LABEL(or_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(or_loop_down_l2)
        LDMDB   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMDB   a1,{v3,v4,v5,lr}        // load target words
        ORR     v3,v3,a3                // OR the four words
        ORR     v4,v4,v1
        ORR     v5,v5,v2
        ORR     lr,lr,ip
        STMDB   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     or_loop_down_l2         // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void xor_loop_down (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be XORed
//       exit
//               xptr ^= yptr for count words
//               a1 - a4, ip destroyed
        EXPORT(xor_loop_down)           // word aligned xor loop down
        DECLARE_FUNCTION(xor_loop_down)
GLABEL(xor_loop_down)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     xor_loop_down_l1        // yup, so branch
        CMP     a4,#2                   // XOR the first 1-3 words
        LDR     a4,[a2,#-4]!            // to align the total to a multiple
        LDR     ip,[a1,#-4]!            // of 4 words
        EOR     ip,ip,a4
        STR     ip,[a1]
        BLT     xor_loop_down_l1        // better to branch than skip instrs.
        LDRGE   a4,[a2,#-4]!
        LDRGE   ip,[a1,#-4]!
        EORGE   ip,ip,a4
        STRGE   ip,[a1]
        LDRGT   a4,[a2,#-4]!
        LDRGT   ip,[a1,#-4]!
        EORGT   ip,ip,a4
        STRGT   ip,[a1]
LABEL(xor_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(xor_loop_down_l2)
        LDMDB   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMDB   a1,{v3,v4,v5,lr}        // load target words
        EOR     v3,v3,a3                // XOR the four words
        EOR     v4,v4,v1
        EOR     v5,v5,v2
        EOR     lr,lr,ip
        STMDB   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     xor_loop_down_l2        // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void and_loop_down (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be ANDed
//       exit
//               xptr &= yptr for count words
//               a1 - a4, ip destroyed
        EXPORT(and_loop_down)           // word aligned and loop down
        DECLARE_FUNCTION(and_loop_down)
GLABEL(and_loop_down)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     and_loop_down_l1        // yup, so branch
        CMP     a4,#2                   // AND the first 1-3 words
        LDR     a4,[a2,#-4]!            // to align the total to a multiple
        LDR     ip,[a1,#-4]!            // of 4 words
        AND     ip,ip,a4
        STR     ip,[a1]
        BLT     and_loop_down_l1        // better to branch than skip instrs.
        LDRGE   a4,[a2,#-4]!
        LDRGE   ip,[a1,#-4]!
        ANDGE   ip,ip,a4
        STRGE   ip,[a1]
        LDRGT   a4,[a2,#-4]!
        LDRGT   ip,[a1,#-4]!
        ANDGT   ip,ip,a4
        STRGT   ip,[a1]
LABEL(and_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(and_loop_down_l2)
        LDMDB   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMDB   a1,{v3,v4,v5,lr}        // load target words
        AND     v3,v3,a3                // AND the four words
        AND     v4,v4,v1
        AND     v5,v5,v2
        AND     lr,lr,ip
        STMDB   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     and_loop_down_l2        // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void eqv_loop_down (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be XORed
//       exit
//               xptr = ~(xptr ^ yptr) for count words
//               a1 - a4, ip destroyed
        EXPORT(eqv_loop_down)           // word aligned eqv loop down
        DECLARE_FUNCTION(eqv_loop_down)
GLABEL(eqv_loop_down)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     eqv_loop_down_l1        // yup, so branch
        CMP     a4,#2                   // EQV the first 1-3 words
        LDR     a4,[a2,#-4]!            // to align the total to a multiple
        LDR     ip,[a1,#-4]!            // of 4 words
        EOR     ip,ip,a4
        MVN     ip,ip
        STR     ip,[a1]
        BLT     eqv_loop_down_l1        // better to branch than skip instrs.
        LDRGE   a4,[a2,#-4]!
        LDRGE   ip,[a1,#-4]!
        EORGE   ip,ip,a4
        MVNGE   ip,ip
        STRGE   ip,[a1]
        BLE     eqv_loop_down_l1        // better to branch than skip instrs.
        LDRGT   a4,[a2,#-4]!
        LDRGT   ip,[a1,#-4]!
        EORGT   ip,ip,a4
        MVNGT   ip,ip
        STRGT   ip,[a1]
LABEL(eqv_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(eqv_loop_down_l2)
        LDMDB   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMDB   a1,{v3,v4,v5,lr}        // load target words
        EOR     v3,v3,a3                // EVQ the four words
        MVN     v3,v3
        EOR     v4,v4,v1
        MVN     v4,v4
        EOR     v5,v5,v2
        MVN     v5,v5
        EOR     lr,lr,ip
        MVN     lr,lr
        STMDB   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     eqv_loop_down_l2        // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void nand_loop_down (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be NANDed
//       exit
//               xptr = ~(xptr & yptr) for count words
//               a1 - a4, ip destroyed
        EXPORT(nand_loop_down)          // word aligned nand loop down
        DECLARE_FUNCTION(nand_loop_down)
GLABEL(nand_loop_down)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     nand_loop_down_l1       // yup, so branch
        CMP     a4,#2                   // NAND the first 1-3 words
        LDR     a4,[a2,#-4]!            // to align the total to a multiple
        LDR     ip,[a1,#-4]!            // of 4 words
        AND     ip,ip,a4
        MVN     ip,ip
        STR     ip,[a1]
        BLT     nand_loop_down_l1       // better to branch than skip instrs.
        LDRGE   a4,[a2,#-4]!
        LDRGE   ip,[a1,#-4]!
        ANDGE   ip,ip,a4
        MVNGE   ip,ip
        STRGE   ip,[a1]
        BLE     nand_loop_down_l1       // better to branch than skip instrs.
        LDRGT   a4,[a2,#-4]!
        LDRGT   ip,[a1,#-4]!
        ANDGT   ip,ip,a4
        MVNGT   ip,ip
        STRGT   ip,[a1]
LABEL(nand_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(nand_loop_down_l2)
        LDMDB   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMDB   a1,{v3,v4,v5,lr}        // load target words
        AND     v3,v3,a3                // NAND the four words
        MVN     v3,v3
        AND     v4,v4,v1
        MVN     v4,v4
        AND     v5,v5,v2
        MVN     v5,v5
        AND     lr,lr,ip
        MVN     lr,lr
        STMDB   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     nand_loop_down_l2       // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void nor_loop_down (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be NORed
//       exit
//               xptr = ~(xptr | yptr) for count words
//               a1 - a4, ip destroyed
        EXPORT(nor_loop_down)           // word aligned nor loop down
        DECLARE_FUNCTION(nor_loop_down)
GLABEL(nor_loop_down)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     nor_loop_down_l1        // yup, so branch
        CMP     a4,#2                   // NOR the first 1-3 words
        LDR     a4,[a2,#-4]!            // to align the total to a multiple
        LDR     ip,[a1,#-4]!            // of 4 words
        ORR     ip,ip,a4
        MVN     ip,ip
        STR     ip,[a1]
        BLT     nor_loop_down_l1        // better to branch than skip instrs.
        LDRGE   a4,[a2,#-4]!
        LDRGE   ip,[a1,#-4]!
        ORRGE   ip,ip,a4
        MVNGE   ip,ip
        STRGE   ip,[a1]
        BLE     nor_loop_down_l1        // better to branch than skip instrs.
        LDRGT   a4,[a2,#-4]!
        LDRGT   ip,[a1,#-4]!
        ORRGT   ip,ip,a4
        MVNGT   ip,ip
        STRGT   ip,[a1]
LABEL(nor_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(nor_loop_down_l2)
        LDMDB   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMDB   a1,{v3,v4,v5,lr}        // load target words
        ORR     v3,v3,a3                // NOR the four words
        MVN     v3,v3
        ORR     v4,v4,v1
        MVN     v4,v4
        ORR     v5,v5,v2
        MVN     v5,v5
        ORR     lr,lr,ip
        MVN     lr,lr
        STMDB   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     nor_loop_down_l2        // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void andc2_loop_down (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be ANDC2ed
//       exit
//               xptr = xptr & ~yptr for count words
//               a1 - a4, ip destroyed
        EXPORT(andc2_loop_down)         // word aligned andc2 loop down
        DECLARE_FUNCTION(andc2_loop_down)
GLABEL(andc2_loop_down)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     andc2_loop_down_l1      // yup, so branch
        CMP     a4,#2                   // ANDC2 the first 1-3 words
        LDR     a4,[a2,#-4]!            // to align the total to a multiple
        LDR     ip,[a1,#-4]!            // of 4 words
        BIC     ip,ip,a4
        STR     ip,[a1]
        BLT     andc2_loop_down_l1      // better to branch than skip instrs.
        LDRGE   a4,[a2,#-4]!
        LDRGE   ip,[a1,#-4]!
        BICGE   ip,ip,a4
        STRGE   ip,[a1]
        LDRGT   a4,[a2,#-4]!
        LDRGT   ip,[a1,#-4]!
        BICGT   ip,ip,a4
        STRGT   ip,[a1]
LABEL(andc2_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(andc2_loop_down_l2)
        LDMDB   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMDB   a1,{v3,v4,v5,lr}        // load target words
        BIC     v3,v3,a3                // ANDC2 the four words
        BIC     v4,v4,v1
        BIC     v5,v5,v2
        BIC     lr,lr,ip
        STMDB   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     andc2_loop_down_l2      // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void orc2_loop_down (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be XORed
//       exit
//               xptr = xptr | ~yptr for count words
//               a1 - a4, ip destroyed
        EXPORT(orc2_loop_down)          // word aligned orc2 loop down
        DECLARE_FUNCTION(orc2_loop_down)
GLABEL(orc2_loop_down)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     orc2_loop_down_l1       // yup, so branch
        CMP     a4,#2                   // ORC2 the first 1-3 words
        LDR     a4,[a2,#-4]!            // to align the total to a multiple
        LDR     ip,[a1,#-4]!            // of 4 words
        MVN     a4,a4
        ORR     ip,ip,a4
        STR     ip,[a1]
        BLT     orc2_loop_down_l1       // better to branch than skip instrs.
        LDRGE   a4,[a2,#-4]!
        LDRGE   ip,[a1,#-4]!
        MVNGE   a4,a4
        ORRGE   ip,ip,a4
        STRGE   ip,[a1]
        BLE     orc2_loop_down_l1       // better to branch than skip instrs.
        LDRGT   a4,[a2,#-4]!
        LDRGT   ip,[a1,#-4]!
        MVNGT   a4,a4
        ORRGT   ip,ip,a4
        STRGT   ip,[a1]
LABEL(orc2_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(orc2_loop_down_l2)
        LDMDB   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMDB   a1,{v3,v4,v5,lr}        // load target words
        MVN     a3,a3                   // ORC2 the four words
        ORR     v3,v3,a3
        MVN     v1,v1
        ORR     v4,v4,v1
        MVN     v2,v2
        ORR     v5,v5,v2
        MVN     ip,ip
        ORR     lr,lr,ip
        STMDB   a1!,{v3,v4,v5,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     orc2_loop_down_l2       // if count still positive then loop
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern void not_loop_down (uintD* xptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = count of words to be NOTed
//       exit
//               xptr = ~xptr for count words
//               a1 - a4, ip destroyed
        EXPORT(not_loop_down)           // word aligned not loop down
        DECLARE_FUNCTION(not_loop_down)
GLABEL(not_loop_down)
        ANDS    a3,a2,#3                // multiple of 4 words ?
        BEQ     not_loop_down_l1        // yup, so branch
        CMP     a3,#2                   // NOT the first 1-3 words
        LDR     a3,[a1,#-4]!            // to align the total to a multiple
        MVN     a3,a3                   // of 4 words
        STR     a3,[a1]
        BLT     not_loop_down_l1        // better to branch than skip instrs.
        LDRGE   a3,[a1,#-4]!
        MVNGE   a3,a3
        STRGE   a3,[a1]
        LDRGT   a3,[a1,#-4]!
        MVNGT   a3,a3
        STRGT   a3,[a1]
LABEL(not_loop_down_l1)
        BICS    a4,a2,#3                // set counter to multiple of 4
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{lr}                // save work regs
LABEL(not_loop_down_l2)
        LDMDB   a1,{a2,a3,ip,lr}        // load 4 words in one go,NO writeback
        MVN     a2,a2                   // NOT the four words
        MVN     a3,a3
        MVN     ip,ip
        MVN     lr,lr
        STMDB   a1!,{a2,a3,ip,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     not_loop_down_l2        // if count still positive then loop
        LDMFD   sp!,{pc}^               // restore work regs and return

// extern void and_test_loop_down (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be AND_TESTed
//       exit
//               a1 = TRUE if any words ANDed together are non-zero else FALSE
//               a2 - a4, ip destroyed
        EXPORT(and_test_loop_down)      // word aligned and_test loop down
        DECLARE_FUNCTION(and_test_loop_down)
GLABEL(and_test_loop_down)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     and_test_loop_down_l1   // yup, so branch
        CMP     a4,#2
        LDR     a4,[a2,#-4]!            // AND_TEST the first 1-3 words
        LDR     ip,[a1,#-4]!            // to align the total to a multiple
        TST     ip,a4                   // of 4 words
        MOVNE   a1,#1                   // return TRUE if AND_TEST ok
        MOVNES  pc,lr
        BCC     and_test_loop_down_l1   // better to branch than skip instrs.
        LDRGE   a4,[a2,#-4]!
        LDRGE   ip,[a1,#-4]!
        TSTGE   ip,a4
        MOVNE   a1,#1
        MOVNES  pc,lr
        ANDS    a4,a3,#3
        CMP     a4,#2
        BLE     and_test_loop_down_l1   // better to branch than skip instrs.
        LDRGT   a4,[a2,#-4]!
        LDRGT   ip,[a1,#-4]!
        TSTGT   ip,a4
        MOVNE   a1,#1
        MOVNES  pc,lr
LABEL(and_test_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // return FALSE
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v6,lr}          // save work regs
        MOV     v6,a1                   // move xptr to v6
        MOV     a1,#1                   // set result to TRUE
LABEL(and_test_loop_down_l2)
        LDMDB   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMDB   v6!,{v3,v4,v5,lr}       // load target words
        TST     v3,a3                   // AND_TEST the four words
        TSTEQ   v4,v1
        TSTEQ   v5,v2
        TSTEQ   lr,ip
        LDMNEFD sp!,{v1-v6,pc}^
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     and_test_loop_down_l2   // if count still positive then loop
        MOV     a1,#0
        LDMFD   sp!,{v1-v6,pc}^         // restore work regs and return

// extern void compare_loop_down (uintD* xptr, uintD* yptr, uintC count);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be COMPAREd
//       exit
//               a1 = +1 if first non-equal word in xptr[] and yptr[]
//                       xptr[i] > yptr[i]
//                    -1 if xptr[i] < yptr[i]
//                     0 otherwise
//               a2 - a4, ip destroyed
        EXPORT(compare_loop_down)       // word aligned compare loop down
        DECLARE_FUNCTION(compare_loop_down)
GLABEL(compare_loop_down)
        ANDS    a4,a3,#3                // multiple of 4 words ?
        BEQ     compare_loop_down_l1    // yup, so branch
        LDR     a4,[a2,#-4]!            // COMPARE the first 1-3 words
        LDR     ip,[a1,#-4]!            // to align the total to a multiple
        CMP     ip,a4                   // of 4 words
        MVNLO   a1,#0                   // x < y -> -1
        MOVHI   a1,#1                   // x > y -> +1
        MOVNES  pc,lr                   // and return result if not equal
        ANDS    a4,a3,#3
        CMP     a4,#2
        BLT     compare_loop_down_l1    // need to branch 'cos PSR used
        LDR     a4,[a2,#-4]!
        LDR     ip,[a1,#-4]!
        CMP     ip,a4
        MVNLO   a1,#0
        MOVHI   a1,#1
        MOVNES  pc,lr
        ANDS    a4,a3,#3
        CMP     a4,#2
        BLE     compare_loop_down_l1    // need to branch 'cos PSR used
        LDR     a4,[a2,#-4]!
        LDR     ip,[a1,#-4]!
        CMP     ip,a4
        MVNLO   a1,#0
        MOVHI   a1,#1
        MOVNES  pc,lr
LABEL(compare_loop_down_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // xptr[] == yptr[] -> 0
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v6,lr}          // save work regs
        MOV     v6,a1                   // move xptr to v6
        MOV     a1,#1                   // set result to +1
LABEL(compare_loop_down_l2)
        LDMDB   a2!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMDB   v6!,{v3,v4,v5,lr}       // load test words
        CMP     lr,ip                   // COMPARE the four words
        CMPEQ   v5,v2
        CMPEQ   v4,v1
        CMPEQ   v3,a3
        MVNLO   a1,#0                   // x < y -> -1 (a1 already holds +1)
        LDMNEFD sp!,{v1-v6,pc}^
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     compare_loop_down_l2    // if count still positive then loop
        MOV     a1,#0
        LDMFD   sp!,{v1-v6,pc}^         // restore work regs and return

// extern uintD addto_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
//       entry
//               a1 = sourceptr
//               a2 = destptr
//               a3 = count of words to be added
//       exit
//               destptr[] = sourceptr[] + destptr[]
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(addto_loop_up)           // word aligned addto loop up
        DECLARE_FUNCTION(addto_loop_up)
GLABEL(addto_loop_up)
                MOV     a4,a3           // set regs for a call
                MOV     a3,a2           // to add_loop_up
                                        // and drop into add_loop_up

// extern uintD add_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
//       entry
//               a1 = sourceptr1
//               a2 = sourceptr2
//               a3 = destptr
//               a4 = count of words to be added
//       exit
//               destptr[] = sourceptr1[] + sourceptr2[]
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(add_loop_up)             // word aligned add loop up
        DECLARE_FUNCTION(add_loop_up)
GLABEL(add_loop_up)
        ANDS    ip,a4,#3                // multiple of 4 words ?
        BEQ     add_loop_up_l1          // yup, so branch
        STMFD   sp!,{v6,lr}
        LDR     v6,[a2],#4              // add the first 1-3 words
        LDR     lr,[a1],#4              // to align the total to a multiple
        ADDS    lr,lr,v6                // of 4 words
        STR     lr,[a3],#4
        TEQ     ip,#1
        BEQ     add_loop_up_l0          // need to branch 'cos PSR used
        LDR     v6,[a2],#4
        LDR     lr,[a1],#4
        ADCS    lr,lr,v6
        STR     lr,[a3],#4
        TEQ     ip,#2
        BEQ     add_loop_up_l0          // need to branch 'cos PSR used
        LDR     v6,[a2],#4
        LDR     lr,[a1],#4
        ADCS    lr,lr,v6
        STR     lr,[a3],#4
LABEL(add_loop_up_l0)                   // at least one add has happened
        BICS    a4,a4,#3                // set counter to multiple of 4
        BNE     add_loop_up_l3          // branch if more adds to do
        ADCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        LDMEQFD sp!,{v6,pc}^            // and return
LABEL(add_loop_up_l1)
        BICS    a4,a4,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // no adds, so C = 0
        MOVEQS  pc,lr                   // if zero then we're done
        CMN     a4,#0                   // clear carry bit
        STMFD   sp!,{v6,lr}
LABEL(add_loop_up_l3)
        STMFD   sp!,{v1-v5}             // save work regs
LABEL(add_loop_up_l2)
        LDMIA   a2!,{v1,v2,v3,ip}       // load 4 words in one go
        LDMIA   a1!,{v4,v5,v6,lr}       // and from source2
        ADCS    v4,v4,v1                // add the four words with carry
        ADCS    v5,v5,v2
        ADCS    v6,v6,v3
        ADCS    lr,lr,ip
        STMIA   a3!,{v4,v5,v6,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4, preserve C
        TEQ     a4,#0                   // are we done ?
        BNE     add_loop_up_l2          // if count non-zero then loop
        ADC     a1,a4,a4                // set result to Carry (a4 is 0)
        LDMFD   sp!,{v1-v6,pc}^         // restore work regs and return

// extern uintD inc_loop_up (uintD* ptr, uintC count);
//       entry
//               a1 = ptr
//               a2 = count of words to be INCed
//       exit
//               a1 = 0 if any words are non-zero after increment else 1
//                      stop incrementing when first word becomes non-zero
//               a2 - a4, ip destroyed
        EXPORT(inc_loop_up)             // word aligned inc loop up
        DECLARE_FUNCTION(inc_loop_up)
GLABEL(inc_loop_up)
        ANDS    a3,a2,#1                // multiple of 2 words ?
        BEQ     inc_loop_up_l1          // yup, so branch
        LDR     a4,[a1]                 // INC the first word
        ADDS    a4,a4,#1                // align the total to a multiple of 2
        STR     a4,[a1],#4
        MOVNE   a1,#0                   // set result to 0
        MOVNES  pc,lr                   // return 0 if non-zero result
LABEL(inc_loop_up_l1)
        BICS    a4,a2,#1                // set counter to multiple of 2
        MOVEQ   a1,#1                   // return 1
        MOVEQS  pc,lr                   // if zero then we're done
        MOV     ip,a1                   // move ptr to ip
        MOV     a1,#0                   // set result to 0
        ANDS    a3,a4,#3
        BEQ     inc_loop_up_l3
        LDMIA   ip,{a2,a3}              // load 2 words in one go
        ADDS    a2,a2,#1                // INC the two words
        ADDEQS  a3,a3,#1                // stopping when first word non-zero
        STMIA   ip!,{a2,a3}             // store 2 results
        MOVNES  pc,lr                   // return 0 if any result non-zero
        SUBS    a4,a4,#2                // decrement counter by 2
        MOVEQ   a1,#1                   // if finished loop then
        MOVEQS  pc,lr                   // return 1
LABEL(inc_loop_up_l3)                   // now a multiple of 4 words
        STMFD   sp!,{v1,lr}             // save work regs
LABEL(inc_loop_up_l2)
        LDMIA   ip,{a2,a3,v1,lr}        // load 4 words in one go
        ADDS    a2,a2,#1                // INC the four words
        ADDEQS  a3,a3,#1                // stopping when first word non-zero
        ADDEQS  v1,v1,#1
        ADDEQS  lr,lr,#1
        STMIA   ip!,{a2,a3,v1,lr}       // store 4 results
        LDMNEFD sp!,{v1,pc}^            // return 0 if any result non-zero
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     inc_loop_up_l2          // if count still positive then loop
        MOV     a1,#1
        LDMFD   sp!,{v1,pc}^            // restore work regs and return 1

// extern uintD sub_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count);
//       entry
//               a1 = sourceptr1
//               a2 = sourceptr2
//               a3 = destptr
//               a4 = count of words to be subtracted
//       exit
//               destptr[] = sourceptr1[] -  sourceptr2[]
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(sub_loop_up)             // word aligned sub loop up
        DECLARE_FUNCTION(sub_loop_up)
GLABEL(sub_loop_up)
        ANDS    ip,a4,#3                // multiple of 4 words ?
        BEQ     sub_loop_up_l1          // yup, so branch
        STMFD   sp!,{v6,lr}
        LDR     v6,[a2],#4              // subtract the first 1-3 words
        LDR     lr,[a1],#4              // to align the total to a multiple
        SUBS    lr,lr,v6                // of 4 words
        STR     lr,[a3],#4
        TEQ     ip,#1
        BNE     sub_loop_up_l0          // branch if more than one subtract
LABEL(sub_loop_up_l4)                   // drop through for better instr. timings
        BICS    a4,a4,#3                // set counter to multiple of 4
        SBCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        LDMEQFD sp!,{v6,pc}^            // and return
        STMFD   sp!,{v1-v5}             // save work regs
        B       sub_loop_up_l2          // branch if more subtracts to do
LABEL(sub_loop_up_l0)
        LDR     v6,[a2],#4
        LDR     lr,[a1],#4
        SBCS    lr,lr,v6
        STR     lr,[a3],#4
        TEQ     ip,#2
        BEQ     sub_loop_up_l4          // need to branch 'cos PSR used
        LDR     v6,[a2],#4
        LDR     lr,[a1],#4
        SBCS    lr,lr,v6
        STR     lr,[a3],#4
        B       sub_loop_up_l4
LABEL(sub_loop_up_l1)
        BICS    a4,a4,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // no subtracts, so C = 0
        MOVEQS  pc,lr                   // if zero then we're done
        CMP     a4,#0                   // set carry bit, since a4 > 0
        STMFD   sp!,{v1-v6,lr}          // save work regs
LABEL(sub_loop_up_l2)
        LDMIA   a2!,{v1,v2,v3,ip}       // load 4 words in one go
        LDMIA   a1!,{v4,v5,v6,lr}       // and from source2
        SBCS    v4,v4,v1                // subtract the four words with carry
        SBCS    v5,v5,v2
        SBCS    v6,v6,v3
        SBCS    lr,lr,ip
        STMIA   a3!,{v4,v5,v6,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4, preserve C
        TEQ     a4,#0                   // are we done ?
        BNE     sub_loop_up_l2          // if count non-zero then loop
        SBC     a1,a4,a4                // set result to Carry (a4 is 0)
        LDMFD   sp!,{v1-v6,pc}^         // restore work regs and return

// extern uintD subx_loop_up (uintD* sourceptr1, uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);
//       entry
//               a1 = sourceptr1
//               a2 = sourceptr2
//               a3 = destptr
//               a4 = count of words to be subtracted
//               [sp] = carry
//       exit
//               destptr[] = sourceptr1[] -  sourceptr2[]
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(subx_loop_up)            // word aligned xsub loop up
        DECLARE_FUNCTION(subx_loop_up)
GLABEL(subx_loop_up)
        LDR     ip,[sp]                 // get starting value of carry
LABEL(subx_loop_up_lsub)
        RSBS    ip,ip,#0                // set carry in PSR
        ANDS    ip,a4,#3                // multiple of 4 words ?
        BEQ     subx_loop_up_l1         // yup, so branch
        STMFD   sp!,{v6,lr}
        LDR     v6,[a2],#4              // subtract the first 1-3 words
        LDR     lr,[a1],#4              // to align the total to a multiple
        SBCS    lr,lr,v6                // of 4 words
        STR     lr,[a3],#4
        TEQ     ip,#1
        BNE     subx_loop_up_l0         // branch if more than one subtract
LABEL(subx_loop_up_l4)                  // drop through for better instr. timings
        BICS    a4,a4,#3                // set counter to multiple of 4
        SBCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        LDMEQFD sp!,{v6,pc}^            // and return
        STMFD   sp!,{v1-v5}             // save work regs
        B       subx_loop_up_l2         // branch if more subtracts to do
LABEL(subx_loop_up_l0)
        LDR     v6,[a2],#4
        LDR     lr,[a1],#4
        SBCS    lr,lr,v6
        STR     lr,[a3],#4
        TEQ     ip,#2
        BEQ     subx_loop_up_l4         // need to branch 'cos PSR used
        LDR     v6,[a2],#4
        LDR     lr,[a1],#4
        SBCS    lr,lr,v6
        STR     lr,[a3],#4
        B       subx_loop_up_l4
LABEL(subx_loop_up_l1)
        BICS    a4,a4,#3                // set counter to multiple of 4
        SBCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{v1-v6,lr}          // save work regs
LABEL(subx_loop_up_l2)
        LDMIA   a2!,{v1,v2,v3,ip}       // load 4 words in one go
        LDMIA   a1!,{v4,v5,v6,lr}       // and from source2
        SBCS    v4,v4,v1                // subtract the four words with carry
        SBCS    v5,v5,v2
        SBCS    v6,v6,v3
        SBCS    lr,lr,ip
        STMIA   a3!,{v4,v5,v6,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4, preserve C
        TEQ     a4,#0                   // are we done ?
        BNE     subx_loop_up_l2         // if count non-zero then loop
        SBC     a1,a4,a4                // set result to Carry (a4 is 0)
        LDMFD   sp!,{v1-v6,pc}^         // restore work regs and return

// extern uintD subfrom_loop_up (uintD* sourceptr, uintD* destptr, uintC count);
//       entry
//               a1 = sourceptr
//               a2 = destptr
//               a3 = count of words to be subtracted
//       exit
//               destptr[] = destptr[] - sourceptr[]
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(subfrom_loop_up)         // word aligned subfrom loop up
        DECLARE_FUNCTION(subfrom_loop_up)
GLABEL(subfrom_loop_up)
        ANDS    ip,a3,#3                // multiple of 4 words ?
        BEQ     subfrom_loop_up_l1      // yup, so branch
        STMFD   sp!,{lr}
        LDR     a4,[a1],#4              // subtract the first 1-3 words
        LDR     lr,[a2]                 // to align the total to a multiple
        SUBS    lr,lr,a4                // of 4 words
        STR     lr,[a2],#4
        TEQ     ip,#1
        BNE     subfrom_loop_up_l0      // branch if more than one subtract
LABEL(subfrom_loop_up_l4)                    // drop through for better instr. timings
        BICS    a4,a3,#3                // set counter to multiple of 4
        SBCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        LDMEQFD sp!,{pc}^               // and return
        STMFD   sp!,{v1-v5}             // save work regs
        B       subfrom_loop_up_l2      // branch if more subtracts to do
LABEL(subfrom_loop_up_l0)
        LDR     a4,[a1],#4
        LDR     lr,[a2]
        SBCS    lr,lr,a4
        STR     lr,[a2],#4
        TEQ     ip,#2
        BEQ     subfrom_loop_up_l4      // need to branch 'cos PSR used
        LDR     a4,[a1],#4
        LDR     lr,[a2]
        SBCS    lr,lr,a4
        STR     lr,[a2],#4
        B       subfrom_loop_up_l4
LABEL(subfrom_loop_up_l1)
        BICS    a4,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,#0                   // no subtracts, so C = 0
        MOVEQS  pc,lr                   // if zero then we're done
        CMP     a4,#0                   // set carry bit, since a4 > 0
        STMFD   sp!,{v1-v5,lr}          // save work regs
LABEL(subfrom_loop_up_l2)
        LDMIA   a1!,{a3,v1,v2,ip}       // load 4 words in one go
        LDMIA   a2,{v3,v4,v5,lr}        // and from destptr
        SBCS    v3,v3,a3                // subtract the four words with carry
        SBCS    v4,v4,v1
        SBCS    v5,v5,v2
        SBCS    lr,lr,ip
        STMIA   a2!,{v3,v4,v5,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4, preserve C
        TEQ     a4,#0                   // are we done ?
        BNE     subfrom_loop_up_l2      // if count non-zero then loop
        SBC     a1,a4,a4                // set result to Carry (a4 is 0)
        LDMFD   sp!,{v1-v5,pc}^         // restore work regs and return

// extern uintD dec_loop_up (uintD* ptr, uintC count);
//       entry
//               a1 = ptr
//               a2 = count of words to be DECed
//       exit
//               a1 = 0 if any words are non-zero before decrement else -1
//                      stop decrementing when first word is non-zero
//               a2 - a4, ip destroyed
        EXPORT(dec_loop_up)             // word aligned dec loop up
        DECLARE_FUNCTION(dec_loop_up)
GLABEL(dec_loop_up)
        ANDS    a3,a2,#1                // multiple of 2 words ?
        BEQ     dec_loop_up_l1          // yup, so branch
        LDR     a4,[a1]                 // DEC the first word
        SUBS    a4,a4,#1                // align the total to a multiple of 2
        STR     a4,[a1],#4
        MOVCS   a1,#0                   // set result to 0
        MOVCSS  pc,lr                   // return 0 if non-zero result
LABEL(dec_loop_up_l1)
        BICS    a4,a2,#1                // set counter to multiple of 2
        MVNEQ   a1,#0                   // return -1
        MOVEQS  pc,lr                   // if zero then we're done
        MOV     ip,a1                   // move ptr to ip
        MOV     a1,#0                   // set result to 0
        ANDS    a3,a4,#3
        BEQ     dec_loop_up_l3
        LDMIA   ip,{a2,a3}              // load 2 words in one go
        SUBS    a2,a2,#1                // DEC the two words
        SUBCCS  a3,a3,#1                // stopping when first word non-zero
        STMIA   ip!,{a2,a3}             // store 2 results
        MOVCSS  pc,lr                   // return 0 if any result non-zero
        SUBS    a4,a4,#2                // decrement counter by 2
        MVNEQ   a1,#0                   // if finished loop then
        MOVEQS  pc,lr                   // return -1
LABEL(dec_loop_up_l3)                   // now a multiple of 4 words
        STMFD   sp!,{v1,lr}             // save work regs
LABEL(dec_loop_up_l2)
        LDMIA   ip,{a2,a3,v1,lr}        // load 4 words in one go
        SUBS    a2,a2,#1                // DEC the four words
        SUBCCS  a3,a3,#1                // stopping when first word non-zero
        SUBCCS  v1,v1,#1
        SUBCCS  lr,lr,#1
        STMIA   ip!,{a2,a3,v1,lr}       // store 4 results
        LDMCSFD sp!,{v1,pc}^            // return 0 if any carry
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     dec_loop_up_l2          // if count still positive then loop
        MVN     a1,#0
        LDMFD   sp!,{v1,pc}^            // restore work regs and return -1

// extern void neg_loop_up (uintD* ptr, uintC count);
//       entry
//               a1 = ptr
//               a2 = count of words. The long integer is to be NEGated
//       exit
//               ptr[] = -ptr[] for count words
//               a1 = last carry
//               a2 - a4, ip destroyed
        EXPORT(neg_loop_up)             // word aligned neg loop up
        DECLARE_FUNCTION(neg_loop_up)
GLABEL(neg_loop_up)
        CMPS    a2,#0                   // count = 0 ?
        MOVEQ   a1,#0                   // yup, so return 0
        MOVEQS  pc,lr
LABEL(neg_loop_up_l1)                   // skip all the zero words first
        LDR     a3,[a1],#4              // compare words against zero
        CMPS    a3,#0                   // upwards in memory
        BNE     neg_loop_up_l2          // non-zero, so negate rest of words
        SUBS    a2,a2,#1                // reduce count of words
        BNE     neg_loop_up_l1          // more ?, so loop
        MOV     a1,#0                   // return 0
        MOVS    pc,lr
LABEL(neg_loop_up_l2)
        RSB     a3,a3,#0                // first non-zero word = -word
        STR     a3,[a1,#-4]
        SUBS    a2,a2,#1
        MVNEQ   a1,#0                   // done ? -> return -1
        MOVEQS  pc,lr
                                        // now NOT rest of the words
        ANDS    a3,a2,#3                // multiple of 4 words ?
        BEQ     neg_loop_up_l3          // yup, so branch
        CMP     a3,#2                   // NOT the first 1-3 words
        LDR     a3,[a1]                 // to align the total to a multiple
        MVN     a3,a3                   // of 4 words
        STR     a3,[a1],#4
        BLT     neg_loop_up_l3          // better to branch than skip instrs.
        LDRGE   a3,[a1]
        MVNGE   a3,a3
        STRGE   a3,[a1],#4
        LDRGT   a3,[a1]
        MVNGT   a3,a3
        STRGT   a3,[a1],#4
LABEL(neg_loop_up_l3)
        BICS    a4,a2,#3                // set counter to multiple of 4
        MVNEQ   a1,#0                   // set result to -1
        MOVEQS  pc,lr                   // if zero then we're done
        STMFD   sp!,{lr}                // save work regs
LABEL(neg_loop_up_l4)
        LDMIA   a1,{a2,a3,ip,lr}        // load 4 words in one go,NO writeback
        MVN     a2,a2                   // NOT the four words
        MVN     a3,a3
        MVN     ip,ip
        MVN     lr,lr
        STMIA   a1!,{a2,a3,ip,lr}       // store 4 results
        SUBS    a4,a4,#4                // decrement counter by 4
        BGT     neg_loop_up_l4          // if count still positive then loop
        MVN     a1,#0                   // set result to -1
        LDMFD   sp!,{pc}^               // restore work regs and return -1

// extern uintD shift1left_loop_up (uintD* ptr, uintC count);
//       entry
//               a1 = ptr
//               a2 = count of words to be shifted left
//       exit
//               a1 = carry out from last shift left
//               a2 - a4, ip destroyed
        EXPORT(shift1left_loop_up)      // word aligned shift1left loop up
        DECLARE_FUNCTION(shift1left_loop_up)
GLABEL(shift1left_loop_up)
        CMN     a1,#0                   // clear carry bit, since a1 > 0
        ANDS    a3,a2,#1                // multiple of 2 words ?
        BEQ     shift1left_loop_up_l1   // yup, so branch
        LDR     a4,[a1]                 // shift left the first word
        ADDS    a4,a4,a4
        STR     a4,[a1],#4
LABEL(shift1left_loop_up_l1)
        BICS    a4,a2,#1                // set counter to multiple of 2
        ADCEQ   a1,a4,a4                // if zero set result to C (a4 is 0)
        MOVEQS  pc,lr                   // and return
        ANDS    a3,a4,#3                // multiple of 4 words ?
        BEQ     shift1left_loop_up_l3   // yup, so branch
        LDMIA   a1,{a2,a3}              // load 2 words in one go
        ADCS    a2,a2,a2                // shift left the two words
        ADCS    a3,a3,a3
        STMIA   a1!,{a2,a3}             // store 2 results
        BICS    a4,a4,#2                // decrement counter by 2
        ADCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        MOVEQS  pc,lr                   // and return
LABEL(shift1left_loop_up_l3)            // now a multiple of 4 words
        STMFD   sp!,{lr}                // save work regs
LABEL(shift1left_loop_up_l2)
        LDMIA   a1,{a2,a3,ip,lr}        // load 4 words in one go
        ADCS    a2,a2,a2                // shift left the four words
        ADCS    a3,a3,a3
        ADCS    ip,ip,ip
        ADCS    lr,lr,lr
        STMIA   a1!,{a2,a3,ip,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4
        TEQ     a4,#0                   // are we done ?
        BNE     shift1left_loop_up_l2   // if count non-zero then loop
        ADC     a1,a4,a4                // set result to Carry (a4 is 0)
        LDMFD   sp!,{pc}^               // restore work regs and return 1

// extern uintD shiftleft_loop_up (uintD* ptr, uintC count, uintC i, uintD carry);
//       entry
//               a1 = ptr
//               a2 = count of words to be shifted left
//               a3 = size of left shift
//               a4 = value to ORR in for first shift
//       exit
//               a1 = shift out from last shift left
//               a2 - a4, ip destroyed
        EXPORT(shiftleft_loop_up)       // word aligned shiftleft loop up
        DECLARE_FUNCTION(shiftleft_loop_up)
GLABEL(shiftleft_loop_up)
        STMFD   sp!,{v6,lr}
        RSB     v6,a3,#32               // size of complementary right shift
        ANDS    ip,a2,#3                // multiple of 4 words ?
        BEQ     shiftleft_loop_up_l1    // yup, so branch
        LDR     lr,[a1]                 // shiftleft the first 1-3 words
        ORR     a4,a4,lr,ASL a3         // to align the total to a multiple
        STR     a4,[a1],#4              // of 4 words
        MOV     a4,lr,LSR v6
        CMP     ip,#2
        BLT     shiftleft_loop_up_l1    // better to branch than skip instrs.
        LDRGE   lr,[a1]
        ORRGE   a4,a4,lr,ASL a3
        STRGE   a4,[a1],#4
        MOVGE   a4,lr,LSR v6
        LDRGT   lr,[a1]
        ORRGT   a4,a4,lr,ASL a3
        STRGT   a4,[a1],#4
        MOVGT   a4,lr,LSR v6
LABEL(shiftleft_loop_up_l1)
        BICS    ip,a2,#3                // set counter to multiple of 4
        MOVEQ   a1,a4                   // if zero then we're done
        LDMEQFD sp!,{v6,pc}^            // so return last shift out
        STMFD   sp!,{v1-v3}             // save work regs
LABEL(shiftleft_loop_up_l2)
        LDMIA   a1,{v1,v2,v3,lr}        // load 4 words in one go
        ORR     a2,a4,v1,ASL a3         // shiftleft the four words
        MOV     a4,v1,LSR v6            // keep carry in a4
        ORR     v1,a4,v2,ASL a3         // and store results down a register
        MOV     a4,v2,LSR v6            // to regs a2,v1-v3
        ORR     v2,a4,v3,ASL a3
        MOV     a4,v3,LSR v6
        ORR     v3,a4,lr,ASL a3
        MOV     a4,lr,LSR v6
        STMIA   a1!,{a2,v1,v2,v3}       // store 4 results
        SUBS    ip,ip,#4                // decrement counter by 4
        BGT     shiftleft_loop_up_l2    // if count still positive then loop
        MOV     a1,a4                   // result = last shift out
        LDMFD   sp!,{v1-v3,v6,pc}^      // restore work regs and return

#endif

// extern uintD shiftleftcopy_loop_up (uintD* sourceptr, uintD* destptr, uintC count, uintC i);
//       entry
//               a1 = sourceptr
//               a2 = destptr
//               a3 = count of words to be shifted left
//               a4 = size of left shift
//       exit
//               a1 = shift out from last shift left
//               a2 - a4, ip destroyed
        EXPORT(shiftleftcopy_loop_up)   // word aligned shiftleftcopy loop up
        DECLARE_FUNCTION(shiftleftcopy_loop_up)
GLABEL(shiftleftcopy_loop_up)
        STMFD   sp!,{v5,v6,lr}
        MOV     v5,#0                   // initial shift carry
        RSB     v6,a4,#32               // size of complementary right shift
        ANDS    ip,a3,#3                // multiple of 4 words ?
        BEQ     shiftleftcopy_loop_up_l1 // yup, so branch
        LDR     lr,[a1],#4              // shiftleft the first 1-3 words
        ORR     v5,v5,lr,ASL a4         // to align the total to a multiple
        STR     v5,[a2],#4              // of 4 words
        MOV     v5,lr,LSR v6
        CMP     ip,#2
        BLT     shiftleftcopy_loop_up_l1 // better to branch than skip instrs.
        LDRGE   lr,[a1],#4
        ORRGE   v5,v5,lr,ASL a4
        STRGE   v5,[a2],#4
        MOVGE   v5,lr,LSR v6
        LDRGT   lr,[a1],#4
        ORRGT   v5,v5,lr,ASL a4
        STRGT   v5,[a2],#4
        MOVGT   v5,lr,LSR v6
LABEL(shiftleftcopy_loop_up_l1)
        BICS    ip,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,v5                   // if zero then we're done
        LDMEQFD sp!,{v5,v6,pc}^         // so return last shift out
        STMFD   sp!,{v1-v3}             // save work regs
LABEL(shiftleftcopy_loop_up_l2)
        LDMIA   a1!,{v1,v2,v3,lr}       // load 4 words in one go
        ORR     a3,v5,v1,ASL a4         // shiftleft the four words
        MOV     v5,v1,LSR v6            // keep carry in v5
        ORR     v1,v5,v2,ASL a4         // and store results down a register
        MOV     v5,v2,LSR v6            // to regs a3,v1-v3
        ORR     v2,v5,v3,ASL a4
        MOV     v5,v3,LSR v6
        ORR     v3,v5,lr,ASL a4
        MOV     v5,lr,LSR v6
        STMIA   a2!,{a3,v1,v2,v3}       // store 4 results
        SUBS    ip,ip,#4                // decrement counter by 4
        BGT     shiftleftcopy_loop_up_l2 // if count still positive then loop
        MOV     a1,v5                   // result = last shift out
        LDMFD   sp!,{v1-v3,v5,v6,pc}^   // restore work regs and return

#if !CL_DS_BIG_ENDIAN_P

// extern uintD shift1right_loop_down (uintD* ptr, uintC count, uintD carry);
//       entry
//               a1 = ptr
//               a2 = count of words to be shifted right
//               a3 = carry
//       exit
//               a1 = carry out from last shift right
//               a2 - a4, ip destroyed
        EXPORT(shift1right_loop_down)   // word aligned shift1right loop down
        DECLARE_FUNCTION(shift1right_loop_down)
GLABEL(shift1right_loop_down)
        MOVS    a3,a3,LSR #1            // set carry
        ANDS    a3,a2,#1                // multiple of 2 words ?
        BEQ     shift1right_loop_down_l1 // yup, so branch
        LDR     a4,[a1,#-4]!            // shift right the first word
        MOVS    a4,a4,RRX
        STR     a4,[a1]
LABEL(shift1right_loop_down_l1)
        BICS    a4,a2,#1                // set counter to multiple of 2
        MOVEQ   a1,a4,RRX               // if zero set result to C (a4 is 0)
        MOVEQS  pc,lr                   // and return
        ANDS    a3,a4,#3                // multiple of 4 words ?
        BEQ     shift1right_loop_down_l3 // yup, so branch
        LDMDB   a1,{a2,a3}              // load 2 words in one go
        MOVS    a3,a3,RRX               // shift right the two words
        MOVS    a2,a2,RRX
        STMDB   a1!,{a2,a3}             // store 2 results
        BICS    a4,a4,#2                // decrement counter by 2
        ADCEQ   a1,a4,a4                // set result to Carry (a4 is 0)
        MOVEQS  pc,lr                   // and return
LABEL(shift1right_loop_down_l3)         // now a multiple of 4 words
        STMFD   sp!,{lr}                // save work regs
LABEL(shift1right_loop_down_l2)
        LDMDB   a1,{a2,a3,ip,lr}        // load 4 words in one go
        MOVS    lr,lr,RRX               // shift right the four words
        MOVS    ip,ip,RRX
        MOVS    a3,a3,RRX
        MOVS    a2,a2,RRX
        STMDB   a1!,{a2,a3,ip,lr}       // store 4 results
        SUB     a4,a4,#4                // decrement counter by 4
        TEQ     a4,#0                   // are we done ?
        BNE     shift1right_loop_down_l2 // if count non-zero then loop
        MOV     a1,a4,RRX               // set result to Carry (a4 is 0)
        LDMFD   sp!,{pc}^               // restore work regs and return 1

// extern uintD shiftright_loop_down (uintD* ptr, uintC count, uintC i);
//       entry
//               a1 = ptr
//               a2 = count of words to be shifted right
//               a3 = size of right shift
//       exit
//               a1 = shift out from last shift right
//               a2 - a4, ip destroyed
        EXPORT(shiftright_loop_down)    // word aligned shiftright loop down
        DECLARE_FUNCTION(shiftright_loop_down)
GLABEL(shiftright_loop_down)
        STMFD   sp!,{v6,lr}
        MOV     a4,#0                   // initial shift carry
        RSB     v6,a3,#32               // size of complementary left shift
LABEL(shiftright_loop_down_l0)
        ANDS    ip,a2,#3                // multiple of 4 words ?
        BEQ     shiftright_loop_down_l1 // yup, so branch
        LDR     lr,[a1,#-4]!            // shiftright the first 1-3 words
        ORR     a4,a4,lr,LSR a3         // to align the total to a multiple
        STR     a4,[a1]                 // of 4 words
        MOV     a4,lr,ASL v6
        CMP     ip,#2
        BLT     shiftright_loop_down_l1 // better to branch than skip instrs.
        LDRGE   lr,[a1,#-4]!
        ORRGE   a4,a4,lr,LSR a3
        STRGE   a4,[a1]
        MOVGE   a4,lr,ASL v6
        LDRGT   lr,[a1,#-4]!
        ORRGT   a4,a4,lr,LSR a3
        STRGT   a4,[a1]
        MOVGT   a4,lr,ASL v6
LABEL(shiftright_loop_down_l1)
        BICS    ip,a2,#3                // set counter to multiple of 4
        MOVEQ   a1,a4                   // if zero then we're done
        LDMEQFD sp!,{v6,pc}^            // so return last shift out
        STMFD   sp!,{v1-v3}             // save work regs
LABEL(shiftright_loop_down_l2)
        LDMDB   a1,{a2,v1,v2,v3}        // load 4 words in one go
        ORR     lr,a4,v3,LSR a3         // shiftright the four words
        MOV     a4,v3,ASL v6            // keep carry in a4
        ORR     v3,a4,v2,LSR a3         // and store results up a register
        MOV     a4,v2,ASL v6            // to regs v1-v3,lr
        ORR     v2,a4,v1,LSR a3
        MOV     a4,v1,ASL v6
        ORR     v1,a4,a2,LSR a3
        MOV     a4,a2,ASL v6
        STMDB   a1!,{v1,v2,v3,lr}       // store 4 results
        SUBS    ip,ip,#4                // decrement counter by 4
        BGT     shiftright_loop_down_l2 // if count still positive then loop
        MOV     a1,a4                   // result = last shift out
        LDMFD   sp!,{v1-v3,v6,pc}^      // restore work regs and return

// extern uintD shiftrightsigned_loop_down (uintD* ptr, uintC count, uintC i);
//       entry
//               a1 = ptr
//               a2 = count of words to be shifted right signed
//               a3 = size of right shift
//       exit
//               a1 = shift out from last shift right
//               a2 - a4, ip destroyed
        EXPORT(shiftrightsigned_loop_down)// word aligned shiftrightsigned loop down
        DECLARE_FUNCTION(shiftrightsigned_loop_down)
GLABEL(shiftrightsigned_loop_down)
        STMFD   sp!,{v6,lr}
        RSB     v6,a3,#32               // size of complementary left shift
        LDR     lr,[a1,#-4]             // setup carry for first shift.
        MOV     a4,lr,ASR #31           // this is the sign extended bits
        AND     a4,a4,a4,LSL v6         // 31->(32-i) of the first word
        B       shiftright_loop_down_l0 // use right shift code now

// extern uintD shiftrightcopy_loop_down (uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry);
//       entry
//               a1 = sourceptr
//               a2 = destptr
//               a3 = count of words to be shifted right
//               a4 = size of right shift
//               [sp] = carry for first shift
//       exit
//               a1 = shift out from last shift right
//               a2 - a4, ip destroyed
        EXPORT(shiftrightcopy_loop_down)// word aligned shiftrightcopy loop down
        DECLARE_FUNCTION(shiftrightcopy_loop_down)
GLABEL(shiftrightcopy_loop_down)
        STMFD   sp!,{v5,v6,lr}
        LDR     v5,[sp,#12]             // initial shift carry
        RSB     v6,a4,#32               // size of complementary left shift
        MOV     v5,v5,ASL v6
LABEL(shiftrightcopy_loop_down_l0)
        ANDS    ip,a3,#3                // multiple of 4 words ?
        BEQ     shiftrightcopy_loop_down_l1 // yup, so branch
        LDR     lr,[a1,#-4]!            // shiftright the first 1-3 words
        ORR     v5,v5,lr,LSR a4         // to align the total to a multiple
        STR     v5,[a2,#-4]!            // of 4 words
        MOV     v5,lr,ASL v6
        CMP     ip,#2
        BLT     shiftrightcopy_loop_down_l1 // better to branch than skip instrs.
        LDRGE   lr,[a1,#-4]!
        ORRGE   v5,v5,lr,LSR a4
        STRGE   v5,[a2,#-4]!
        MOVGE   v5,lr,ASL v6
        LDRGT   lr,[a1,#-4]!
        ORRGT   v5,v5,lr,LSR a4
        STRGT   v5,[a2,#-4]!
        MOVGT   v5,lr,ASL v6
LABEL(shiftrightcopy_loop_down_l1)
        BICS    ip,a3,#3                // set counter to multiple of 4
        MOVEQ   a1,v5                   // if zero then we're done
        LDMEQFD sp!,{v5,v6,pc}^         // so return last shift out
        STMFD   sp!,{v1-v3}             // save work regs
LABEL(shiftrightcopy_loop_down_l2)
        LDMDB   a1!,{a3,v1,v2,v3}       // load 4 words in one go
        ORR     lr,v5,v3,LSR a4         // shiftright the four words
        MOV     v5,v3,ASL v6            // keep carry in v5
        ORR     v3,v5,v2,LSR a4         // and store results up a register
        MOV     v5,v2,ASL v6            // to regs v1-v3,lr
        ORR     v2,v5,v1,LSR a4
        MOV     v5,v1,ASL v6
        ORR     v1,v5,a3,LSR a4
        MOV     v5,a3,ASL v6
        STMDB   a2!,{v1,v2,v3,lr}       // store 4 results
        SUBS    ip,ip,#4                // decrement counter by 4
        BGT     shiftrightcopy_loop_down_l2 // if count still positive then loop
        MOV     a1,v5                   // result = last shift out
        LDMFD   sp!,{v1-v3,v5,v6,pc}^   // restore work regs and return

#ifndef HAVE_umull
// mulu32_64_vregs
//       entry
//               a1 = x
//               ip = y
//       exit
//               v1 = low32(x*y)
//               ip = high32(x*y)
//               v2,v3,v4 destroyed
LABEL(mulu32_64_vregs)
        MOV     v1,a1,LSR #16           // temp := top half of x
        MOV     v2,ip,LSR #16           // hi := top half of y
        BIC     v3,a1,v1,LSL #16        // x  := bottom half of x
        BIC     ip,ip,v2,LSL #16        // y  := bottom half of y
        MUL     v4,v3,ip                // low section of result
        MUL     ip,v1,ip                // ) middle sections
        MUL     v3,v2,v3                // )   of result
        MUL     v2,v1,v2                // high section of result
        ADDS    ip,ip,v3                // add middle sections
                                        // (can't use mla as we need carry)
        ADDCS   v2,v2,#0x10000          // carry from above add
        ADDS    v1,v4,ip,LSL #16        // x is now bottom 32 bits of result
        ADC     ip,v2,ip,LSR #16        // hi is top 32 bits
        MOVS    pc,lr
#endif

// extern uintD mulusmall_loop_up (uintD digit, uintD* ptr, uintC len, uintD newdigit);
//       entry
//               a1 = digit
//               a2 = ptr
//               a3 = count of words to be multiplied up
//               a4 = new digit = carry
//       exit
//               a1 = final carry of multiply
//               a2 - a4, ip destroyed
        EXPORT(mulusmall_loop_up)
        DECLARE_FUNCTION(mulusmall_loop_up)
GLABEL(mulusmall_loop_up)
        CMP     a3,#0
        MOVEQ   a1,a4
        MOVEQS  pc,lr
#ifdef HAVE_umull
        STMFD   sp!,{v1,lr}
LABEL(mulusmall_loop_up_l1)
        LDR     ip,[a2]
        UMULL   v1,ip,a1,ip             // muluD(digit,*--ptr,hi=,lo=)
        ADDS    v1,v1,a4                // lo += carry
        ADC     a4,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a2],#4              // *ptr++ = lo
        SUBS    a3,a3,#1                // len--
        BNE     mulusmall_loop_up_l1    // until len==0
        MOV     a1,a4                   // return carry
        LDMFD   sp!,{v1,pc}^
#else
        STMFD   sp!,{v1-v2,lr}
LABEL(mulusmall_loop_up_l1)
        LDR     ip,[a2]

//      BL      mulu32_64_vregs         // muluD(digit,*ptr,hi=,lo=)
// replaced by multiplication of a small x = a1 and a big y = ip :
        MOV     v1,ip,LSR #16           // top half of y
        BIC     ip,ip,v1,LSL #16        // bottom half of y
        MUL     v2,a1,v1                // middle section of result
        MUL     v1,a1,ip                // low section of result
        MOV     ip,#0                   // high section of result
        ADDS    v1,v1,v2,LSL #16        // bottom 32 bits of result
        ADC     ip,ip,v2,LSR #16        // top 32 bits of result

        ADDS    v1,v1,a4                // lo += carry
        ADC     a4,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a2],#4              // *ptr++ = lo
        SUBS    a3,a3,#1                // len--
        BNE     mulusmall_loop_up_l1    // until len==0
        MOV     a1,a4                   // return carry
        LDMFD   sp!,{v1-v2,pc}^
#endif

// extern void mulu_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
//       entry
//               a1 = digit
//               a2 = sourceptr
//               a3 = destptr
//               a4 = count of words to be multiplied up
//       exit
//               a1 - a4, ip destroyed
        EXPORT(mulu_loop_up)
        DECLARE_FUNCTION(mulu_loop_up)
GLABEL(mulu_loop_up)
#ifdef HAVE_umull
        STMFD   sp!,{v1,v5,lr}
        MOV     v5,#0
LABEL(mulu_loop_up_l1)
        LDR     ip,[a2],#4
        UMULL   v1,ip,a1,ip             // muluD(digit,*sourceptr++,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a3],#4              // *destptr++ = lo
        SUBS    a4,a4,#1                // len--
        BNE     mulu_loop_up_l1         // until len==0
        STR     v5,[a3],#4              // *destptr++ = carry
        LDMFD   sp!,{v1,v5,pc}^
#else
        STMFD   sp!,{v1-v5,lr}
        MOV     v5,#0
LABEL(mulu_loop_up_l1)
        LDR     ip,[a2],#4
        BL      mulu32_64_vregs         // muluD(digit,*sourceptr++,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a3],#4              // *destptr++ = lo
        SUBS    a4,a4,#1                // len--
        BNE     mulu_loop_up_l1         // until len==0
        STR     v5,[a3],#4              // *destptr++ = carry
        LDMFD   sp!,{v1-v5,pc}^
#endif

// extern void muluadd_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
//       entry
//               a1 = digit
//               a2 = sourceptr
//               a3 = destptr
//               a4 = count of words to be multiplied added up
//       exit
//               a1 - a4, ip destroyed
        EXPORT(muluadd_loop_up)
        DECLARE_FUNCTION(muluadd_loop_up)
GLABEL(muluadd_loop_up)
#ifdef HAVE_umull
        STMFD   sp!,{v1,v5,lr}
        MOV     v5,#0
LABEL(muluadd_loop_up_l1)
        LDR     ip,[a2],#4
        UMULL   v1,ip,a1,ip             // muluD(digit,*sourceptr++,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADCCS   ip,ip,#0                // if (lo<carry) { hi += 1 };
        LDR     v5,[a3]                 // carry = *destptr
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a3],#4              // *destptr++ = lo
        SUBS    a4,a4,#1                // len--
        BNE     muluadd_loop_up_l1      // until len==0
        MOV     a1,v5                   // return carry
        LDMFD   sp!,{v1,v5,pc}^
#else
        STMFD   sp!,{v1-v5,lr}
        MOV     v5,#0
LABEL(muluadd_loop_up_l1)
        LDR     ip,[a2],#4
        BL      mulu32_64_vregs         // muluD(digit,*sourceptr++,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADCCS   ip,ip,#0                // if (lo<carry) { hi += 1 };
        LDR     v5,[a3]                 // carry = *destptr
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 }; carry=hi
        STR     v1,[a3],#4              // *destptr++ = lo
        SUBS    a4,a4,#1                // len--
        BNE     muluadd_loop_up_l1      // until len==0
        MOV     a1,v5                   // return carry
        LDMFD   sp!,{v1-v5,pc}^
#endif

// extern void mulusub_loop_up (uintD digit, uintD* sourceptr, uintD* destptr, uintC len);
//       entry
//               a1 = digit
//               a2 = sourceptr
//               a3 = destptr
//               a4 = count of words to be multiplied subtracted up
//       exit
//               a1 - a4, ip destroyed
        EXPORT(mulusub_loop_up)
        DECLARE_FUNCTION(mulusub_loop_up)
GLABEL(mulusub_loop_up)
#ifdef HAVE_umull
        STMFD   sp!,{v1,v5,lr}
        MOV     v5,#0
LABEL(mulusub_loop_up_l1)
        LDR     ip,[a2],#4
        UMULL   v1,ip,a1,ip             // muluD(digit,*sourceptr++,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 };
        LDR     ip,[a3]                 // carry = *destptr
        SUBS    ip,ip,v1
        STR     ip,[a3],#4              // *destptr++ = carry - lo
        ADDCC   v5,v5,#1                // if (carry<lo) { hi += 1 }; carry=hi
        SUBS    a4,a4,#1                // len--
        BNE     mulusub_loop_up_l1      // until len==0
        MOV     a1,v5                   // return carry
        LDMFD   sp!,{v1,v5,pc}^
#else
        STMFD   sp!,{v1-v5,lr}
        MOV     v5,#0
LABEL(mulusub_loop_up_l1)
        LDR     ip,[a2],#4
        BL      mulu32_64_vregs         // muluD(digit,*sourceptr++,hi=,lo=)
        ADDS    v1,v1,v5                // lo += carry
        ADC     v5,ip,#0                // if (lo<carry) { hi += 1 };
        LDR     ip,[a3]                 // carry = *destptr
        SUBS    ip,ip,v1
        STR     ip,[a3],#4              // *destptr++ = carry - lo
        ADDCC   v5,v5,#1                // if (carry<lo) { hi += 1 }; carry=hi
        SUBS    a4,a4,#1                // len--
        BNE     mulusub_loop_up_l1      // until len==0
        MOV     a1,v5                   // return carry
        LDMFD   sp!,{v1-v5,pc}^
#endif

#endif

// extern void shiftxor_loop_up (uintD* xptr, const uintD* yptr, uintC count, uintC i);
//       entry
//               a1 = xptr
//               a2 = yptr
//               a3 = count of words to be shifted left
//               a4 = size of left shift
//       exit
//               a1 - a4, ip destroyed
        EXPORT(shiftxor_loop_up)        // word aligned shiftxor loop up
        DECLARE_FUNCTION(shiftxor_loop_up)
GLABEL(shiftxor_loop_up)
        STMFD   sp!,{v5,v6,lr}
        RSB     lr,a4,#32               // size of complementary right shift
        LDR     ip,[a1]                 // get first *xptr
        ANDS    v6,a3,#3                // multiple of 4 words ?
        BEQ     shiftxor_loop_up_l1     // yup, so branch
        LDR     v5,[a2],#4              // get *yptr
        EOR     ip,ip,v5,ASL a4         // combine with modified *xptr
        STR     ip,[a1],#4              // save new *xptr
        LDR     ip,[a1]                 // get next *xptr
        EOR     ip,ip,v5,LSR lr         // combine with *xptr
        CMP     v6,#2
        BLT     shiftxor_loop_up_l1     // better to branch than skip instrs.
        LDR     v5,[a2],#4              // get *yptr
        EOR     ip,ip,v5,ASL a4         // combine with modified *xptr
        STR     ip,[a1],#4              // save new *xptr
        LDR     ip,[a1]                 // get next *xptr
        EOR     ip,ip,v5,LSR lr         // combine with *xptr
        LDRGT   v5,[a2],#4              // get *yptr
        EORGT   ip,ip,v5,ASL a4         // combine with modified *xptr
        STRGT   ip,[a1],#4              // save new *xptr
        LDRGT   ip,[a1]                 // get next *xptr
        EORGT   ip,ip,v5,LSR lr         // combine with *xptr
LABEL(shiftxor_loop_up_l1)
        BICS    a3,a3,#3                // set counter to multiple of 4
        STREQ   ip,[a1]
        LDMEQFD sp!,{v5,v6,pc}^         // return if done
        STMFD   sp!,{v1-v4}             // save work regs
LABEL(shiftxor_loop_up_l2)
        LDMIA   a2!,{v3,v4,v5,v6}       // load 4 words yptr[0..3] in one go
        EOR     v1,ip,v3,ASL a4         // combine with modified *xptr
        LDR     v2,[a1,#4]
        EOR     v2,v2,v3,LSR lr
        EOR     v2,v2,v4,ASL a4         // combine with modified *xptr
        LDR     v3,[a1,#8]
        EOR     v3,v3,v4,LSR lr
        EOR     v3,v3,v5,ASL a4         // combine with modified *xptr
        LDR     v4,[a1,#12]
        EOR     v4,v4,v5,LSR lr
        EOR     v4,v4,v6,ASL a4         // combine with modified *xptr
        STMIA   a1!,{v1,v2,v3,v4}       // store 4 words xptr[0..3] in one go
        LDR     ip,[a1]
        EOR     ip,ip,v6,LSR lr
        SUBS    a3,a3,#4                // decrement counter by 4
        BGT     shiftxor_loop_up_l2
        STR     ip,[a1]
        LDMFD   sp!,{v1-v6,pc}^         // restore work regs and return

        END
