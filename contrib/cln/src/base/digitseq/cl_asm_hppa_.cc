// Externe Routinen zu ARILEV1.D
// Prozessor: HPPA, wegen XMPYU nur auf HPPA 1.1 (etwa HP9000/720)
// Compiler: GNU-C oder HP-C
// Parameter-Übergabe: in Registern %arg0,%arg1,%arg2, Rückgabewert in %ret0.
// Einstellungen: intCsize=32, intDsize=32.

// Großenteils abgeschrieben von hppa.s aus der PARI/GP-Distribution.

                .SHORTDATA
                .IMPORT $global$,DATA


                .CODE
                .EXPORT length32
// Liefert integer-size (>=1, <=32) des Arguments /=0.
length32        .PROC
                .CALLINFO
                .ENTER          // Input in %arg0, Output in %ret0
                // y = 1;
                LDI             1,%ret0
                // if (x & (bit(31-15)*(bit(16)-1)) == 0)
                EXTRU,<>        %arg0,15,16,%r0
                SHD,TR          %arg0,%r0,16,%arg0      // x = x<<(32-16); else
                ADDI            16,%ret0,%ret0          // y = y+16;
                // if (x & (bit(31-7)*(bit(8)-1)) == 0)
                EXTRU,<>        %arg0,7,8,%r0
                SHD,TR          %arg0,%r0,24,%arg0      // x = x<<(32-24); else
                ADDI            8,%ret0,%ret0           // y = y+8;
                // if (x & (bit(31-3)*(bit(4)-1)) == 0)
                EXTRU,<>        %arg0,3,4,%r0
                SHD,TR          %arg0,%r0,28,%arg0      // x = x<<(32-28); else
                ADDI            4,%ret0,%ret0           // y = y+4;
                // if (x & (bit(31-1)*(bit(2)-1)) == 0)
                EXTRU,<>        %arg0,1,2,%r0
                SHD,TR          %arg0,%r0,30,%arg0      // x = x<<(32-30); else
                ADDI            2,%ret0,%ret0           // y = y+2;
                // if (x & (bit(31-0)*(bit(1)-1)) != 0)
                EXTRU,=         %arg0,0,1,%r0
                ADDI            1,%ret0,%ret0           // y = y+1;
                .LEAVE
                .PROCEND


#ifndef __GNUC__ /* mit GNU-C machen wir mulu32() als Macro, der inline multipliziert */

                .SHORTDATA
                .EXPORT mulu32_high
                .ALIGN 8
mulu32_high     .WORD           // 8 Byte Platz
                .WORD

                .CODE
                .EXPORT mulu32_
// extern struct { uint32 lo; uint32 hi; } mulu32_ (uint32 arg1, uint32 arg2);
// 2^32*hi+lo := arg1*arg2.
mulu32_         .PROC
                .CALLINFO
                .ENTER  // Input in %arg0,%arg1, Output in %ret0,mulu32_high
                LDIL    L'mulu32_high-$global$,%r1
                LDO     R'mulu32_high-$global$(%r1),%r1
                                                // %r1 = &x
                STW     %arg0,0(%r1)            // x abspeichern
                FLDWS   0(%r1),%fr4             // und in den Coprozessor laden
                STW     %arg1,0(%r1)            // y abspeichern
                FLDWS   0(%r1),%fr5             // und in den Coprozessor laden
                XMPYU   %fr4,%fr5,%fr6          // beides multiplizieren
                FSTDS   %fr6,0(%r1)             // Ergebnis (64 Bit) abspeichern
                LDWS    4(%r1),%ret0            // low 32 Bit als Ergebnis
                .LEAVE
                .PROCEND

#endif


                .END
