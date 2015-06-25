#include "atlas_asm.h"
/*
 * This file does a 1x12 unrolled r1_sse with these params:
 *    CL=16, ORDER=clmajor
 */
#ifndef ATL_GAS_x8664
   #error "This kernel requires x86-64 assembly!"
#endif
/*
 * Integer register assignment
 */
#define M       %rdi
#define N       %rsi
#define pA0     %r8
#define lda     %rax
#define pX      %rdx
#define pY      %r9
#define II      %rbx
#define pX0     %r11
#define Mr      %rcx
#define incAn   %r15
#define incII   $16*1
#define lda3    %r12
#define lda5    %r13
#define lda7    %rbp
#define lda9    %r14
#define lda11   %r10
#define incAXm $-64
/*
 * SSE register assignment
 */
#define rA0     %xmm0
#define rX0     %xmm1
#define ryt     %xmm2
#define rY0     %xmm3
#define rY1     %xmm4
#define rY2     %xmm5
#define rY3     %xmm6
#define rY4     %xmm7
#define rY5     %xmm8
#define rY6     %xmm9
#define rY7     %xmm10
#define rY8     %xmm11
#define rY9     %xmm12
#define rY10     %xmm13
#define rY11     %xmm14

/*
 * macros
 */
#ifndef MOVA
   #define MOVA movups
#endif
#define movapd movaps
#define movupd movups
#define xorpd xorps
#define addpd addps
#define mulpd mulps
#define addsd addss
#define mulsd mulss
#define movsd movss
#define haddpd haddps
/*
 * Define macros controlling prefetch
 */
#ifndef PFDIST
   #define PFDIST 256
#endif
#ifndef PFADIST
   #define PFADIST PFDIST
#endif
#ifndef PFYDIST
   #define PFYDIST 64
#endif
#ifndef PFXDIST
   #define PFXDIST 64
#endif
#ifndef PFIY
   #define PFIY prefetchnta
#endif
#ifndef PFIA
   #ifdef ATL_3DNow
      #define PFIA prefetchw
   #else
      #define PFIA prefetcht0
   #endif
#endif
#if PFADIST == 0                /* flag for no prefetch */
   #define prefA(mem)
#else
   #define prefA(mem) PFIA mem
#endif
#if PFYDIST == 0                /* flag for no prefetch */
   #define prefY(mem)
#else
   #define prefY(mem) PFIY mem
#endif
#if PFXDIST == 0                /* flag for no prefetch */
   #define prefX(mem)
#else
   #define prefX(mem) PFIX mem
#endif
.text
/*
 *                      %rdi        %rsi           %rdx          %rcx
 * void ATL_UGERK(ATL_CINT M, ATL_CINT N, const TYPE *X, const TYPE *Y,
 *                    %r8      %r9
 *                TYPE *A, ATL_CINT lda)
 */
.text
.global ATL_asmdecor(ATL_UGERK)
ALIGN64
ATL_asmdecor(ATL_UGERK):

/*
 * Save callee-saved iregs
 */
   movq %rbp, -8(%rsp)
   movq %rbx, -16(%rsp)
   movq %r12, -24(%rsp)
   movq %r13, -32(%rsp)
   movq %r14, -40(%rsp)
   movq %r15, -48(%rsp)
/*
 * Compute M = (M/MU)*MU, Mr = M - (M/MU)*MU
 * NOTE: Mr is %rcx reg, so we can use jcx to go to cleanup loop
 */
   mov  %r9, lda        /* move lda to assigned register, rax */
   mov  %rcx, pY        /* move pY to assigned register, r9 */
   mov  M, Mr           /* Mr = M */
   shr $4, M            /* M = M / MU */
   shl $4, M            /* M = (M/MU)*MU */
   sub M, Mr            /* Mr = M - (M/MU)*MU */
/*
 * Setup constants
 */
   mov lda, incAn       /* incAn = lda */
   sub M, incAn         /* incAn = lda - (M/MU)*MU */
   sub Mr, incAn        /* incAn = lda - M */
   shl $2, incAn        /* incAn = (lda-M)*sizeof */
   shl $2, lda          /* lda *= sizeof */
   sub $-128, pA0       /* code compaction by using signed 1-byte offsets */
   sub $-128, pX        /* code compaction by using signed 1-byte offsets */
   mov pX, pX0          /* save for restore after M loops */
   lea (lda, lda,2), lda3       /* lda3 = 3*lda */
   lea (lda, lda,4), lda5       /* lda5 = 5*lda */
   lea (lda3,lda,4), lda7       /* lda7 = 7*lda */
   lea (lda5,lda,4), lda9       /* lda9 = 9*lda */
   lea (lda,lda5,2), lda11      /* lda11 = 11*lda */
   add lda11, incAn             /* incAn = (12*lda-M)*sizeof */
   mov M, II

   ALIGN32
   LOOPN:
      movaps (pY), rY3          /* rY3 = {Y3, Y2, Y1, Y0} */
      pshufd $0x00, rY3, rY0    /* rY0 = {Y0, Y0, Y0, Y0} */
      movaps 16(pY), rY7        /* rY7 = {Y3, Y2, Y1, Y0} */
      pshufd $0x55, rY3, rY1    /* rY1 = {Y1, Y1, Y1, Y1} */
      movaps 32(pY), rY11       /* rY11= {Y3, Y2, Y1, Y0} */
      pshufd $0xAA, rY3, rY2    /* rY2 = {Y2, Y2, Y2, Y2} */
      pshufd $0xFF, rY3, rY3    /* rY3 = {Y3, Y3, Y3, Y3} */

      pshufd $0x00, rY7, rY4    /* rY0 = {Y0, Y0, Y0, Y0} */
      pshufd $0x55, rY7, rY5    /* rY1 = {Y1, Y1, Y1, Y1} */
      pshufd $0xAA, rY7, rY6    /* rY2 = {Y2, Y2, Y2, Y2} */
      pshufd $0xFF, rY7, rY7    /* rY3 = {Y3, Y3, Y3, Y3} */

      pshufd $0x00, rY11, rY8   /* rY8 = {Y0, Y0, Y0, Y0} */
      pshufd $0x55, rY11, rY9   /* rY9 = {Y1, Y1, Y1, Y1} */
      pshufd $0xAA, rY11, rY10  /* rY10= {Y2, Y2, Y2, Y2} */
      pshufd $0xFF, rY11, rY11  /* rY11= {Y3, Y3, Y3, Y3} */

      LOOPM:
         movapd 0-128(pX), rX0
         movapd rY0, ryt
         MOVA   0-128(pA0), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0)
         prefA(PFADIST+0(pA0))
         movapd rY1, ryt
         MOVA   0-128(pA0,lda), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0,lda)
         prefA(PFADIST+0(pA0,lda))
         movapd rY2, ryt
         MOVA   0-128(pA0,lda,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0,lda,2)
         prefA(PFADIST+0(pA0,lda,2))
         movapd rY3, ryt
         MOVA   0-128(pA0,lda3), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0,lda3)
         prefA(PFADIST+0(pA0,lda3))
         movapd rY4, ryt
         MOVA   0-128(pA0,lda,4), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0,lda,4)
         prefA(PFADIST+0(pA0,lda,4))
         movapd rY5, ryt
         MOVA   0-128(pA0,lda5), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0,lda5)
         prefA(PFADIST+0(pA0,lda5))
         movapd rY6, ryt
         MOVA   0-128(pA0,lda3,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0,lda3,2)
         prefA(PFADIST+0(pA0,lda3,2))
         movapd rY7, ryt
         MOVA   0-128(pA0,lda7), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0,lda7)
         prefA(PFADIST+0(pA0,lda7))
         movapd rY8, ryt
         MOVA   0-128(pA0,lda,8), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0,lda,8)
         prefA(PFADIST+0(pA0,lda,8))
         movapd rY9, ryt
         MOVA   0-128(pA0,lda9), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0,lda9)
         prefA(PFADIST+0(pA0,lda9))
         movapd rY10, ryt
         MOVA   0-128(pA0,lda5,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0,lda5,2)
         prefA(PFADIST+0(pA0,lda5,2))
         movapd rY11, ryt
         MOVA   0-128(pA0,lda11), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 0-128(pA0,lda11)
         prefA(PFADIST+0(pA0,lda11))

         movapd 16-128(pX), rX0
         movapd rY0, ryt
         MOVA   16-128(pA0), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0)
         movapd rY1, ryt
         MOVA   16-128(pA0,lda), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0,lda)
         movapd rY2, ryt
         MOVA   16-128(pA0,lda,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0,lda,2)
         movapd rY3, ryt
         MOVA   16-128(pA0,lda3), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0,lda3)
         movapd rY4, ryt
         MOVA   16-128(pA0,lda,4), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0,lda,4)
         movapd rY5, ryt
         MOVA   16-128(pA0,lda5), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0,lda5)
         movapd rY6, ryt
         MOVA   16-128(pA0,lda3,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0,lda3,2)
         movapd rY7, ryt
         MOVA   16-128(pA0,lda7), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0,lda7)
         movapd rY8, ryt
         MOVA   16-128(pA0,lda,8), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0,lda,8)
         movapd rY9, ryt
         MOVA   16-128(pA0,lda9), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0,lda9)
         movapd rY10, ryt
         MOVA   16-128(pA0,lda5,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0,lda5,2)
         movapd rY11, ryt
         MOVA   16-128(pA0,lda11), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 16-128(pA0,lda11)

         movapd 32-128(pX), rX0
         movapd rY0, ryt
         MOVA   32-128(pA0), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0)
         movapd rY1, ryt
         MOVA   32-128(pA0,lda), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0,lda)
         movapd rY2, ryt
         MOVA   32-128(pA0,lda,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0,lda,2)
         movapd rY3, ryt
         MOVA   32-128(pA0,lda3), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0,lda3)
         movapd rY4, ryt
         MOVA   32-128(pA0,lda,4), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0,lda,4)
         movapd rY5, ryt
         MOVA   32-128(pA0,lda5), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0,lda5)
         movapd rY6, ryt
         MOVA   32-128(pA0,lda3,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0,lda3,2)
         movapd rY7, ryt
         MOVA   32-128(pA0,lda7), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0,lda7)
         movapd rY8, ryt
         MOVA   32-128(pA0,lda,8), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0,lda,8)
         movapd rY9, ryt
         MOVA   32-128(pA0,lda9), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0,lda9)
         movapd rY10, ryt
         MOVA   32-128(pA0,lda5,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0,lda5,2)
         movapd rY11, ryt
         MOVA   32-128(pA0,lda11), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 32-128(pA0,lda11)

         movapd 48-128(pX), rX0
         movapd rY0, ryt
         MOVA   48-128(pA0), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0)
         movapd rY1, ryt
         MOVA   48-128(pA0,lda), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0,lda)
         movapd rY2, ryt
         MOVA   48-128(pA0,lda,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0,lda,2)
         movapd rY3, ryt
         MOVA   48-128(pA0,lda3), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0,lda3)
         movapd rY4, ryt
         MOVA   48-128(pA0,lda,4), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0,lda,4)
         movapd rY5, ryt
         MOVA   48-128(pA0,lda5), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0,lda5)
         movapd rY6, ryt
         MOVA   48-128(pA0,lda3,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0,lda3,2)
         movapd rY7, ryt
         MOVA   48-128(pA0,lda7), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0,lda7)
         movapd rY8, ryt
         MOVA   48-128(pA0,lda,8), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0,lda,8)
         movapd rY9, ryt
         MOVA   48-128(pA0,lda9), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0,lda9)
         movapd rY10, ryt
         MOVA   48-128(pA0,lda5,2), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0,lda5,2)
         movapd rY11, ryt
         MOVA   48-128(pA0,lda11), rA0
         mulpd rX0, ryt
         addpd ryt, rA0
         MOVA   rA0, 48-128(pA0,lda11)

         sub incAXm, pX
         sub incAXm, pA0
      sub incII, II
      jnz LOOPM

      cmp $0, Mr
      jz  MCLEANED

      mov Mr, II
      LOOPMCU:
         movsd -128(pX), rX0
         movsd rX0, rA0
         mulsd rY0, rA0
         addsd -128(pA0), rA0
         movsd rA0, -128(pA0)
         movsd rX0, rA0
         mulsd rY1, rA0
         addsd -128(pA0,lda), rA0
         movsd rA0, -128(pA0,lda)
         movsd rX0, rA0
         mulsd rY2, rA0
         addsd -128(pA0,lda,2), rA0
         movsd rA0, -128(pA0,lda,2)
         movsd rX0, rA0
         mulsd rY3, rA0
         addsd -128(pA0,lda3), rA0
         movsd rA0, -128(pA0,lda3)
         movsd rX0, rA0
         mulsd rY4, rA0
         addsd -128(pA0,lda,4), rA0
         movsd rA0, -128(pA0,lda,4)
         movsd rX0, rA0
         mulsd rY5, rA0
         addsd -128(pA0,lda5), rA0
         movsd rA0, -128(pA0,lda5)
         movsd rX0, rA0
         mulsd rY6, rA0
         addsd -128(pA0,lda3,2), rA0
         movsd rA0, -128(pA0,lda3,2)
         movsd rX0, rA0
         mulsd rY7, rA0
         addsd -128(pA0,lda7), rA0
         movsd rA0, -128(pA0,lda7)
         movsd rX0, rA0
         mulsd rY8, rA0
         addsd -128(pA0,lda,8), rA0
         movsd rA0, -128(pA0,lda,8)
         movsd rX0, rA0
         mulsd rY9, rA0
         addsd -128(pA0,lda9), rA0
         movsd rA0, -128(pA0,lda9)
         movsd rX0, rA0
         mulsd rY10, rA0
         addsd -128(pA0,lda5,2), rA0
         movsd rA0, -128(pA0,lda5,2)
         movsd rX0, rA0
         mulsd rY11, rA0
         addsd -128(pA0,lda11), rA0
         movsd rA0, -128(pA0,lda11)
         add $4, pX
         add $4, pA0
      dec II
      jnz LOOPMCU

MCLEANED:
      prefY(12*4+PFYDIST(pY))
      add $12*4, pY
      add incAn, pA0
      mov pX0, pX
      mov M, II
   sub $12, N
   jnz LOOPN
/*
 * EPILOGUE: restore registers and return
 */
   movq -8(%rsp), %rbp
   movq -16(%rsp), %rbx
   movq -24(%rsp), %r12
   movq -32(%rsp), %r13
   movq -40(%rsp), %r14
   movq -48(%rsp), %r15
   ret
