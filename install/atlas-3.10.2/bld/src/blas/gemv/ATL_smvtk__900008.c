#include "atlas_asm.h"
/*
 * This file does a 1x12 unrolled mvt_sse with these params:
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
#define pA0     %rdx
#define lda     %rax
#define pX      %r8
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
#define rY0     %xmm2
#define rY1     %xmm3
#define rY2     %xmm4
#define rY3     %xmm5
#define rY4     %xmm6
#define rY5     %xmm7
#define rY6     %xmm8
#define rY7     %xmm9
#define rY8     %xmm10
#define rY9     %xmm11
#define rY10     %xmm12
#define rY11     %xmm13

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
   #ifdef ATL_3DNow
      #define PFIY prefetchw
   #else
      #define PFIY prefetchnta
   #endif
#endif
#ifndef PFIX
   #define PFIX prefetcht0
#endif
#ifndef PFIA
   #define PFIA prefetchnta
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
/*
 *                      %rdi        %rsi           %rdx          %rcx
 * void ATL_UGEMV(ATL_CINT M, ATL_CINT N, const TYPE *A, ATL_CINT lda,
 *                          %r8      %r9
 *                const TYPE *X, TYPE *Y)
 */
.text
.text
.global ATL_asmdecor(ATL_UGEMV)
ALIGN64
ATL_asmdecor(ATL_UGEMV):

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
   mov  %rcx, lda       /* move lda to assigned register, rax */
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
      #ifdef BETA0
         xorpd rY0, rY0
         xorpd rY1, rY1
         xorpd rY2, rY2
         xorpd rY3, rY3
         xorpd rY4, rY4
         xorpd rY5, rY5
         xorpd rY6, rY6
         xorpd rY7, rY7
         xorpd rY8, rY8
         xorpd rY9, rY9
         xorpd rY10, rY10
         xorpd rY11, rY11
      #else
         movsd 0(pY), rY0
         movsd 4(pY), rY1
         movsd 8(pY), rY2
         movsd 12(pY), rY3
         movsd 16(pY), rY4
         movsd 20(pY), rY5
         movsd 24(pY), rY6
         movsd 28(pY), rY7
         movsd 32(pY), rY8
         movsd 36(pY), rY9
         movsd 40(pY), rY10
         movsd 44(pY), rY11
      #endif

      LOOPM:
         movapd 0-128(pX), rX0
         MOVA   0-128(pA0), rA0
         mulpd rX0, rA0
         addpd rA0, rY0
         prefA(PFADIST+0(pA0))

         MOVA   0-128(pA0,lda), rA0
         mulpd rX0, rA0
         addpd rA0, rY1
         prefA(PFADIST+0(pA0,lda))
         MOVA   0-128(pA0,lda,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY2
         prefA(PFADIST+0(pA0,lda,2))
         MOVA   0-128(pA0,lda3), rA0
         mulpd rX0, rA0
         addpd rA0, rY3
         prefA(PFADIST+0(pA0,lda3))
         MOVA   0-128(pA0,lda,4), rA0
         mulpd rX0, rA0
         addpd rA0, rY4
         prefA(PFADIST+0(pA0,lda,4))
         MOVA   0-128(pA0,lda5), rA0
         mulpd rX0, rA0
         addpd rA0, rY5
         prefA(PFADIST+0(pA0,lda5))
         MOVA   0-128(pA0,lda3,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY6
         prefA(PFADIST+0(pA0,lda3,2))
         MOVA   0-128(pA0,lda7), rA0
         mulpd rX0, rA0
         addpd rA0, rY7
         prefA(PFADIST+0(pA0,lda7))
         MOVA   0-128(pA0,lda,8), rA0
         mulpd rX0, rA0
         addpd rA0, rY8
         prefA(PFADIST+0(pA0,lda,8))
         MOVA   0-128(pA0,lda9), rA0
         mulpd rX0, rA0
         addpd rA0, rY9
         prefA(PFADIST+0(pA0,lda9))
         MOVA   0-128(pA0,lda5,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY10
         prefA(PFADIST+0(pA0,lda5,2))
         MOVA   0-128(pA0,lda11), rA0
         mulpd rX0, rA0
         addpd rA0, rY11
         prefA(PFADIST+0(pA0,lda11))

         movapd 16-128(pX), rX0
         MOVA   16-128(pA0), rA0
         mulpd rX0, rA0
         addpd rA0, rY0

         MOVA   16-128(pA0,lda), rA0
         mulpd rX0, rA0
         addpd rA0, rY1
         MOVA   16-128(pA0,lda,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY2
         MOVA   16-128(pA0,lda3), rA0
         mulpd rX0, rA0
         addpd rA0, rY3
         MOVA   16-128(pA0,lda,4), rA0
         mulpd rX0, rA0
         addpd rA0, rY4
         MOVA   16-128(pA0,lda5), rA0
         mulpd rX0, rA0
         addpd rA0, rY5
         MOVA   16-128(pA0,lda3,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY6
         MOVA   16-128(pA0,lda7), rA0
         mulpd rX0, rA0
         addpd rA0, rY7
         MOVA   16-128(pA0,lda,8), rA0
         mulpd rX0, rA0
         addpd rA0, rY8
         MOVA   16-128(pA0,lda9), rA0
         mulpd rX0, rA0
         addpd rA0, rY9
         MOVA   16-128(pA0,lda5,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY10
         MOVA   16-128(pA0,lda11), rA0
         mulpd rX0, rA0
         addpd rA0, rY11

         movapd 32-128(pX), rX0
         MOVA   32-128(pA0), rA0
         mulpd rX0, rA0
         addpd rA0, rY0

         MOVA   32-128(pA0,lda), rA0
         mulpd rX0, rA0
         addpd rA0, rY1
         MOVA   32-128(pA0,lda,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY2
         MOVA   32-128(pA0,lda3), rA0
         mulpd rX0, rA0
         addpd rA0, rY3
         MOVA   32-128(pA0,lda,4), rA0
         mulpd rX0, rA0
         addpd rA0, rY4
         MOVA   32-128(pA0,lda5), rA0
         mulpd rX0, rA0
         addpd rA0, rY5
         MOVA   32-128(pA0,lda3,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY6
         MOVA   32-128(pA0,lda7), rA0
         mulpd rX0, rA0
         addpd rA0, rY7
         MOVA   32-128(pA0,lda,8), rA0
         mulpd rX0, rA0
         addpd rA0, rY8
         MOVA   32-128(pA0,lda9), rA0
         mulpd rX0, rA0
         addpd rA0, rY9
         MOVA   32-128(pA0,lda5,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY10
         MOVA   32-128(pA0,lda11), rA0
         mulpd rX0, rA0
         addpd rA0, rY11

         movapd 48-128(pX), rX0
         MOVA   48-128(pA0), rA0
         mulpd rX0, rA0
         addpd rA0, rY0

         MOVA   48-128(pA0,lda), rA0
         mulpd rX0, rA0
         addpd rA0, rY1
         MOVA   48-128(pA0,lda,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY2
         MOVA   48-128(pA0,lda3), rA0
         mulpd rX0, rA0
         addpd rA0, rY3
         MOVA   48-128(pA0,lda,4), rA0
         mulpd rX0, rA0
         addpd rA0, rY4
         MOVA   48-128(pA0,lda5), rA0
         mulpd rX0, rA0
         addpd rA0, rY5
         MOVA   48-128(pA0,lda3,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY6
         MOVA   48-128(pA0,lda7), rA0
         mulpd rX0, rA0
         addpd rA0, rY7
         MOVA   48-128(pA0,lda,8), rA0
         mulpd rX0, rA0
         addpd rA0, rY8
         MOVA   48-128(pA0,lda9), rA0
         mulpd rX0, rA0
         addpd rA0, rY9
         MOVA   48-128(pA0,lda5,2), rA0
         mulpd rX0, rA0
         addpd rA0, rY10
         MOVA   48-128(pA0,lda11), rA0
         mulpd rX0, rA0
         addpd rA0, rY11

         sub incAXm, pX
         sub incAXm, pA0
      sub incII, II
      jnz LOOPM

      cmp $0, Mr
      jz  MCLEANED

      mov Mr, II
      LOOPMCU:
         movsd -128(pX), rX0
         movsd -128(pA0), rA0
         mulsd rX0, rA0
         addsd rA0, rY0
         movsd -128(pA0,lda), rA0
         mulsd rX0, rA0
         addsd rA0, rY1
         movsd -128(pA0,lda,2), rA0
         mulsd rX0, rA0
         addsd rA0, rY2
         movsd -128(pA0,lda3), rA0
         mulsd rX0, rA0
         addsd rA0, rY3
         movsd -128(pA0,lda,4), rA0
         mulsd rX0, rA0
         addsd rA0, rY4
         movsd -128(pA0,lda5), rA0
         mulsd rX0, rA0
         addsd rA0, rY5
         movsd -128(pA0,lda3,2), rA0
         mulsd rX0, rA0
         addsd rA0, rY6
         movsd -128(pA0,lda7), rA0
         mulsd rX0, rA0
         addsd rA0, rY7
         movsd -128(pA0,lda,8), rA0
         mulsd rX0, rA0
         addsd rA0, rY8
         movsd -128(pA0,lda9), rA0
         mulsd rX0, rA0
         addsd rA0, rY9
         movsd -128(pA0,lda5,2), rA0
         mulsd rX0, rA0
         addsd rA0, rY10
         movsd -128(pA0,lda11), rA0
         mulsd rX0, rA0
         addsd rA0, rY11
         add $4, pX
         add $4, pA0
      dec II
      jnz LOOPMCU

MCLEANED:
      haddps rY1, rY0
      haddps rY3, rY2
      haddps rY2, rY0
      movaps rY0, 0(pY)
      haddps rY5, rY4
      haddps rY7, rY6
      haddps rY6, rY4
      movaps rY4, 16(pY)
      haddps rY9, rY8
      haddps rY11, rY10
      haddps rY10, rY8
      movaps rY8, 32(pY)
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
