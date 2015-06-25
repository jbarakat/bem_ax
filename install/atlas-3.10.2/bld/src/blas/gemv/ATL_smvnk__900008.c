#include "atlas_asm.h"
/*
 * This file does a 1x12 unrolled mvn_sse with these params:
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
#define pY0     %r11
#define Mr      %rcx
#define incAn   %r15
#define incII   $16*1
#define lda3    %r12
#define lda5    %r13
#define lda7    %rbp
#define lda9    %r14
#define lda11   %r10
#define incAYm $-64
/*
 * SSE register assignment
 */
#define rA0     %xmm0
#define rY      %xmm1
#define rX0     %xmm2
#define rX1     %xmm3
#define rX2     %xmm4
#define rX3     %xmm5
#define rX4     %xmm6
#define rX5     %xmm7
#define rX6     %xmm8
#define rX7     %xmm9
#define rX8     %xmm10
#define rX9     %xmm11
#define rX10     %xmm12
#define rX11     %xmm13
/*
 * macros
 */
#ifndef MOVA
   #define MOVA movaps
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
      #define PFIY prefetcht0
   #endif
#endif
#ifndef PFIX
   #define PFIX prefetchnta
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
   sub $-128, pY        /* code compaction by using signed 1-byte offsets */
   mov pY, pY0          /* save for restore after M loops */
   lea (lda, lda,2), lda3       /* lda3 = 3*lda */
   lea (lda, lda,4), lda5       /* lda5 = 5*lda */
   lea (lda3,lda,4), lda7       /* lda7 = 7*lda */
   lea (lda5,lda,4), lda9       /* lda9 = 9*lda */
   lea (lda,lda5,2), lda11      /* lda11 = 11*lda */
   add lda11, incAn             /* incAn = (12*lda-M)*sizeof */
   mov M, II
/*
 * Zero Y if beta = 0;  Has error if there is Mr and/or M isn't mul of veclen
 */
   #ifdef BETA0
       add Mr, II
      shr $2, II
      xorpd rY, rY
      LOOPZERO:
         movapd rY, -128(pY)
         add $16, pY
      dec II
      jnz LOOPZERO
      lea (M, Mr), II
      bt $1, II
      jnc DONE_ZERO_2
      movlps rY, -128(pY)
      add $8, pY
DONE_ZERO_2:
      bt $0, II
      jnc DONE_ZERO_CLEAN
      movsd rY, -128(pY)
DONE_ZERO_CLEAN:
      mov pY0, pY
      mov M, II
   #endif

   ALIGN32
   LOOPN:
      movaps (pX), rX3          /* rX3 = {X3, X2, X1, X0} */
      pshufd $0x00, rX3, rX0    /* rX0 = {X0, X0, X0, X0} */
      movaps 16(pX), rX7        /* rX7 = {X3, X2, X1, X0} */
      pshufd $0x55, rX3, rX1    /* rX1 = {X1, X1, X1, X1} */
      movaps 32(pX), rX11       /* rX11= {X3, X2, X1, X0} */
      pshufd $0xAA, rX3, rX2    /* rX2 = {X2, X2, X2, X2} */
      pshufd $0xFF, rX3, rX3    /* rX3 = {X3, X3, X3, X3} */

      pshufd $0x00, rX7, rX4    /* rX0 = {X0, X0, X0, X0} */
      pshufd $0x55, rX7, rX5    /* rX1 = {X1, X1, X1, X1} */
      pshufd $0xAA, rX7, rX6    /* rX2 = {X2, X2, X2, X2} */
      pshufd $0xFF, rX7, rX7    /* rX3 = {X3, X3, X3, X3} */

      pshufd $0x00, rX11, rX8   /* rX8 = {X0, X0, X0, X0} */
      pshufd $0x55, rX11, rX9   /* rX9 = {X1, X1, X1, X1} */
      pshufd $0xAA, rX11, rX10  /* rX10= {X2, X2, X2, X2} */
      pshufd $0xFF, rX11, rX11  /* rX11= {X3, X3, X3, X3} */

      LOOPM:
         MOVA   0-128(pA0), rY
         mulpd rX0, rY
         addpd 0-128(pY), rY
         prefA(PFADIST+0(pA0))

         MOVA   0-128(pA0,lda), rA0
         mulpd rX1, rA0
         addpd rA0, rY
         prefA(PFADIST+0(pA0,lda))
         MOVA   0-128(pA0,lda,2), rA0
         mulpd rX2, rA0
         addpd rA0, rY
         prefA(PFADIST+0(pA0,lda,2))
         MOVA   0-128(pA0,lda3), rA0
         mulpd rX3, rA0
         addpd rA0, rY
         prefA(PFADIST+0(pA0,lda3))
         MOVA   0-128(pA0,lda,4), rA0
         mulpd rX4, rA0
         addpd rA0, rY
         prefA(PFADIST+0(pA0,lda,4))
         MOVA   0-128(pA0,lda5), rA0
         mulpd rX5, rA0
         addpd rA0, rY
         prefA(PFADIST+0(pA0,lda5))
         MOVA   0-128(pA0,lda3,2), rA0
         mulpd rX6, rA0
         addpd rA0, rY
         prefA(PFADIST+0(pA0,lda3,2))
         MOVA   0-128(pA0,lda7), rA0
         mulpd rX7, rA0
         addpd rA0, rY
         prefA(PFADIST+0(pA0,lda7))
         MOVA   0-128(pA0,lda,8), rA0
         mulpd rX8, rA0
         addpd rA0, rY
         prefA(PFADIST+0(pA0,lda,8))
         MOVA   0-128(pA0,lda9), rA0
         mulpd rX9, rA0
         addpd rA0, rY
         prefA(PFADIST+0(pA0,lda9))
         MOVA   0-128(pA0,lda5,2), rA0
         mulpd rX10, rA0
         addpd rA0, rY
         prefA(PFADIST+0(pA0,lda5,2))
         MOVA   0-128(pA0,lda11), rA0
         mulpd rX11, rA0
         addpd rA0, rY
         prefA(PFADIST+0(pA0,lda11))
         movapd rY, 0-128(pY)

         MOVA   16-128(pA0), rY
         mulpd rX0, rY
         addpd 16-128(pY), rY

         MOVA   16-128(pA0,lda), rA0
         mulpd rX1, rA0
         addpd rA0, rY
         MOVA   16-128(pA0,lda,2), rA0
         mulpd rX2, rA0
         addpd rA0, rY
         MOVA   16-128(pA0,lda3), rA0
         mulpd rX3, rA0
         addpd rA0, rY
         MOVA   16-128(pA0,lda,4), rA0
         mulpd rX4, rA0
         addpd rA0, rY
         MOVA   16-128(pA0,lda5), rA0
         mulpd rX5, rA0
         addpd rA0, rY
         MOVA   16-128(pA0,lda3,2), rA0
         mulpd rX6, rA0
         addpd rA0, rY
         MOVA   16-128(pA0,lda7), rA0
         mulpd rX7, rA0
         addpd rA0, rY
         MOVA   16-128(pA0,lda,8), rA0
         mulpd rX8, rA0
         addpd rA0, rY
         MOVA   16-128(pA0,lda9), rA0
         mulpd rX9, rA0
         addpd rA0, rY
         MOVA   16-128(pA0,lda5,2), rA0
         mulpd rX10, rA0
         addpd rA0, rY
         MOVA   16-128(pA0,lda11), rA0
         mulpd rX11, rA0
         addpd rA0, rY
         movapd rY, 16-128(pY)

         MOVA   32-128(pA0), rY
         mulpd rX0, rY
         addpd 32-128(pY), rY

         MOVA   32-128(pA0,lda), rA0
         mulpd rX1, rA0
         addpd rA0, rY
         MOVA   32-128(pA0,lda,2), rA0
         mulpd rX2, rA0
         addpd rA0, rY
         MOVA   32-128(pA0,lda3), rA0
         mulpd rX3, rA0
         addpd rA0, rY
         MOVA   32-128(pA0,lda,4), rA0
         mulpd rX4, rA0
         addpd rA0, rY
         MOVA   32-128(pA0,lda5), rA0
         mulpd rX5, rA0
         addpd rA0, rY
         MOVA   32-128(pA0,lda3,2), rA0
         mulpd rX6, rA0
         addpd rA0, rY
         MOVA   32-128(pA0,lda7), rA0
         mulpd rX7, rA0
         addpd rA0, rY
         MOVA   32-128(pA0,lda,8), rA0
         mulpd rX8, rA0
         addpd rA0, rY
         MOVA   32-128(pA0,lda9), rA0
         mulpd rX9, rA0
         addpd rA0, rY
         MOVA   32-128(pA0,lda5,2), rA0
         mulpd rX10, rA0
         addpd rA0, rY
         MOVA   32-128(pA0,lda11), rA0
         mulpd rX11, rA0
         addpd rA0, rY
         movapd rY, 32-128(pY)

         MOVA   48-128(pA0), rY
         mulpd rX0, rY
         addpd 48-128(pY), rY

         MOVA   48-128(pA0,lda), rA0
         mulpd rX1, rA0
         addpd rA0, rY
         MOVA   48-128(pA0,lda,2), rA0
         mulpd rX2, rA0
         addpd rA0, rY
         MOVA   48-128(pA0,lda3), rA0
         mulpd rX3, rA0
         addpd rA0, rY
         MOVA   48-128(pA0,lda,4), rA0
         mulpd rX4, rA0
         addpd rA0, rY
         MOVA   48-128(pA0,lda5), rA0
         mulpd rX5, rA0
         addpd rA0, rY
         MOVA   48-128(pA0,lda3,2), rA0
         mulpd rX6, rA0
         addpd rA0, rY
         MOVA   48-128(pA0,lda7), rA0
         mulpd rX7, rA0
         addpd rA0, rY
         MOVA   48-128(pA0,lda,8), rA0
         mulpd rX8, rA0
         addpd rA0, rY
         MOVA   48-128(pA0,lda9), rA0
         mulpd rX9, rA0
         addpd rA0, rY
         MOVA   48-128(pA0,lda5,2), rA0
         mulpd rX10, rA0
         addpd rA0, rY
         MOVA   48-128(pA0,lda11), rA0
         mulpd rX11, rA0
         addpd rA0, rY
         movapd rY, 48-128(pY)

         sub incAYm, pY
         sub incAYm, pA0
      sub incII, II
      jnz LOOPM

      cmp $0, Mr
      jz  MCLEANED

      mov Mr, II
      LOOPMCU:
         movsd -128(pY), rY
         movsd -128(pA0), rA0
         mulsd rX0, rA0
         addsd rA0, rY
         movsd  -128(pA0,lda), rA0
         mulsd rX1, rA0
         addsd rA0, rY
         movsd  -128(pA0,lda,2), rA0
         mulsd rX2, rA0
         addsd rA0, rY
         movsd  -128(pA0,lda3), rA0
         mulsd rX3, rA0
         addsd rA0, rY
         movsd  -128(pA0,lda,4), rA0
         mulsd rX4, rA0
         addsd rA0, rY
         movsd  -128(pA0,lda5), rA0
         mulsd rX5, rA0
         addsd rA0, rY
         movsd  -128(pA0,lda3,2), rA0
         mulsd rX6, rA0
         addsd rA0, rY
         movsd  -128(pA0,lda7), rA0
         mulsd rX7, rA0
         addsd rA0, rY
         movsd  -128(pA0,lda,8), rA0
         mulsd rX8, rA0
         addsd rA0, rY
         movsd  -128(pA0,lda9), rA0
         mulsd rX9, rA0
         addsd rA0, rY
         movsd  -128(pA0,lda5,2), rA0
         mulsd rX10, rA0
         addsd rA0, rY
         movsd  -128(pA0,lda11), rA0
         mulsd rX11, rA0
         addsd rA0, rY
         movsd rY, -128(pY)
         add $4, pY
         add $4, pA0
      dec II
      jnz LOOPMCU

MCLEANED:
      prefX(12*4+PFXDIST(pX))
      add $12*4, pX
      add incAn, pA0
      mov pY0, pY
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
