#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_sJIK48x48x48NN0x0x0_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=N, TB=N, MB=48, NB=48, KB=48, 
 * lda=0, ldb=0, ldc=0, mu=12, nu=1, ku=48, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const float *stM = A + 48;
   const float *stN = B + (((ldb) << 5)+((ldb) << 4));
   const int incAk = (lda);
   const int incAm = 12 - (((lda) << 5)+((lda) << 4)), incAn = -48;
   #define incBk 48
   const int incBm = -48, incBn = (ldb);
   #define incCm 12
   const int incCn = (ldc) - 48;
   float *pC0=C;
   const float *pA0=A;
   const float *pB0=B;
   register int k;
   register float rA0, rA1, rA2, rA3, rA4, rA5, rA6, rA7, rA8, rA9, rA10, rA11;
   register float rB0;
   register float rC0_0, rC1_0, rC2_0, rC3_0, rC4_0, rC5_0, rC6_0, rC7_0, rC8_0, rC9_0, rC10_0, rC11_0;
   do /* N-loop */
   {
      do /* M-loop */
      {
         rA0 = beta;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         rC1_0 = pC0[1];
         rC1_0 *= rA0;
         rC2_0 = pC0[2];
         rC2_0 *= rA0;
         rC3_0 = pC0[3];
         rC3_0 *= rA0;
         rC4_0 = pC0[4];
         rC4_0 *= rA0;
         rC5_0 = pC0[5];
         rC5_0 *= rA0;
         rC6_0 = pC0[6];
         rC6_0 *= rA0;
         rC7_0 = pC0[7];
         rC7_0 *= rA0;
         rC8_0 = pC0[8];
         rC8_0 *= rA0;
         rC9_0 = pC0[9];
         rC9_0 *= rA0;
         rC10_0 = pC0[10];
         rC10_0 *= rA0;
         rC11_0 = pC0[11];
         rC11_0 *= rA0;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[1];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[2];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[3];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[4];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[5];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[6];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[7];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[8];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[9];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[10];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[11];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[12];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[13];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[14];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[15];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[16];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[17];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[18];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[19];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[20];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[21];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[22];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[23];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[24];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[25];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[26];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[27];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[28];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[29];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[30];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[31];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[32];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[33];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[34];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[35];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[36];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[37];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[38];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[39];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[40];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[41];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[42];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[43];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[44];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[45];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[46];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[47];
         rA1 = pA0[1];
         rA2 = pA0[2];
         rA3 = pA0[3];
         rA4 = pA0[4];
         rA5 = pA0[5];
         rA6 = pA0[6];
         rA7 = pA0[7];
         rA8 = pA0[8];
         rA9 = pA0[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA0[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA0[11];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC4_0 += rA4 * rB0;
         rC5_0 += rA5 * rB0;
         rC6_0 += rA6 * rB0;
         rC7_0 += rA7 * rB0;
         rC8_0 += rA8 * rB0;
         rC9_0 += rA9 * rB0;
         rC10_0 += rA10 * rB0;
         rC11_0 += rA11 * rB0;
         pA0 += incAk;
         pB0 += incBk;
         *pC0 = rC0_0;
         pC0[1] = rC1_0;
         pC0[2] = rC2_0;
         pC0[3] = rC3_0;
         pC0[4] = rC4_0;
         pC0[5] = rC5_0;
         pC0[6] = rC6_0;
         pC0[7] = rC7_0;
         pC0[8] = rC8_0;
         pC0[9] = rC9_0;
         pC0[10] = rC10_0;
         pC0[11] = rC11_0;
         pC0 += incCm;
         pA0 += incAm;
         pB0 += incBm;
      }
      while(pA0 != stM);
      pC0 += incCn;
      pA0 += incAn;
      pB0 += incBn;
   }
   while(pB0 != stN);
}
#ifdef incAm
   #undef incAm
#endif
#ifdef incAn
   #undef incAn
#endif
#ifdef incAk
   #undef incAk
#endif
#ifdef incBm
   #undef incBm
#endif
#ifdef incBn
   #undef incBn
#endif
#ifdef incBk
   #undef incBk
#endif
#ifdef incCm
   #undef incCm
#endif
#ifdef incCn
   #undef incCn
#endif
#ifdef incCk
   #undef incCk
#endif
#ifdef Mb
   #undef Mb
#endif
#ifdef Nb
   #undef Nb
#endif
#ifdef Kb
   #undef Kb
#endif
