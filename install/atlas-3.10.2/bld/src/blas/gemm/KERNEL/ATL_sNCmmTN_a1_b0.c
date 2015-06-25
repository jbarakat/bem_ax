#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_sJIK48x48x48TN0x0x0_a1_b0
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=48, NB=48, KB=48, 
 * lda=0, ldb=0, ldc=0, mu=12, nu=1, ku=48, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const float *stM = A + (((lda) << 5)+((lda) << 4));
   const float *stN = B + (((ldb) << 5)+((ldb) << 4));
   #define incAk 48
   const int incAm = ((((lda) << 3)+((lda) << 2)) - 48), incAn = -(((lda) << 5)+((lda) << 4));
   #define incBk 48
   const int incBm = -48, incBn = (ldb);
   #define incCm 12
   const int incCn = (ldc) - 48;
   float *pC0=C;
   const float *pA0=A, *pA1=pA0+(lda), *pA2=pA1+(lda), *pA3=pA2+(lda), *pA4=pA3+(lda), *pA5=pA4+(lda), *pA6=pA5+(lda), *pA7=pA6+(lda), *pA8=pA7+(lda), *pA9=pA8+(lda), *pA10=pA9+(lda), *pA11=pA10+(lda);
   const float *pB0=B;
   register int k;
   register float rA0, rA1, rA2, rA3, rA4, rA5, rA6, rA7, rA8, rA9, rA10, rA11;
   register float rB0;
   register float rC0_0, rC1_0, rC2_0, rC3_0, rC4_0, rC5_0, rC6_0, rC7_0, rC8_0, rC9_0, rC10_0, rC11_0;
   do /* N-loop */
   {
      do /* M-loop */
      {
/*
 *       Feeble prefetch of C
 */
         rC0_0 = *pC0;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = *pA1;
         rA2 = *pA2;
         rA3 = *pA3;
         rA4 = *pA4;
         rA5 = *pA5;
         rA6 = *pA6;
         rA7 = *pA7;
         rA8 = *pA8;
         rA9 = *pA9;
         rC0_0 = rA0 * rB0;
         rA10 = *pA10;
         rC1_0 = rA1 * rB0;
         rA11 = *pA11;
         rC2_0 = rA2 * rB0;
         rC3_0 = rA3 * rB0;
         rC4_0 = rA4 * rB0;
         rC5_0 = rA5 * rB0;
         rC6_0 = rA6 * rB0;
         rC7_0 = rA7 * rB0;
         rC8_0 = rA8 * rB0;
         rC9_0 = rA9 * rB0;
         rC10_0 = rA10 * rB0;
         rC11_0 = rA11 * rB0;
         rA0 = pA0[1];
         rB0 = pB0[1];
         rA1 = pA1[1];
         rA2 = pA2[1];
         rA3 = pA3[1];
         rA4 = pA4[1];
         rA5 = pA5[1];
         rA6 = pA6[1];
         rA7 = pA7[1];
         rA8 = pA8[1];
         rA9 = pA9[1];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[1];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[1];
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
         rA0 = pA0[2];
         rB0 = pB0[2];
         rA1 = pA1[2];
         rA2 = pA2[2];
         rA3 = pA3[2];
         rA4 = pA4[2];
         rA5 = pA5[2];
         rA6 = pA6[2];
         rA7 = pA7[2];
         rA8 = pA8[2];
         rA9 = pA9[2];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[2];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[2];
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
         rA0 = pA0[3];
         rB0 = pB0[3];
         rA1 = pA1[3];
         rA2 = pA2[3];
         rA3 = pA3[3];
         rA4 = pA4[3];
         rA5 = pA5[3];
         rA6 = pA6[3];
         rA7 = pA7[3];
         rA8 = pA8[3];
         rA9 = pA9[3];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[3];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[3];
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
         rA0 = pA0[4];
         rB0 = pB0[4];
         rA1 = pA1[4];
         rA2 = pA2[4];
         rA3 = pA3[4];
         rA4 = pA4[4];
         rA5 = pA5[4];
         rA6 = pA6[4];
         rA7 = pA7[4];
         rA8 = pA8[4];
         rA9 = pA9[4];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[4];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[4];
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
         rA0 = pA0[5];
         rB0 = pB0[5];
         rA1 = pA1[5];
         rA2 = pA2[5];
         rA3 = pA3[5];
         rA4 = pA4[5];
         rA5 = pA5[5];
         rA6 = pA6[5];
         rA7 = pA7[5];
         rA8 = pA8[5];
         rA9 = pA9[5];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[5];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[5];
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
         rA0 = pA0[6];
         rB0 = pB0[6];
         rA1 = pA1[6];
         rA2 = pA2[6];
         rA3 = pA3[6];
         rA4 = pA4[6];
         rA5 = pA5[6];
         rA6 = pA6[6];
         rA7 = pA7[6];
         rA8 = pA8[6];
         rA9 = pA9[6];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[6];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[6];
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
         rA0 = pA0[7];
         rB0 = pB0[7];
         rA1 = pA1[7];
         rA2 = pA2[7];
         rA3 = pA3[7];
         rA4 = pA4[7];
         rA5 = pA5[7];
         rA6 = pA6[7];
         rA7 = pA7[7];
         rA8 = pA8[7];
         rA9 = pA9[7];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[7];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[7];
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
         rA0 = pA0[8];
         rB0 = pB0[8];
         rA1 = pA1[8];
         rA2 = pA2[8];
         rA3 = pA3[8];
         rA4 = pA4[8];
         rA5 = pA5[8];
         rA6 = pA6[8];
         rA7 = pA7[8];
         rA8 = pA8[8];
         rA9 = pA9[8];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[8];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[8];
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
         rA0 = pA0[9];
         rB0 = pB0[9];
         rA1 = pA1[9];
         rA2 = pA2[9];
         rA3 = pA3[9];
         rA4 = pA4[9];
         rA5 = pA5[9];
         rA6 = pA6[9];
         rA7 = pA7[9];
         rA8 = pA8[9];
         rA9 = pA9[9];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[9];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[9];
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
         rA0 = pA0[10];
         rB0 = pB0[10];
         rA1 = pA1[10];
         rA2 = pA2[10];
         rA3 = pA3[10];
         rA4 = pA4[10];
         rA5 = pA5[10];
         rA6 = pA6[10];
         rA7 = pA7[10];
         rA8 = pA8[10];
         rA9 = pA9[10];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[10];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[10];
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
         rA0 = pA0[11];
         rB0 = pB0[11];
         rA1 = pA1[11];
         rA2 = pA2[11];
         rA3 = pA3[11];
         rA4 = pA4[11];
         rA5 = pA5[11];
         rA6 = pA6[11];
         rA7 = pA7[11];
         rA8 = pA8[11];
         rA9 = pA9[11];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[11];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[11];
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
         rA0 = pA0[12];
         rB0 = pB0[12];
         rA1 = pA1[12];
         rA2 = pA2[12];
         rA3 = pA3[12];
         rA4 = pA4[12];
         rA5 = pA5[12];
         rA6 = pA6[12];
         rA7 = pA7[12];
         rA8 = pA8[12];
         rA9 = pA9[12];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[12];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[12];
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
         rA0 = pA0[13];
         rB0 = pB0[13];
         rA1 = pA1[13];
         rA2 = pA2[13];
         rA3 = pA3[13];
         rA4 = pA4[13];
         rA5 = pA5[13];
         rA6 = pA6[13];
         rA7 = pA7[13];
         rA8 = pA8[13];
         rA9 = pA9[13];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[13];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[13];
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
         rA0 = pA0[14];
         rB0 = pB0[14];
         rA1 = pA1[14];
         rA2 = pA2[14];
         rA3 = pA3[14];
         rA4 = pA4[14];
         rA5 = pA5[14];
         rA6 = pA6[14];
         rA7 = pA7[14];
         rA8 = pA8[14];
         rA9 = pA9[14];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[14];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[14];
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
         rA0 = pA0[15];
         rB0 = pB0[15];
         rA1 = pA1[15];
         rA2 = pA2[15];
         rA3 = pA3[15];
         rA4 = pA4[15];
         rA5 = pA5[15];
         rA6 = pA6[15];
         rA7 = pA7[15];
         rA8 = pA8[15];
         rA9 = pA9[15];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[15];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[15];
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
         rA0 = pA0[16];
         rB0 = pB0[16];
         rA1 = pA1[16];
         rA2 = pA2[16];
         rA3 = pA3[16];
         rA4 = pA4[16];
         rA5 = pA5[16];
         rA6 = pA6[16];
         rA7 = pA7[16];
         rA8 = pA8[16];
         rA9 = pA9[16];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[16];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[16];
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
         rA0 = pA0[17];
         rB0 = pB0[17];
         rA1 = pA1[17];
         rA2 = pA2[17];
         rA3 = pA3[17];
         rA4 = pA4[17];
         rA5 = pA5[17];
         rA6 = pA6[17];
         rA7 = pA7[17];
         rA8 = pA8[17];
         rA9 = pA9[17];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[17];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[17];
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
         rA0 = pA0[18];
         rB0 = pB0[18];
         rA1 = pA1[18];
         rA2 = pA2[18];
         rA3 = pA3[18];
         rA4 = pA4[18];
         rA5 = pA5[18];
         rA6 = pA6[18];
         rA7 = pA7[18];
         rA8 = pA8[18];
         rA9 = pA9[18];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[18];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[18];
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
         rA0 = pA0[19];
         rB0 = pB0[19];
         rA1 = pA1[19];
         rA2 = pA2[19];
         rA3 = pA3[19];
         rA4 = pA4[19];
         rA5 = pA5[19];
         rA6 = pA6[19];
         rA7 = pA7[19];
         rA8 = pA8[19];
         rA9 = pA9[19];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[19];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[19];
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
         rA0 = pA0[20];
         rB0 = pB0[20];
         rA1 = pA1[20];
         rA2 = pA2[20];
         rA3 = pA3[20];
         rA4 = pA4[20];
         rA5 = pA5[20];
         rA6 = pA6[20];
         rA7 = pA7[20];
         rA8 = pA8[20];
         rA9 = pA9[20];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[20];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[20];
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
         rA0 = pA0[21];
         rB0 = pB0[21];
         rA1 = pA1[21];
         rA2 = pA2[21];
         rA3 = pA3[21];
         rA4 = pA4[21];
         rA5 = pA5[21];
         rA6 = pA6[21];
         rA7 = pA7[21];
         rA8 = pA8[21];
         rA9 = pA9[21];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[21];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[21];
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
         rA0 = pA0[22];
         rB0 = pB0[22];
         rA1 = pA1[22];
         rA2 = pA2[22];
         rA3 = pA3[22];
         rA4 = pA4[22];
         rA5 = pA5[22];
         rA6 = pA6[22];
         rA7 = pA7[22];
         rA8 = pA8[22];
         rA9 = pA9[22];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[22];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[22];
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
         rA0 = pA0[23];
         rB0 = pB0[23];
         rA1 = pA1[23];
         rA2 = pA2[23];
         rA3 = pA3[23];
         rA4 = pA4[23];
         rA5 = pA5[23];
         rA6 = pA6[23];
         rA7 = pA7[23];
         rA8 = pA8[23];
         rA9 = pA9[23];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[23];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[23];
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
         rA0 = pA0[24];
         rB0 = pB0[24];
         rA1 = pA1[24];
         rA2 = pA2[24];
         rA3 = pA3[24];
         rA4 = pA4[24];
         rA5 = pA5[24];
         rA6 = pA6[24];
         rA7 = pA7[24];
         rA8 = pA8[24];
         rA9 = pA9[24];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[24];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[24];
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
         rA0 = pA0[25];
         rB0 = pB0[25];
         rA1 = pA1[25];
         rA2 = pA2[25];
         rA3 = pA3[25];
         rA4 = pA4[25];
         rA5 = pA5[25];
         rA6 = pA6[25];
         rA7 = pA7[25];
         rA8 = pA8[25];
         rA9 = pA9[25];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[25];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[25];
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
         rA0 = pA0[26];
         rB0 = pB0[26];
         rA1 = pA1[26];
         rA2 = pA2[26];
         rA3 = pA3[26];
         rA4 = pA4[26];
         rA5 = pA5[26];
         rA6 = pA6[26];
         rA7 = pA7[26];
         rA8 = pA8[26];
         rA9 = pA9[26];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[26];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[26];
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
         rA0 = pA0[27];
         rB0 = pB0[27];
         rA1 = pA1[27];
         rA2 = pA2[27];
         rA3 = pA3[27];
         rA4 = pA4[27];
         rA5 = pA5[27];
         rA6 = pA6[27];
         rA7 = pA7[27];
         rA8 = pA8[27];
         rA9 = pA9[27];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[27];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[27];
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
         rA0 = pA0[28];
         rB0 = pB0[28];
         rA1 = pA1[28];
         rA2 = pA2[28];
         rA3 = pA3[28];
         rA4 = pA4[28];
         rA5 = pA5[28];
         rA6 = pA6[28];
         rA7 = pA7[28];
         rA8 = pA8[28];
         rA9 = pA9[28];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[28];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[28];
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
         rA0 = pA0[29];
         rB0 = pB0[29];
         rA1 = pA1[29];
         rA2 = pA2[29];
         rA3 = pA3[29];
         rA4 = pA4[29];
         rA5 = pA5[29];
         rA6 = pA6[29];
         rA7 = pA7[29];
         rA8 = pA8[29];
         rA9 = pA9[29];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[29];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[29];
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
         rA0 = pA0[30];
         rB0 = pB0[30];
         rA1 = pA1[30];
         rA2 = pA2[30];
         rA3 = pA3[30];
         rA4 = pA4[30];
         rA5 = pA5[30];
         rA6 = pA6[30];
         rA7 = pA7[30];
         rA8 = pA8[30];
         rA9 = pA9[30];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[30];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[30];
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
         rA0 = pA0[31];
         rB0 = pB0[31];
         rA1 = pA1[31];
         rA2 = pA2[31];
         rA3 = pA3[31];
         rA4 = pA4[31];
         rA5 = pA5[31];
         rA6 = pA6[31];
         rA7 = pA7[31];
         rA8 = pA8[31];
         rA9 = pA9[31];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[31];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[31];
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
         rA0 = pA0[32];
         rB0 = pB0[32];
         rA1 = pA1[32];
         rA2 = pA2[32];
         rA3 = pA3[32];
         rA4 = pA4[32];
         rA5 = pA5[32];
         rA6 = pA6[32];
         rA7 = pA7[32];
         rA8 = pA8[32];
         rA9 = pA9[32];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[32];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[32];
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
         rA0 = pA0[33];
         rB0 = pB0[33];
         rA1 = pA1[33];
         rA2 = pA2[33];
         rA3 = pA3[33];
         rA4 = pA4[33];
         rA5 = pA5[33];
         rA6 = pA6[33];
         rA7 = pA7[33];
         rA8 = pA8[33];
         rA9 = pA9[33];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[33];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[33];
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
         rA0 = pA0[34];
         rB0 = pB0[34];
         rA1 = pA1[34];
         rA2 = pA2[34];
         rA3 = pA3[34];
         rA4 = pA4[34];
         rA5 = pA5[34];
         rA6 = pA6[34];
         rA7 = pA7[34];
         rA8 = pA8[34];
         rA9 = pA9[34];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[34];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[34];
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
         rA0 = pA0[35];
         rB0 = pB0[35];
         rA1 = pA1[35];
         rA2 = pA2[35];
         rA3 = pA3[35];
         rA4 = pA4[35];
         rA5 = pA5[35];
         rA6 = pA6[35];
         rA7 = pA7[35];
         rA8 = pA8[35];
         rA9 = pA9[35];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[35];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[35];
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
         rA0 = pA0[36];
         rB0 = pB0[36];
         rA1 = pA1[36];
         rA2 = pA2[36];
         rA3 = pA3[36];
         rA4 = pA4[36];
         rA5 = pA5[36];
         rA6 = pA6[36];
         rA7 = pA7[36];
         rA8 = pA8[36];
         rA9 = pA9[36];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[36];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[36];
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
         rA0 = pA0[37];
         rB0 = pB0[37];
         rA1 = pA1[37];
         rA2 = pA2[37];
         rA3 = pA3[37];
         rA4 = pA4[37];
         rA5 = pA5[37];
         rA6 = pA6[37];
         rA7 = pA7[37];
         rA8 = pA8[37];
         rA9 = pA9[37];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[37];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[37];
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
         rA0 = pA0[38];
         rB0 = pB0[38];
         rA1 = pA1[38];
         rA2 = pA2[38];
         rA3 = pA3[38];
         rA4 = pA4[38];
         rA5 = pA5[38];
         rA6 = pA6[38];
         rA7 = pA7[38];
         rA8 = pA8[38];
         rA9 = pA9[38];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[38];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[38];
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
         rA0 = pA0[39];
         rB0 = pB0[39];
         rA1 = pA1[39];
         rA2 = pA2[39];
         rA3 = pA3[39];
         rA4 = pA4[39];
         rA5 = pA5[39];
         rA6 = pA6[39];
         rA7 = pA7[39];
         rA8 = pA8[39];
         rA9 = pA9[39];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[39];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[39];
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
         rA0 = pA0[40];
         rB0 = pB0[40];
         rA1 = pA1[40];
         rA2 = pA2[40];
         rA3 = pA3[40];
         rA4 = pA4[40];
         rA5 = pA5[40];
         rA6 = pA6[40];
         rA7 = pA7[40];
         rA8 = pA8[40];
         rA9 = pA9[40];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[40];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[40];
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
         rA0 = pA0[41];
         rB0 = pB0[41];
         rA1 = pA1[41];
         rA2 = pA2[41];
         rA3 = pA3[41];
         rA4 = pA4[41];
         rA5 = pA5[41];
         rA6 = pA6[41];
         rA7 = pA7[41];
         rA8 = pA8[41];
         rA9 = pA9[41];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[41];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[41];
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
         rA0 = pA0[42];
         rB0 = pB0[42];
         rA1 = pA1[42];
         rA2 = pA2[42];
         rA3 = pA3[42];
         rA4 = pA4[42];
         rA5 = pA5[42];
         rA6 = pA6[42];
         rA7 = pA7[42];
         rA8 = pA8[42];
         rA9 = pA9[42];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[42];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[42];
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
         rA0 = pA0[43];
         rB0 = pB0[43];
         rA1 = pA1[43];
         rA2 = pA2[43];
         rA3 = pA3[43];
         rA4 = pA4[43];
         rA5 = pA5[43];
         rA6 = pA6[43];
         rA7 = pA7[43];
         rA8 = pA8[43];
         rA9 = pA9[43];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[43];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[43];
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
         rA0 = pA0[44];
         rB0 = pB0[44];
         rA1 = pA1[44];
         rA2 = pA2[44];
         rA3 = pA3[44];
         rA4 = pA4[44];
         rA5 = pA5[44];
         rA6 = pA6[44];
         rA7 = pA7[44];
         rA8 = pA8[44];
         rA9 = pA9[44];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[44];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[44];
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
         rA0 = pA0[45];
         rB0 = pB0[45];
         rA1 = pA1[45];
         rA2 = pA2[45];
         rA3 = pA3[45];
         rA4 = pA4[45];
         rA5 = pA5[45];
         rA6 = pA6[45];
         rA7 = pA7[45];
         rA8 = pA8[45];
         rA9 = pA9[45];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[45];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[45];
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
         rA0 = pA0[46];
         rB0 = pB0[46];
         rA1 = pA1[46];
         rA2 = pA2[46];
         rA3 = pA3[46];
         rA4 = pA4[46];
         rA5 = pA5[46];
         rA6 = pA6[46];
         rA7 = pA7[46];
         rA8 = pA8[46];
         rA9 = pA9[46];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[46];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[46];
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
         rA0 = pA0[47];
         rB0 = pB0[47];
         rA1 = pA1[47];
         rA2 = pA2[47];
         rA3 = pA3[47];
         rA4 = pA4[47];
         rA5 = pA5[47];
         rA6 = pA6[47];
         rA7 = pA7[47];
         rA8 = pA8[47];
         rA9 = pA9[47];
         rC0_0 += rA0 * rB0;
         rA10 = pA10[47];
         rC1_0 += rA1 * rB0;
         rA11 = pA11[47];
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
         pA1 += incAk;
         pA2 += incAk;
         pA3 += incAk;
         pA4 += incAk;
         pA5 += incAk;
         pA6 += incAk;
         pA7 += incAk;
         pA8 += incAk;
         pA9 += incAk;
         pA10 += incAk;
         pA11 += incAk;
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
         pA1 += incAm;
         pA2 += incAm;
         pA3 += incAm;
         pA4 += incAm;
         pA5 += incAm;
         pA6 += incAm;
         pA7 += incAm;
         pA8 += incAm;
         pA9 += incAm;
         pA10 += incAm;
         pA11 += incAm;
         pB0 += incBm;
      }
      while(pA0 != stM);
      pC0 += incCn;
      pA0 += incAn;
      pA1 += incAn;
      pA2 += incAn;
      pA3 += incAn;
      pA4 += incAn;
      pA5 += incAn;
      pA6 += incAn;
      pA7 += incAn;
      pA8 += incAn;
      pA9 += incAn;
      pA10 += incAn;
      pA11 += incAn;
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
