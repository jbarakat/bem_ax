/* Calling CGEQRF and CUNGQR to compute Q without workspace querying */

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <cblas.h>

int main()
{
   lapack_complex_float *a,*tau,*r,one,zero;
   int info,m,n,lda;
   int i,j;
   float err=0.0;
   m = 10;   n = 5;   lda = m;
   one = lapack_make_complex_float(1.0,0.0);
   zero= lapack_make_complex_float(0.0,0.0);
   a = (lapack_complex_float*) calloc(m*n,sizeof(lapack_complex_float));
   r = (lapack_complex_float*) calloc(n*n,sizeof(lapack_complex_float));
   tau = (lapack_complex_float*) calloc(m,sizeof(lapack_complex_float));
   for(j=0;j<n;j++)
      for(i=0;i<m;i++)
         a[i+j*m] = lapack_make_complex_float(i+1,j+1);
   info = LAPACKE_cgeqrf(LAPACK_COL_MAJOR,m,n,a,lda,tau);
   info = LAPACKE_cungqr(LAPACK_COL_MAJOR,m,n,n,a,lda,tau);
   for(j=0;j<n;j++)
      for(i=0;i<n;i++)
         r[i+j*n]=(i==j)?-one:zero;
   cblas_cgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,
               n,n,m,&one,a,lda,a,lda,&one,r,n );
   for(i=0;i<n;i++)
      for(j=0;j<n;j++)
	 if(cabs(r[i+j*n])>err)
	   err = cabs(r[i+j*n]);
   printf("error=%e\n",err);
   free(tau);
   free(r);
   free(a);
   return(info);
}
