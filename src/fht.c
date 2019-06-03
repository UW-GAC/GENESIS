#include <R.h>
#include <Rinternals.h>

// (a,b) -> (a+b,a-b) without overflow
void rotate(double x[], int a, int b)
{
    static double t;
    t = x[a];
    x[a] = x[a] + x[b];
    x[b] = t - x[b];
}

// Integer log2
int ilog2(int x)
{
    int l2 = 0;
    for (; x; x >>=1) ++l2;
    return l2;
}

/**
 * Fast Walsh-Hadamard transform
 */

void fwht(double data[], int *size )
{
    const int l2 = ilog2(*size) - 1;

    for (int i = 0; i < l2; ++i)
    {
        for (int j = 0; j < (1 << l2); j += 1 << (i+1))
        for (int k = 0; k < (1 << i ); ++k)
	  rotate(data, j + k, j + k + (1<<i) );
    }
}

/**
 *  Matrix version
 */
 
void mfwht(double data[], int *nrow, int *ncol){

  for(int i=0; i < *ncol; i++){ fwht(data+(i * *nrow), nrow);}

}

SEXP big_mfwht(SEXP data){

  SEXP answer = PROTECT(Rf_duplicate(data));
  
  SEXP dims = PROTECT(getAttrib(data, Rf_install("dim")));
  int nrow = INTEGER(dims)[0];
  int ncol = INTEGER(dims)[1];
  UNPROTECT(1); /* dims */
  
  for(int i=0; i < ncol; i++){ fwht( REAL(answer)+(i * nrow), &nrow);}

  UNPROTECT(1); /* answer */
  
  return answer;
}
