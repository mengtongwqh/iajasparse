#include <stdio.h>
#include <assert.h>
#include "include/mindeg.h"
#include "include/colamd.h"


int mindeg            /* return (1) if OK, (0) otherwise */
(
    size_t* ia_in,          /* row pointers of A */
    size_t* ja_in,         /* col indices of A */
    size_t n_in,            /* number of rows and columns of A */
    size_t nnz_in,
    size_t* order_out        /* output permutation, size n */
) {

   int     stats[COLAMD_STATS]; /* output statistics and error codes */
   double* knobs;
   int     test, i;

   /* perform casting of input data */
   int n = (int) n_in;
   int nnz = (int) nnz_in;
   int* ia = (int*) malloc((n+1)*sizeof(int));
   int* ja = (int*) malloc(nnz*sizeof(nnz));
   int* order = (int*) malloc(n*sizeof(int));
   int* perm =  (int*) malloc((n+1)*sizeof(int));

   for (int i = 0; i < n+1; ++i) ia[i] = (int) ia_in[i];
   for (int i = 0; i < nnz; ++i) ja[i] = (int) ja_in[i];

   knobs = NULL;  /* use default params */

   assert(perm != NULL);
   /* note this has to be of size n+1 */

   test = symamd (n, ja, ia, perm, knobs,  stats , calloc, free);
 
   if (test == 0) {
      symamd_report(stats);
   }
   assert( test != 0);
  
   for (i=0; i<n; i++) {
     order[i] = perm[i];
   }

   for (int i = 0; i < n; ++i) 
       order_out[i] = (size_t) order[i];

   free(perm); perm = NULL;
   free(ia); ia = NULL;
   free(ja); ja = NULL;
   free(order); order = NULL;

   return(test);
}
