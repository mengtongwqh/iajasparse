// [> CS 472/672 Assignment #3 - mindeg.h file <]

// int mindeg(int *ia,int *ja,int N,int *order);
/*
   Minimum degree ordering

   input:

      N = number of rows
      ia, ja = ia-ja sparse matrix data structure

   output:

      order[new_order] = old_order

   returns 1 if all ok
           0 if failure
*/
int mindeg            /* return (1) if OK, (0) otherwise */
(
    unsigned int* ia_in,          /* row pointers of A */
    unsigned int* ja_in,         /* col indices of A */
    unsigned int n_in,            /* number of rows and columns of A */
    unsigned int nnz_in,
    unsigned int* order_out        /* output permutation, size n */
);
