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
    size_t* ia_in,          /* row pointers of A */
    size_t* ja_in,         /* col indices of A */
    size_t n_in,            /* number of rows and columns of A */
    size_t nnz_in,
    size_t* order_out        /* output permutation, size n */
);
