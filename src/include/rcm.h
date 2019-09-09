/* CS 472/672 Assignment #3 - rcm.h file */

void rcm(unsigned int *ia, unsigned int *ja, 
        unsigned int N, unsigned int NNZ, unsigned int *order);
/* 
   Reverse Cuthill-McKee ordering

   input:

      N = number of rows
      ia, ja = ia-ja sparse matrix data structure

   output:

      order[new_order] = old_order
*/

/*
  These are work routines used by rcm().  You should not have to call these.
*/
void levels(int *ia,int *ja,int n,int iroot,int *nlvl,int *ial,int *jal);
void rcm1(int *ia,int *ja,int n,int iroot,int nlvl,int *ial,int *jal,int *order);
void pseudo(int *ia,int *ja,int n,int *ial,int *jal,int *iroot,int *nlvl);

