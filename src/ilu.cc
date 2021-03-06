
#include <iaja/global_defs.h>
#include <iaja/ilu.h>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <string>


// headers for reordering schemes
extern "C" {
    #include "include/mindeg.h"
    #include "include/rcm.h"
}


/* ============================================ *
 *              SPARSEILUROW                    *
 * ============================================ */
SparseILURow::SparseILURow(const SparseILURow& rhs)
    : SparseVector(rhs), level_of_fill(rhs.level_of_fill), diag(rhs.diag) {}

SparseILURow::SparseILURow(SparseILURow&& rhs)
    : SparseVector<double>(std::move(rhs)),
    level_of_fill(std::move(rhs.level_of_fill)),
    diag(rhs.diag) {}

/* ============================================ *
 *                SPARSEILU                     *
 * ============================================ */

/***** Constructors and Destructors *****/

SparseILU::SparseILU(SparseMatrix<double>& mtrx)
    : A(mtrx), n(mtrx.nrow()),
    order_new2old(mtrx.nrow()),
    order_old2new(mtrx.nrow()) {
    assert(mtrx.nrow() == mtrx.ncol());
    // reserve enough space for rows to prevent realloc
    rows.reserve(n);
}


/***** Operator Overloading ******/

// print the ILU-factored matrix
std::ostream& operator<<(std::ostream& os, const SparseILU& ilu) {
    for (const auto i : ilu.rows) {
        os << i << "\n"; // calling the SparseVector:: operator<<
    }
    return os;
}


/***** Public Methods ******/
// compute reordering
void SparseILU::reorder(const std::string& reorder_method) {

    // reordering
    if ( reorder_method == "natural" ) {
        for (unsigned int i = 0; i < n; ++i) order_new2old[i] = i;
    } else if (reorder_method == "mindeg") {
        mindeg(A.ia, A.ja, n, A.nonzeros(), order_new2old);
    } else if (reorder_method == "rcm") {
        rcm(A.ia, A.ja, n, A.nonzeros(), order_new2old);
    } else {
        std::cerr << reorder_method
            << " is not one of the implemented reordering methods"
            << std::endl;
        exit(EXIT_FAILURE);
    }

    // inverse ordering
    for ( unsigned int i = 0; i < n; ++i )
        order_old2new[order_new2old[i]] = i;
}

// brute-force symbolic factorization
void SparseILU::analyse(const unsigned int max_level_of_fill) {

    // clean the vector row if not empty
    if (!rows.empty()) {
        rows.clear(); rows.reserve(n);
    }

    // allocate space for work vectors
    // list: implied linked list of nonzero locations of this row
    std::vector<unsigned int> row_linked_list(n, 0);
    std::vector<unsigned int> row_level_of_fill(n, 0);

    // loop thru each row of the matrix
    for (unsigned int i = 0; i < n; ++i) {

        // load nonzeros of this row into linked list
        unsigned int iold = order_new2old[i];
        unsigned int nnzrow = A.ia[iold+1] - A.ia[iold];

        // for each entry, convert from old to new ordering
        std::vector<unsigned int> jaf(nnzrow);
        for (unsigned int jj = 0, j = A.ia[iold]; j < A.ia[iold+1]; ++j, ++jj)
            jaf[jj] = order_old2new[A.ja[j]];
        // sort the entries ascending
        if (nnzrow > 0) std::sort(jaf.begin(), jaf.end());


        // ILU level-of-fill symbolic factorization
        if ( max_level_of_fill == 0 ) {

            /* ILU(0) */
            // factorized matrix has the same sparsity pattern
            // as the original matrix
            rows.push_back(SparseILURow(n, nnzrow));
            for (unsigned int k = 0; k < nnzrow; ++k) {
                rows[i].ja[k] = jaf[k];
            }

        } else {
            /* ILU(maxLevelOfFill) */
            // merge row i with each row in the linked list
            // now we know the number of nonzero entries, allocate ILU

            // reset row_level_of_fill to record LevelOfFill for current row entries
            std::fill(row_level_of_fill.begin(), row_level_of_fill.end(), n+1);

            // construct linked list
            unsigned int klist, list_begin;
            klist = list_begin = jaf[0]; // linked list begin
            for ( unsigned int j = 1;
                    j < nnzrow;
                    ++j, klist = row_linked_list[klist] ) {
                row_linked_list[klist] = jaf[j];
            }
            row_linked_list[klist] = n + 1; // linked list end

            // set the level of original entries to 0
            for ( klist = list_begin;
                    klist != n+1;
                    klist = row_linked_list[klist] ) {
                row_level_of_fill[klist] = 0;
            }

            // merging induced fills
            merge_linked_list(i, max_level_of_fill, list_begin,
                    row_linked_list, row_level_of_fill, nnzrow);

            // allocate sparseILURow for this row
            rows.push_back(SparseILURow(n, nnzrow));
            // skip thru linked list and extract level-of-fill and jaf
            for (unsigned int j = 0, klist = list_begin;
                    klist != n+1;
                    klist = row_linked_list[klist], ++j) {
                rows[i].ja[j] = klist;
                rows[i].level_of_fill[j] = row_level_of_fill[klist];
            }
        } // if max_level_of_fill

        // locate diagonal entry
        unsigned int jdiag = 0;
        for ( ; rows[i].ja[jdiag] < i && jdiag < nnzrow; ++jdiag ) {}
        if ( rows[i].ja[jdiag] == i ) {
            rows[i].diag = jdiag;
        } else {
            std::cerr << "ILU analyse: Diagonal or row " << i << " is zero" << std::endl;
            exit(EXIT_FAILURE);
        }
    } // i-loop thru all rows
}

void SparseILU::merge_linked_list(
        const unsigned int rowid,
        const unsigned int max_level_of_fill,
        const unsigned int list_begin,
        std::vector<unsigned int>& row_linked_list,
        std::vector<unsigned int>& row_level_of_fill,
        unsigned int& nnzrow
        ) {

    // loop over each column index in the implied linked list
    for ( unsigned int i = list_begin; i < rowid; i = row_linked_list[i] ) {

        // only consider the entries where current level-of-fill
        // is less than max level-of-fill
        if ( row_level_of_fill[i] < max_level_of_fill ) {

            unsigned int klist = i;

            // loop thru each column in the ith-row
            for ( unsigned int jj = rows[i].diag; jj < rows[i].nnz; ++jj ) {

                unsigned int j = rows[i].ja[jj];
                unsigned int levelij = std::min(
                        rows[i].level_of_fill[jj] + row_level_of_fill[i] + 1,
                        row_level_of_fill[j] );

                // if the fill has been considered before
                // no need to add it again to the linked list
                if ( levelij <= max_level_of_fill && row_level_of_fill[j] > max_level_of_fill ) {
                    // otherwise add this fill to the linked list
                    while ( row_linked_list[klist] < j )
                        // jump forward to the correct insert position
                        klist = row_linked_list[klist];
                    unsigned int knext = row_linked_list[klist];
                    // modify linked list pointers
                    row_linked_list[klist] = j;
                    row_linked_list[j]     = knext;
                    // increment nonzero count for this row
                    ++nnzrow;
                }

                row_level_of_fill[j] = levelij;
            } // jj-loop
        }
    } // i-loop
}

void SparseILU::factor() {

    assert(order_new2old && order_old2new);

    // allocate temp vector
    FullVector<double> row_temp(n);

    // loop through each row
    for ( unsigned int i = 0; i < n; ++i ) {

        unsigned int iold = order_new2old[i];

        // scatter the entries of ith row under new ordering
        for ( unsigned int j = A.ia[iold]; j < A.ia[iold+1]; ++j )
            row_temp[ order_old2new[ A.ja[j]] ] = A.a[j];

        // scan work vector columnwise, determine multiplier
        // only elements left of the diagonal are considered
        for ( unsigned int j = 0; j < rows[i].diag; ++j ) {

            unsigned int idx_col = rows[i].ja[j];

            // compute multiplier for this entry
            double mult = row_temp[idx_col] / rows[idx_col].a[rows[idx_col].diag];
            row_temp[idx_col] = mult;

            // subtract the mult*(j-th row) from temp vector entries
            // to the right of this entry
            for ( unsigned int jj = rows[idx_col].diag+1; jj < rows[idx_col].nnz; ++jj ) {
                int idx_row = rows[idx_col].ja[jj];
                row_temp[idx_row] -= mult * rows[idx_col].a[jj];
            }
        } // j-loop for L of previously factored rows

        // gather nonzero entries in row_temp into packed form
        for (unsigned int j = 0; j < rows[i].nnz; ++j) {
            unsigned int idx = rows[i].ja[j];
            rows[i].a[j] = row_temp[idx];
            row_temp[idx] = 0.0;
        }
    } // i(row)-loop
}


void SparseILU::solve(const FullVector<double>& b, FullVector<double>& x) const {

    assert(x.length() == n && b.length() == n);

    // alloc temp vector
    FullVector<double> y(n);

    // load b (in old order) into y (in new order)
    for (unsigned int i = 0; i < n; ++i)
        y[i] = b[ order_new2old[i] ];

    // solve Ly = b
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < rows[i].diag; ++j) {
            y[i] -=  rows[i].a[j] * y[ rows[i].ja[j] ];
        }
    }

    // solve U(y_new) = y
    for (int i = n-1; i >= 0; --i) {
        for (unsigned int j = rows[i].diag+1; j < rows[i].nnz; ++j) {
            y[i] -= rows[i].a[j] * y[ rows[i].ja[j] ];
        }
        y[i] /= rows[i].a[ rows[i].diag ];
    }

    // revert y (new order) to x (old order)
    for (unsigned int i = 0; i < n; ++i)
        x[ order_new2old[i] ] = y[i];
}


void SparseILU::print_level_of_fill(std::ostream& os) const {

    const unsigned int width = PRINT_WIDTH_UNSIGNED_INT;

    for (auto i = rows.cbegin(); i != rows.cend(); ++i) {
        for (unsigned int j = 0, jj = 0; j < n; ++j) {
            if ( jj < i->nnz && i->ja[jj] == j ) {
                os << std::setw(width) << std::right << i->level_of_fill[jj];
                ++jj;
            } else {
                os <<std::setw(width) << std::right << "";
            }
        }
        os << "\n";
    }
}

