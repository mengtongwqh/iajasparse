#include <iaja/iaja_config.h>
#include <iaja/ilu.h>

#include <algorithm>
#include <cassert>
#include <iomanip>

IAJA_NAMESPACE_OPEN

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
        for (size_type i = 0; i < n; ++i) order_new2old[i] = i;
    } else if (reorder_method == "mindeg") {
        // mindeg(A.ia, A.ja, n, A.nonzeros(), order_new2old);
    } else if (reorder_method == "rcm") {
        // rcm(A.ia, A.ja, n, A.nonzeros(), order_new2old);
    } else {
        std::cerr << reorder_method
            << " is not one of the implemented reordering methods"
            << std::endl;
        exit(EXIT_FAILURE);
    }

    // inverse ordering
    for ( size_type i = 0; i < n; ++i )
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
    std::vector<size_type> row_linked_list(n, 0);
    std::vector<unsigned int> row_level_of_fill(n, 0);

    // loop thru each row of the matrix
    for (size_type i = 0; i < n; ++i) {

        // load nonzeros of this row into linked list
        size_type iold = order_new2old[i];
        size_type nnzrow = A.ia[iold+1] - A.ia[iold];

        // for each entry, convert from old to new ordering
        std::vector<size_type> jaf(nnzrow);
        for (size_type jj = 0, j = A.ia[iold]; j < A.ia[iold+1]; ++j, ++jj)
            jaf[jj] = order_old2new[A.ja[j]];
        // sort the entries ascending
        if (nnzrow > 0) std::sort(jaf.begin(), jaf.end());


        // ILU level-of-fill symbolic factorization
        if ( max_level_of_fill == 0 ) {

            /* ILU(0) */
            // factorized matrix has the same sparsity pattern
            // as the original matrix
            rows.push_back(SparseILURow(n, nnzrow));
            for (size_type k = 0; k < nnzrow; ++k) {
                rows[i].ja[k] = jaf[k];
            }

        } else {
            /* ILU(maxLevelOfFill) */
            // merge row i with each row in the linked list
            // now we know the number of nonzero entries, allocate ILU

            // reset row_level_of_fill to record LevelOfFill for current row entries
            std::fill(row_level_of_fill.begin(), row_level_of_fill.end(), n+1);

            // construct linked list
            size_type klist, list_begin;
            klist = list_begin = jaf[0]; // linked list begin
            for ( size_type j = 1;
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
            for (size_type j = 0, klist = list_begin;
                    klist != n+1;
                    klist = row_linked_list[klist], ++j) {
                rows[i].ja[j] = klist;
                rows[i].level_of_fill[j] = row_level_of_fill[klist];
            }
        } // if max_level_of_fill

        // locate diagonal entry
        size_type jdiag = 0;
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
        const size_type rowid,
        const unsigned int max_level_of_fill,
        const size_type list_begin,
        std::vector<size_type>& row_linked_list,
        std::vector<unsigned int>& row_level_of_fill,
        size_type& nnzrow
        ) {

    // loop over each column index in the implied linked list
    for ( size_type i = list_begin; i < rowid; i = row_linked_list[i] ) {

        // only consider the entries where current level-of-fill
        // is less than max level-of-fill
        if ( row_level_of_fill[i] < max_level_of_fill ) {

            size_type klist = i;

            // loop thru each column in the ith-row
            for ( size_type jj = rows[i].diag; jj < rows[i].nnz; ++jj ) {

                size_type j = rows[i].ja[jj];
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
                    size_type knext = row_linked_list[klist];
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

#ifdef PROFILING
    clock_t begin = clock();
#endif

    // allocate temp vector
    FullVector<double> row_temp(n);

    // loop through each row
    for ( size_type i = 0; i < n; ++i ) {

        size_type iold = order_new2old[i];

        // scatter the entries of ith row under new ordering
        for ( size_type j = A.ia[iold]; j < A.ia[iold+1]; ++j )
            row_temp[ order_old2new[ A.ja[j]] ] = A.a[j];

        // scan work vector columnwise, determine multiplier
        // only elements left of the diagonal are considered
        for ( size_type j = 0; j < rows[i].diag; ++j ) {

            size_type idx_col = rows[i].ja[j];

            // compute multiplier for this entry
            double mult = row_temp[idx_col] / rows[idx_col].a[rows[idx_col].diag];
            row_temp[idx_col] = mult;

            // subtract the mult*(j-th row) from temp vector entries
            // to the right of this entry
            for ( size_type jj = rows[idx_col].diag+1; jj < rows[idx_col].nnz; ++jj ) {
                size_type idx_row = rows[idx_col].ja[jj];
                row_temp[idx_row] -= mult * rows[idx_col].a[jj];
            }
        } // j-loop for L of previously factored rows

        // gather nonzero entries in row_temp into packed form
        for (size_type j = 0; j < rows[i].nnz; ++j) {
            size_type idx = rows[i].ja[j];
            rows[i].a[j] = row_temp[idx];
            row_temp[idx] = 0.0;
        }
    } // i(row)-loop

#ifdef PROFILING
    clock_t end = clock();
    std::cout << "Timing for ILU numerical factorization: " << std::scientific <<
        ( static_cast<double>(end - begin)/static_cast<double>(CLOCKS_PER_SEC) )
        << "\n";
#endif
}


void SparseILU::solve(const FullVector<double>& b, FullVector<double>& x) const {

    assert(x.length() == n && b.length() == n);

    // alloc temp vector
    FullVector<double> y(n);

    // load b (in old order) into y (in new order)
    for (size_type i = 0; i < n; ++i)
        y[i] = b[ order_new2old[i] ];

    // solve Ly = b
    for (size_type i = 0; i < n; ++i) {
        for (size_type j = 0; j < rows[i].diag; ++j) {
            y[i] -=  rows[i].a[j] * y[ rows[i].ja[j] ];
        }
    }

    // solve U(y_new) = y
    for (size_type i = n; i >= 1; --i) {
        for (size_type j = rows[i-1].diag+1; j < rows[i-1].nnz; ++j) {
            y[i-1] -= rows[i-1].a[j] * y[ rows[i-1].ja[j] ];
        }
        y[i-1] /= rows[i-1].a[ rows[i-1].diag ];
    }

    // revert y (new order) to x (old order)
    for (size_type i = 0; i < n; ++i)
        x[ order_new2old[i] ] = y[i];
}


void SparseILU::print_level_of_fill(std::ostream& os) const {

    const unsigned int width = PRINT_WIDTH_UNSIGNED_INT;

    for (auto i = rows.cbegin(); i != rows.cend(); ++i) {
        for (size_type j = 0, jj = 0; j < n; ++j) {
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

IAJA_NAMESPACE_CLOSE
