#include <iaja/iaja_config.h>
#include <iaja/ilu.h>
#include <iaja/util.h>

#include <algorithm>
#include <cassert>
#include <iomanip>

// #ifdef PROFILE_ALL
    // #define PROFILING
// #else
    // #undef PROFILING
// #endif

IAJA_NAMESPACE_OPEN

/* ============================================ *
 *              SPARSEFACTORROW                    *
 * ============================================ */

IncompleteFactor::SparseFactorRow::SparseFactorRow(size_type n,
        FullVector<size_type>&& ja, FullVector<unsigned int>&& lof):
    SparseVector<FloatType>(n, std::move(ja)),
    level_of_fill(std::move(lof)), diag(0) {}

IncompleteFactor::SparseFactorRow::SparseFactorRow(SparseFactorRow&& rhs):
    SparseVector<FloatType>(std::move(rhs)),
    level_of_fill(std::move(rhs.level_of_fill)),
    diag(rhs.diag) {rhs.diag = 0;}

IncompleteFactor::SparseFactorRow&
IncompleteFactor::SparseFactorRow:: operator= (
        SparseFactorRow&& rhs) {
    if (this != &rhs) {
        SparseVector<FloatType>:: operator= (std::move(rhs));
        level_of_fill = std::move(rhs.level_of_fill);
        diag = rhs.diag;
        rhs.diag = 0;
    }
    return *this;
}


/* ============================================ *
 *              INCOMPLETEFACTOR                *
 * ============================================ */

/* ------------------------- *
 *     Ctor, Dtor, Assign    *
 * ------------------------- */

IncompleteFactor::IncompleteFactor(SparseMatrixIaja<FloatType>& mtrx,
        const std::string& reorder_method)
    : A(mtrx), n(mtrx.nrow()),
    order_new2old(mtrx.nrow()),
    order_old2new(mtrx.nrow()) {
    assert(mtrx.nrow() == mtrx.ncol());
    // reserve enough space for rows to prevent realloc
    rows.reserve(n);
    // init reordering
    reorder(reorder_method);
}

IncompleteFactor::IncompleteFactor(IncompleteFactor&& rhs):
    A(rhs.A), n(rhs.n),
    order_new2old(std::move(rhs.order_new2old)),
    order_old2new(std::move(rhs.order_old2new)),
    rows(std::move(rows)) {}

IncompleteFactor& IncompleteFactor:: operator = (IncompleteFactor&& rhs) {
    if (this != &rhs) {
        A = rhs.A; n = rhs.n;
        order_new2old = std::move(rhs.order_new2old);
        order_old2new = std::move(rhs.order_old2new);
        rows = std::move(rhs.rows);
    }
    return *this;
}


/* ---------- *
 *    I/O     *
 * ---------- */

std::ostream& operator<<(std::ostream& os, const IncompleteFactor& ifac) {
    for (const auto& i : ifac.rows) {
        // calling the SparseVector:: operator<<
        os << i << "\n";
    }
    return os;
}

void IncompleteFactor::print_level_of_fill(std::ostream& os) const {
    const unsigned int width = 2;
    for (auto i = rows.cbegin(); i != rows.cend(); ++i) {
        for (size_type j = 0, jj = 0; j < n; ++j) {
            if ( jj < i->nnonzero() && i->get_ja(jj) == j ) {
                os << std::setw(width) << std::right << i->level_of_fill[++jj];
            } else {
                os <<std::setw(width) << std::right << " ";
            }
        }
        os << "\n";
    }
}

/* ------------------------- *
 *      Type Conversion      *
 * ------------------------- */

IncompleteFactor:: operator SparseMatrixIaja<FloatType>() const {

    std::vector<size_type> iaf, jaf;
    std::vector<FloatType> af;
    iaf.push_back(0);

    for (size_type i = 0, col_counter = 0; i < n; ++i) {
        for (size_type j = 0; j < rows[i].nnonzero(); ++j) {
            jaf.push_back(rows[i].get_ja(j));
            af.push_back(rows[i][j]);
        }
        iaf.push_back(col_counter+=rows[i].nnonzero());
    }
    return SparseMatrixIaja<FloatType>(n, iaf, jaf, af);
}


/* ---------------------- *
 *  Numerical Operations  *
 * ---------------------- */

void IncompleteFactor::reorder(const std::string& reorder_method) {

    // --------------
    //   reordering
    // --------------

    if ( reorder_method == "natural" ) {
        for (size_type i = 0; i < n; ++i) order_new2old[i] = i;
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


void IncompleteFactor::analyse(const unsigned int max_level_of_fill) {

    // ------------------------------------
    // brute-force symbolic factorization
    // ------------------------------------

#ifdef PROFILING
    Timer timer;
#endif

    // clean the vector row if not empty
    if (!rows.empty()) {
        rows.clear(); rows.reserve(n);
    }

    // allocate space for work vectors
    // list: implied linked list of nonzero locations of this row
    std::vector<size_type> row_linked_list(n, 0);
    std::vector<unsigned int> row_level_of_fill(n, n+1);

    // loop thru each row of the matrix
    for (size_type i = 0; i < n; ++i) {

        // load nonzeros of this row into linked list
        size_type iold = order_new2old[i];
        size_type nnzrow = A.get_ia(iold+1) - A.get_ia(iold);

        // for each entry, convert from old to new ordering
        FullVector<size_type> jaf(nnzrow);
        for (size_type jj = 0, j = A.get_ia(iold); j < A.get_ia(iold+1); ++j, ++jj)
            jaf[jj] = order_old2new[A.get_ja(j)];
        // sort the entries ascending
        if (nnzrow > 0) std::sort(jaf.begin(), jaf.end());


        // ILU level-of-fill symbolic factorization
        if ( max_level_of_fill == 0 ) {

            /* ILU(0) */
            // factorized matrix has the same sparsity pattern
            // as the original matrix
            FullVector<size_type> ja_this_row(nnzrow);
            FullVector<unsigned int> lof_this_row(nnzrow);
            for (size_type k = 0; k < nnzrow; ++k) {
                ja_this_row[k] = jaf[k];
                lof_this_row[k] = 0;
            }
            rows.emplace_back(n, std::move(ja_this_row), std::move(lof_this_row));

        } else {
            /* ILU(maxLevelOfFill) */

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
            for ( size_type klist = list_begin;
                    klist != n+1;
                    klist = row_linked_list[klist] ) {
                row_level_of_fill[klist] = 0;
            }

            // merge row i with each row in the linked list
            merge_linked_list(i, max_level_of_fill, list_begin,
                    row_linked_list, row_level_of_fill, nnzrow);

            // now we know the number of nonzero entries, allocate SparseFactorRow
            FullVector<size_type> ja_this_row(nnzrow);
            FullVector<unsigned int> lof_this_row(nnzrow);
            // skip thru linked list and extract level-of-fill and jaf
            for (size_type klist = list_begin, j = 0;
                    klist != n+1;
                    klist = row_linked_list[klist], ++j) {
                ja_this_row[j] = klist;
                lof_this_row[j] = row_level_of_fill[klist];
                // reset row_level_of_fill to record lof from next row
                row_level_of_fill[klist] = n+1;
            }
            rows.emplace_back(n, std::move(ja_this_row), std::move(lof_this_row));
        } // if max_level_of_fill == 0

        // locate diagonal entry
        size_type jdiag = 0;
        for ( ; rows[i].get_ja(jdiag) < i && jdiag < nnzrow; ++jdiag ) {}
        if ( rows.back().get_ja(jdiag) == i ) {
            rows.back().diag = jdiag;
        } else {
            std::cerr << "ILU analyse: Diagonal or row " << i << " is zero" << std::endl;
            exit(EXIT_FAILURE);
        }
    } // i-loop thru all rows

#ifdef PROFILING
    std::cout << "Timing for ILU symbolic factorization (analyse): " <<
        std::scientific << timer.elapsed() << "\n";
#endif
}

void IncompleteFactor::merge_linked_list(
        const size_type rowid,
        const unsigned int max_level_of_fill,
        const size_type list_begin,
        std::vector<size_type>& row_linked_list,
        std::vector<unsigned int>& row_level_of_fill,
        size_type& nnzrow
        ) {

    // ------------------------------------
    // helper function for symbolic_factor
    // merge induced fills into linked list
    // ------------------------------------

    // loop over each column index in the implied linked list
    for ( size_type i = list_begin; i < rowid; i = row_linked_list[i] ) {

        // only consider the entries where current level-of-fill
        // is less than max level-of-fill
        if ( row_level_of_fill[i] < max_level_of_fill ) {

            size_type klist = i;

            // loop thru each column in the ith-row
            for ( size_type jj = rows[i].diag+1; jj < rows[i].nnonzero(); ++jj ) {

                size_type j = rows[i].get_ja(jj);
                unsigned int levelij = std::min(
                        rows[i].level_of_fill[jj] + row_level_of_fill[i] + 1,
                        row_level_of_fill[j] );

                // if the fill has been considered before
                // no need to add it again to the linked list
                if ( levelij <= max_level_of_fill &&
                        row_level_of_fill[j] > max_level_of_fill ) {
                    // otherwise add this fill to the linked list
                    while ( row_linked_list[klist] < j ) {
                        // jump forward to the correct insert position
                        klist = row_linked_list[klist];
                    }
                    assert(klist < j);
                    assert(row_linked_list[klist] > j);
                    size_type knext = row_linked_list[klist];
                    // insert this entry into the linked list
                    row_linked_list[klist] = j;
                    row_linked_list[j]     = knext;
                    klist = j;
                    // increment nonzero count for this row
                    ++nnzrow;
                }
                row_level_of_fill[j] = levelij;
            } // jj-loop
        }
    } // i-loop

#ifdef DEBUG
        // make sure the list is in order
        size_type j = 0;
        for ( size_type i = list_begin;
                i != n+1;
                i = row_linked_list[i], ++j) {
            assert(i < row_linked_list[i]);
        }
        assert(j == nnzrow);
#endif
}



/* ============================================ *
 *                SPARSEILU                     *
 * ============================================ */

/* ------------------------- *
 *     Ctor, Dtor, Assign    *
 * ------------------------- */
SparseILU::SparseILU(SparseMatrixIaja<FloatType>& mtrx,
        const std::string& reorder_method) :
    IncompleteFactor(mtrx, reorder_method) {}

SparseILU::SparseILU(SparseILU&& rhs):
    IncompleteFactor(std::move(rhs)) {}

SparseILU& SparseILU:: operator = (SparseILU&& rhs) {
    if (this != & rhs) {
        IncompleteFactor:: operator = (std::move(rhs));
    }
    return *this;
}


/* ------------------------- *
 *     Numerical Methods     *
 * ------------------------- */

void SparseILU::factor() {

    // ------------------------------------
    // Numeric factorization
    // ------------------------------------

#ifdef PROFILING
    Timer timer;
#endif

    // allocate temp vector
    FullVector<FloatType> row_temp(n, 0.0);

    // loop through each row
    for ( size_type i = 0; i < n; ++i ) {

        size_type iold = order_new2old[i];

        // scatter the entries of ith row under new ordering
        for ( size_type j = A.get_ia(iold); j < A.get_ia(iold+1); ++j )
            row_temp[ order_old2new[ A.get_ja(j)] ] = A[j];

        // scan work vector columnwise, determine multiplier
        // only elements left of the diagonal are considered
        for ( size_type j = 0; j < rows[i].diag; ++j ) {

            size_type idx_col = rows[i].get_ja(j);

            // compute multiplier for this entry
            FloatType mult = row_temp[idx_col] / rows[idx_col][rows[idx_col].diag];
            row_temp[idx_col] = mult;

            // subtract the mult*(j-th row) from temp vector entries
            // to the right of this entry
            for ( size_type jj = rows[idx_col].diag+1; jj < rows[idx_col].nnonzero(); ++jj ) {
                size_type idx_row = rows[idx_col].get_ja(jj);
                row_temp[idx_row] -= mult * rows[idx_col][jj];
            }
        } // j-loop for L of previously factored rows

        // gather nonzero entries in row_temp into packed form
        for (size_type j = 0; j < rows[i].nnonzero(); ++j) {
            size_type idx = rows[i].get_ja(j);
            rows[i][j] = row_temp[idx];
            row_temp[idx] = 0.0;
        }
    } // i(row)-loop

#ifdef PROFILING
    std::cout << "Timing for ILU numeric factorization (factor): " <<
        std::scientific << timer.elapsed() << "\n";
#endif
}



void SparseILU::solve(const FullVector<FloatType>& b, FullVector<FloatType>& x) const {

    // ------------------------------------
    // Forward and back solve
    // ------------------------------------
#ifdef PROFILING
    Timer timer;
#endif

    assert(x.length() == n && b.length() == n);

    // alloc temp vector
    FullVector<FloatType> y(n);

    // load b (in old order) into y (in new order)
    for (size_type i = 0; i < n; ++i)
        y[i] = b[ order_new2old[i] ];

    // solve Ly = b
    for (size_type i = 0; i < n; ++i) {
        for (size_type j = 0; j < rows[i].diag; ++j) {
            y[i] -=  rows[i][j] * y[ rows[i].get_ja(j) ];
        }
    }

    // solve U(y_new) = y
    for (size_type i = n; i >= 1; --i) {
        for (size_type j = rows[i-1].diag+1; j < rows[i-1].nnonzero(); ++j) {
            y[i-1] -= rows[i-1][j] * y[ rows[i-1].get_ja(j) ];
        }
        y[i-1] /= rows[i-1][rows[i-1].diag];
    }

    // revert y (new order) to x (old order)
    for (size_type i = 0; i < n; ++i)
        x[ order_new2old[i] ] = y[i];

#ifdef PROFILING
    std::cout << "Timing for ILU forward/back solve (solve): " <<
        std::scientific << timer.elapsed() << "\n";
#endif
}



/* ============================================ *
 *                SPARSEICHOL                   *
 * ============================================ */

/* ------------------------- *
 *   Ctor, Dtor, Assign      *
 * ------------------------- */

SparseIchol::SparseIchol(SparseMatrixIaja<FloatType>& mtrx,
        const std::string& reorder_method):
    IncompleteFactor(mtrx, reorder_method) {}

SparseIchol::SparseIchol(SparseIchol&& rhs):
    IncompleteFactor(std::move(rhs)) {}

SparseIchol& SparseIchol:: operator = (SparseIchol&& rhs) {
    if (this != &rhs)
        IncompleteFactor:: operator = (std::move(rhs));
    return *this;
}

/* ------------------------- *
 *   Numerical Operation     *
 * ------------------------- */

void SparseIchol::factor() {

    FullVector<FloatType> row_temp(n, 0.0);
    FullVector<size_type> upper_offset(n);
    for (size_type i = 0; i < n; ++i)
        upper_offset[i] = rows[i].diag + 1;

    for ( size_type i = 0; i < n; ++i ) {

        size_type iold = order_new2old[i];

        // scatter the entries of ith row under new ordering
        for ( size_type j = A.get_ia(iold); j < A.get_ia(iold+1); ++j )
            // TODO only load the lower-trig part?
            row_temp[ order_old2new[ A.get_ja(j) ] ] = A[j];

        // compute off-diagonal entries for this row
        for ( size_type j = 0; j < rows[i].diag; ++j ) {

            size_type idx_col = rows[i].get_ja(j);

            // off-diagonal entries
            for ( size_type k = 0; k < rows[idx_col].diag; ++k ) {
                size_type idx_row = rows[j].get_ja(k);
                row_temp[idx_col] -= row_temp[idx_row]*rows[idx_col][k];
            }
            // TODO throw exception if zero diagonal
            assert(fabs(rows[idx_col][rows[idx_col].diag]) > FLOATING_POINT_NEARLY_ZERO);
            row_temp[idx_col] /= rows[idx_col][rows[idx_col].diag];
        }

        // compute diagonal entry for this row
        double diag_entry = row_temp[rows[i].get_ja(rows[i].diag)];
        for (size_type j = 0; j < rows[i].diag; ++j) diag_entry -= rows[i][j];
        row_temp[rows[i].get_ja(rows[i].diag)] = sqrt(diag_entry);

        // gather nonzero entries in row_temp into packed form
        for (size_type j = 0; j <= rows[i].diag; ++j) {
            // set entries left of diagonal
            FloatType  this_entry = row_temp[rows[i].get_ja(j)];
            rows[i][j] = this_entry;
            // mirror reflection w.r.t diagonal
            size_type idx_row = rows[i].get_ja(j);
            size_type idx_col = upper_offset[idx_row]++;
            rows[idx_row][idx_col] = this_entry;
        }

        // reset row_temp to zero;
        for ( size_type j= A.get_ia(iold); j< A.get_ia(iold+1); ++j )
            row_temp[A.get_ja(j)] = 0.0;
    } // i(row)-loop
}


void SparseIchol::solve(const FullVector<FloatType>& b,
        FullVector<FloatType>& x) const {

    assert(x.length() == n && b.length() == n);

    // allocate temp vector
    FullVector<FloatType> y(n);

    // load b (in old order) into y (in new order)
    for (size_type i = 0; i < n; ++i)
        y[i] = b[order_new2old[i]];

    // solve Ly = b
    for (size_type i = 0; i < n; ++i) {
        for (size_type j = 0; j < rows[i].diag; ++j) {
            y[i] -= rows[i][j] * y[rows[i].get_ja(j)];
        }
        y[i] /= rows[i][rows[i].diag];
    }


    // solve L^{T}(y_new) = y
    for (size_type i = 0; i < n; ++i) {
        for (size_type j = rows[n-i-1].diag+1; j < rows[n-i-1].nnonzero(); ++j) {
            y[n-i-1] -= y[rows[n-i-1].get_ja(j)] * rows[n-i-1][j];
        }
        y[n-i-1] /= rows[n-i-1][rows[n-i-1].diag];
    }

    // revert y (new order) to x (older)
    for (size_type i = 0; i < n; ++i) {
        x[ order_new2old[i] ] = y[i];
    }
}

IAJA_NAMESPACE_CLOSE
