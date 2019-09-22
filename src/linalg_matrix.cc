#include <iaja/iaja_config.h>
#include <iaja/linalg_matrix.h>

IAJA_NAMESPACE_OPEN

/* ============================================ *
 *             MATRIX SPARSITY                  *
 * ============================================ */

/* ---------------------- *
 *   Ctor, Dtor, Assign   *
 * ---------------------- */

SparsityIaja::SparsityIaja():
    nr(0), nc(0), nnz(0), ia(), ja() {}

SparsityIaja::SparsityIaja(size_type nc,
        const FullVector<size_type>& ia, const FullVector<size_type>& ja):
     nr(ia.length()-1), nc(nc), nnz(ja.length()), ia(ia), ja(ja) {}

SparsityIaja::SparsityIaja(SparsityIaja&& rhs):
    nr(rhs.nr), nc(rhs.nc), nnz(rhs.nnz),
    ia(std::move(rhs.ia)), ja(std::move(rhs.ja)) {
    rhs.nr = rhs.nc = rhs.nnz = 0;
}

SparsityIaja& SparsityIaja:: operator= (SparsityIaja&& rhs) {
    if ( this != &rhs ) {
        nr = rhs.nr; nc = rhs.nc; nnz = rhs.nnz;
        ia = std::move(rhs.ia); ja = std::move(rhs.ja);
        rhs.nr = rhs.nc = rhs.nnz = 0;
    }
    return *this;
}


/* -------- *
 *   I/O    *
 * -------- */

std::ostream& operator<<(std::ostream& os, const SparsityIaja& sp) {

    auto width = SPARSITY_PRINT_WIDTH;

    for (auto i = sp.ia.cbegin(), jj = sp.ja.cbegin(); i != sp.ia.cend(); ++i) {

        for (SparsityIaja::size_type j = 0 ; j < sp.nc; ++j) {
            if ( jj != sp.ja.cend() && *jj == j ) {
                ++jj;
                os << std::setw(width) << std::left << SPARSITY_FILL_PATTERN;
            } else {
                os << std::setw(width) << std::left << SPARSITY_EMPTY_PATTERN;
            }
        }

        os << "\n";
    }

    return os;
}

IAJA_NAMESPACE_CLOSE
