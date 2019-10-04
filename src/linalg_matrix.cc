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

SparsityIaja::SparsityIaja(size_type nr, size_type nc, size_type nnz, 
        const size_type* ia, const size_type* ja):
    nr(nr), nc(nc), nnz(nnz), ia(nr+1, ia), ja(nnz, ja) {}

SparsityIaja::SparsityIaja(size_type nr, size_type nc, size_type nnz,
        const long int* ia_in, const long int* ja_in):
    nr(nr), nc(nc), nnz(nnz), ia(nr+1), ja(nnz) {

        size_type j = 0;
        for (auto i = ia.begin(); i != ia.end(); ++i, ++j)
            *i = ia_in[j] >= 0 ? ia_in[j] : nnz-ia_in[j]-1;
        j = 0;
        for (auto i = ja.begin(); i != ja.end(); ++i, ++j)
            *i = ja_in[j] >= 0 ? ja_in[j] : nr-ja_in[j]-1;
    }

SparsityIaja::SparsityIaja(size_type nc,
        const FullVector<size_type>& ia,
        const FullVector<size_type>& ja):
     nr(ia.length()-1), nc(nc), nnz(ja.length()), ia(ia), ja(ja) {}

SparsityIaja::SparsityIaja(size_type nc,
        FullVector<size_type>&& ia,
        FullVector<size_type>&& ja):
    nr(ia.length()-1), nc(nc), nnz(ja.length()),
    ia(std::move(ia)), ja(std::move(ja)) {}

SparsityIaja::SparsityIaja(SparsityIaja&& rhs):
    nr(rhs.nr), nc(rhs.nc), nnz(rhs.nnz),
    ia(std::move(rhs.ia)), ja(std::move(rhs.ja)) {
    rhs.nr = rhs.nc = rhs.nnz = 0;
}

SparsityIaja::SparsityIaja(size_type nc,
        const std::vector<size_type>& ia,
        const std::vector<size_type>& ja):
    nr(ia.size()-1), nc(nc), nnz(ja.size()), ia(ia), ja(ja) {}

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


/* -------------- *
 *   Dimensions   *
 * -------------- */
void SparsityIaja::compress_storage() {
    ia.remove_trailing_zeros();
    nr = ia.length() - 1;
    size_type true_leng = *(ia.cend()-1);
    if ( nnz > true_leng ) {
        nnz = true_leng;
        ja.resize(true_leng);
    }
}


IAJA_NAMESPACE_CLOSE
