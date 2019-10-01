
#include <iaja/linalg_matrix.h>
#include <iaja/iaja_config.h>
#include <iostream>
using iaja::FullVector;
using iaja::SizeType;
using iaja::SparsityIaja;
using std::cout;
using std::endl;

int main() {

    SizeType ia[] = {0, 4, 6, 9, 13, 16, 20};
    SizeType ja[] = {0, 2, 3, 5, 1, 3, 0, 2, 4, 0, 1, 3, 5, 2, 4, 5, 0, 3, 4, 5};
    SizeType nr, nc, nnz;
    nr = nc = sizeof(ia)/sizeof(SizeType) - 1;
    nnz = sizeof(ja)/sizeof(SizeType);
    FullVector<SizeType> iav(nr+1, ia);
    FullVector<SizeType> jav(nnz, ja);
    SparsityIaja spp_lhs(nc, iav, jav);
    SparsityIaja spp_rhs(spp_lhs);
    cout << spp_rhs << endl;
}
