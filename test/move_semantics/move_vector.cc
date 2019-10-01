#include <iaja/linalg_vector.h>
#include <iaja/iaja_config.h>
#include <iostream>
using iaja::FullVector;
using iaja::SizeType;
using iaja::SparseVector;
using std::cout;
using std::endl;

int main() {

    double a[] = {1, 2, 3, 4, 5};
    SizeType na = sizeof(a)/sizeof(double);
    FullVector<double> rhs(na, a);
    FullVector<double> lhs(std::move(rhs));
    // lhs = std::move(rhs);
    // lhs(std::move(rhs));
    cout << "Testing move semantics for FullVector" << endl;
    cout << " ----- RHS -----" << endl;
    cout << "LENGTH: " << rhs.length() << rhs << endl;
    cout << " ----- LHS -----" << endl;
    cout << "LENGTH: " << lhs.length() << lhs << endl;

    SizeType ia[] = {1, 2, 3, 5,  7};
    SizeType nnz = na;
    SizeType n = 10;
    SparseVector<double> rhs_sp(n, nnz, ia, a);
    SparseVector<double> lhs_sp(rhs_sp);
    lhs_sp.print_compressed(cout) << endl;
    cout << rhs_sp.length() << endl;

}
