#ifndef _IAJA_CONFIG_H_
#define _IAJA_CONFIG_H_
#include <cstddef>

namespace iaja {
// parameters controlling print styles
const unsigned int PRINT_WIDTH_DOUBLE = 15;
const unsigned int PRINT_WIDTH_UNSIGNED_INT = 5;
const unsigned int PRINT_PRECISION_DOUBLE = 6;

const char  SPARSITY_FILL_PATTERN = 'X';
const char  SPARSITY_EMPTY_PATTERN = ' ';
const unsigned int SPARSITY_PRINT_WIDTH = 2;

const double FLOATING_POINT_NEARLY_ZERO = 1e-10;
using SizeType = std::size_t;
using FloatType = double;
} // namespace iaja

// DEINITIONS
#define PROFILING
#define VERBOSE
// #define SHOW_MOVE_COPY
// #define SHOW_CALLED_FCN

// macros to start and end iaja namespace
#define IAJA_NAMESPACE_OPEN namespace iaja {
#define IAJA_NAMESPACE_CLOSE }

#endif
