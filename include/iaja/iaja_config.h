#ifndef _IAJA_CONFIG_H_
#define _IAJA_CONFIG_H_
#include <cstddef>

// parameters controlling print styles
const unsigned int PRINT_WIDTH_DOUBLE = 15;
const unsigned int PRINT_WIDTH_UNSIGNED_INT = 5;
using SizeType = std::size_t;

#define PROFILING

// macros to start and end iaja namespace
#define IAJA_NAMESPACE_OPEN namespace iaja {
#define IAJA_NAMESPACE_CLOSE }

#endif
