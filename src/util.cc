#include <iaja/util.h>

#include <string>
#include <sstream>

IAJA_NAMESPACE_OPEN

std::string TimerToFile::extract_function_name(const std::string pretty_fcn) {

    std::istringstream ss(pretty_fcn);
    std::string fcn_name;

    while (ss >> fcn_name) {
        std::string::size_type start = fcn_name.find('(');
        if ( start != std::string::npos ) {
            std::string::size_type end = fcn_name.find(')');
            fcn_name.erase(start, end);
            return fcn_name;
        }
    }
    return "Unknown_fcn";
}

IAJA_NAMESPACE_CLOSE
