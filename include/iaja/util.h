#include <iaja/iaja_config.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

IAJA_NAMESPACE_OPEN

#ifdef PROFILING
    #define TIMER_BEGIN TimerToFile timer(__PRETTY_FUNCTION__);
    #define TIMER_END   timer.write_to_file();
#else 
    #define TIMER_BEGIN
    #define TIMER_END
#endif


class Timer {
 /* ----------------------- *
  * High resolution timing  *
  * ----------------------- */
 public:
  Timer():timer_begin(clock_t::now()) {}

  void reset() {timer_begin = clock_t::now();}

  double elapsed() const {
      return std::chrono::duration_cast<second_t>
          (clock_t::now() - timer_begin).count();
  }

 protected:
  using clock_t = std::chrono::high_resolution_clock;
  using second_t = std::chrono::duration< double, std::ratio<1> >;
  std::chrono::time_point<clock_t> timer_begin;
};


class TimerToFile : public Timer {
 /* -------------------------------------- *
  * Output Timing to a file                *
  * whose name is the name of the function *
  * from which the timer is called         *
  * -------------------------------------- */
 public:
  explicit TimerToFile(std::string func_name):
      Timer(),
      output_file(extract_function_name(func_name)+".timer", std::ios_base::app) {}

  // if you forget to end the timer by write_to_file()
  // the timer will record the execution
  // time up to when the timer goes out of scope
  ~TimerToFile() { write_to_file(); }

  void write_to_file() {
      output_file << std::scientific
      << std::setprecision(PRINT_PRECISION_DOUBLE)
      << elapsed() << "\n"; output_file.close();
  }

 private:
  std::ofstream output_file;
  std::string extract_function_name(const std::string pretty_fcn);
};

IAJA_NAMESPACE_CLOSE
