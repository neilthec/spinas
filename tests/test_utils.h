// Small standard C++ test utility shared by the SPINAS validations.
#ifndef SPINAS_TEST_UTILS_H
#define SPINAS_TEST_UTILS_H

#include <cmath>
#include <complex>
#include <exception>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>

namespace spinas_test {

typedef void (*test_function)();

class registration {
public:
  registration(const char* suite, const char* name, test_function function,
               const char* file, int line);
};

void check(bool passed, const std::string& description,
           const char* file, int line);
void message(const std::string& text);
int run_all(int argc, char** argv);

template<class T>
typename std::enable_if<std::is_arithmetic<T>::value, bool>::type
finite(const T& value) {
  return std::isfinite(static_cast<long double>(value));
}

template<class T>
typename std::enable_if<!std::is_arithmetic<T>::value, bool>::type
finite(const T&) { return true; }

template<class T>
bool finite(const std::complex<T>& value) {
  return std::isfinite(value.real()) && std::isfinite(value.imag());
}

template<class T>
std::string text(const T& value) {
  std::ostringstream output;
  output << std::setprecision(std::numeric_limits<long double>::max_digits10)
         << value;
  return output.str();
}

template<class A, class B>
void check_equal(const A& actual, const B& expected, const char* actual_name,
                 const char* expected_name, const char* file, int line) {
  const bool passed = finite(actual) && finite(expected) && actual == expected;
  if (passed) check(true, "", file, line);
  else check(false, std::string(actual_name) + " == " + expected_name +
    "; actual: " + text(actual) + "; expected: " + text(expected), file, line);
}

template<class A, class B>
void check_small(const A& value, const B& tolerance, const char* value_name,
                 const char* tolerance_name, const char* file, int line) {
  using std::abs;
  const bool passed = finite(value) && finite(tolerance) && tolerance >= 0 &&
    abs(value) <= tolerance;
  if (passed) check(true, "", file, line);
  else check(false, std::string("abs(") + value_name + ") <= " + tolerance_name +
    "; actual: " + text(value) + "; tolerance: " + text(tolerance), file, line);
}

template<class Exception, class Function>
void check_throw(Function function, const char* expression, const char* type,
                 const char* file, int line) {
  try { function(); }
  catch (const Exception&) { check(true, "", file, line); return; }
  catch (const std::exception& error) {
    check(false, std::string(expression) + " should throw " + type +
      "; different exception: " + error.what(), file, line);
    return;
  }
  catch (...) {
    check(false, std::string(expression) + " should throw " + type +
      "; different nonstandard exception", file, line);
    return;
  }
  check(false, std::string(expression) + " should throw " + type +
    "; no exception was thrown", file, line);
}

} // namespace spinas_test

#define SPINAS_TEST_SUITE(name) namespace name { \
  static const char* spinas_test_suite_name() { return #name; }
#define SPINAS_TEST_SUITE_END() }
#define SPINAS_TEST_CASE(name) \
  static void name(); \
  static ::spinas_test::registration name##_registration( \
    spinas_test_suite_name(), #name, &name, __FILE__, __LINE__); \
  static void name()
#define SPINAS_CHECK(expression) \
  ::spinas_test::check(static_cast<bool>(expression), #expression, __FILE__, __LINE__)
#define SPINAS_CHECK_EQUAL(actual, expected) \
  ::spinas_test::check_equal(actual, expected, #actual, #expected, __FILE__, __LINE__)
#define SPINAS_CHECK_SMALL(value, tolerance) \
  ::spinas_test::check_small(value, tolerance, #value, #tolerance, __FILE__, __LINE__)
#define SPINAS_CHECK_THROW(expression, exception_type) \
  ::spinas_test::check_throw<exception_type>([&]() { (void)(expression); }, \
    #expression, #exception_type, __FILE__, __LINE__)
#define SPINAS_TEST_MESSAGE(expression) do { \
  std::ostringstream spinas_test_message; \
  spinas_test_message << expression; \
  ::spinas_test::message(spinas_test_message.str()); \
} while (false)

#endif
