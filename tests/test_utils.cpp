#include "test_utils.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

// The existing SPINAS random-momentum helpers use this generator.
namespace spinas { extern std::mt19937 gen; }

namespace spinas_test {
namespace {

struct test_case {
  std::string name;
  test_function function;
  const char* file;
  int line;
};

std::vector<test_case>& tests() {
  static std::vector<test_case> entries;
  return entries;
}

std::size_t checks = 0;
std::size_t failures = 0;
std::string current_test;

std::uint32_t read_seed(const std::string& value) {
  if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
    throw std::invalid_argument("The seed must be an unsigned 32-bit integer.");
  const unsigned long long result = std::stoull(value);
  if (result > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument("The seed must be an unsigned 32-bit integer.");
  return static_cast<std::uint32_t>(result);
}

void usage(const char* program) {
  std::cout << "Usage: " << program
    << " [--seed UINT32] [--filter NAME] [--list]\n"
    "The default seed comes from the current time and is printed for replay.\n"
    "A filter selects suite/test names containing NAME.\n";
}

} // namespace

registration::registration(const char* suite, const char* name,
                           test_function function, const char* file, int line) {
  tests().push_back(test_case{std::string(suite) + "/" + name, function, file, line});
}

void check(bool passed, const std::string& description,
           const char* file, int line) {
  ++checks;
  if (passed) return;
  ++failures;
  std::cerr << "FAIL " << current_test << "\n  " << file << ':' << line
            << ": " << description << '\n';
}

void message(const std::string& text) { std::cout << text << '\n'; }

int run_all(int argc, char** argv) {
  const std::uint64_t now = static_cast<std::uint64_t>(
    std::chrono::high_resolution_clock::now().time_since_epoch().count());
  std::uint32_t seed = static_cast<std::uint32_t>(now ^ (now >> 32));
  std::string filter;
  bool list = false;
  try {
    for (int argument = 1; argument < argc; ++argument) {
      const std::string option = argv[argument];
      if (option == "--help" || option == "-h") { usage(argv[0]); return 0; }
      if (option == "--list") { list = true; continue; }
      if (option == "--seed" || option == "--filter") {
        if (++argument >= argc)
          throw std::invalid_argument("Missing value after " + option);
        if (option == "--seed") seed = read_seed(argv[argument]);
        else filter = argv[argument];
      } else if (option.compare(0, 7, "--seed=") == 0) {
        seed = read_seed(option.substr(7));
      } else if (option.compare(0, 9, "--filter=") == 0) {
        filter = option.substr(9);
      } else throw std::invalid_argument("Unknown test option: " + option);
    }
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    usage(argv[0]);
    return 2;
  }

  std::sort(tests().begin(), tests().end(),
    [](const test_case& a, const test_case& b) { return a.name < b.name; });
  checks = failures = 0;
  std::size_t selected = 0, failed_cases = 0;
  spinas::gen.seed(seed);
  if (!list) std::cout << "SPINAS test seed: " << seed << '\n';
  for (const test_case& test : tests()) {
    if (!filter.empty() && test.name.find(filter) == std::string::npos) continue;
    ++selected;
    if (list) { std::cout << test.name << '\n'; continue; }
    current_test = test.name;
    const std::size_t before = failures;
    try { test.function(); }
    catch (const std::exception& error) {
      check(false, std::string("Unhandled exception: ") + error.what(),
        test.file, test.line);
    }
    catch (...) { check(false, "Unhandled nonstandard exception", test.file, test.line); }
    if (failures != before) ++failed_cases;
    std::cout << (failures == before ? "PASS " : "FAIL ") << test.name << '\n';
  }
  if (selected == 0) {
    std::cerr << "No tests matched the requested selection.\n";
    return 2;
  }
  if (list) return 0;
  std::cout << "SPINAS tests: " << selected - failed_cases << '/' << selected
            << " cases passed; " << checks - failures << '/' << checks
            << " checks passed.\n";
  return failures == 0 ? 0 : 1;
}

} // namespace spinas_test

#ifndef SPINAS_TEST_NO_MAIN
int main(int argc, char** argv) { return spinas_test::run_all(argc, argv); }
#endif
