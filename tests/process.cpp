// Direct checks of the spin-1 symmetrization used by process amplitudes.
#include "test_utils.h"
#include "spinas.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <vector>

using namespace spinas;

SPINAS_TEST_SUITE(process_tests)

namespace {

typedef std::vector<int> spin_assignment;
typedef std::set<spin_assignment> assignment_set;

// Build the Cartesian product directly, without using the implementation's
// loop count, bit positions, or spin-index helpers as an oracle.
assignment_set expected_assignments(const std::vector<int>& states) {
  assignment_set assignments;
  assignments.insert(spin_assignment());
  for (int state : states) {
    assignment_set next;
    for (const spin_assignment& previous : assignments) {
      spin_assignment first = previous;
      first.push_back(state < 0 ? -1 : 1);
      first.push_back(state > 0 ? 1 : -1);
      next.insert(first);
      if (state == 0) {
        spin_assignment second = previous;
        second.push_back(-1);
        second.push_back(1);
        next.insert(second);
      }
    }
    assignments.swap(next);
  }
  return assignments;
}

spin_assignment combine(const std::vector<int>& first,
                        const std::vector<int>& second) {
  spin_assignment result;
  for (std::size_t leg = 0; leg < first.size(); ++leg) {
    result.push_back(first[leg]);
    result.push_back(second[leg]);
  }
  return result;
}

int fixed_loop_count(process& amplitude, const std::vector<int>& s) {
  switch (s.size()) {
    case 1: return amplitude.get_num_spin_loops(s[0]);
    case 2: return amplitude.get_num_spin_loops(s[0], s[1]);
    case 3: return amplitude.get_num_spin_loops(s[0], s[1], s[2]);
    case 4: return amplitude.get_num_spin_loops(s[0], s[1], s[2], s[3]);
    case 5: return amplitude.get_num_spin_loops(s[0], s[1], s[2], s[3], s[4]);
    default: return amplitude.get_num_spin_loops(s[0], s[1], s[2], s[3], s[4], s[5]);
  }
}

ldouble fixed_normalization(process& amplitude, const std::vector<int>& s) {
  switch (s.size()) {
    case 1: return amplitude.get_spin_normalization(s[0]);
    case 2: return amplitude.get_spin_normalization(s[0], s[1]);
    case 3: return amplitude.get_spin_normalization(s[0], s[1], s[2]);
    case 4: return amplitude.get_spin_normalization(s[0], s[1], s[2], s[3]);
    case 5: return amplitude.get_spin_normalization(s[0], s[1], s[2], s[3], s[4]);
    default: return amplitude.get_spin_normalization(s[0], s[1], s[2], s[3], s[4], s[5]);
  }
}

void fixed_spin_indices(process& amplitude, const std::vector<int>& s,
                        std::vector<int>& a, std::vector<int>& b, int iterator) {
  switch (s.size()) {
    case 1:
      amplitude.get_spinor_spins(s[0], a[0], b[0], iterator);
      break;
    case 2:
      amplitude.get_spinor_spins(s[0], a[0], b[0], s[1], a[1], b[1], iterator);
      break;
    case 3:
      amplitude.get_spinor_spins(s[0], a[0], b[0], s[1], a[1], b[1],
                                s[2], a[2], b[2], iterator);
      break;
    case 4:
      amplitude.get_spinor_spins(s[0], a[0], b[0], s[1], a[1], b[1],
                                s[2], a[2], b[2], s[3], a[3], b[3], iterator);
      break;
    case 5:
      amplitude.get_spinor_spins(s[0], a[0], b[0], s[1], a[1], b[1],
                                s[2], a[2], b[2], s[3], a[3], b[3],
                                s[4], a[4], b[4], iterator);
      break;
    default:
      amplitude.get_spinor_spins(s[0], a[0], b[0], s[1], a[1], b[1],
                                s[2], a[2], b[2], s[3], a[3], b[3],
                                s[4], a[4], b[4], s[5], a[5], b[5], iterator);
      break;
  }
}

} // namespace

SPINAS_TEST_CASE(spin_one_symmetrization_through_six_legs) {
  process amplitude;
  // The existing normalization uses sqrt(2.0), so allow ordinary double
  // rounding even on machines with a wider long double representation.
  const ldouble tolerance = 64 * std::max(
    ldouble(std::numeric_limits<double>::epsilon()),
    std::numeric_limits<ldouble>::epsilon());

  int configurations = 1;
  for (int legs = 1; legs <= 6; ++legs) {
    configurations *= 3;
    for (int configuration = 0; configuration < configurations; ++configuration) {
      std::vector<int> states(legs), first(legs), second(legs);
      std::vector<int> fixed_first(legs), fixed_second(legs);
      std::vector<const int*> state_pointers(legs);
      std::vector<int*> first_pointers(legs), second_pointers(legs);
      int state_code = configuration;
      for (int leg = 0; leg < legs; ++leg) {
        states[leg] = 2 * (state_code % 3 - 1);
        state_code /= 3;
        state_pointers[leg] = &states[leg];
        first_pointers[leg] = &first[leg];
        second_pointers[leg] = &second[leg];
      }

      const assignment_set expected = expected_assignments(states);
      const int count = amplitude.get_num_spin_loops(states.data(), legs);
      SPINAS_CHECK_EQUAL(count, static_cast<int>(expected.size()));
      SPINAS_CHECK_EQUAL(fixed_loop_count(amplitude, states), count);
      const ldouble normalization = amplitude.get_spin_normalization(states.data(), legs);
      const ldouble expected_normalization = 1 / std::sqrt(ldouble(expected.size()));
      SPINAS_CHECK_SMALL(normalization - expected_normalization, tolerance);
      SPINAS_CHECK_SMALL(normalization * normalization * count - 1, tolerance);
      SPINAS_CHECK_EQUAL(fixed_normalization(amplitude, states), normalization);

      assignment_set actual;
      for (int iterator = 0; iterator < count; ++iterator) {
        // Sentinels reveal any output index that was not assigned.
        std::fill(first.begin(), first.end(), 7);
        std::fill(second.begin(), second.end(), 7);
        std::fill(fixed_first.begin(), fixed_first.end(), 7);
        std::fill(fixed_second.begin(), fixed_second.end(), 7);
        amplitude.get_spinor_spins(state_pointers.data(), first_pointers.data(),
                                   second_pointers.data(), legs, iterator);
        const spin_assignment assignment = combine(first, second);
        SPINAS_CHECK(expected.count(assignment) == 1);
        SPINAS_CHECK(actual.insert(assignment).second);
        fixed_spin_indices(amplitude, states, fixed_first, fixed_second, iterator);
        SPINAS_CHECK(combine(fixed_first, fixed_second) == assignment);
      }
      SPINAS_CHECK(actual == expected);
    }
  }
}

SPINAS_TEST_SUITE_END()
