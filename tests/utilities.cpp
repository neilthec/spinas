#include "test_utils.h"
#include <cmath>
#include <limits>

#include "types.h"
#include "utilities.h"

using namespace spinas;

SPINAS_TEST_SUITE(utilities_tests)

SPINAS_TEST_CASE(tiny_and_signed_rotations) {
  const ldouble tolerance = 128*std::numeric_limits<ldouble>::epsilon();
  const ldouble axis[3] = {0,0,1};
  for(const ldouble angle : {ldouble(1e-12L),ldouble(-1e-12L),ldouble(-0.4L),ldouble(4)}){
    ldouble momentum[4] = {7,1,0,0};
    rotate_momentum(momentum,axis,angle);
    SPINAS_CHECK_EQUAL(momentum[0],ldouble(7));
    SPINAS_CHECK_SMALL(std::abs(momentum[1]-std::cos(angle)),tolerance);
    SPINAS_CHECK_SMALL(std::abs(momentum[2]/std::sin(angle)-ldouble(1)),tolerance);
  }
  const ldouble diagonalAxis[3] = {std::sqrt(ldouble(0.5)),0,std::sqrt(ldouble(0.5))};
  const ldouble angle = ldouble(1e-12L);
  ldouble momentum[4] = {7,1,0,0};
  rotate_momentum(momentum,diagonalAxis,angle);
  SPINAS_CHECK_SMALL(std::abs(momentum[3]/(angle*angle/4)-ldouble(1)),tolerance);
}

SPINAS_TEST_CASE(zero_and_tiny_boosts_preserve_generic_momenta) {
  const ldouble zero[3] = {0,0,0};
  for(const int sign : {-1,1}){
    ldouble momentum[4] = {sign*ldouble(7),1,-2,3};
    boost_momentum(momentum,zero);
    SPINAS_CHECK_EQUAL(momentum[0],sign*ldouble(7));
    SPINAS_CHECK_EQUAL(momentum[1],ldouble(1));
    SPINAS_CHECK_EQUAL(momentum[2],ldouble(-2));
    SPINAS_CHECK_EQUAL(momentum[3],ldouble(3));
  }
  const ldouble small = ldouble(1e-12L), velocity[3] = {small,small,0};
  ldouble momentum[4] = {0,0,1,0};
  boost_momentum(momentum,velocity);
  const ldouble tolerance = 128*std::numeric_limits<ldouble>::epsilon();
  SPINAS_CHECK_SMALL(std::abs(momentum[1]/(small*small/2)-ldouble(1)),tolerance);
  SPINAS_CHECK_SMALL(std::abs(momentum[0]/(-small)-ldouble(1)),tolerance);
}

SPINAS_TEST_CASE(boosting_an_already_fast_particle) {
  const ldouble energy = 1e6L, mass = 1;
  const ldouble magnitude = std::sqrt(energy-mass)*std::sqrt(energy+mass);
  for(const ldouble speed : {ldouble(0.9L),std::nextafter(ldouble(1),ldouble(0))}){
    const ldouble velocity[3] = {0,0,speed};
    const ldouble gamma = 1/std::sqrt(std::fma(-speed,speed,ldouble(1)));
    const ldouble expectedEnergy = gamma*std::fma(-speed,magnitude,energy);
    const ldouble expectedMomentum = gamma*std::fma(-speed,energy,magnitude);
    for(const int sign : {-1,1}){
      ldouble momentum[4] = {sign*energy,0,0,sign*magnitude};
      boost_momentum(momentum,velocity);
      const ldouble tolerance = 128*std::numeric_limits<ldouble>::epsilon()*energy;
      SPINAS_CHECK_SMALL(std::abs(momentum[0]-sign*expectedEnergy),tolerance);
      SPINAS_CHECK_SMALL(std::abs(momentum[3]-sign*expectedMomentum),tolerance);
      SPINAS_CHECK_EQUAL(momentum[1],ldouble(0));
      SPINAS_CHECK_EQUAL(momentum[2],ldouble(0));
    }
  }
}

SPINAS_TEST_SUITE_END()
