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

static ldouble minkowski_product(const ldouble a[4], const ldouble b[4]) {
  return a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];
}

SPINAS_TEST_CASE(lorentz_invariants_and_inverse_transformations) {
  const ldouble axis[3] = {ldouble(2)/7,ldouble(-3)/7,ldouble(6)/7};
  const ldouble angle = ldouble(0.63L);
  const ldouble velocity[3] = {ldouble(0.25L),ldouble(-0.3L),ldouble(0.15L)};
  const ldouble inverse[3] = {-velocity[0],-velocity[1],-velocity[2]};
  const ldouble tolerance = 4096*std::numeric_limits<ldouble>::epsilon();
  for(int sample=0;sample<25;sample++){
    ldouble originalA[4], originalB[4];
    choose_random_momentum(originalA,-5,5);
    choose_random_momentum(originalB,-5,5);
    if(sample%2==0) originalB[0] = -originalB[0];
    const ldouble aa = minkowski_product(originalA,originalA);
    const ldouble bb = minkowski_product(originalB,originalB);
    const ldouble ab = minkowski_product(originalA,originalB);
    // Check rotations, boosts, and their noncommuting composition separately.
    for(int transformation=0;transformation<3;transformation++){
      ldouble a[4], b[4];
      for(int component=0;component<4;component++){
        a[component] = originalA[component];
        b[component] = originalB[component];
      }
      if(transformation!=1){
        rotate_momentum(a,axis,angle);
        rotate_momentum(b,axis,angle);
      }
      if(transformation!=0){
        boost_momentum(a,velocity);
        boost_momentum(b,velocity);
      }
      SPINAS_CHECK_SMALL(std::abs(minkowski_product(a,a)-aa),tolerance*(1+std::abs(aa)));
      SPINAS_CHECK_SMALL(std::abs(minkowski_product(b,b)-bb),tolerance*(1+std::abs(bb)));
      SPINAS_CHECK_SMALL(std::abs(minkowski_product(a,b)-ab),tolerance*(1+std::abs(ab)));
      if(transformation!=0){
        boost_momentum(a,inverse);
        boost_momentum(b,inverse);
      }
      if(transformation!=1){
        rotate_momentum(a,axis,-angle);
        rotate_momentum(b,axis,-angle);
      }
      for(int component=0;component<4;component++){
        SPINAS_CHECK_SMALL(std::abs(a[component]-originalA[component]),tolerance*(1+std::abs(originalA[component])));
        SPINAS_CHECK_SMALL(std::abs(b[component]-originalB[component]),tolerance*(1+std::abs(originalB[component])));
      }
    }
  }
}

SPINAS_TEST_SUITE_END()
