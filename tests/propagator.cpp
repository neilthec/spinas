


//File:  SPINAS/boost/propagator.cpp

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "types.h"
#include "aliases.h"
#include "utilities.h"
#include "propagator.h"


using namespace spinas;

BOOST_AUTO_TEST_SUITE(propagator_tests)

BOOST_AUTO_TEST_CASE(test_denominator) {
  BOOST_TEST_MESSAGE("Testing propagator: ");
  BOOST_TEST_MESSAGE("\t* Constructor");
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
  ldouble m = 3.87;
  ldouble width = 0.734;
  ldouble p1[4] = {13, -4, 7, -9};
  propagator prop = propagator(m, 0);
  cdouble expected = p1[0]*p1[0];
  for (int i=1;i<4;i++) expected -= p1[i]*p1[i];
  expected -= m*m;
  
  BOOST_CHECK_SMALL(std::abs(prop.denominator(p1) - expected), epsilon);
  
  prop.set_width(width);
  expected -= cdouble(0, 1) * m * width;
  BOOST_CHECK_SMALL(std::abs(prop.denominator(p1) - expected), epsilon);
  
  m = 9;
  prop.set_mass(m);
  prop.set_width(0);
  expected = p1[0]*p1[0];
  for (int i=1;i<4;i++) expected -= p1[i]*p1[i];
  expected -= m*m;
  BOOST_CHECK_SMALL(std::abs(prop.denominator(p1) - expected), epsilon);
}

BOOST_AUTO_TEST_SUITE_END()

