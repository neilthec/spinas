



//File:  SPINAS/boost/sproduct.cpp

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "types.h"
#include "utilities.h"
#include "cmatrix.h"
#include "cvector.h"
#include "particle.h"
#include "psproduct.h"


using namespace spinas;


BOOST_AUTO_TEST_SUITE(psproduct_tests)

//p1|1}}=0
BOOST_AUTO_TEST_CASE(psproduct_base_tests) {
  BOOST_TEST_MESSAGE("Testing psproduct:");
  BOOST_TEST_MESSAGE("\t* p1|1}}=0");
  constexpr ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 10000;
  ldouble m1;
  ldouble mom1[4];
  particle p1;
  psproduct ps11c = psproduct(&p1, &p1);
  for(int i=0;i<100;i++){
    m1 = choose_random_momentum(mom1,-50,50);
    p1.set_mass(m1);
    p1.set_momentum(mom1);
    ps11c.update();
    for(int j=-2;j<=2;j+=2){
      BOOST_CHECK_SMALL(std::abs(ps11c.v(j)), epsilon);
    }	
  }
}








BOOST_AUTO_TEST_SUITE_END()




