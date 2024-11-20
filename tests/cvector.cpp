


//File:  SPINAS/boost/cvector.cpp

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "types.h"
#include "aliases.h"
#include "utilities.h"
#include "cmatrix.h"
#include "cvector.h"


constexpr int NTESTS = 100;

using namespace spinas;

BOOST_AUTO_TEST_SUITE(cvector_tests)

BOOST_AUTO_TEST_CASE(constructor_test) {
  BOOST_TEST_MESSAGE("Testing cvector: ");
  BOOST_TEST_MESSAGE("\t* Constructor and Equality");
  cdouble vec0, vec1;
  cvector svec1, svec2;
  ldouble ld1, ld2;

  for(int i=0;i<NTESTS;i++){
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    svec1 = cvector(vec0,vec1);
    
    BOOST_CHECK_EQUAL(svec1.get(0), vec0);
    BOOST_CHECK_EQUAL(svec1.get(1), vec1);
    
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    svec2 = cvector(vec0, vec1);
    
    BOOST_CHECK(svec1 == svec1);
    BOOST_CHECK(!(svec1 == svec2));
    BOOST_CHECK(!(svec1 != svec1));
    BOOST_CHECK(svec1 != svec2);
  }
}

BOOST_AUTO_TEST_CASE(conjugate_test) {
  BOOST_TEST_MESSAGE("\t* Conjugation");
  cdouble vec0, vec1;
  cvector svec;

  for(int i=0;i<NTESTS;i++){
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    svec = cvector(vec0, vec1);
    svec = svec.get_conjugate();

BOOST_CHECK_EQUAL(svec.get(0), std::conj(vec0));
BOOST_CHECK_EQUAL(svec.get(1), std::conj(vec1));
  }
}

BOOST_AUTO_TEST_CASE(multiplication_test) {
  BOOST_TEST_MESSAGE("\t* Multiplication");
  cdouble vec0, vec1, vec2, vec3, sp;
  cvector svec, svec2, svec3;
  ldouble p[4];
  cmatrix mat;

    

  for(int i=0;i<NTESTS;i++){
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    svec = cvector(vec0, vec1);
    vec2 = choose_random_cdouble(-50,50);
    vec3 = choose_random_cdouble(-50,50);
    svec2 = cvector(vec2, vec3);
    sp = svec * svec2;
    
    BOOST_CHECK_EQUAL(sp, vec0 * vec2 + vec1 * vec3);
    
    choose_random_momentum(p, -50, 50);
    mat = cmatrix(p, true);
    svec3 = mat * svec;
    
    BOOST_CHECK_EQUAL(svec3.get(0), mat.get(0, 0) * svec.get(0) + mat.get(0, 1) * svec.get(1));
    BOOST_CHECK_EQUAL(svec3.get(1), mat.get(1, 0) * svec.get(0) + mat.get(1, 1) * svec.get(1));
    
    svec3 = svec * mat;
    
    BOOST_CHECK_EQUAL(svec3.get(0), svec.get(0) * mat.get(0, 0) + svec.get(1) * mat.get(1, 0));
    BOOST_CHECK_EQUAL(svec3.get(1), svec.get(0) * mat.get(0, 1) + svec.get(1) * mat.get(1, 1));
  }
}

BOOST_AUTO_TEST_SUITE_END()


