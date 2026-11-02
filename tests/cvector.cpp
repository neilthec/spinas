


//File:  SPINAS/tests/cvector.cpp

#include "test_utils.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "types.h"
//#include "aliases.h"
#include "utilities.h"
#include "cmatrix.h"
#include "cvector.h"


constexpr int NTESTS = 100;

using namespace spinas;

SPINAS_TEST_SUITE(cvector_tests)

SPINAS_TEST_CASE(constructor_test) {
  SPINAS_TEST_MESSAGE("Testing cvector: ");
  SPINAS_TEST_MESSAGE("\t* Constructor and Equality");
  cdouble vec0, vec1;
  cvector svec1, svec2;
  ldouble ld1, ld2;

  for(int i=0;i<NTESTS;i++){
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    svec1 = cvector(vec0,vec1);
    
    SPINAS_CHECK_EQUAL(svec1.get(0), vec0);
    SPINAS_CHECK_EQUAL(svec1.get(1), vec1);
    
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    svec2 = cvector(vec0, vec1);
    
    SPINAS_CHECK(svec1 == svec1);
    SPINAS_CHECK(!(svec1 == svec2));
    SPINAS_CHECK(!(svec1 != svec1));
    SPINAS_CHECK(svec1 != svec2);
  }
}

SPINAS_TEST_CASE(conjugate_test) {
  SPINAS_TEST_MESSAGE("\t* Conjugation");
  cdouble vec0, vec1;
  cvector svec;

  for(int i=0;i<NTESTS;i++){
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    svec = cvector(vec0, vec1);
    svec = svec.get_conjugate();

SPINAS_CHECK_EQUAL(svec.get(0), std::conj(vec0));
SPINAS_CHECK_EQUAL(svec.get(1), std::conj(vec1));
  }
}

SPINAS_TEST_CASE(multiplication_test) {
  SPINAS_TEST_MESSAGE("\t* Multiplication");
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
    
    SPINAS_CHECK_EQUAL(sp, vec0 * vec2 + vec1 * vec3);
    
    choose_random_momentum(p, -50, 50);
    mat = cmatrix(p, true);
    svec3 = mat * svec;
    
    SPINAS_CHECK_EQUAL(svec3.get(0), mat.get(0, 0) * svec.get(0) + mat.get(0, 1) * svec.get(1));
    SPINAS_CHECK_EQUAL(svec3.get(1), mat.get(1, 0) * svec.get(0) + mat.get(1, 1) * svec.get(1));
    
    svec3 = svec * mat;
    
    SPINAS_CHECK_EQUAL(svec3.get(0), svec.get(0) * mat.get(0, 0) + svec.get(1) * mat.get(1, 0));
    SPINAS_CHECK_EQUAL(svec3.get(1), svec.get(0) * mat.get(0, 1) + svec.get(1) * mat.get(1, 1));
  }
}

SPINAS_TEST_SUITE_END()


