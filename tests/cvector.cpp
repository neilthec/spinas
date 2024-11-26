


//File:  SPINAS/boost/cvector.cpp

#include <boost/test/unit_test.hpp>
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

BOOST_AUTO_TEST_SUITE(cvector_tests)

BOOST_AUTO_TEST_CASE(constructor_test) {
  BOOST_TEST_MESSAGE("Testing cvector: ");
  BOOST_TEST_MESSAGE("\t* Constructor and Equality");
  cdouble vec0, vec1, vec2;
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



  for(int i=0;i<NTESTS;i++){
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    vec2 = choose_random_cdouble(-50,50);
    svec1 = cvector(vec0,vec1,vec2);
    
    BOOST_CHECK_EQUAL(svec1.get(0), vec0);
    BOOST_CHECK_EQUAL(svec1.get(1), vec1);
    BOOST_CHECK_EQUAL(svec1.get(2), vec2);
    
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    vec2 = choose_random_cdouble(-50,50);
    svec2 = cvector(vec0, vec1, vec2);
    
    BOOST_CHECK(svec1 == svec1);
    BOOST_CHECK(!(svec1 == svec2));
    BOOST_CHECK(!(svec1 != svec1));
    BOOST_CHECK(svec1 != svec2);
  }
}

BOOST_AUTO_TEST_CASE(conjugate_test) {
  BOOST_TEST_MESSAGE("\t* Conjugation");
  cdouble vec0, vec1, vec2;
  cvector svec;

  for(int i=0;i<NTESTS;i++){
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    svec = cvector(vec0, vec1);
    svec = svec.get_conjugate();

    BOOST_CHECK_EQUAL(svec.get(0), std::conj(vec0));
    BOOST_CHECK_EQUAL(svec.get(1), std::conj(vec1));
  }

  for(int i=0;i<NTESTS;i++){
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    vec2 = choose_random_cdouble(-50,50);
    svec = cvector(vec0, vec1, vec2);
    svec = svec.get_conjugate();

    BOOST_CHECK_EQUAL(svec.get(0), std::conj(vec0));
    BOOST_CHECK_EQUAL(svec.get(1), std::conj(vec1));
    BOOST_CHECK_EQUAL(svec.get(2), std::conj(vec2));
  }
}

BOOST_AUTO_TEST_CASE(multiplication_test) {
  BOOST_TEST_MESSAGE("\t* Multiplication");
  cdouble vec0, vec1, vec2, vec3, vec4, vec5, sp;
  cvector svec, svec2, svec3;
  ldouble p[4];
  cmatrix mat;

    

  for(int i=0;i<NTESTS;i++){
    vec0 = choose_random_cdouble(-50,50);
    vec1 = choose_random_cdouble(-50,50);
    vec2 = choose_random_cdouble(-50,50);
    svec = cvector(vec0, vec1, vec2);
    vec3 = choose_random_cdouble(-50,50);
    vec4 = choose_random_cdouble(-50,50);
    vec5 = choose_random_cdouble(-50,50);
    svec2 = cvector(vec3, vec4, vec5);
    sp = svec * svec2;
    
    BOOST_CHECK_EQUAL(sp, vec0 * vec3 + vec1 * vec4 + vec2 * vec5);
    
    mat = cmatrix(choose_random_cdouble(-50,50), choose_random_cdouble(-50,50), choose_random_cdouble(-50,50), 
      choose_random_cdouble(-50,50), choose_random_cdouble(-50,50), choose_random_cdouble(-50,50),
      choose_random_cdouble(-50,50), choose_random_cdouble(-50,50), choose_random_cdouble(-50,50));
    svec3 = mat * svec;
    
    BOOST_CHECK_EQUAL(svec3.get(0), mat.get(0, 0) * svec.get(0) + mat.get(0, 1) * svec.get(1) + mat.get(0, 2) * svec.get(2));
    BOOST_CHECK_EQUAL(svec3.get(1), mat.get(1, 0) * svec.get(0) + mat.get(1, 1) * svec.get(1) + mat.get(1, 2) * svec.get(2));
    BOOST_CHECK_EQUAL(svec3.get(2), mat.get(2, 0) * svec.get(0) + mat.get(2, 1) * svec.get(1) + mat.get(2, 2) * svec.get(2));
    
    svec3 = svec * mat;
    
    BOOST_CHECK_EQUAL(svec3.get(0), svec.get(0) * mat.get(0, 0) + svec.get(1) * mat.get(1, 0) + svec.get(2) * mat.get(2, 0));
    BOOST_CHECK_EQUAL(svec3.get(1), svec.get(0) * mat.get(0, 1) + svec.get(1) * mat.get(1, 1) + svec.get(2) * mat.get(2, 1));
    BOOST_CHECK_EQUAL(svec3.get(2), svec.get(0) * mat.get(0, 2) + svec.get(1) * mat.get(1, 2) + svec.get(2) * mat.get(2, 2));
  }
}

BOOST_AUTO_TEST_SUITE_END()


