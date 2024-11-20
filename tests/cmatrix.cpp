


//File:  SPINAS/boost/cmatrix.cpp

#define BOOST_TEST_MODULE spinas_tests
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <complex>

#include "types.h"
#include "aliases.h"
#include "utilities.h"
#include "cmatrix.h"

constexpr int NTESTS = 100;

using namespace spinas;

struct GlobalFixture {
  GlobalFixture() {
  }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);



BOOST_AUTO_TEST_SUITE(cmatrix_tests)


BOOST_AUTO_TEST_CASE(constructor) {
  BOOST_TEST_MESSAGE("Testing cmatrix:");
  BOOST_TEST_MESSAGE("\t* Constructor and Equality");
  ldouble p1[4], p2[4];
  cmatrix m1, m2;

  //Test equality in random matrices
  for(int i=0;i<NTESTS;i++){
    choose_random_momentum(p1, -50, 50);
    choose_random_momentum(p2, -50, 50);
    //We are taking a very small chance that p1 == p2 could occur.
    m1 = cmatrix(p1,true);
    m2 = cmatrix(p2,true);
    //Check equality
    BOOST_CHECK(m1 == m1);
    BOOST_CHECK(!(m1 == m2));
    BOOST_CHECK(!(m1 != m1));
    BOOST_CHECK(m1 != m2);
  }
}

BOOST_AUTO_TEST_CASE(determinant) {
  BOOST_TEST_MESSAGE("\t* Determinant");
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
  ldouble p[4], mass2;
  cmatrix m_true, m_false;

  //Test Random Momenta
  for(int i=0;i<NTESTS;i++){
    choose_random_momentum(p, -50, 50);
    m_true = cmatrix(p,true);
    m_false = cmatrix(p,false);
    mass2 = p[0]*p[0];
    for (int i=1;i<4;i++) mass2 -= p[i]*p[i];
    BOOST_CHECK_SMALL(std::abs(det(m_true) - mass2), epsilon);
    BOOST_CHECK_SMALL(std::abs(det(m_false) - mass2), epsilon);
  }
}

BOOST_AUTO_TEST_CASE(addition) {
  BOOST_TEST_MESSAGE("\t* Addition");
  ldouble p1[4], p2[4], p3[4];
  cmatrix m1, m2, m3, expected;

  for(int i=0;i<NTESTS;i++){
    choose_random_momentum(p1, -50, 50);
    m1 = cmatrix(p1,true);
    choose_random_momentum(p2, -50, 50);
    m2 = cmatrix(p2,true);
    m3 = m1 + m2;
    
    for(int j=0;j<4;j++) p3[j] = p1[j] + p2[j];
    expected = cmatrix(p3,true);
    BOOST_CHECK_EQUAL(m3, expected);
    
    m3 = m1;
    m3+= m2;
    BOOST_CHECK_EQUAL(m3, expected);
  }
}

BOOST_AUTO_TEST_CASE(subtraction) {
  BOOST_TEST_MESSAGE("\t* Subtraction");
  ldouble p1[4], p2[4], p3[4];
  cmatrix m1, m2, m3, expected;
  
  for(int i=0;i<NTESTS;i++){
    choose_random_momentum(p1, -50, 50);
    m1 = cmatrix(p1,true);
    choose_random_momentum(p2, -50, 50);
    m2 = cmatrix(p2,true);
    m3 = m1 - m2;
    
    for(int j=0;j<4;j++) p3[j] = p1[j] - p2[j];
    expected = cmatrix(p3,true);
    BOOST_CHECK_EQUAL(m3, expected);
    
    m3 = m1;
    m3-= m2;
    BOOST_CHECK_EQUAL(m3, expected);
  }
}

BOOST_AUTO_TEST_CASE(multiplication) {
  BOOST_TEST_MESSAGE("\t* Multiplication");
  ldouble p[4], p2[4];
  cmatrix m1, m2, m3, expected, expected_cn, expected_ln, expected_mult;
  ldouble mass2, ln1, ln2;
  cdouble cn;
  int in;

  for(int i=0;i<NTESTS;i++){
    choose_random_momentum(p, -50, 50);
    mass2 = p[0]*p[0];
    for(int j=1;j<4;j++) mass2 -= p[j]*p[j];
    m1 = cmatrix(p,true);
    m2 = cmatrix(p,false);
    expected = cmatrix(cdouble(mass2, 0), cdouble(0, 0), cdouble(0, 0), cdouble(mass2, 0));

    
    BOOST_CHECK_EQUAL(m1 * m2, expected);
    BOOST_CHECK_EQUAL(m2 * m1, expected);
    
    m3 = m1;
    m3 *= m2;
    BOOST_CHECK_EQUAL(m3, expected);
    m3 = m2;
    m3 *= m1;
    BOOST_CHECK_EQUAL(m3, expected);

    cn = choose_random_cdouble(-50, 50);
    expected_cn = cmatrix(cdouble(p[0] - p[3], 0) * cn, cdouble(-p[1], p[2]) * cn, cdouble(-p[1], -p[2]) * cn, cdouble(p[0] + p[3], 0) * cn);
    BOOST_CHECK_EQUAL(m1 * cn, expected_cn);
    BOOST_CHECK_EQUAL(cn * m1, expected_cn);
    
    ln1 = choose_random_ldouble(-50, 50);
    expected_ln = cmatrix(cdouble(p[0] - p[3], 0) * ln1, cdouble(-p[1], p[2]) * ln1, cdouble(-p[1], -p[2]) * ln1, cdouble(p[0] + p[3], 0) * ln1);
    BOOST_CHECK_EQUAL(m1 * ln1, expected_ln);
    BOOST_CHECK_EQUAL(ln1 * m1, expected_ln);

    in = choose_random_int(0, 50);
    for(int j=0;j<4;j++) p2[j] = in*p[j];
    expected_mult = cmatrix(p2, true);
    BOOST_CHECK_EQUAL(m1 * in, expected_mult);
    BOOST_CHECK_EQUAL(in * m1, expected_mult);
  }
}

BOOST_AUTO_TEST_CASE(division) {
  BOOST_TEST_MESSAGE("\t* Division");
  ldouble p[4];
  cmatrix m1, expected_cn, expected_ln;
  cdouble cn;
  ldouble ln1, ln2;

  for(int i=0;i<NTESTS;i++){
    choose_random_momentum(p, -50, 50);
    m1 = cmatrix(p, true);
    
    cn = choose_random_cdouble(-50, 50);
    expected_cn = cmatrix(cdouble(p[0] - p[3], 0) / cn, cdouble(-p[1], p[2]) / cn, cdouble(-p[1], -p[2]) / cn, cdouble(p[0] + p[3], 0) / cn);
    BOOST_CHECK_EQUAL(m1 / cn, expected_cn);
    
    ln1 = choose_random_ldouble(-50, 50);
    expected_ln = cmatrix(cdouble(p[0] - p[3], 0) / ln1, cdouble(-p[1], p[2]) / ln1, cdouble(-p[1], -p[2]) / ln1, cdouble(p[0] + p[3], 0) / ln1);
    BOOST_CHECK_EQUAL(m1 / ln1, expected_ln);
  }
}

BOOST_AUTO_TEST_SUITE_END()
