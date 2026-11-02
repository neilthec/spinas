


//File:  SPINAS/tests/cmatrix.cpp

#include "test_utils.h"
#include <cmath>
#include <complex>
#include <stdexcept>

#include "types.h"
//#include "aliases.h"
#include "utilities.h"
#include "cmatrix.h"

constexpr int NTESTS = 100;

using namespace spinas;

SPINAS_TEST_SUITE(cmatrix_tests)


SPINAS_TEST_CASE(constructor) {
  SPINAS_TEST_MESSAGE("Testing cmatrix:");
  SPINAS_TEST_MESSAGE("\t* Constructor and Equality");
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
    SPINAS_CHECK(m1 == m1);
    SPINAS_CHECK(!(m1 == m2));
    SPINAS_CHECK(!(m1 != m1));
    SPINAS_CHECK(m1 != m2);
  }
}

SPINAS_TEST_CASE(determinant) {
  SPINAS_TEST_MESSAGE("\t* Determinant");
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
    SPINAS_CHECK_SMALL(std::abs(det(m_true) - mass2), epsilon);
    SPINAS_CHECK_SMALL(std::abs(det(m_false) - mass2), epsilon);
  }
}

SPINAS_TEST_CASE(addition) {
  SPINAS_TEST_MESSAGE("\t* Addition");
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
    SPINAS_CHECK_EQUAL(m3, expected);
    
    m3 = m1;
    m3+= m2;
    SPINAS_CHECK_EQUAL(m3, expected);
  }
}

SPINAS_TEST_CASE(subtraction) {
  SPINAS_TEST_MESSAGE("\t* Subtraction");
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
    SPINAS_CHECK_EQUAL(m3, expected);
    
    m3 = m1;
    m3-= m2;
    SPINAS_CHECK_EQUAL(m3, expected);
  }
}

SPINAS_TEST_CASE(multiplication) {
  SPINAS_TEST_MESSAGE("\t* Multiplication");
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

    
    SPINAS_CHECK_EQUAL(m1 * m2, expected);
    SPINAS_CHECK_EQUAL(m2 * m1, expected);
    
    m3 = m1;
    m3 *= m2;
    SPINAS_CHECK_EQUAL(m3, expected);
    m3 = m2;
    m3 *= m1;
    SPINAS_CHECK_EQUAL(m3, expected);

    cn = choose_random_cdouble(-50, 50);
    expected_cn = cmatrix(cdouble(p[0] - p[3], 0) * cn, cdouble(-p[1], p[2]) * cn, cdouble(-p[1], -p[2]) * cn, cdouble(p[0] + p[3], 0) * cn);
    SPINAS_CHECK_EQUAL(m1 * cn, expected_cn);
    SPINAS_CHECK_EQUAL(cn * m1, expected_cn);
    
    ln1 = choose_random_ldouble(-50, 50);
    expected_ln = cmatrix(cdouble(p[0] - p[3], 0) * ln1, cdouble(-p[1], p[2]) * ln1, cdouble(-p[1], -p[2]) * ln1, cdouble(p[0] + p[3], 0) * ln1);
    SPINAS_CHECK_EQUAL(m1 * ln1, expected_ln);
    SPINAS_CHECK_EQUAL(ln1 * m1, expected_ln);

    in = choose_random_int(0, 50);
    for(int j=0;j<4;j++) p2[j] = in*p[j];
    expected_mult = cmatrix(p2, true);
    SPINAS_CHECK_EQUAL(m1 * in, expected_mult);
    SPINAS_CHECK_EQUAL(in * m1, expected_mult);
  }
}

SPINAS_TEST_CASE(division) {
  SPINAS_TEST_MESSAGE("\t* Division");
  ldouble p[4];
  cmatrix m1, expected_cn, expected_ln;
  cdouble cn;
  ldouble ln1, ln2;

  for(int i=0;i<NTESTS;i++){
    choose_random_momentum(p, -50, 50);
    m1 = cmatrix(p, true);
    
    cn = choose_random_cdouble(-50, 50);
    expected_cn = cmatrix(cdouble(p[0] - p[3], 0) / cn, cdouble(-p[1], p[2]) / cn, cdouble(-p[1], -p[2]) / cn, cdouble(p[0] + p[3], 0) / cn);
    SPINAS_CHECK_EQUAL(m1 / cn, expected_cn);
    
    ln1 = choose_random_ldouble(-50, 50);
    expected_ln = cmatrix(cdouble(p[0] - p[3], 0) / ln1, cdouble(-p[1], p[2]) / ln1, cdouble(-p[1], -p[2]) / ln1, cdouble(p[0] + p[3], 0) / ln1);
    SPINAS_CHECK_EQUAL(m1 / ln1, expected_ln);
  }
}

SPINAS_TEST_CASE(general_complex_algebra_and_aliasing) {
  // Momentum matrices are Hermitian. Use independent entries here so that
  // exchanging a row and a column cannot accidentally give the same answer.
  const cdouble a[2][2] = {{cdouble(1,2),cdouble(-3,1)},
                           {cdouble(4,-2),cdouble(2,5)}};
  const cdouble b[2][2] = {{cdouble(-2,3),cdouble(1,-4)},
                           {cdouble(3,2),cdouble(-1,1)}};
  cmatrix first(a[0][0],a[0][1],a[1][0],a[1][1]);
  const cmatrix second(b[0][0],b[0][1],b[1][0],b[1][1]);
  const cmatrix product = first*second;
  cmatrix square = first;
  square *= square;
  cmatrix doubled = first, cleared = first;
  doubled += doubled;
  cleared -= cleared;
  for(int row=0;row<2;row++)
    for(int column=0;column<2;column++){
      const cdouble expected = a[row][0]*b[0][column]+a[row][1]*b[1][column];
      const cdouble expectedSquare = a[row][0]*a[0][column]+a[row][1]*a[1][column];
      SPINAS_CHECK_EQUAL(product.get(row,column),expected);
      SPINAS_CHECK_EQUAL(square.get(row,column),expectedSquare);
      SPINAS_CHECK_EQUAL((first+second).get(row,column),a[row][column]+b[row][column]);
      SPINAS_CHECK_EQUAL((first-second).get(row,column),a[row][column]-b[row][column]);
      SPINAS_CHECK_EQUAL((-first).get(row,column),-a[row][column]);
      SPINAS_CHECK_EQUAL(doubled.get(row,column),cdouble(2,0)*a[row][column]);
      SPINAS_CHECK_EQUAL(cleared.get(row,column),cdouble(0,0));
    }

  const cdouble factor(2,-1);
  cmatrix scaled = first, divided = first;
  scaled *= factor;
  divided /= factor;
  cmatrix realScaled = first, realDivided = first, intScaled = first, intDivided = first;
  realScaled *= ldouble(-2); realDivided /= ldouble(-2);
  intScaled *= 3; intDivided /= 3;
  for(int row=0;row<2;row++)
    for(int column=0;column<2;column++){
      SPINAS_CHECK_EQUAL(scaled.get(row,column),a[row][column]*factor);
      SPINAS_CHECK_EQUAL(divided.get(row,column),a[row][column]/factor);
      SPINAS_CHECK_EQUAL(realScaled.get(row,column),a[row][column]*ldouble(-2));
      SPINAS_CHECK_EQUAL(realDivided.get(row,column),a[row][column]/ldouble(-2));
      SPINAS_CHECK_EQUAL(intScaled.get(row,column),a[row][column]*ldouble(3));
      SPINAS_CHECK_EQUAL(intDivided.get(row,column),a[row][column]/ldouble(3));
    }
}

SPINAS_TEST_CASE(invalid_indices_and_complex_determinant) {
  const cmatrix matrix(cdouble(1,1),0,0,1);
  SPINAS_CHECK_THROW(matrix.get(-1,0),std::out_of_range);
  SPINAS_CHECK_THROW(matrix.get(0,2),std::out_of_range);
  // get_det intentionally returns a real value and rejects a complex one.
  SPINAS_CHECK_THROW(matrix.get_det(),std::domain_error);
}

SPINAS_TEST_CASE(nonfinite_entries_never_compare_equal) {
  const ldouble invalid[] = {std::numeric_limits<ldouble>::quiet_NaN(),
    std::numeric_limits<ldouble>::infinity(),-std::numeric_limits<ldouble>::infinity()};
  const cmatrix zero;
  for(const ldouble value : invalid)
    for(int position=0;position<4;position++)
      for(int component=0;component<2;component++){
        cdouble entries[4] = {};
        entries[position] = component==0 ? cdouble(value,0) : cdouble(0,value);
        const cmatrix matrix(entries[0],entries[1],entries[2],entries[3]);
        SPINAS_CHECK(!(matrix==zero));
        SPINAS_CHECK(!(zero==matrix));
        SPINAS_CHECK(!(matrix==matrix));
        SPINAS_CHECK(matrix!=zero);
        SPINAS_CHECK(matrix!=matrix);
      }
}

SPINAS_TEST_SUITE_END()
