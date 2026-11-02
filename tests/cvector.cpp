


//File:  SPINAS/tests/cvector.cpp

#include "test_utils.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>
#include <stdexcept>

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

SPINAS_TEST_CASE(general_complex_outer_products_and_scalar_operations) {
  const cdouble a[2] = {cdouble(1,2),cdouble(-3,4)};
  const cdouble b[2] = {cdouble(2,-1),cdouble(5,3)};
  const cvector first(a), second(b);
  const cmatrix product = outer(first,second);
  for(int row=0;row<2;row++)
    for(int column=0;column<2;column++)
      SPINAS_CHECK_EQUAL(product.get(row,column),a[row]*b[column]);
  // These spinor contractions are bilinear: outer() must not conjugate.
  SPINAS_CHECK_EQUAL(first*second,a[0]*b[0]+a[1]*b[1]);
  const cdouble factor(2,-3);
  for(int component=0;component<2;component++){
    SPINAS_CHECK_EQUAL((-first).get(component),-a[component]);
    SPINAS_CHECK_EQUAL((first*factor).get(component),a[component]*factor);
    SPINAS_CHECK_EQUAL((factor*first).get(component),factor*a[component]);
    SPINAS_CHECK_EQUAL((first*ldouble(-2)).get(component),a[component]*ldouble(-2));
    SPINAS_CHECK_EQUAL((ldouble(-2)*first).get(component),ldouble(-2)*a[component]);
    SPINAS_CHECK_EQUAL((first*3).get(component),a[component]*ldouble(3));
    SPINAS_CHECK_EQUAL((3*first).get(component),ldouble(3)*a[component]);
    SPINAS_CHECK_EQUAL(first.get_conjugate().get_conjugate().get(component),a[component]);
  }
  const cmatrix matrix(cdouble(2,1),cdouble(-1,3),cdouble(4,-2),cdouble(3,5));
  const cvector rowProduct = first*matrix, columnProduct = matrix*first;
  for(int component=0;component<2;component++){
    SPINAS_CHECK_EQUAL(rowProduct.get(component),a[0]*matrix.get(0,component)+a[1]*matrix.get(1,component));
    SPINAS_CHECK_EQUAL(columnProduct.get(component),matrix.get(component,0)*a[0]+matrix.get(component,1)*a[1]);
  }
  SPINAS_CHECK_THROW(first.get(-1),std::out_of_range);
  SPINAS_CHECK_THROW(first.get(2),std::out_of_range);
}

SPINAS_TEST_CASE(nonfinite_entries_never_compare_equal) {
  const ldouble invalid[] = {std::numeric_limits<ldouble>::quiet_NaN(),
    std::numeric_limits<ldouble>::infinity(),-std::numeric_limits<ldouble>::infinity()};
  const cvector zero;
  for(const ldouble value : invalid)
    for(int position=0;position<2;position++)
      for(int component=0;component<2;component++){
        cdouble entries[2] = {};
        entries[position] = component==0 ? cdouble(value,0) : cdouble(0,value);
        const cvector vector(entries);
        SPINAS_CHECK(!(vector==zero));
        SPINAS_CHECK(!(zero==vector));
        SPINAS_CHECK(!(vector==vector));
        SPINAS_CHECK(vector!=zero);
        SPINAS_CHECK(vector!=vector);
      }
}

SPINAS_TEST_SUITE_END()


