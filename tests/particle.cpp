


//File:  SPINAS/boost/particle.cpp

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>
#include <limits>

#include "types.h"
#include "aliases.h"
#include "utilities.h"
#include "cmatrix.h"
#include "cvector.h"
#include "particle.h"


constexpr int NTESTS = 100;

using namespace spinas;


bool test_cdouble_equality(const cdouble& z1, const cdouble& z2) {
  const ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
  cdouble z3 = z1 - z2;
  if (std::abs(z3) > epsilon) return false;
  return true;
}


BOOST_AUTO_TEST_SUITE(particle_tests)

BOOST_AUTO_TEST_CASE(constructor_test) {
  BOOST_TEST_MESSAGE("Testing particle: ");
  BOOST_TEST_MESSAGE("\t* Constructor");
  ldouble mom1[4], mom2[4], mass2;
  particle p1;

  for(int i=0;i<NTESTS;i++){
    //Choose Mom1
    choose_random_momentum(mom1, -50, 50);
    mass2 = mom1[0]*mom1[0];
    for(int j=1;j<4;j++) mass2 -= mom1[j]*mom1[j];
    p1 = particle(mom1, std::sqrt(mass2));
    BOOST_CHECK(p1.test_angles());
    
    //Choose Mom2 with same mass
    choose_random_momentum(mom2, -50, 50);
    mom2[0] = mass2;
    for(int j=1;j<4;j++) mom2[0] += mom2[j]*mom2[j];
    mom2[0] = std::sqrt(mom2[0]);
    p1.set_momentum(mom2);
    BOOST_CHECK(p1.test_angles());
  }
}

BOOST_AUTO_TEST_CASE(momentum_matrices_test) {
  BOOST_TEST_MESSAGE("\t* Momentum Matrices");
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
  int i = 0;
  cdouble zero = cdouble(0, 0), one = cdouble(1, 0);
  ldouble mom1[4] = {13, -4, 7, -9}, mom3[4] = {6, 1, -2, std::sqrt(8)};
  particle p1 = particle(mom3, std::sqrt(23));
  p1.set_momentum(mom1);
  BOOST_CHECK_EQUAL(p1.umat(), cmatrix(mom1, true));
  BOOST_CHECK_EQUAL(p1.lmat(), cmatrix(mom1, false));
  
  // Determinant
  cdouble p1detm23 = p1.umat().get_det() - cdouble(23, 0);
  BOOST_CHECK_SMALL(std::abs(p1detm23), epsilon);
  BOOST_CHECK_EQUAL(p1.lmat(), cmatrix(mom1, false));
  
  p1detm23 = p1.lmat().get_det() - cdouble(23, 0);
  BOOST_CHECK_SMALL(std::abs(p1detm23), epsilon);
  
  // Anticommutation: p1*p2 + p2*p1 = 2p1.p2
  ldouble mom2[4] = {5, 1, -2, 3};
  particle p2 = particle(std::sqrt(25 - 1 - 4 - 9));
  p2.set_momentum(mom2);
  cmatrix id = cmatrix(one, zero, zero, one);
  cmatrix anticom = p1.umat() * p2.lmat() + p2.umat() * p1.lmat();
  BOOST_CHECK_EQUAL(anticom, 2 * p1.dot(p2) * id);
  anticom = p1.lmat() * p2.umat() + p2.lmat() * p1.umat();
  BOOST_CHECK_EQUAL(anticom,  2 * p1.dot(p2) * id);
}

BOOST_AUTO_TEST_CASE(helicity_spinors_test) {
  BOOST_TEST_MESSAGE("\t* Helicity Spinors");
  ldouble mom1[4] = {std::sqrt(16 + 49 + 81), -4, 7, -9};
  particle p1 = particle(mom1, 0);
  
  //|1> == [1|*
  BOOST_CHECK_EQUAL(p1.rangle(), p1.lsquare().get_conjugate());
  
  //<1| == |1]*
  BOOST_CHECK_EQUAL(p1.langle(), p1.rsquare().get_conjugate());
  
  // p1|1> == 0
  BOOST_CHECK(p1.umat() * p1.rangle() == cvector(cdouble(0, 0), cdouble(0, 0)));
  
  // [1|p1 == 0
  BOOST_CHECK(p1.lsquare() * p1.umat() == cvector(cdouble(0, 0), cdouble(0, 0)));
  
  // <11> == 0
  BOOST_CHECK(p1.langle() * p1.rangle() == cdouble(0, 0));
  
  // [11] == 0
  BOOST_CHECK(p1.lsquare() * p1.rsquare() == cdouble(0, 0));
  
  // |1>[1| = p1
  BOOST_CHECK(outer(p1.rangle(), p1.lsquare()) == p1.lmat());
  
  // |1]<1| = p1
  BOOST_CHECK(outer(p1.rsquare(), p1.langle()) == p1.umat());
  
  // p1|1] == 0
  BOOST_CHECK(p1.lmat() * p1.rsquare() == cvector(cdouble(0, 0), cdouble(0, 0)));
  
  // <1|p1 == 0
  BOOST_CHECK(p1.langle() * p1.lmat() == cvector(cdouble(0, 0), cdouble(0, 0)));
}


BOOST_AUTO_TEST_CASE(spin_spinors) {
  BOOST_TEST_MESSAGE("\t* Spin Spinors");
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
  ldouble mom1[4] = {std::sqrt(100 + 16 + 49 + 81), -4, 7, -9};
  particle p1 = particle(mom1, 10);
  ldouble ten = 10;
  cmatrix epsU = cmatrix(cdouble(0, 0), cdouble(1, 0), cdouble(-1, 0), cdouble(0, 0));
  cmatrix epsL = cmatrix(cdouble(0,0),cdouble(-1,0),cdouble(1,0),cdouble(0,0));
  cmatrix delta = cmatrix(cdouble(1, 0), cdouble(0, 0), cdouble(0, 0), cdouble(1, 0));
  
  // Conjugation
  for (int j = -1; j < 2; j = j + 2){
    //|1>^I == ([1|_I)*
    BOOST_CHECK_EQUAL(p1.rangle(j), p1.lsquare(j, false).get_conjugate());

    //<1|^I == (|1]_I)*
    BOOST_CHECK_EQUAL(p1.langle(j), p1.rsquare(j, false).get_conjugate());

    //|1>_I == -([1|^I)*
    BOOST_CHECK_EQUAL(p1.rangle(j, false), -p1.lsquare(j).get_conjugate());

    //<1|_I == -(|1]^I)*
    BOOST_CHECK_EQUAL(p1.langle(j, false), -p1.rsquare(j).get_conjugate());
  }
  
  

  // Spinor Products
  for(int i=-1;i<2;i=i+2)
    for(int j=-1;j<2;j=j+2){
      int epsi = 0, epsj=0;
      if(i==1) epsi=1;
      if(j==1) epsj=1;

      //<11>^{IJ} == -m eps^{IJ}
      BOOST_CHECK_SMALL(std::abs(p1.langle(i)*p1.rangle(j) - (-ten*epsU.get(epsi,epsj))), epsilon);
      
      //<11>^I_J == -m delta^I_J
      BOOST_CHECK_SMALL(std::abs(p1.langle(i)*p1.rangle(j,false) - (-ten*delta.get(epsi,epsj))), epsilon);

      //<11>_I^J == m delta_I^J
      BOOST_CHECK_SMALL(std::abs(p1.langle(i,false)*p1.rangle(j) - ten*delta.get(epsi,epsj)), epsilon);

      //<11>_{IJ} == m eps_{IJ}
      BOOST_CHECK_SMALL(std::abs(p1.langle(i,false)*p1.rangle(j,false) - ten*epsL.get(epsi,epsj)), epsilon);

      //[11]_{IJ} == -m eps_{IJ}
      BOOST_CHECK_SMALL(std::abs(p1.lsquare(i,false)*p1.rsquare(j,false) - (-ten*epsL.get(epsi,epsj))), epsilon);

      //[11]_I^J == -m delta_I^J
      BOOST_CHECK_SMALL(std::abs(p1.lsquare(i,false)*p1.rsquare(j) - (-ten*delta.get(epsi,epsj))), epsilon);

      //[11]^I_J == m delta^I_J
      BOOST_CHECK_SMALL(std::abs(p1.lsquare(i)*p1.rsquare(j,false) - ten*delta.get(epsi,epsj)), epsilon);

      //[11]^{IJ} == m eps^{IJ}
      BOOST_CHECK_SMALL(std::abs(p1.lsquare(i)*p1.rsquare(j) - ten*epsU.get(epsi,epsj)), epsilon);
      
    }

  
  //Outer Products
  //|1>^I[1|_I = p1
  cmatrix outerProd = outer(p1.rangle(-1), p1.lsquare(-1, false)) +
    outer(p1.rangle(1), p1.lsquare(1, false));
  BOOST_CHECK(outerProd == p1.lmat());
  
  //|1>_I[1|^I = -p1
  outerProd = outer(p1.rangle(-1, false), p1.lsquare(-1)) +
    outer(p1.rangle(1, false), p1.lsquare(1));
  BOOST_CHECK(outerProd == -p1.lmat());
  
  //|1]_I<1|^I = p1
  outerProd = outer(p1.rsquare(-1, false), p1.langle(-1)) +
    outer(p1.rsquare(1, false), p1.langle(1));
  BOOST_CHECK(outerProd == p1.umat());
  
  //|1]^I<1|_I = -p1
  outerProd = outer(p1.rsquare(-1), p1.langle(-1, false)) +
    outer(p1.rsquare(1), p1.langle(1, false));
  BOOST_CHECK(outerProd == -p1.umat());

  // Momentum times Spinor Gives Mass times Spinor
  for(int j=-1;j<2;j=j+2){
    //p1|1>^I=-m|1]^I
    BOOST_CHECK_EQUAL(p1.umat()*p1.rangle(j), -ten*p1.rsquare(j));

    //<1|^I*p1=[1|^I*m
    BOOST_CHECK_EQUAL(p1.langle(j)*p1.lmat(), ten*p1.lsquare(j));

    //p1|1>_I=-m|1]_I
    BOOST_CHECK_EQUAL(p1.umat()*p1.rangle(j,false), -ten*p1.rsquare(j,false));

    //<1|_I*p1=[1|_I*m
    BOOST_CHECK_EQUAL(p1.langle(j,false)*p1.lmat(), ten*p1.lsquare(j,false));

    //p1|1]^I=-m|1>^I
    BOOST_CHECK_EQUAL(p1.lmat()*p1.rsquare(j), -ten*p1.rangle(j));

    //[1|^I*p1=<1|^I*m
    BOOST_CHECK_EQUAL(p1.lsquare(j)*p1.umat(), ten*p1.langle(j));

    //p1|1]_I=-m|1>_I
    BOOST_CHECK_EQUAL(p1.lmat()*p1.rsquare(j,false), -ten*p1.rangle(j,false));

    //[1|_I*p1=<1|_I*m
    BOOST_CHECK_EQUAL(p1.lsquare(j,false)*p1.umat(), ten*p1.langle(j,false));
  }

  
  
}



BOOST_AUTO_TEST_SUITE_END()

