//File:  SPINAS/boost/schain.cpp

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
#include "sproduct.h"
#include "schain.h"
#include "supleproduct.h"


using namespace spinas;


BOOST_AUTO_TEST_SUITE(supleproduct_tests)

//[123]=-[132] & <123>=-<132>
BOOST_AUTO_TEST_CASE(stp123stp132_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 10000000;// 
  BOOST_TEST_MESSAGE("Testing supleproduct:");
  BOOST_TEST_MESSAGE("\t* [123]=-[132] & <123>=-<132>");
  ldouble m1,m2,m3;
  ldouble mom1[4], mom2[4], mom3[4];
  for(int i=0;i<100;i++)
  for(int j1=0;j1<2;j1++){
    m1=0;
    if(j1==0) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    particle p1=particle(mom1,m1);
    schain c1s = schain(&p1,SQUARE,3);
    schain c1a = schain(&p1,ANGLE,3);
    for(int j2=0;j2<2;j2++){
      m2=0;
      if(j2==0) choose_random_massless_momentum(mom2,-50,50);
      else m2 = choose_random_momentum(mom2,-50,50);
      particle p2=particle(mom2,m2);
      schain c2s = schain(&p2,SQUARE,3);
      schain c2a = schain(&p2,ANGLE,3);
      for(int j3=0;j3<2;j3++){
        m3=0;
        if(j3==0) choose_random_massless_momentum(mom3,-50,50);
        else m3 = choose_random_momentum(mom3,-50,50);
        particle p3=particle(mom3,m3);
        schain c3s = schain(&p3,SQUARE,3);
        schain c3a = schain(&p3,ANGLE,3);

        //All the particles are set up.  Now we can test the supleproducts.
        supleproduct s123 = supleproduct(&c1s,&c2s,&c3s);
        supleproduct s132 = supleproduct(&c1s,&c3s,&c2s);
        supleproduct a123 = supleproduct(&c1a,&c2a,&c3a);
        supleproduct a132 = supleproduct(&c1a,&c3a,&c2a);

        if(j1==0&&j2==0&&j3==0){
          BOOST_CHECK_SMALL(std::abs(s123.v()+s132.v()),epsilon);
          BOOST_CHECK_SMALL(std::abs(a123.v()+a132.v()),epsilon);
        }
        else if((j1==0&&j2==0)||(j1==0&&j3==0)||(j2==0&&j3==0)){
          for(int ds1=-2;ds1<=2;ds1+=2){
            BOOST_CHECK_SMALL(std::abs(s123.v(ds1)+s132.v(ds1)),epsilon);
            BOOST_CHECK_SMALL(std::abs(a123.v(ds1)+a132.v(ds1)),epsilon);
          }
        }
        else if(j1==0){
          for(int ds2=-2;ds2<=2;ds2+=2)
            for(int ds3=-2;ds3<=2;ds3+=2){
              BOOST_CHECK_SMALL(std::abs(s123.v(ds2,ds3)+s132.v(ds3,ds2)),epsilon);
              BOOST_CHECK_SMALL(std::abs(a123.v(ds2,ds3)+a132.v(ds3,ds2)),epsilon);
            }
        }
        else if(j2==0||j3==0){
          for(int ds1=-2;ds1<=2;ds1+=2)
            for(int ds2=-2;ds2<=2;ds2+=2){
              BOOST_CHECK_SMALL(std::abs(s123.v(ds1,ds2)+s132.v(ds1,ds2)),epsilon);
              BOOST_CHECK_SMALL(std::abs(a123.v(ds1,ds2)+a132.v(ds1,ds2)),epsilon);
            }
        }
        else{
          for(int ds1=-2;ds1<=2;ds1+=2)
            for(int ds2=-2;ds2<=2;ds2+=2)
              for(int ds3=-2;ds3<=2;ds3+=2){
                BOOST_CHECK_SMALL(std::abs(s123.v(ds1,ds2,ds3)+s132.v(ds1,ds3,ds2)),epsilon);
                BOOST_CHECK_SMALL(std::abs(a123.v(ds1,ds2,ds3)+a132.v(ds1,ds3,ds2)),epsilon);
              }
        }
      }
    }
  }
}

//[111]=0 & <111>=0
BOOST_AUTO_TEST_CASE(stp111_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;// 
  BOOST_TEST_MESSAGE("Testing supleproduct:");
  BOOST_TEST_MESSAGE("\t* [111]=0 & <111>=0");
  ldouble m1=0;
  ldouble mom1[4];
  for(int i=0;i<100;i++){
    choose_random_massless_momentum(mom1,-50,50);
    particle p1=particle(mom1,m1);
    schain c1s = schain(&p1,SQUARE,3);
    schain c1a = schain(&p1,ANGLE,3);
    //All the particles are set up.  Now we can test the supleproducts.
    supleproduct s111 = supleproduct(&c1s,&c1s,&c1s);
    supleproduct a111 = supleproduct(&c1a,&c1a,&c1a);
    BOOST_CHECK_SMALL(std::abs(s111.v()),epsilon);
    BOOST_CHECK_SMALL(std::abs(a111.v()),epsilon);
  }
}

//[111]=m1^3epsilon & <111>=m1^3epsilon
BOOST_AUTO_TEST_CASE(stp111m3_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;// 
  BOOST_TEST_MESSAGE("Testing supleproduct:");
  BOOST_TEST_MESSAGE("\t* [111]=m1^3*epsilon & <111>=m1^3*epsilon");
  ldouble m1;
  ldouble mom1[4];
  for(int i=0;i<100;i++){
    m1 = choose_random_momentum(mom1,-50,50);
    particle p1=particle(mom1,m1);
    schain c1s = schain(&p1,SQUARE,3);
    schain c1a = schain(&p1,ANGLE,3);
    //All the particles are set up.  Now we can test the supleproducts.
    supleproduct s111 = supleproduct(&c1s,&c1s,&c1s);
    supleproduct a111 = supleproduct(&c1a,&c1a,&c1a);
      BOOST_CHECK_SMALL(std::abs(s111.v(-2,0,2)-m1*m1*m1),epsilon);
      BOOST_CHECK_SMALL(std::abs(a111.v(-2,0,2)-m1*m1*m1),epsilon);
      BOOST_CHECK_SMALL(std::abs(s111.v(0,2,-2)-m1*m1*m1),epsilon);
      BOOST_CHECK_SMALL(std::abs(a111.v(0,2,-2)-m1*m1*m1),epsilon);
      BOOST_CHECK_SMALL(std::abs(s111.v(2,-2,0)-m1*m1*m1),epsilon);
      BOOST_CHECK_SMALL(std::abs(a111.v(2,-2,0)-m1*m1*m1),epsilon);

      BOOST_CHECK_SMALL(std::abs(s111.v(-2,2,0)+m1*m1*m1),epsilon);
      BOOST_CHECK_SMALL(std::abs(a111.v(-2,2,0)+m1*m1*m1),epsilon);
      BOOST_CHECK_SMALL(std::abs(s111.v(0,-2,2)+m1*m1*m1),epsilon);
      BOOST_CHECK_SMALL(std::abs(a111.v(0,-2,2)+m1*m1*m1),epsilon);
      BOOST_CHECK_SMALL(std::abs(s111.v(2,0,-2)+m1*m1*m1),epsilon);
      BOOST_CHECK_SMALL(std::abs(a111.v(2,0,-2)+m1*m1*m1),epsilon);

      BOOST_CHECK_SMALL(std::abs(s111.v(2,2,0)),epsilon);
      BOOST_CHECK_SMALL(std::abs(a111.v(2,2,0)),epsilon);
  }
}

  
BOOST_AUTO_TEST_SUITE_END()