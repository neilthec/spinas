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
#include "starproduct.h"


using namespace spinas;


BOOST_AUTO_TEST_SUITE(starproduct_tests)

//[[12|p3|34>>=0 for massless p3
BOOST_AUTO_TEST_CASE(star12334_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 10000000000000;// 
  BOOST_TEST_MESSAGE("Testing starproduct:");
  BOOST_TEST_MESSAGE("\t* [[12|p3|34>>=0 for m3=0");
  ldouble m1,m2,m3,m4;
  ldouble mom1[4], mom2[4], mom3[4], mom4[4];
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
      m3=0;
      choose_random_massless_momentum(mom3,-50,50);
      particle p3=particle(mom3,m3);
      schain c3s = schain(&p3,SQUARE,3);
      schain c3a = schain(&p3,ANGLE,3);
      for(int j4=0;j4<2;j4++){
        m4=0;
        if(j4==0) choose_random_massless_momentum(mom4,-50,50);
        else m4 = choose_random_momentum(mom4,-50,50);
        particle p4=particle(mom4,m4);
        schain c4s = schain(&p4,SQUARE,3);
        schain c4a = schain(&p4,ANGLE,3);

        //All the particles are set up.  Now we can test the starproducts.
        starproduct a12334s = starproduct(&c1a,&c2a,&p3,&c3s,&c4s);
        starproduct a12343s = starproduct(&c1a,&c2a,&p3,&c4s,&c3s);
        starproduct a13324s = starproduct(&c1a,&c3a,&p3,&c2s,&c4s);
        starproduct a31324s = starproduct(&c3a,&c1a,&p3,&c2s,&c4s);
        starproduct s12334a = starproduct(&c1s,&c2s,&p3,&c3a,&c4a);
        starproduct s12343a = starproduct(&c1s,&c2s,&p3,&c4a,&c3a);
        starproduct s13324a = starproduct(&c1s,&c3s,&p3,&c2a,&c4a);
        starproduct s31324a = starproduct(&c3s,&c1s,&p3,&c2a,&c4a);

        if(j1==0&&j2==0&&j4==0){
          BOOST_CHECK_SMALL(std::abs(a12334s.v()),epsilon);
          BOOST_CHECK_SMALL(std::abs(a12343s.v()),epsilon);
          BOOST_CHECK_SMALL(std::abs(a13324s.v()),epsilon);
          BOOST_CHECK_SMALL(std::abs(a31324s.v()),epsilon);
          BOOST_CHECK_SMALL(std::abs(s12334a.v()),epsilon);
          BOOST_CHECK_SMALL(std::abs(s12343a.v()),epsilon);
          BOOST_CHECK_SMALL(std::abs(s13324a.v()),epsilon);
          BOOST_CHECK_SMALL(std::abs(s31324a.v()),epsilon);
        }
        else if((j1==0&&j2==0)||(j1==0&&j4==0)||(j2==0&&j4==0)){
          for(int ds1=-2;ds1<=2;ds1+=2){
            BOOST_CHECK_SMALL(std::abs(a12334s.v(ds1)),epsilon);
            BOOST_CHECK_SMALL(std::abs(a12343s.v(ds1)),epsilon);
            BOOST_CHECK_SMALL(std::abs(a13324s.v(ds1)),epsilon);
            BOOST_CHECK_SMALL(std::abs(a31324s.v(ds1)),epsilon);
            BOOST_CHECK_SMALL(std::abs(s12334a.v(ds1)),epsilon);
            BOOST_CHECK_SMALL(std::abs(s12343a.v(ds1)),epsilon);
            BOOST_CHECK_SMALL(std::abs(s13324a.v(ds1)),epsilon);
            BOOST_CHECK_SMALL(std::abs(s31324a.v(ds1)),epsilon);
          }
        }
        else if(j1==0||j2==0||j4==0){
          for(int ds1=-2;ds1<=2;ds1+=2)
            for(int ds2=-2;ds2<=2;ds2+=2){
              BOOST_CHECK_SMALL(std::abs(a12334s.v(ds1,ds2)),epsilon);
              BOOST_CHECK_SMALL(std::abs(a12343s.v(ds1,ds2)),epsilon);
              BOOST_CHECK_SMALL(std::abs(a13324s.v(ds1,ds2)),epsilon);
              BOOST_CHECK_SMALL(std::abs(a31324s.v(ds1,ds2)),epsilon);
              BOOST_CHECK_SMALL(std::abs(s12334a.v(ds1,ds2)),epsilon);
              BOOST_CHECK_SMALL(std::abs(s12343a.v(ds1,ds2)),epsilon);
              BOOST_CHECK_SMALL(std::abs(s13324a.v(ds1,ds2)),epsilon);
              BOOST_CHECK_SMALL(std::abs(s31324a.v(ds1,ds2)),epsilon);
            }
        }
        else{
          for(int ds1=-2;ds1<=2;ds1+=2)
            for(int ds2=-2;ds2<=2;ds2+=2)
              for(int ds4=-2;ds4<=2;ds4+=2){
                BOOST_CHECK_SMALL(std::abs(a12334s.v(ds1,ds2,ds4)),epsilon);
                BOOST_CHECK_SMALL(std::abs(a12343s.v(ds1,ds2,ds4)),epsilon);
                BOOST_CHECK_SMALL(std::abs(a13324s.v(ds1,ds2,ds4)),epsilon);
                BOOST_CHECK_SMALL(std::abs(a31324s.v(ds1,ds2,ds4)),epsilon);
                BOOST_CHECK_SMALL(std::abs(s12334a.v(ds1,ds2,ds4)),epsilon);
                BOOST_CHECK_SMALL(std::abs(s12343a.v(ds1,ds2,ds4)),epsilon);
                BOOST_CHECK_SMALL(std::abs(s13324a.v(ds1,ds2,ds4)),epsilon);
                BOOST_CHECK_SMALL(std::abs(s31324a.v(ds1,ds2,ds4)),epsilon);
              }
        }
      }
    }
  }
}




  
BOOST_AUTO_TEST_SUITE_END()