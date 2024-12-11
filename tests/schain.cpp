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


using namespace spinas;


BOOST_AUTO_TEST_SUITE(schain_tests)

//[12]=[12] & <12>=<12>
BOOST_AUTO_TEST_CASE(s12s_a12a_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1;// 
  BOOST_TEST_MESSAGE("Testing schain:");
  BOOST_TEST_MESSAGE("\t* [12]=[12] & <12>=<12>");
  ldouble m1,m2;
  ldouble mom1[4], mom2[4];
  int ni;
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m1=0;
    m2=0;
    if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom2,-50,50);
    else m2 = choose_random_momentum(mom2,-50,50);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    sproduct s12s = sproduct(SQUARE,&p1,&p2,3);	
    sproduct a12a = sproduct(ANGLE,&p1,&p2,3);
    schain c2s = schain(&p2,SQUARE,3);
    schain c2a = schain(&p2,ANGLE,3);
    if(o==0)
      for(int ds1=-2;ds1<=2;ds1+=2)
        for(int ds2=-2;ds2<=2;ds2+=2){
          BOOST_CHECK_SMALL(std::abs(s12s.v(ds1,ds2)-p1.lsquare(ds1,3)*c2s.v(ds2)),epsilon);
          BOOST_CHECK_SMALL(std::abs(a12a.v(ds1,ds2)-p1.langle(ds1,3)*c2a.v(ds2)),epsilon);
        }
    else if(o==1)
      for(int ds2=-2;ds2<=2;ds2+=2){
        BOOST_CHECK_SMALL(std::abs(s12s.v(ds2)-p1.lsquare(3)*c2s.v(ds2)),epsilon);
        BOOST_CHECK_SMALL(std::abs(a12a.v(ds2)-p1.langle(3)*c2a.v(ds2)),epsilon);
      }
    else if(o==2)
      for(int ds1=-2;ds1<=2;ds1+=2){
        BOOST_CHECK_SMALL(std::abs(s12s.v(ds1)-p1.lsquare(ds1,3)*c2s.v()),epsilon);
        BOOST_CHECK_SMALL(std::abs(a12a.v(ds1)-p1.langle(ds1,3)*c2a.v()),epsilon);
      }
    else if(o==3){
      BOOST_CHECK_SMALL(std::abs(s12s.v()-p1.lsquare(3)*c2s.v()),epsilon);
      BOOST_CHECK_SMALL(std::abs(a12a.v()-p1.langle(3)*c2a.v()),epsilon);
    }
  }
}



//[1|p3|2>=[1|p3|2> & <1|p3|2]=<1|p3|2]
BOOST_AUTO_TEST_CASE(s132a_a132s_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000000;// 
  BOOST_TEST_MESSAGE("Testing schain:");
  BOOST_TEST_MESSAGE("\t* [1|p3|2>=[1|p3|2> & <1|p3|2]=<1|p3|2]");
  ldouble m1,m2, m3;
  ldouble mom1[4], mom2[4], mom3[4];
  int ni;
  for(int i=0;i<100;i++)
  for(int j=0;j<2;j++){
    m3=0;
    if(j==1) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50); 
    particle p3=particle(mom3,m3);
    for(int o=0;o<4;o++){
      m1=0;
      m2=0;
      if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
      else m1 = choose_random_momentum(mom1,-50,50);
      if(o==2||o==3) choose_random_massless_momentum(mom2,-50,50);
      else m2 = choose_random_momentum(mom2,-50,50);
      particle p1=particle(mom1,m1);
      particle p2=particle(mom2,m2);
      sproduct s132a = sproduct(SQUARE,&p1,&p3,&p2,3);	
      sproduct a132s = sproduct(ANGLE,&p1,&p3,&p2,3);
      schain c32s = schain(&p3,&p2,SQUARE,3);
      schain c32a = schain(&p3,&p2,ANGLE,3);
      if(o==0)
        for(int ds1=-2;ds1<=2;ds1+=2)
          for(int ds2=-2;ds2<=2;ds2+=2){
            BOOST_CHECK_SMALL(std::abs(s132a.v(ds1,ds2)-p1.lsquare(ds1,3)*c32a.v(ds2)),epsilon);
            BOOST_CHECK_SMALL(std::abs(a132s.v(ds1,ds2)-p1.langle(ds1,3)*c32s.v(ds2)),epsilon);
          }
      else if(o==1)
        for(int ds2=-2;ds2<=2;ds2+=2){
          BOOST_CHECK_SMALL(std::abs(s132a.v(ds2)-p1.lsquare(3)*c32a.v(ds2)),epsilon);
          BOOST_CHECK_SMALL(std::abs(a132s.v(ds2)-p1.langle(3)*c32s.v(ds2)),epsilon);
        }
      else if(o==2)
        for(int ds1=-2;ds1<=2;ds1+=2){
          BOOST_CHECK_SMALL(std::abs(s132a.v(ds1)-p1.lsquare(ds1,3)*c32a.v()),epsilon);
          BOOST_CHECK_SMALL(std::abs(a132s.v(ds1)-p1.langle(ds1,3)*c32s.v()),epsilon);
        }
      else if(o==3){
        BOOST_CHECK_SMALL(std::abs(s132a.v()-p1.lsquare(3)*c32a.v()),epsilon);
        BOOST_CHECK_SMALL(std::abs(a132s.v()-p1.langle(3)*c32s.v()),epsilon);
      }
    }
  }
}
  



//[1|p3p4|2]=[1|p3p4|2] & <1|p3p4|2>=<1|p3p4|2>
BOOST_AUTO_TEST_CASE(s1342s_a1342a_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 10000000000000;// 
  BOOST_TEST_MESSAGE("Testing schain:");
  BOOST_TEST_MESSAGE("\t* [1|p3p4|2]=[1|p3p4|2] & <1|p3p4|2>=<1|p3p4|2>");
  ldouble m1,m2, m3,m4;
  ldouble mom1[4], mom2[4], mom3[4], mom4[4];
  for(int i=0;i<100;i++)
  for(int j=0;j<4;j++){
    m3=0;
    m4=0;
    if(j==1||j==3) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50); 
    if(j==2||j==3) choose_random_massless_momentum(mom4,-50,50);
    else m4 = choose_random_momentum(mom4,-50,50);
    particle p3=particle(mom3,m3);
    particle p4=particle(mom4,m4);
    for(int o=0;o<4;o++){
      m1=0;
      m2=0;
      if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
      else m1 = choose_random_momentum(mom1,-50,50);
      if(o==2||o==3) choose_random_massless_momentum(mom2,-50,50);
      else m2 = choose_random_momentum(mom2,-50,50);
      particle p1=particle(mom1,m1);
      particle p2=particle(mom2,m2);
      sproduct s1342s = sproduct(SQUARE,&p1,&p3,&p4,&p2,3);	
      sproduct a1342a = sproduct(ANGLE,&p1,&p3,&p4,&p2,3);
      schain c342s = schain(&p3,&p4,&p2,SQUARE,3);
      schain c342a = schain(&p3,&p4,&p2,ANGLE,3);
      if(o==0)
        for(int ds1=-2;ds1<=2;ds1+=2)
          for(int ds2=-2;ds2<=2;ds2+=2){
            BOOST_CHECK_SMALL(std::abs(s1342s.v(ds1,ds2)-p1.lsquare(ds1,3)*c342s.v(ds2)),epsilon);
            BOOST_CHECK_SMALL(std::abs(a1342a.v(ds1,ds2)-p1.langle(ds1,3)*c342a.v(ds2)),epsilon);
          }
      else if(o==1)
        for(int ds2=-2;ds2<=2;ds2+=2){
          BOOST_CHECK_SMALL(std::abs(s1342s.v(ds2)-p1.lsquare(3)*c342s.v(ds2)),epsilon);
          BOOST_CHECK_SMALL(std::abs(a1342a.v(ds2)-p1.langle(3)*c342a.v(ds2)),epsilon);
        }
      else if(o==2)
        for(int ds1=-2;ds1<=2;ds1+=2){
          BOOST_CHECK_SMALL(std::abs(s1342s.v(ds1)-p1.lsquare(ds1,3)*c342s.v()),epsilon);
          BOOST_CHECK_SMALL(std::abs(a1342a.v(ds1)-p1.langle(ds1,3)*c342a.v()),epsilon);
        }
      else if(o==3){
        BOOST_CHECK_SMALL(std::abs(s1342s.v()-p1.lsquare(3)*c342s.v()),epsilon);
        BOOST_CHECK_SMALL(std::abs(a1342a.v()-p1.langle(3)*c342a.v()),epsilon);
      }
    }
  }
}
  
BOOST_AUTO_TEST_SUITE_END()


