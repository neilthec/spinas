



//File:  SPINAS/boost/sproduct.cpp

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
#include "particle.h"
#include "sproduct.h"


using namespace spinas;

void test_sproduct(sproduct* sp1, sproduct* sp2, const int& ni, const int& np, particle* p1, particle* p2, particle* p3, particle* p4, particle* p5, const ldouble& expected, const char* spstring, const char* resstring, const ldouble& epsFactor);

void test_sproduct_spinsum(sproduct* sp1, sproduct* sp2, sproduct* product, const ldouble& m1, const ldouble& m3, const ldouble& factor);


BOOST_AUTO_TEST_SUITE(sproduct_tests)

//[12]<21>=2p1.p2
BOOST_AUTO_TEST_CASE(s12sa21a_tests) {
  BOOST_TEST_MESSAGE("Testing sproduct:");
  BOOST_TEST_MESSAGE("\t* [12]<21>=2p1.p2");
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
    ni=2;
    if(o==1||o==2) ni=1;
    if(o==3) ni=0;    
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    sproduct s12s = sproduct(SQUARE,&p1,&p2);	
    ldouble expected = 2.*p1.dot(p2);//[12]
    test_sproduct(&s12s, &s12s, ni, 2, &p1, &p2, &p2, &p2, &p2, expected, "[12][12]*", "-2p1.p2",10);//
  }
}

//[12][21]=-2m1m2
//Not zero when m1=m2=0.  Don't test
BOOST_AUTO_TEST_CASE(s12ss21s_tests) {
  BOOST_TEST_MESSAGE("\t* [12][21]=-2m1m2 and <12><21>=-2m1m2");
  ldouble m1,m2;
  ldouble mom1[4], mom2[4];
  int ni;
  for(int i=0;i<100;i++)
  for(int o=0;o<3;o++){
    m1=0;
    m2=0;
    if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom2,-50,50);
    else m2 = choose_random_momentum(mom2,-50,50);
    ni=2;
    if(o==1||o==2) ni=1;
    if(o==3) ni=0;    
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    sproduct s12s = sproduct(SQUARE,&p1,&p2);
    sproduct a12a = sproduct(ANGLE,&p1,&p2);
    ldouble expected = -2.*m1*m2;
test_sproduct(&s12s, &a12a, ni, 2, &p1, &p2, &p2, &p2, &p2, expected, "[12]<12>*", "-2m1m2",1000);//
test_sproduct(&a12a, &s12s, ni, 2, &p1, &p2, &p2, &p2, &p2, expected, "<12>[12]*", "-2m1m2",1000);//
  }
}

//-[12]<231] = 2m1p2.p3
//If m1=0, this formula does not apply.  So, we must skip it.
BOOST_AUTO_TEST_CASE(s12sa231s_tests) {
  BOOST_TEST_MESSAGE("\t* -[12]<231] = 2m1p2.p3");
  ldouble m1,m2,m3;
  ldouble mom1[4], mom2[4], mom3[4];
  int ni;
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m2=0;m3=0;
    m1 = choose_random_momentum(mom1,-50,50);
    if(o==1||o==3) choose_random_massless_momentum(mom2,-50,50);
    else m2 = choose_random_momentum(mom2,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50);
    ni=2;
    if(o==1||o==3) ni=1;
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    sproduct s12s = sproduct(SQUARE,&p1,&p2);//[12]	
    sproduct a132s = sproduct(ANGLE,&p1,&p3,&p2);//<132]* = -<231]
    ldouble expected = 2.*m1*p2.dot(p3);
test_sproduct(&s12s, &a132s, ni, 3, &p1, &p2, &p3, &p3, &p3, expected, "[12]<132]*", "2m1p2.p3",100000);//
  }
}

//[12][2341] = -2m1m2p3.p4
//If m1=0 or m2=0, this formula does not apply.  So, we must skip it.
BOOST_AUTO_TEST_CASE(s12ss2341s_tests) {
  BOOST_TEST_MESSAGE("\t* [12][2341] = -2m1m2p3.p4");
  ldouble m1,m2,m3,m4;
  ldouble mom1[4], mom2[4], mom3[4], mom4[4];
  int ni=2;
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m3=0;m4=0;
    m1=choose_random_momentum(mom1,-50,50);
    m2=choose_random_momentum(mom2,-50,50);
    if(o==1||o==3) choose_random_massless_momentum(mom3,-50,50);
    else m3=choose_random_momentum(mom3,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom4,-50,50);
    else m4=choose_random_momentum(mom4,-50,50);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    particle p4=particle(mom4,m4);
    sproduct s12s = sproduct(SQUARE,&p1,&p2);//[12]	
    sproduct a1432a = sproduct(ANGLE,&p1,&p4,&p3,&p2);//<1432>* = [2341]	
    ldouble expected = -2.*m1*m2*p3.dot(p4);
    test_sproduct(&s12s, &a1432a, ni, 4, &p1, &p2, &p3, &p4, &p4, expected, "[12]<1432>*", "-2m1m2p3.p4",500000);//
  }
}

//[12]<2341> = 2p1.p2p3.p4-2p1.p3p2.p4+2p1.p4p2.p3
//We are only checking the real part since we don't have an epsilon tensor to check the imaginary part.
BOOST_AUTO_TEST_CASE(s12sa2341a_tests) {
  BOOST_TEST_MESSAGE("\t* [12]<2341> = 2p1.p2p3.p4-2p1.p3p2.p4+2p1.p4p2.p3");
  ldouble m1,m2,m3,m4;
  ldouble mom1[4], mom2[4], mom3[4], mom4[4];
  int ni;
  for(int i=0;i<100;i++)
  for(int o=0;o<16;o++){
    m1=0;m2=0;m3=0;m4=0;
    if(o==1||o==5||o==6||o==7||o==11||o==12||o==13||o==15) choose_random_massless_momentum(mom1,-50,50);
    else m1=choose_random_momentum(mom1,-50,50);
    if(o==2||o==5||o==8||o==9||o==11||o==12||o==14||o==15) choose_random_massless_momentum(mom2,-50,50);
    else m2=choose_random_momentum(mom2,-50,50);
    if(o==3||o==6||o==8||o==10||o==11||o==13||o==14||o==15) choose_random_massless_momentum(mom3,-50,50);
    else m3=choose_random_momentum(mom3,-50,50);
    if(o==4||o==7||o==9||o==10||o==12||o==13||o==14||o==15) choose_random_massless_momentum(mom4,-50,50);
    else m4=choose_random_momentum(mom4,-50,50);
    ni=2;
    if(o==1||o==2||o==6||o==7||o==8||o==9||o==13||o==14) ni=1;
    if(o==5||o==11||o==12||o==15) ni=0;
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    particle p4=particle(mom4,m4);
    sproduct s12s = sproduct(SQUARE,&p1,&p2);//[12]
    sproduct s1432s = sproduct(SQUARE,&p1,&p4,&p3,&p2);//[1432]
    ldouble expected = 2.*p1.dot(p2)*p3.dot(p4)-2.*p1.dot(p3)*p2.dot(p4)+2.*p1.dot(p4)*p2.dot(p3);
    test_sproduct(&s12s, &s1432s, ni, 4, &p1, &p2, &p3, &p4, &p4, expected, "[12][1432]*", "2p1.p2p3.p4-2p1.p3p2.p4+2p1.p4p2.p3",100000);//
  }
}

//[142><132]* = -[142><231] = -2m1m2p3.p4
BOOST_AUTO_TEST_CASE(s142aa132sC_tests) {
  BOOST_TEST_MESSAGE("\t* [142><132]* = -[142><231] = -2m1m2p3.p4");
  ldouble m1,m2,m3,m4;
  ldouble mom1[4], mom2[4], mom3[4], mom4[4];
  int ni=2;
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m3=0;m4=0;
    m1=choose_random_momentum(mom1,-50,50);
    m2=choose_random_momentum(mom2,-50,50);
    if(o==1||o==3) choose_random_massless_momentum(mom3,-50,50);
    else m3=choose_random_momentum(mom3,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom4,-50,50);
    else m4=choose_random_momentum(mom4,-50,50);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    particle p4=particle(mom4,m4);
    sproduct sp1 = sproduct(SQUARE,&p1,&p4,&p2);//[142>
    sproduct sp2 = sproduct(ANGLE,&p1,&p3,&p2);//<132]
    ldouble expected = -2.*m1*m2*p3.dot(p4);
    test_sproduct(&sp1, &sp2, ni, 4, &p1, &p2, &p3, &p4, &p4, expected, "[142><132]*", "-2m1m2p3.p4",1000000);//
  }
}

//[142>[132>* = -[142>[231> = 2p1.p4p2.p3-2p1.p2p3.p4+2p1.p3p2.p4
BOOST_AUTO_TEST_CASE(s142as132aC_tests) {
  BOOST_TEST_MESSAGE("\t* [142>[132>* = -[142>[231> = 2p1.p4p2.p3-2p1.p2p3.p4+2p1.p3p2.p4");
  ldouble m1,m2,m3,m4;
  ldouble mom1[4], mom2[4], mom3[4], mom4[4];
  int ni;
  for(int i=0;i<100;i++)
  for(int o=0;o<16;o++){
    m1=0;m2=0;m3=0;m4=0;
    if(o==1||o==5||o==6||o==7||o==11||o==12||o==13||o==15) choose_random_massless_momentum(mom1,-50,50);
    else m1=choose_random_momentum(mom1,-50,50);
    if(o==2||o==5||o==8||o==9||o==11||o==12||o==14||o==15) choose_random_massless_momentum(mom2,-50,50);
    else m2=choose_random_momentum(mom2,-50,50);
    if(o==3||o==6||o==8||o==10||o==11||o==13||o==14||o==15) choose_random_massless_momentum(mom3,-50,50);
    else m3=choose_random_momentum(mom3,-50,50);
    if(o==4||o==7||o==9||o==10||o==12||o==13||o==14||o==15) choose_random_massless_momentum(mom4,-50,50);
    else m4=choose_random_momentum(mom4,-50,50);
    ni=2;
    if(o==1||o==2||o==6||o==7||o==8||o==9||o==13||o==14) ni=1;
    if(o==5||o==11||o==12||o==15) ni=0;
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    particle p4=particle(mom4,m4);
    sproduct sp1 = sproduct(SQUARE,&p1,&p4,&p2);//[142>
    sproduct sp2 = sproduct(SQUARE,&p1,&p3,&p2);//[132>
    ldouble expected = 2.*p1.dot(p4)*p2.dot(p3)-2.*p1.dot(p2)*p3.dot(p4)+2.*p1.dot(p3)*p2.dot(p4);
    test_sproduct(&sp1, &sp2, ni, 4, &p1, &p2, &p3, &p4, &p4, expected, "[142>[132>*", "2p1.p4p2.p3-2p1.p2p3.p4+2p1.p3p2.p4",100000);//
  }
}

//[142>[1352]* = [142><2531> = -2m2(p1.p4p3.p5-p1.p5p3.p4+p1.p3p4.p5)
BOOST_AUTO_TEST_CASE(s142as1352sC_tests) {
  BOOST_TEST_MESSAGE("\t* [142>[1352]* = [142><2531> = -2m2(p1.p4p3.p5-p1.p5p3.p4+p1.p3p4.p5)");
  ldouble m1,m2,m3,m4,m5;
  ldouble mom1[4], mom2[4], mom3[4], mom4[4], mom5[4];
  int ni;
  for(int i=0;i<100;i++)
  for(int o=0;o<14;o++){
    m1=0;m3=0;m4=0;m5=0;
    if(o==1||o==5||o==6||o==7||o==10||o==11||o==12||o==13) choose_random_massless_momentum(mom1,-50,50);
    else m1=choose_random_momentum(mom1,-50,50);
    m2=choose_random_momentum(mom2,-50,50);
    if(o==2||o==5||o==8||o==9||o==10||o==11||o==13) choose_random_massless_momentum(mom3,-50,50);
    else m3=choose_random_momentum(mom3,-50,50);
    if(o==3||o==6||o==8||o==10||o==12||o==13) choose_random_massless_momentum(mom4,-50,50);
    else m4=choose_random_momentum(mom4,-50,50);
    if(o==4||o==7||o==9||o==11||o==12||o==13) choose_random_massless_momentum(mom5,-50,50);
    else m5=choose_random_momentum(mom5,-50,50);
    ni=2;
    if(o==1||o==5||o==6||o==7||o==10||o==11||o==12||o==13) ni=1;
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    particle p4=particle(mom4,m4);
    particle p5=particle(mom5,m5);
    sproduct sp1 = sproduct(SQUARE,&p1,&p4,&p2);//[142>
    sproduct sp2 = sproduct(SQUARE,&p1,&p3,&p5,&p2);//[1352]
    ldouble expected = -2.*m2*(p1.dot(p4)*p3.dot(p5)-p1.dot(p5)*p3.dot(p4)+p1.dot(p3)*p4.dot(p5));
    test_sproduct(&sp1, &sp2, ni, 5, &p1, &p2, &p3, &p4, &p5, expected, "[142>[1352]*", "-2m2(p1.p4p3.p5-p1.p5p3.p4+p1.p3p4.p5)",200000000);//
  }
}



//[12_I]<2^I3>=[123>
BOOST_AUTO_TEST_CASE(s123a_lu_tests) {
  BOOST_TEST_MESSAGE("\t* [12_I]<2^I3>=[123>");
  ldouble m1, m2, m3;
  ldouble mom1[4], mom2[4], mom3[4];
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m1=0;
    m3=0;
    if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    m2 = choose_random_momentum(mom2,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    sproduct s12s = sproduct(SQUARE,&p1,&p2,LOWER);
    sproduct a23a = sproduct(ANGLE,&p2,&p3);
    sproduct s12sa23a = sproduct(SQUARE,&p1,&p2,&p3);//[123>
    test_sproduct_spinsum(&s12s, &a23a, &s12sa23a, m1, m3, 1);
  }
}


//[12^I]<2_I3>=-[123>
BOOST_AUTO_TEST_CASE(s123a_ul_tests) {
  BOOST_TEST_MESSAGE("\t* [12^I]<2_I3>=-[123>");
  ldouble m1, m2, m3;
  ldouble mom1[4], mom2[4], mom3[4];
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m1=0;
    m3=0;
    if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    m2 = choose_random_momentum(mom2,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    sproduct s12s = sproduct(SQUARE,&p1,&p2);
    sproduct a23a = sproduct(ANGLE,&p2, LOWER, &p3);
    sproduct s12sa23a = sproduct(SQUARE,&p1,&p2,&p3);//[123>
    test_sproduct_spinsum(&s12s, &a23a, &s12sa23a, m1, m3, -1);
  }
}



//<12^I>[2_I3]=<123]
BOOST_AUTO_TEST_CASE(a123s_ul_tests) {
  BOOST_TEST_MESSAGE("\t* <12^I>[2_I3]=<123]");
  ldouble m1, m2, m3;
  ldouble mom1[4], mom2[4], mom3[4];
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m1=0;
    m3=0;
    if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    m2 = choose_random_momentum(mom2,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    sproduct a12a = sproduct(ANGLE,&p1,&p2);
    sproduct s23s = sproduct(SQUARE,&p2, LOWER, &p3);
    sproduct a12as23s = sproduct(ANGLE,&p1,&p2,&p3);//<123]
    test_sproduct_spinsum(&a12a, &s23s, &a12as23s, m1, m3, 1);
  }
}



//<12_I>[2^I3]=-<123]
BOOST_AUTO_TEST_CASE(a123s_lu_tests) {
  BOOST_TEST_MESSAGE("\t* <12_I>[2^I3]=-<123]");
  ldouble m1, m2, m3;
  ldouble mom1[4], mom2[4], mom3[4];
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m1=0;
    m3=0;
    if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    m2 = choose_random_momentum(mom2,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    sproduct a12a = sproduct(ANGLE,&p1,&p2, LOWER);
    sproduct s23s = sproduct(SQUARE,&p2, &p3);
    sproduct a12as23s = sproduct(ANGLE,&p1,&p2,&p3);//<123]
    test_sproduct_spinsum(&a12a, &s23s, &a12as23s, m1, m3, -1);
  }
}


//[12_I][2^I3]=m2[13]
BOOST_AUTO_TEST_CASE(m2s13s_lu_tests) {
  BOOST_TEST_MESSAGE("\t* [12_I][2^I3]=m2[13]");
  ldouble m1, m2, m3;
  ldouble mom1[4], mom2[4], mom3[4];
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m1=0;
    m3=0;
    if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    m2 = choose_random_momentum(mom2,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    sproduct s12s = sproduct(SQUARE,&p1,&p2,LOWER);
    sproduct s23s = sproduct(SQUARE,&p2,&p3);
    sproduct s13s = sproduct(SQUARE,&p1,&p3);//[13]
    test_sproduct_spinsum(&s12s, &s23s, &s13s, m1, m3, m2);
  }
}

//[12^I][2_I3]=-m2[13]
BOOST_AUTO_TEST_CASE(m2s13s_ul_tests) {
  BOOST_TEST_MESSAGE("\t* [12^I][2_I3]=-m2[13]");
  ldouble m1, m2, m3;
  ldouble mom1[4], mom2[4], mom3[4];
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m1=0;
    m3=0;
    if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    m2 = choose_random_momentum(mom2,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    sproduct s12s = sproduct(SQUARE,&p1,&p2);
    sproduct s23s = sproduct(SQUARE,&p2, LOWER,&p3);
    sproduct s13s = sproduct(SQUARE,&p1,&p3);//[13]
    test_sproduct_spinsum(&s12s, &s23s, &s13s, m1, m3, -m2);
  }
}

//<12^I><2_I3>=m2<13>
BOOST_AUTO_TEST_CASE(m2a13a_ul_tests) {
  BOOST_TEST_MESSAGE("\t* <12^I><2_I3>=m2<13>");
  ldouble m1, m2, m3;
  ldouble mom1[4], mom2[4], mom3[4];
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m1=0;
    m3=0;
    if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    m2 = choose_random_momentum(mom2,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    sproduct a12a = sproduct(ANGLE,&p1,&p2);
    sproduct a23a = sproduct(ANGLE,&p2, LOWER,&p3);
    sproduct a13a = sproduct(ANGLE,&p1,&p3);//<13>
    test_sproduct_spinsum(&a12a, &a23a, &a13a, m1, m3, m2);
  }
}

//<12_I><2^I3>=-m2<13>
BOOST_AUTO_TEST_CASE(m2a13a_lu_tests) {
  BOOST_TEST_MESSAGE("\t* <12_I><2^I3>=-m2<13>");
  ldouble m1, m2, m3;
  ldouble mom1[4], mom2[4], mom3[4];
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m1=0;
    m3=0;
    if(o==1||o==3) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    m2 = choose_random_momentum(mom2,-50,50);
    if(o==2||o==3) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    sproduct a12a = sproduct(ANGLE,&p1,&p2, LOWER);
    sproduct a23a = sproduct(ANGLE,&p2,&p3);
    sproduct a13a = sproduct(ANGLE,&p1,&p3);//<13>
    test_sproduct_spinsum(&a12a, &a23a, &a13a, m1, m3, -m2);
  }
}



//[1|p4|2>=m2[12]-m1<12>-[1|p3|2>
BOOST_AUTO_TEST_CASE(s142a_m2s12s_m1a12a_s132a_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 10000000000;//
  BOOST_TEST_MESSAGE("\t* [1|p4|2>=m2[12]-m1<12>-[1|p3|2>");
  ldouble m1, m2, m3, m4;
  ldouble mom1[4], mom2[4], mom3[4], mom4[4];
  for(int i=0;i<100;i++)
  for(int o=0;o<8;o++){
    m1 = 0;
    if(o>3) choose_random_massless_momentum(mom1,-50,50);
    else m1 = choose_random_momentum(mom1,-50,50);
    m2 = 0;
    if(o==2||o==3||o==6||o==7) choose_random_massless_momentum(mom2,-50,50);
    else m2 = choose_random_momentum(mom2,-50,50);
    m3 = 0;
    if(o==1||o==3||o==5||o==7) choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50);
    for(int j=0;j<4;j++)
      mom4[j] = - mom1[j] - mom2[j] - mom3[j];
    m4 = mom4[0]*mom4[0];
    for(int j=1;j<4;j++)
      m4 -= mom4[j]*mom4[j];
    m4 = std::sqrt(m4);
    particle p1=particle(mom1,m1);
    particle p2=particle(mom2,m2);
    particle p3=particle(mom3,m3);
    particle p4=particle(mom4,m4);
    sproduct s142a = sproduct(SQUARE,&p1,&p4,&p2);
    sproduct s132a = sproduct(SQUARE,&p1,&p3,&p2);
    sproduct s12s = sproduct(SQUARE,&p1,&p2);
    sproduct a12a = sproduct(ANGLE,&p1,&p2);
    if(o==0||o==1)
      for(int ds1=-1;ds1<=1;ds1+=2)
	for(int ds2=-1;ds2<=1;ds2+=2)
	  BOOST_CHECK_SMALL(std::abs(s142a.v(ds1,ds2)-m2*s12s.v(ds1,ds2)+m1*a12a.v(ds1,ds2)+s132a.v(ds1,ds2)),epsilon);
    else if(o==2||o==3)
      for(int ds1=-1;ds1<=1;ds1+=2)
	BOOST_CHECK_SMALL(std::abs(s142a.v(ds1)+m1*a12a.v(ds1)+s132a.v(ds1)),epsilon);
    else if(o==4||o==5)
      for(int ds2=-1;ds2<=1;ds2+=2)
	BOOST_CHECK_SMALL(std::abs(s142a.v(ds2)-m2*s12s.v(ds2)+s132a.v(ds2)),epsilon);
    else if(o==6||o==7)
      BOOST_CHECK_SMALL(std::abs(s142a.v()+s132a.v()),epsilon);
  }
}



//[1|p4|2>+[2|p4|1>=-[1|p3|2>-[2|p3|1> (m1=m2)
BOOST_AUTO_TEST_CASE(s142a_s241a_s132a_s231a_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000000;//
  BOOST_TEST_MESSAGE("\t* [1|p4|2>+[2|p4|1>=-[1|p3|2>-[2|p3|1> (m1=m2)");
  ldouble m12, m3, m4;
  ldouble mom1[4], mom2[4], mom3[4], mom4[4];
  for(int i=0;i<100;i++)
  for(int o=0;o<4;o++){
    m12 = 0;
    if(o==2||o==3) {
      choose_random_massless_momentum(mom1,-50,50);
      choose_random_massless_momentum(mom2,-50,50);
    }
    else{
      m12 = choose_random_momentum(mom1,-50,50);
      choose_random_momentum(mom2,-50,50);
      mom2[0] = std::sqrt(mom2[1]*mom2[1]+mom2[2]*mom2[2]+mom2[3]*mom2[3]+m12*m12);
    }
    m3 = 0;
    if(o==1||o==3)  choose_random_massless_momentum(mom3,-50,50);
    else m3 = choose_random_momentum(mom3,-50,50);
    for(int j=0;j<4;j++)
      mom4[j] = - mom1[j] - mom2[j] - mom3[j];
    m4 = mom4[0]*mom4[0];
    for(int j=1;j<4;j++)
      m4 -= mom4[j]*mom4[j];
    m4 = std::sqrt(m4);
    particle p1=particle(mom1,m12);
    particle p2=particle(mom2,m12);
    particle p3=particle(mom3,m3);
    particle p4=particle(mom4,m4);
    sproduct s142a = sproduct(SQUARE,&p1,&p4,&p2);
    sproduct s241a = sproduct(SQUARE,&p2,&p4,&p1);
    sproduct s132a = sproduct(SQUARE,&p1,&p3,&p2);
    sproduct s231a = sproduct(SQUARE,&p2,&p3,&p1);
    sproduct s12s = sproduct(SQUARE,&p1,&p2);
    sproduct a12a = sproduct(ANGLE,&p1,&p2);
    if(o<2)
      for(int ds1=-1;ds1<=1;ds1+=2)
	for(int ds2=-1;ds2<=1;ds2+=2)
	  BOOST_CHECK_SMALL(std::abs(s142a.v(ds1,ds2)+s241a.v(ds2,ds1)+s132a.v(ds1,ds2)+s231a.v(ds2,ds1)),epsilon);
    else
      BOOST_CHECK_SMALL(std::abs(s142a.v()+s241a.v()+s132a.v()+s231a.v()),epsilon);
  }
}





BOOST_AUTO_TEST_SUITE_END()


void test_sproduct_spinsum(sproduct* sp1, sproduct* sp2, sproduct* product, const ldouble& m1, const ldouble& m3, const ldouble& factor){
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000000;//
  if(m1==0 && m3==0){
    cdouble spinorProd(0,0);
    for(int ds2=-1;ds2<=1;ds2+=2)
      spinorProd += sp1->v(ds2)*sp2->v(ds2);
    BOOST_CHECK_SMALL(std::abs( spinorProd - factor*product->v() ), epsilon);
  }
  else if(m1==0)
    for(int ds3=-1;ds3<=1;ds3+=2){
      cdouble spinorProd(0,0);
      for(int ds2=-1;ds2<=1;ds2+=2)
	spinorProd += sp1->v(ds2)*sp2->v(ds2,ds3);	
      BOOST_CHECK_SMALL(std::abs( spinorProd - factor*product->v(ds3) ), epsilon);
    }
  else if(m3==0)
    for(int ds1=-1;ds1<=1;ds1+=2){
      cdouble spinorProd(0,0);
      for(int ds2=-1;ds2<=1;ds2+=2)
	spinorProd += sp1->v(ds1,ds2)*sp2->v(ds2);
      BOOST_CHECK_SMALL(std::abs( spinorProd - factor*product->v(ds1) ), epsilon);
    }
  else
    for(int ds1=-1;ds1<=1;ds1+=2)
      for(int ds3=-1;ds3<=1;ds3+=2){
	cdouble spinorProd(0,0);
	for(int ds2=-1;ds2<=1;ds2+=2)
	  spinorProd += sp1->v(ds1,ds2)*sp2->v(ds2,ds3);	
	BOOST_CHECK_SMALL(std::abs( spinorProd - factor*product->v(ds1,ds3) ), epsilon);
      }
}


//ni is the number of spin indices, np is the number of particles.
void test_sproduct(sproduct* sp1, sproduct* sp2, const int& ni, const int& np, particle* p1, particle* p2, particle* p3, particle* p4, particle* p5, const ldouble& expected, const char* spstring, const char* resstring, const ldouble& epsFactor){
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 100000000*epsFactor;
  int i=0;
  ldouble mom1[4], mom2[4], mom3[4], mom4[4], mom5[4];
  ldouble mom1o[4], mom2o[4], mom3o[4], mom4o[4], mom5o[4];
  for(int j=0;j<4;j++){
    mom1[j] = p1->get_momentum(j);
    mom1o[j] = mom1[j];
    mom2[j] = p2->get_momentum(j);
    mom2o[j] = mom2[j];
    mom3[j] = p3->get_momentum(j);
    mom3o[j] = mom3[j];
    mom4[j] = p4->get_momentum(j);
    mom4o[j] = mom4[j];
    mom5[j] = p5->get_momentum(j);
    mom5o[j] = mom5[j];
  }
  sp1->update();
  sp2->update();
  for(int l=0;l<10;l++){
    ldouble spinProd = 0;
    if(ni==2)
      for(int j=0;j<2;j++)
	for(int k=0;k<2;k++)
	  spinProd += std::real(sp1->v(2*j-1,2*k-1)*std::conj(sp2->v(2*j-1,2*k-1)));
    else if(ni==1)
      for(int k=0;k<2;k++)
	spinProd += std::real(sp1->v(2*k-1)*std::conj(sp2->v(2*k-1)));
    else
      spinProd += std::real(sp1->v()*std::conj(sp2->v()));
    
    //Do the Boost Check
    BOOST_CHECK_SMALL(std::abs(spinProd-expected), epsilon);
    
    //Reset Momentum before Rotations and Boosts
    for(int k=0;k<4;k++){
      mom1[k] = mom1o[k];
      mom2[k] = mom2o[k];
      mom3[k] = mom3o[k];
      mom4[k] = mom4o[k];
      mom5[k] = mom5o[k];
    }
    //Random rotation
    ldouble u[3], um, angle, v[3], vm;
    for(int k=0;k<3;k++) u[k] = choose_random_ldouble(0,1000);
    um = std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
    for(int k=0;k<3;k++) u[k] /= um;
    angle = choose_random_ldouble(0,2.*3.141592653589793238462643383279502884);
    rotate_momentum(mom1,u,angle);
    rotate_momentum(mom2,u,angle);
    rotate_momentum(mom3,u,angle);
    rotate_momentum(mom4,u,angle);
    rotate_momentum(mom5,u,angle);
    //Random Boost
    vm = 2;
    while(vm>=0.99){
      for(int k=0;k<3;k++) v[k] = choose_random_ldouble(0,1);
      vm = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    }
    boost_momentum(mom1,v);
    boost_momentum(mom2,v);
    boost_momentum(mom3,v);
    boost_momentum(mom4,v);
    boost_momentum(mom5,v);
    //Random Rotation
    for(int k=0;k<3;k++) u[k] = choose_random_ldouble(0,1000);
    um = std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
    for(int k=0;k<3;k++) u[k] /= um;
    angle = choose_random_ldouble(0,2.*3.141592653589793238462643383279502884);
    rotate_momentum(mom1,u,angle);
    rotate_momentum(mom2,u,angle);
    rotate_momentum(mom3,u,angle);
    rotate_momentum(mom4,u,angle);
    rotate_momentum(mom5,u,angle);
    p1->set_momentum(mom1);
    p2->set_momentum(mom2);
    p3->set_momentum(mom3);
    p4->set_momentum(mom4);
    p5->set_momentum(mom5);
    sp1->update();
    sp2->update();
  }
}


