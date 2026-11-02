



//File:  SPINAS/tests/sproduct.cpp

#include "test_utils.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>
#include <algorithm>
#include <stdexcept>

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


SPINAS_TEST_SUITE(sproduct_tests)

//[12]<21>=2p1.p2
SPINAS_TEST_CASE(s12sa21a_tests) {
  SPINAS_TEST_MESSAGE("Testing sproduct:");
  SPINAS_TEST_MESSAGE("\t* [12]<21>=2p1.p2");
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
SPINAS_TEST_CASE(s12ss21s_tests) {
  SPINAS_TEST_MESSAGE("\t* [12][21]=-2m1m2 and <12><21>=-2m1m2");
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
SPINAS_TEST_CASE(s12sa231s_tests) {
  SPINAS_TEST_MESSAGE("\t* -[12]<231] = 2m1p2.p3");
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
SPINAS_TEST_CASE(s12ss2341s_tests) {
  SPINAS_TEST_MESSAGE("\t* [12][2341] = -2m1m2p3.p4");
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
SPINAS_TEST_CASE(s12sa2341a_tests) {
  SPINAS_TEST_MESSAGE("\t* [12]<2341> = 2p1.p2p3.p4-2p1.p3p2.p4+2p1.p4p2.p3");
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
SPINAS_TEST_CASE(s142aa132sC_tests) {
  SPINAS_TEST_MESSAGE("\t* [142><132]* = -[142><231] = -2m1m2p3.p4");
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
SPINAS_TEST_CASE(s142as132aC_tests) {
  SPINAS_TEST_MESSAGE("\t* [142>[132>* = -[142>[231> = 2p1.p4p2.p3-2p1.p2p3.p4+2p1.p3p2.p4");
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
SPINAS_TEST_CASE(s142as1352sC_tests) {
  SPINAS_TEST_MESSAGE("\t* [142>[1352]* = [142><2531> = -2m2(p1.p4p3.p5-p1.p5p3.p4+p1.p3p4.p5)");
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
SPINAS_TEST_CASE(s123a_lu_tests) {
  SPINAS_TEST_MESSAGE("\t* [12_I]<2^I3>=[123>");
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
SPINAS_TEST_CASE(s123a_ul_tests) {
  SPINAS_TEST_MESSAGE("\t* [12^I]<2_I3>=-[123>");
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
SPINAS_TEST_CASE(a123s_ul_tests) {
  SPINAS_TEST_MESSAGE("\t* <12^I>[2_I3]=<123]");
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
SPINAS_TEST_CASE(a123s_lu_tests) {
  SPINAS_TEST_MESSAGE("\t* <12_I>[2^I3]=-<123]");
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
SPINAS_TEST_CASE(m2s13s_lu_tests) {
  SPINAS_TEST_MESSAGE("\t* [12_I][2^I3]=m2[13]");
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
SPINAS_TEST_CASE(m2s13s_ul_tests) {
  SPINAS_TEST_MESSAGE("\t* [12^I][2_I3]=-m2[13]");
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
SPINAS_TEST_CASE(m2a13a_ul_tests) {
  SPINAS_TEST_MESSAGE("\t* <12^I><2_I3>=m2<13>");
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
SPINAS_TEST_CASE(m2a13a_lu_tests) {
  SPINAS_TEST_MESSAGE("\t* <12_I><2^I3>=-m2<13>");
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
SPINAS_TEST_CASE(s142a_m2s12s_m1a12a_s132a_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 10000000000;//
  SPINAS_TEST_MESSAGE("\t* [1|p4|2>=m2[12]-m1<12>-[1|p3|2>");
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
	  SPINAS_CHECK_SMALL(std::abs(s142a.v(ds1,ds2)-m2*s12s.v(ds1,ds2)+m1*a12a.v(ds1,ds2)+s132a.v(ds1,ds2)),epsilon);
    else if(o==2||o==3)
      for(int ds1=-1;ds1<=1;ds1+=2)
	SPINAS_CHECK_SMALL(std::abs(s142a.v(ds1)+m1*a12a.v(ds1)+s132a.v(ds1)),epsilon);
    else if(o==4||o==5)
      for(int ds2=-1;ds2<=1;ds2+=2)
	SPINAS_CHECK_SMALL(std::abs(s142a.v(ds2)-m2*s12s.v(ds2)+s132a.v(ds2)),epsilon);
    else if(o==6||o==7)
      SPINAS_CHECK_SMALL(std::abs(s142a.v()+s132a.v()),epsilon);
  }
}



//[1|p4|2>+[2|p4|1>=-[1|p3|2>-[2|p3|1> (m1=m2)
SPINAS_TEST_CASE(s142a_s241a_s132a_s231a_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000000;//
  SPINAS_TEST_MESSAGE("\t* [1|p4|2>+[2|p4|1>=-[1|p3|2>-[2|p3|1> (m1=m2)");
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
	  SPINAS_CHECK_SMALL(std::abs(s142a.v(ds1,ds2)+s241a.v(ds2,ds1)+s132a.v(ds1,ds2)+s231a.v(ds2,ds1)),epsilon);
    else
      SPINAS_CHECK_SMALL(std::abs(s142a.v()+s241a.v()+s132a.v()+s231a.v()),epsilon);
  }
}





// Direct component contractions check the complex value of each spin state,
// without using sproduct or cmatrix multiplication as the reference.
static void set_chain_particle(particle& value, ldouble mass,
                               ldouble x, ldouble y, ldouble z) {
  const ldouble momentum[4] = {std::sqrt(mass*mass+x*x+y*y+z*z),x,y,z};
  value.set_mass(mass);
  value.set_momentum(momentum);
}

static sproduct make_chain(bool angle, particle* left, particle** middle,
                           int count, particle* right,
                           bool leftUpper = true, bool rightUpper = true) {
  switch(count) {
  case 0: return sproduct(angle,left,leftUpper,right,rightUpper);
  case 1: return sproduct(angle,left,leftUpper,middle[0],right,rightUpper);
  case 2: return sproduct(angle,left,leftUpper,middle[0],middle[1],right,rightUpper);
  case 3: return sproduct(angle,left,leftUpper,middle[0],middle[1],middle[2],right,rightUpper);
  case 4: return sproduct(angle,left,leftUpper,middle[0],middle[1],middle[2],middle[3],right,rightUpper);
  case 5: return sproduct(angle,left,leftUpper,middle[0],middle[1],middle[2],middle[3],middle[4],right,rightUpper);
  case 6: return sproduct(angle,left,leftUpper,middle[0],middle[1],middle[2],middle[3],middle[4],middle[5],right,rightUpper);
  default: throw std::invalid_argument("Unsupported chain length in test.");
  }
}

static cdouble explicit_chain(bool angle, particle& left, particle** middle,
                              int count, particle& right, int spinLeft,
                              int spinRight, bool leftUpper, bool rightUpper) {
  const cvector start = left.mass()!=0 ?
    (angle ? left.langle(spinLeft,leftUpper) : left.lsquare(spinLeft,leftUpper)) :
    (angle ? left.langle() : left.lsquare());
  cdouble row[2] = {start.get(0),start.get(1)};
  for(int insertion=0; insertion<count; ++insertion) {
    const particle& value = *middle[insertion];
    const ldouble energy=value.get_momentum(0), x=value.get_momentum(1),
      y=value.get_momentum(2), z=value.get_momentum(3);
    const ldouble sign = (angle != (insertion%2!=0)) ? 1 : -1;
    const cdouble matrix[2][2] = {{energy+sign*z,cdouble(sign*x,-sign*y)},
                                 {cdouble(sign*x,sign*y),energy-sign*z}};
    cdouble next[2] = {0,0};
    for(int column=0; column<2; ++column)
      for(int index=0; index<2; ++index)
        next[column] += row[index]*matrix[index][column];
    row[0]=next[0]; row[1]=next[1];
  }
  const bool rightAngle = angle != (count%2!=0);
  const cvector finish = right.mass()!=0 ?
    (rightAngle ? right.rangle(spinRight,rightUpper) : right.rsquare(spinRight,rightUpper)) :
    (rightAngle ? right.rangle() : right.rsquare());
  return row[0]*finish.get(0)+row[1]*finish.get(1);
}

static cdouble chain_value(sproduct& value, particle& left, particle& right,
                           int spinLeft, int spinRight) {
  if(left.mass()!=0 && right.mass()!=0) return value.v(spinLeft,spinRight);
  if(left.mass()!=0) return value.v(spinLeft);
  if(right.mass()!=0) return value.v(spinRight);
  return value.v();
}

static void check_chain(sproduct& value, bool angle, particle& left,
                        particle** middle, int count, particle& right,
                        bool leftUpper=true, bool rightUpper=true) {
  for(int spinLeft=-1; spinLeft<=1; spinLeft+=2)
    for(int spinRight=-1; spinRight<=1; spinRight+=2) {
      const cdouble expected=explicit_chain(angle,left,middle,count,right,
        spinLeft,spinRight,leftUpper,rightUpper);
      const ldouble tolerance=4096*std::numeric_limits<ldouble>::epsilon()*
        std::max(ldouble(1),std::abs(expected));
      const cdouble actual=chain_value(value,left,right,spinLeft,spinRight);
      SPINAS_CHECK_SMALL(std::abs(actual-expected),tolerance);
      // Repeat after all computations are cached, rather than only testing a
      // newly constructed product.
      SPINAS_CHECK_EQUAL(chain_value(value,left,right,spinLeft,spinRight),actual);
    }
}

SPINAS_TEST_CASE(all_chain_lengths_and_endpoint_indices) {
  particle middleValues[6];
  particle* middle[6];
  for(int i=0; i<6; ++i) {
    set_chain_particle(middleValues[i],ldouble(1)+ldouble(i)/5,
      ldouble(0.2)+ldouble(i)/10,-ldouble(0.3),ldouble(0.4)-ldouble(i)/20);
    middle[i]=&middleValues[i];
  }
  for(int massMask=0; massMask<4; ++massMask) {
    particle left,right;
    set_chain_particle(left,massMask&1 ? ldouble(1.3) : 0,0.3,-0.7,0.4);
    set_chain_particle(right,massMask&2 ? ldouble(0.8) : 0,-0.6,0.2,0.5);
    for(bool angle : {false,true})
      for(bool leftUpper : {false,true})
        for(bool rightUpper : {false,true}) {
          if((!leftUpper && left.mass()==0) || (!rightUpper && right.mass()==0)) continue;
          for(int count=0; count<=6; ++count) {
            sproduct value=make_chain(angle,&left,middle,count,&right,leftUpper,rightUpper);
            check_chain(value,angle,left,middle,count,right,leftUpper,rightUpper);
          }
        }
  }
}

SPINAS_TEST_CASE(adjacent_equal_momenta_reduce_to_mass_squared) {
  particle left,right,internal;
  set_chain_particle(left,1.3,0.3,-0.7,0.4);
  set_chain_particle(right,0.8,-0.6,0.2,0.5);
  set_chain_particle(internal,2.5,0.4,-0.3,0.6);
  particle* middle[6]={&internal,&internal,&internal,&internal,&internal,&internal};
  for(bool angle : {false,true}) {
    sproduct bracket=make_chain(angle,&left,middle,0,&right);
    ldouble factor=1;
    for(int count=2; count<=6; count+=2) {
      factor *= internal.mass()*internal.mass();
      sproduct value=make_chain(angle,&left,middle,count,&right);
      for(int spinLeft : {-1,1}) for(int spinRight : {-1,1}) {
        const cdouble expected=factor*bracket.v(spinLeft,spinRight);
        SPINAS_CHECK_SMALL(std::abs(value.v(spinLeft,spinRight)-expected),
          4096*std::numeric_limits<ldouble>::epsilon()*std::max(ldouble(1),std::abs(expected)));
      }
    }
  }
}

SPINAS_TEST_CASE(update_rebuilds_cached_spins_momenta_and_mass_categories) {
  particle left,right,internal[3];
  set_chain_particle(left,1.3,0.3,-0.7,0.4);
  set_chain_particle(right,0.8,-0.6,0.2,0.5);
  particle* middle[3]={&internal[0],&internal[1],&internal[2]};
  for(int i=0; i<3; ++i) set_chain_particle(internal[i],1+i,0.2,0.3+i/10.,-0.4);
  sproduct angle=make_chain(ANGLE,&left,middle,3,&right);
  sproduct square=make_chain(SQUARE,&left,middle,3,&right);
  for(int stage=0; stage<7; ++stage) {
    if(stage==1) set_chain_particle(left,1.3,-0.2,0.1,0.9);
    if(stage==2) set_chain_particle(internal[1],0.7,-0.5,0.6,0.2);
    if(stage==3) set_chain_particle(left,0,0.3,-0.7,0.4);
    if(stage==4) set_chain_particle(right,0,-0.6,0.2,0.5);
    if(stage==5) set_chain_particle(left,2,0.1,0.4,0.7);
    if(stage==6) set_chain_particle(right,1.5,-0.3,0.7,0.2);
    angle.update(); square.update();
    check_chain(angle,ANGLE,left,middle,3,right);
    check_chain(square,SQUARE,left,middle,3,right);
  }
}

SPINAS_TEST_CASE(incompatible_spin_arity_and_massless_lower_indices_throw) {
  particle massive1,massive2,massless1,massless2;
  set_chain_particle(massive1,1,0.3,-0.7,0.4);
  set_chain_particle(massive2,2,-0.6,0.2,0.5);
  set_chain_particle(massless1,0,0.2,0.5,-0.3);
  set_chain_particle(massless2,0,-0.4,0.1,0.6);
  for(bool angle : {false,true}) {
    sproduct both(angle,&massive1,&massive2), neither(angle,&massless1,&massless2);
    sproduct left(angle,&massive1,&massless1), right(angle,&massless1,&massive1);
    SPINAS_CHECK_THROW(both.v(),std::runtime_error);
    SPINAS_CHECK_THROW(both.v(1),std::runtime_error);
    SPINAS_CHECK_THROW(neither.v(1),std::runtime_error);
    SPINAS_CHECK_THROW(neither.v(1,1),std::runtime_error);
    SPINAS_CHECK_THROW(left.v(),std::runtime_error);
    SPINAS_CHECK_THROW(left.v(1,1),std::runtime_error);
    SPINAS_CHECK_THROW(right.v(),std::runtime_error);
    SPINAS_CHECK_THROW(right.v(1,1),std::runtime_error);
    SPINAS_CHECK_THROW(sproduct(angle,&massless1,LOWER,&massive1),std::runtime_error);
    SPINAS_CHECK_THROW(sproduct(angle,&massive1,&massless1,LOWER),std::runtime_error);
    sproduct lowered(angle,&massive1,LOWER,&massive2);
    (void)lowered.v(1,-1);
    set_chain_particle(massive1,0,0.3,-0.7,0.4);
    SPINAS_CHECK_THROW(lowered.update(),std::runtime_error);
    set_chain_particle(massive1,1,0.3,-0.7,0.4);
  }
}

SPINAS_TEST_SUITE_END()


void test_sproduct_spinsum(sproduct* sp1, sproduct* sp2, sproduct* product, const ldouble& m1, const ldouble& m3, const ldouble& factor){
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000000;//
  if(m1==0 && m3==0){
    cdouble spinorProd(0,0);
    for(int ds2=-1;ds2<=1;ds2+=2)
      spinorProd += sp1->v(ds2)*sp2->v(ds2);
    SPINAS_CHECK_SMALL(std::abs( spinorProd - factor*product->v() ), epsilon);
  }
  else if(m1==0)
    for(int ds3=-1;ds3<=1;ds3+=2){
      cdouble spinorProd(0,0);
      for(int ds2=-1;ds2<=1;ds2+=2)
	spinorProd += sp1->v(ds2)*sp2->v(ds2,ds3);	
      SPINAS_CHECK_SMALL(std::abs( spinorProd - factor*product->v(ds3) ), epsilon);
    }
  else if(m3==0)
    for(int ds1=-1;ds1<=1;ds1+=2){
      cdouble spinorProd(0,0);
      for(int ds2=-1;ds2<=1;ds2+=2)
	spinorProd += sp1->v(ds1,ds2)*sp2->v(ds2);
      SPINAS_CHECK_SMALL(std::abs( spinorProd - factor*product->v(ds1) ), epsilon);
    }
  else
    for(int ds1=-1;ds1<=1;ds1+=2)
      for(int ds3=-1;ds3<=1;ds3+=2){
	cdouble spinorProd(0,0);
	for(int ds2=-1;ds2<=1;ds2+=2)
	  spinorProd += sp1->v(ds1,ds2)*sp2->v(ds2,ds3);	
	SPINAS_CHECK_SMALL(std::abs( spinorProd - factor*product->v(ds1,ds3) ), epsilon);
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
    
    // Check the spin sum before applying the next Lorentz transformation.
    SPINAS_CHECK_SMALL(std::abs(spinProd-expected), epsilon);
    
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


