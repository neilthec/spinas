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

//[[123]]=-[[132]] & <<123>>=-<<132>>
BOOST_AUTO_TEST_CASE(stp123stp132_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 10000000;// 
  BOOST_TEST_MESSAGE("Testing supleproduct:");
  BOOST_TEST_MESSAGE("\t* [[123]]=-[[132]] & <<123>>=-<<132>>");
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

//[[111]]=0 & <<111>>=0
BOOST_AUTO_TEST_CASE(stp111_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;// 
  BOOST_TEST_MESSAGE("Testing supleproduct:");
  BOOST_TEST_MESSAGE("\t* [[111]]=0 & <<111>>=0");
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

//[[111]]=m1^3epsilon & <<111>>=m1^3epsilon
BOOST_AUTO_TEST_CASE(stp111m3_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;// 
  BOOST_TEST_MESSAGE("Testing supleproduct:");
  BOOST_TEST_MESSAGE("\t* [[111]]=m1^3*epsilon & <<111>>=m1^3*epsilon");
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

//[[ijk]][[klm]]=mk^2[[il]][[jm]]-mk^2[[im]][[jl]]
BOOST_AUTO_TEST_CASE(stp123stp456_tests) {
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 10000000000000;// 
  BOOST_TEST_MESSAGE("Testing supleproduct:");
  BOOST_TEST_MESSAGE("\t* [[ijk]][[klm]]=mk^2[[il]][[jm]]-mk^2[[im]][[jl]]");
  ldouble mi,mj,mk,ml,mm;
  ldouble momi[4], momj[4], momk[4], moml[4], momm[4];
  for(int n=0;n<100;n++){
    mk = choose_random_momentum(momk,-50,50);
    particle pk=particle(momk,mk);
    schain cksU = schain(&pk,SQUARE,3);
    schain cksL = schain(&pk,LOWER,SQUARE,3);
    schain ckaU = schain(&pk,ANGLE,3);
    schain ckaL = schain(&pk,LOWER,ANGLE,3);
    for(int i=0;i<2;i++){
      mi=0;
      if(i==0) choose_random_massless_momentum(momi,-50,50);
      else mi = choose_random_momentum(momi,-50,50);
      particle pi=particle(momi,mi);
      schain cis = schain(&pi,SQUARE,3);
      schain cia = schain(&pi,ANGLE,3);
      for(int j=0;j<2;j++){
        mj=0;
        if(j==0) choose_random_massless_momentum(momj,-50,50);
        else mj = choose_random_momentum(momj,-50,50);
        particle pj=particle(momj,mj);
        schain cjs = schain(&pj,SQUARE,3);
        schain cja = schain(&pj,ANGLE,3);
        for(int l=0;l<2;l++){
          ml=0;
          if(l==0) choose_random_massless_momentum(moml,-50,50);
          else ml = choose_random_momentum(moml,-50,50);
          particle pl=particle(moml,ml);
          schain cls = schain(&pl,SQUARE,3);
          schain cla = schain(&pl,ANGLE,3);
          for(int m=0;m<2;m++){
            mm=0;
            if(m==0) choose_random_massless_momentum(momm,-50,50);
            else mm = choose_random_momentum(momm,-50,50);
            particle pm=particle(momm,mm);
            schain cms = schain(&pm,SQUARE,3);
            schain cma = schain(&pm,ANGLE,3);
            //All the particles are set up.  Now we can test the supleproducts.
            supleproduct sijks = supleproduct(&cis,&cjs,&cksU);
            supleproduct sklms = supleproduct(&cksL,&cls,&cms);
            supleproduct sijka = supleproduct(&cia,&cja,&ckaU);
            supleproduct sklma = supleproduct(&ckaL,&cla,&cma);
            sproduct sils = sproduct(SQUARE,&pi,&pl,3);
            sproduct sjms = sproduct(SQUARE,&pj,&pm,3);
            sproduct sims = sproduct(SQUARE,&pi,&pm,3);
            sproduct sjls = sproduct(SQUARE,&pj,&pl,3);
            sproduct aila = sproduct(ANGLE,&pi,&pl,3);
            sproduct ajma = sproduct(ANGLE,&pj,&pm,3);
            sproduct aima = sproduct(ANGLE,&pi,&pm,3);
            sproduct ajla = sproduct(ANGLE,&pj,&pl,3);
           
            //Now for the tests
            if(i==0&&j==0&&l==0&&m==0){
              cdouble ress = cdouble(0,0), resa = cdouble(0,0);
              for(int dsk=-2;dsk<=2;dsk+=2){
                ress += sijks.v(dsk)*sklms.v(-dsk);
                resa += sijka.v(dsk)*sklma.v(-dsk);
              }
              BOOST_CHECK_SMALL(std::abs(ress-mk*mk*sils.v()*sjms.v()+mk*mk*sims.v()*sjls.v()),epsilon);
              BOOST_CHECK_SMALL(std::abs(resa-mk*mk*aila.v()*ajma.v()+mk*mk*aima.v()*ajla.v()),epsilon);
            }
          }
        }
      }
    }
  }
}

  
BOOST_AUTO_TEST_SUITE_END()