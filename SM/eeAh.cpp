
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2023 Neil Christensen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

//File:  SPINAS/SM/eeAh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eeAh.h"

namespace spinas {

  eeAh::eeAh(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& vev):
    e(echarge), me(masse), mh(massh), v(vev), prop(masse,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(me);
    p2=particle(me);
    p3=particle(0);
    p4=particle(mh);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s342a = sproduct(SQUARE,&p3,&p4,&p2);
    a342s = sproduct(ANGLE,&p3,&p4,&p2);
    s341a = sproduct(SQUARE,&p3,&p4,&p1);
    a341s = sproduct(ANGLE,&p3,&p4,&p1);
    s3243s = sproduct(SQUARE,&p3,&p2,&p4,&p3);
    a3243a = sproduct(ANGLE,&p3,&p2,&p4,&p3);
  }
  void eeAh::set_masses(const ldouble& masse, const ldouble& massh){
    me=masse;
    mh=massh;
    p1.set_mass(me);
    p2.set_mass(me);
    p4.set_mass(mh);
    prop.set_mass(me);
  }
  void eeAh::set_v(const ldouble& vev){
    v = vev;
  }
  void eeAh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s13s.update();
    a13a.update();
    s23s.update();
    a23a.update();
    s12s.update();
    a12a.update();
    s342a.update();
    a342s.update();
    s341a.update();
    a341s.update();
    s3243s.update();
    a3243a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT=prop.den(propTP);
    pDenU=prop.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eeAh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    /*- mh^2[13][23] - me[13][342> - me[23][341> - <12>[3243]
      - mh^2<13><23> - me<13><342] - me<23><341] - [12]<3243>
      Becomes after a sign change due to p3 and p4 being outgoing:*/
    if(ds3>0){
      //- mh^2[13][23] + me[13][342> + me[23][341> + <12>[3243]
      return sqrt2*e*me/v*(- mh*mh*s13s.v(ds1)*s23s.v(ds2) + me*s13s.v(ds1)*s342a.v(ds2) + me*s23s.v(ds2)*s341a.v(ds1) + a12a.v(ds1,ds2)*s3243s.v())/pDenT/pDenU;
    }
    else if(ds3<0){
      //- mh^2<13><23> + me<13><342] + me<23><341] + [12]<3243>
      return sqrt2*e*me/v*(- mh*mh*a13a.v(ds1)*a23a.v(ds2) + me*a13a.v(ds1)*a342s.v(ds2) + me*a23a.v(ds2)*a341s.v(ds1) + s12s.v(ds1,ds2)*a3243a.v())/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eeAh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=4){
	  M = amp(j1,j2,j3);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2^2=1/4
    return amp2/4.0;
  }
  



  //  Tests
  int test_eeAh(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> A , h       :";
    {//amp^2
      int i=0;
      // me=0.0005, mmu=0.105, pspatial=250
      ldouble me=0.0005, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble vev=2.*MW*SW/EE;
      eeAh eeAhAmp = eeAh(EE,me,mh,vev);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.944407096414361E-11,6.831700609227872E-12,4.333250100735940E-12,3.282765227838440E-12,2.717988414450522E-12,2.377174820163977E-12,2.160452329437888E-12,2.022183380354527E-12,1.939434188319940E-12,1.900548289807435E-12,1.900548289807435E-12,1.939434188319939E-12,2.022183380354527E-12,2.160452329437888E-12,2.377174820163978E-12,2.717988414450523E-12,3.282765227838439E-12,4.333250100735940E-12,6.831700609227876E-12,1.944407096414363E-11};
      i += eeAhAmp.test_2to2_amp2([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH);
      i += eeAhAmp.test_2to2_amp2_rotations([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH);
      i += eeAhAmp.test_2to2_amp2_boosts([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH);
      i += eeAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH);
      //Close to Higgs threshold
      pspatial = 125.1;
      ldouble dataCH2[20] = {3.211432817747364E-11,1.128341260465880E-11,7.156907423747919E-12,5.421899563515941E-12,4.489099638641374E-12,3.926203132245237E-12,3.568258687149836E-12,3.339890131177967E-12,3.203219435276269E-12,3.138994484194547E-12,3.138994484194547E-12,3.203219435276269E-12,3.339890131177967E-12,3.568258687149836E-12,3.926203132245237E-12,4.489099638641374E-12,5.421899563515939E-12,7.156907423747918E-12,1.128341260465880E-11,3.211432817747359E-11};
      i += eeAhAmp.test_2to2_amp2([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH2);
      i += eeAhAmp.test_2to2_amp2_rotations([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH2);
      i += eeAhAmp.test_2to2_amp2_boosts([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH2);
      i += eeAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH2);
      //me~mh and close to threshold
      me = 125.1;
      mh = 125;
      pspatial = 95;
      eeAhAmp.set_masses(me,mh);
      ldouble dataCH4[20] = {1.440032241806676E-01,1.150427820970875E-01,9.540296918656953E-02,8.163935252752440E-02,7.179919005350162E-02,6.471588278824776E-02,5.966681578694494E-02,5.619627748107786E-02,5.401951957811184E-02,5.296974106006085E-02,5.296974106006085E-02,5.401951957811182E-02,5.619627748107783E-02,5.966681578694494E-02,6.471588278824776E-02,7.179919005350162E-02,8.163935252752438E-02,9.540296918656953E-02,1.150427820970875E-01,1.440032241806675E-01};
      i += eeAhAmp.test_2to2_amp2([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH4);
      i += eeAhAmp.test_2to2_amp2_rotations([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH4);
      i += eeAhAmp.test_2to2_amp2_boosts([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH4);
      i += eeAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH4);
      //Invert me and mh and close to threshold
      me = 125;
      mh = 0.0005;
      pspatial = 125.1;
      eeAhAmp.set_masses(me,mh);
      ldouble dataCH3[20] = {1.723455244977663E-01,1.271661161741428E-01,1.004181179164352E-01,8.337808956648010E-02,7.200299581975704E-02,6.421279368516000E-02,5.885648858329141E-02,5.526736497806094E-02,5.305405667123540E-02,5.199695637421633E-02,5.199695637421633E-02,5.305405667123538E-02,5.526736497806092E-02,5.885648858329140E-02,6.421279368515999E-02,7.200299581975703E-02,8.337808956648009E-02,1.004181179164352E-01,1.271661161741428E-01,1.723455244977661E-01};
      i += eeAhAmp.test_2to2_amp2([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH3);
      i += eeAhAmp.test_2to2_amp2_rotations([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH3);
      i += eeAhAmp.test_2to2_amp2_boosts([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH3);
      i += eeAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAhAmp.amp2(); }, me,me,0,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
