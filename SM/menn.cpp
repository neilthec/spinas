
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

//File:  SPINAS/SM/menn.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/menn.h"

namespace spinas {
  //Constructors
  menn::menn(const ldouble& echarge, const ldouble& masse, const ldouble& massmu, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), me(masse), mm(massmu), MW(massW), WW(widthW), SW(sinW), prop(massW,widthW),
    p1(particle(mm)), p2(particle(me)),
    p3(particle(0)), p4(particle(0)),
    //[24], ⟨13⟩, [14], ⟨23⟩
    s24s(sproduct(SQUARE,&p2,&p4)),
    a13a(sproduct(ANGLE,&p1,&p3)),
    s14s(sproduct(SQUARE,&p1,&p4)),
    a23a(sproduct(ANGLE,&p2,&p3))
  {}
  void menn::set_masses(const ldouble& masse, const ldouble& massmu, const ldouble& massW){
    me=masse;
    mm=massmu;
    MW=massW;
    p1.set_mass(mm);
    p2.set_mass(me);
    p3.set_mass(0);
    p4.set_mass(0);
    prop.set_mass(massW);
  }
  void menn::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s24s.update();
    a13a.update();
    s14s.update();
    a23a.update();
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]-mom4[j];
    pDenU = prop.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble menn::amp(const int& ds1, const int& ds2){
    //All ingoing: -(2MW^2 [24]⟨13⟩ + mμme[14]⟨23⟩)/ (2MW^2 (u−MW^2))
    //Sign changes due to p3 and p4 being outgoing.
    // (2MW^2 [24]⟨13⟩ + mμme[14]⟨23⟩)/ (2MW^2 (u−MW^2))
    return e*e/2.0/SW/SW*(2.0*MW*MW*s24s.v(ds2)*a13a.v(ds1) + mm*me*s14s.v(ds1)*a23a.v(ds2))/(MW*MW*pDenU);
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble menn::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2){
	M = amp(j1,j2);
	amp2 += std::pow(std::abs(M),2);
      }
    //Average over initial spins 1/2*1/2 = 1/4
    return amp2/4.0;
  }

  



  //  Tests
  int test_menn(){
    int n=0;//Number of fails
    std::cout<<"\t* M , e  -> ne, Nm      :";
    {//amp^2
      int i=0;
      // me=0.0005, mmu=0.106, pspatial=250
      ldouble me=0.0005, mmu=0.106, MW=80.385, WW=0.0;//WW=2.10;
      ldouble sinW=0.474;
      menn mennAmp = menn(0.31333,me,mmu,MW,WW,sinW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.978379120129875E-05,2.969860636913393E-04,9.190786416844540E-04,2.019332477854471E-03,3.767910537270882E-03,6.403200858818891E-03,1.026490446345585E-02,1.584671840415083E-02,2.388276913702579E-02,3.549440113192476E-02,5.244970641426275E-02,7.764471991824608E-02,1.160479212386511E-01,1.766878219287820E-01,2.772166113244994E-01,4.556382683643243E-01,8.054003517176767E-01,1.606117735285715E+00,4.015985608767513E+00,1.755146128790467E+01};
      i += mennAmp.test_2to2_amp2([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH);
      i += mennAmp.test_2to2_amp2_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH);
      i += mennAmp.test_2to2_amp2_boosts([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH);
      i += mennAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH);
      // me=0.0005, mmu=0.106, pspatial=0.05
      pspatial = 0.05;
      ldouble dataCH2[20] = {5.574841028950011E-15,1.791232924954341E-14,3.184787843851914E-14,4.738148983640861E-14,6.451316468374422E-14,8.324290422105935E-14,1.035707096888883E-13,1.254965823277664E-13,1.490205233782298E-13,1.741425340808159E-13,2.008626156760627E-13,2.291807694045093E-13,2.590969965066960E-13,2.906112982231637E-13,3.237236757944545E-13,3.584341304611115E-13,3.947426634636784E-13,4.326492760427002E-13,4.721539694387228E-13,5.132567448922932E-13};
      i += mennAmp.test_2to2_amp2([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH2);
      i += mennAmp.test_2to2_amp2_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH2);
      i += mennAmp.test_2to2_amp2_boosts([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH2);
      i += mennAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH2);
      // me=0.106, mmu=0.0005, pspatial=0.005
      me=0.106;
      mmu=0.0005;
      pspatial = 250;
      mennAmp.set_masses(me,mmu,MW);
      ldouble dataCH3[20] = {2.978379120129873E-05,2.969860636913392E-04,9.190786416844539E-04,2.019332477854471E-03,3.767910537270882E-03,6.403200858818891E-03,1.026490446345585E-02,1.584671840415083E-02,2.388276913702579E-02,3.549440113192476E-02,5.244970641426275E-02,7.764471991824608E-02,1.160479212386511E-01,1.766878219287819E-01,2.772166113244994E-01,4.556382683643243E-01,8.054003517176767E-01,1.606117735285715E+00,4.015985608767513E+00,1.755146128790467E+01};
      i += mennAmp.test_2to2_amp2([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH3);
      i += mennAmp.test_2to2_amp2_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH3);
      i += mennAmp.test_2to2_amp2_boosts([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH3);
      i += mennAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH3);
      // me=0.106, mmu=0.0005, pspatial=0.05
      me=0.106;
      mmu=0.0005;
      pspatial = 0.05;
      mennAmp.set_masses(me,mmu,MW);
      ldouble dataCH4[20] = {5.574841028950017E-15,1.791232924954341E-14,3.184787843851915E-14,4.738148983640862E-14,6.451316468374424E-14,8.324290422105936E-14,1.035707096888883E-13,1.254965823277664E-13,1.490205233782298E-13,1.741425340808159E-13,2.008626156760627E-13,2.291807694045094E-13,2.590969965066960E-13,2.906112982231637E-13,3.237236757944546E-13,3.584341304611115E-13,3.947426634636784E-13,4.326492760427002E-13,4.721539694387229E-13,5.132567448922932E-13};
      i += mennAmp.test_2to2_amp2([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH4);
      i += mennAmp.test_2to2_amp2_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH4);
      i += mennAmp.test_2to2_amp2_boosts([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH4);
      i += mennAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH4);
      // me=0.105, mmu=0.106, pspatial=250
      me=0.105;
      mmu=0.106;
      pspatial=250;
      mennAmp.set_masses(me,mmu,MW);
      ldouble dataCH5[20] = {2.978384384380531E-05,2.969862390647535E-04,9.190789682917295E-04,2.019332992163374E-03,3.767911286528968E-03,6.403201905116417E-03,1.026490588972179E-02,1.584672032322839E-02,2.388277170581320E-02,3.549440457416060E-02,5.244971105871386E-02,7.764472626577304E-02,1.160479300868339E-01,1.766878346202236E-01,2.772166302871948E-01,4.556382984318614E-01,8.054004039105621E-01,1.606117840504570E+00,4.015985891009827E+00,1.755146286864769E+01};
      i += mennAmp.test_2to2_amp2([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH5);
      i += mennAmp.test_2to2_amp2_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH5);
      i += mennAmp.test_2to2_amp2_boosts([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH5);
      i += mennAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH5);
      // me=0.105, mmu=0.106, pspatial=0.05
      pspatial=0.05;
      ldouble dataCH6[20] = {2.988853213248676E-13,3.436064037073394E-13,3.914440133451065E-13,4.423981536166547E-13,4.964688279004736E-13,5.536560395750563E-13,6.139597920188997E-13,6.773800886105044E-13,7.439169327283747E-13,8.135703277510181E-13,8.863402770569463E-13,9.622267840246746E-13,1.041229852032722E-12,1.123349484459610E-12,1.208585684683866E-12,1.296938456084019E-12,1.388407802038603E-12,1.482993725926155E-12,1.580696231125216E-12,1.681515321014330E-12};
      i += mennAmp.test_2to2_amp2([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH6);
      i += mennAmp.test_2to2_amp2_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH6);
      i += mennAmp.test_2to2_amp2_boosts([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH6);
      i += mennAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH6);
      // me=105, mmu=106, MW=0.1, pspatial=500
      me=105;
      mmu=106;
      MW=0.1;
      pspatial=500;
      mennAmp.set_masses(me,mmu,MW);
      ldouble dataCH7[20] = {1.511513540736549E+10,1.513318718353640E+10,1.515331443327512E+10,1.517589690849148E+10,1.520141301879668E+10,1.523047409357314E+10,1.526387398719470E+10,1.530266268942818E+10,1.534825868161382E+10,1.540262611754188E+10,1.546856510021413E+10,1.555020928895327E+10,1.565392691299386E+10,1.579006623606122E+10,1.597663938288907E+10,1.624802474588694E+10,1.667899941381540E+10,1.746863610173603E+10,1.938038535279377E+10,3.038177336148563E+10};
      i += mennAmp.test_2to2_amp2([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH7);
      i += mennAmp.test_2to2_amp2_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH7);
      i += mennAmp.test_2to2_amp2_boosts([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH7);
      i += mennAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH7);
      // me=105, mmu=106, MW=0.1, pspatial=5
      pspatial=5;
      ldouble dataCH8[20] = {5.413153965989237E+10,5.458446914728065E+10,5.504732731653917E+10,5.552043146558590E+10,5.600411218172909E+10,5.649871403027652E+10,5.700459628571535E+10,5.752213370852575E+10,5.805171737094188E+10,5.859375553524999E+10,5.914867458851257E+10,5.971692003793858E+10,6.029895757147826E+10,6.089527418861777E+10,6.150637940678146E+10,6.213280654922648E+10,6.277511412083735E+10,6.343388727880482E+10,6.410973940580651E+10,6.480331379400745E+10};
      i += mennAmp.test_2to2_amp2([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH8);
      i += mennAmp.test_2to2_amp2_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH8);
      i += mennAmp.test_2to2_amp2_boosts([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH8);
      i += mennAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH8);
      // me=0.105, mmu=0.106, MW=0.11, pspatial=250
      me=0.105;
      mmu=0.106;
      MW=0.11;
      pspatial=250;
      mennAmp.set_masses(me,mmu,MW);
      ldouble dataCH9[20] = {1.012842765556309E-02,1.041085926583142E-02,1.107122319110112E-02,1.224488951158999E-02,1.412047360174428E-02,1.696494620082339E-02,2.116312555726878E-02,2.728155429912584E-02,3.617519553965640E-02,4.917242516617579E-02,6.841012340304178E-02,9.747309169596709E-02,1.426935387274610E-01,2.160060500814319E-01,3.418731558769506E-01,5.764312798818740E-01,1.070976863593772E+00,2.349097828058597E+00,7.271072978395182E+00,7.261436443985059E+01};
      i += mennAmp.test_2to2_amp2([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH9);
      i += mennAmp.test_2to2_amp2_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH9);
      i += mennAmp.test_2to2_amp2_boosts([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH9);
      i += mennAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH9);
      // me=0.105, mmu=0.106, MW=0.11, pspatial=0.005
      pspatial=0.005;
      ldouble dataCH10[20] = {8.303602651083473E-02,8.404665544647953E-02,8.507408974500479E-02,8.611865698042978E-02,8.718069240104087E-02,8.826053914184739E-02,8.935854844385668E-02,9.047507988041842E-02,9.161050159089781E-02,9.276519052194866E-02,9.393953267666819E-02,9.513392337192716E-02,9.634876750418181E-02,9.758447982408611E-02,9.884148522023715E-02,1.001202190124002E-01,1.014211272545746E-01,1.027446670482780E-01,1.040913068664413E-01,1.054615268883256E-01};
      i += mennAmp.test_2to2_amp2([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH10);
      i += mennAmp.test_2to2_amp2_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH10);
      i += mennAmp.test_2to2_amp2_boosts([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH10);
      i += mennAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mennAmp.amp2(); }, mmu,me,0,0,pspatial,dataCH10);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
