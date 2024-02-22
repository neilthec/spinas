
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

//File:  SPINAS/SM/eehh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eehh.h"

namespace spinas {
  //Constructors
  eehh::eehh(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW):
    e(echarge), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW),
    prope(me,0), proph(mh,wh),
    p1(particle(me)), p2(particle(me)),
    p3(particle(mh)), p4(particle(mh)),
    s12s(sproduct(SQUARE,&p1,&p2)),
    a12a(sproduct(ANGLE,&p1,&p2)),
    s132a(sproduct(SQUARE,&p1,&p3,&p2)),
    s231a(sproduct(SQUARE,&p2,&p3,&p1))
  {
    prehS = 3*e*e*me*mh*mh/(4*MW*MW*SW*SW);
    prehTU = e*e*me*me/(4.0*MW*MW*SW*SW);
  }
  void eehh::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    prope.set_mass(me);
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    p1.set_mass(me);
    p2.set_mass(me);
    p3.set_mass(mh);
    p4.set_mass(mh);
    prehS = 3.0*e*e*me*mh*mh/(4.0*MW*MW*SW*SW);
    prehTU = e*e*me*me/(4.0*MW*MW*SW*SW);
  }
  void eehh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s12s.update();
    a12a.update();
    s132a.update();
    s231a.update();
    //Propagator Momentum
    ldouble propPS[4], propPT[4], propPU[4];
    for(int j=0;j<4;j++){
      propPS[j] = mom1[j]+mom2[j];
      propPT[j] = mom1[j]-mom3[j];
      propPU[j] = mom1[j]-mom4[j];
    }
    pDenSh = proph.denominator(propPS);
    pDenTe = prope.denominator(propPT);
    pDenUe = prope.denominator(propPU);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eehh::amp(const int& ds1, const int& ds2){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Higgs
    //prehS = 3*e*e*me*mh*mh/(4*MW*MW*SW*SW);
    //S Channel
    //all ingoing: - prehS*([12]+<12>)/(s-mh^2)
    amplitude += - prehS*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))/pDenSh;

    //T Channel
    //prehTU = e*e*me*me/(4*MW*MW*SW*SW);
    //all ingoing:  - prehTU * ( 2me([12]+<12>) + [132> + [231> )/(t-me^2)
    //34 outgoing:  - prehTU * ( 2me([12]+<12>) - [132> - [231> )/(t-me^2)
    amplitude += - prehTU*( two*me*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2)) - s132a.v(ds1,ds2) - s231a.v(ds2,ds1) )/pDenTe;

    //U Channel
    //all ingoing:  - prehTU * ( 2me([12]+<12>) + [142> + [241> )/(u-me^2)
    //all ingoing:  - prehTU * ( 2me([12]+<12>) - [132> - [231> )/(u-me^2)
    //34 outgoing:  - prehTU * ( 2me([12]+<12>) + [132> + [231> )/(u-me^2)
    amplitude += -prehTU*( two*me*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2)) + s132a.v(ds1,ds2) + s231a.v(ds2,ds1) )/pDenUe;

    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eehh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2){
	M = amp(j1,j2);
	amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Symmetry Factor 1/2
    return amp2/8.0;
  }

  



  //  Tests
  int test_eehh(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> h , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### me=0.0005, mh=125, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, mh=125, wh=0, MW=80.385, SW=0.474;//Set width to 0 for comparison with Feynman diagrams.
      eehh eehhAmp = eehh(0.31333,me,mh,wh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.786300494527163E-13,1.786300494285207E-13,1.786300494201812E-13,1.786300494162168E-13,1.786300494139650E-13,1.786300494125639E-13,1.786300494116561E-13,1.786300494110700E-13,1.786300494107166E-13,1.786300494105498E-13,1.786300494105498E-13,1.786300494107166E-13,1.786300494110700E-13,1.786300494116561E-13,1.786300494125639E-13,1.786300494139650E-13,1.786300494162168E-13,1.786300494201812E-13,1.786300494285207E-13,1.786300494527163E-13};
      i += eehhAmp.test_2to2_amp2([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH);
      i += eehhAmp.test_2to2_amp2_rotations([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH);
      i += eehhAmp.test_2to2_amp2_boosts([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH);
      i += eehhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH);
      //std::cout<<"\n########### me=0.0005, mh=125, MW=80.385, pspatial=126\n";
      pspatial = 126;
      ldouble dataCH2[20] = {1.087296572585901E-12,1.087296572587164E-12,1.087296572588259E-12,1.087296572589196E-12,1.087296572589984E-12,1.087296572590630E-12,1.087296572591140E-12,1.087296572591519E-12,1.087296572591770E-12,1.087296572591895E-12,1.087296572591895E-12,1.087296572591770E-12,1.087296572591519E-12,1.087296572591140E-12,1.087296572590630E-12,1.087296572589984E-12,1.087296572589196E-12,1.087296572588259E-12,1.087296572587164E-12,1.087296572585901E-12};
      i += eehhAmp.test_2to2_amp2([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH2);
      i += eehhAmp.test_2to2_amp2_rotations([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH2);
      i += eehhAmp.test_2to2_amp2_boosts([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH2);
      i += eehhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH2);
      //std::cout<<"\n########### me=0.1, mh=0.15, MW=80.385, pspatial=0.2\n";
      me=0.10;
      mh=0.15;
      pspatial = 0.2;
      eehhAmp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {6.375457705183125E-14,3.898994174573297E-14,2.553100804952224E-14,1.725378399826743E-14,1.179585714987176E-14,8.064156016079406E-15,5.490759467066099E-15,3.758848438571753E-15,2.686518398551556E-15,2.172959659802411E-15,2.172959659802407E-15,2.686518398551551E-15,3.758848438571750E-15,5.490759467066096E-15,8.064156016079404E-15,1.179585714987175E-14,1.725378399826741E-14,2.553100804952223E-14,3.898994174573298E-14,6.375457705183121E-14};
      i += eehhAmp.test_2to2_amp2([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH4);
      i += eehhAmp.test_2to2_amp2_rotations([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH4);
      i += eehhAmp.test_2to2_amp2_boosts([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH4);
      i += eehhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH4);
      //std::cout<<"\n########### me=0.1, mh=0.15, MW=0.11, pspatial=0.2\n";
      me=0.10;
      mh=0.15;
      MW=0.11;
      pspatial = 0.2;
      eehhAmp.set_masses(me,mh,MW);
      ldouble dataCH5[20] = {1.818195978042583E-02,1.111941425139924E-02,7.281104870835434E-03,4.920550354551248E-03,3.364022007396509E-03,2.299790338633001E-03,1.565891774531357E-03,1.071973720024310E-03,7.661594152232771E-04,6.196992744049647E-04,6.196992744049635E-04,7.661594152232759E-04,1.071973720024309E-03,1.565891774531356E-03,2.299790338633000E-03,3.364022007396507E-03,4.920550354551244E-03,7.281104870835432E-03,1.111941425139924E-02,1.818195978042582E-02};
      i += eehhAmp.test_2to2_amp2([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH5);
      i += eehhAmp.test_2to2_amp2_rotations([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH5);
      i += eehhAmp.test_2to2_amp2_boosts([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH5);
      i += eehhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH5);
      //std::cout<<"\n########### me=0.1, mh=0.15, MW=0.006, pspatial=0.2\n";
      me=0.10;
      mh=0.15;
      MW=0.006;
      pspatial = 0.2;
      eehhAmp.set_masses(me,mh,MW);
      ldouble dataCH6[20] = {2.054028342169865E+03,1.256167778200125E+03,8.225513612183764E+02,5.558779146680927E+02,3.800358503880577E+02,2.598088761414025E+02,1.768998570286543E+02,1.211015990345364E+02,8.655354936947531E+01,7.000784781298680E+01,7.000784781298665E+01,8.655354936947516E+01,1.211015990345363E+02,1.768998570286542E+02,2.598088761414025E+02,3.800358503880575E+02,5.558779146680923E+02,8.225513612183761E+02,1.256167778200125E+03,2.054028342169864E+03};
      i += eehhAmp.test_2to2_amp2([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH6);
      i += eehhAmp.test_2to2_amp2_rotations([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH6);
      i += eehhAmp.test_2to2_amp2_boosts([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH6);
      i += eehhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eehhAmp.amp2(); }, me,me,mh,mh,pspatial,dataCH6);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
