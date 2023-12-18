
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

//File:  SPINAS/SM/hgdd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/hgdd.h"

namespace spinas {

  hgdd::hgdd(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd, const ldouble& massh, const ldouble& massW, const ldouble& sinW):
    e(echarge), gs(gscharge), md(massd), mh(massh), prop(massd,0), MW(massW), SW(sinW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(mh);
    p2=particle(0);
    p3=particle(md);
    p4=particle(md);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
    a213s = sproduct(ANGLE,&p2,&p1,&p3);
    s214a = sproduct(SQUARE,&p2,&p1,&p4);
    a214s = sproduct(ANGLE,&p2,&p1,&p4);
    s2312s = sproduct(SQUARE,&p2,&p3,&p1,&p2);
    a2312a = sproduct(ANGLE,&p2,&p3,&p1,&p2);
    pre = sqrt2*e*gs*md/(2.0*MW*SW);
  }
  void hgdd::set_masses(const ldouble& massd, const ldouble& massh, const ldouble& massW){
    md=massd;
    mh=massh;
    p1.set_mass(mh);
    p3.set_mass(md);
    p4.set_mass(md);
    prop.set_mass(md);
    pre = sqrt2*e*gs*md/(2.0*MW*SW);
  }
  void hgdd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s24s.update();
    a24a.update();
    s23s.update();
    a23a.update();
    s34s.update();
    a34a.update();
    s213a.update();
    a213s.update();
    s214a.update();
    a214s.update();
    s2312s.update();
    a2312a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT=prop.denominator(propTP);
    pDenU=prop.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble hgdd::amp(const int& ds2, const int& ds3, const int& ds4){//Double Spin
    //pre = sqrt2*e*gs*md/(2.0*MW*SW);
    if(ds2>0){
      //dgDh all in:
      //( + mh^2[12][23] - md[12][243> + md[23][241> - <13>[2342] )/su
      //hgDd: 1<->4:
      //( - mh^2[24][23] + md[24][213> + md[23][214> + <34>[2312] )/tu
      //34 out:
      //- ( + mh^2[24][23] + md[24][213> + md[23][214> + <34>[2312] )/tu
      return - pre*(mh*mh*s24s.v(ds4)*s23s.v(ds3) + md*s24s.v(ds4)*s213a.v(ds3) + md*s23s.v(ds3)*s214a.v(ds4) + a34a.v(ds3,ds4)*s2312s.v())/pDenT/pDenU;
    }
    else if(ds2<0){
      //dgDh all in:
      //( + mh^2<12><23> - md<12><243] + md<23><241] - [13]<2342> )/su
      //hgDd: 1<->4
      //( - mh^2<24><23> + md<24><213] + md<23><214] + [34]<2312> )/tu
      //34 out:
      //- ( + mh^2<24><23> + md<24><213] + md<23><214] + [34]<2312> )/tu
      return - pre*(mh*mh*a24a.v(ds4)*a23a.v(ds3) + md*a24a.v(ds4)*a213s.v(ds3) + md*a23a.v(ds3)*a214s.v(ds4) + s34s.v(ds3,ds4)*a2312a.v())/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble hgdd::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=4)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j2,j3,j4);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(TaTa) = 4 //Both channels have the same color factor
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/8
    return amp2/16.0;
  }
  
  //set_momenta(...) must be called before amp2_Aplus().
  //Positive Helicity Photon Only
  ldouble hgdd::amp2_gplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(2,j3,j4);
	amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(TaTa) = 4
      }
    //Average over initial colors 1/8
    return amp2/8.0;
  }
  



  //  Tests
  int test_hgdd(){
    int n=0;//Number of fails
    std::cout<<"\t* h , g  -> d , D       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# md=0.0075, mh=125, pspatial=250\n";
      ldouble md=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble gs=1.238;
      hgdd hgddAmp = hgdd(EE,gs,md,mh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {6.726806891516043E-08,2.363472705633797E-08,1.499116975049329E-08,1.135694678628791E-08,9.403063470714620E-09,8.223995952194059E-09,7.474229940220273E-09,6.995879224501135E-09,6.709602837058897E-09,6.575074459493757E-09,6.575074459493756E-09,6.709602837058893E-09,6.995879224501131E-09,7.474229940220273E-09,8.223995952194059E-09,9.403063470714619E-09,1.135694678628791E-08,1.499116975049328E-08,2.363472705633795E-08,6.726806891516018E-08};
      i += hgddAmp.test_2to2_amp2([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH);
      i += hgddAmp.test_2to2_amp2_rotations([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH);
      i += hgddAmp.test_2to2_amp2_boosts([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH);
      i += hgddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += hgddAmp.test_2to2_amp2([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH);
      i += hgddAmp.test_2to2_amp2_rotations([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH);
      i += hgddAmp.test_2to2_amp2_boosts([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH);
      i += hgddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH);
      //std::cout<<"\n# md=0.0042, mh=125, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {8.964297848141621E-08,3.149618231790152E-08,1.997757858539258E-08,1.513452925894535E-08,1.253073929119191E-08,1.095948672197510E-08,9.960331241096254E-09,9.322870043902889E-09,8.941371526766171E-09,8.762095907773109E-09,8.762095907773109E-09,8.941371526766170E-09,9.322870043902884E-09,9.960331241096253E-09,1.095948672197510E-08,1.253073929119191E-08,1.513452925894534E-08,1.997757858539255E-08,3.149618231790149E-08,8.964297848141598E-08};
      i += hgddAmp.test_2to2_amp2([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH2);
      i += hgddAmp.test_2to2_amp2_rotations([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH2);
      i += hgddAmp.test_2to2_amp2_boosts([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH2);
      i += hgddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH2);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += hgddAmp.test_2to2_amp2([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH2);
      i += hgddAmp.test_2to2_amp2_rotations([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH2);
      i += hgddAmp.test_2to2_amp2_boosts([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH2);
      i += hgddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH2);
      //std::cout<<"\n# md=125.1, mh=125, pspatial=95\n";
      md = 125.1;
      mh = 125;
      pspatial = 95;
      hgddAmp.set_masses(md,mh,MW);
      ldouble dataCH4[20] = {1.637131935476594E+00,1.621840087462764E+00,1.608361646035582E+00,1.596655403884674E+00,1.586685849449741E+00,1.578422943131582E+00,1.571841929706476E+00,1.566923184933916E+00,1.563652094729969E+00,1.562018965634304E+00,1.562018965634304E+00,1.563652094729969E+00,1.566923184933916E+00,1.571841929706476E+00,1.578422943131582E+00,1.586685849449741E+00,1.596655403884674E+00,1.608361646035582E+00,1.621840087462764E+00,1.637131935476594E+00};
      i += hgddAmp.test_2to2_amp2([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH4);
      i += hgddAmp.test_2to2_amp2_rotations([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH4);
      i += hgddAmp.test_2to2_amp2_boosts([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH4);
      i += hgddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH4);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += hgddAmp.test_2to2_amp2([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH4);
      i += hgddAmp.test_2to2_amp2_rotations([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH4);
      i += hgddAmp.test_2to2_amp2_boosts([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH4);
      i += hgddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH4);
      //std::cout<<"\n# md=125, mh=0.0005, pspatial=125.1\n";
      md = 125;
      mh = 0.0005;
      pspatial = 125.1;
      hgddAmp.set_masses(md,mh,MW);
      ldouble dataCH3[20] = {1.621250195430904E+00,1.619851112875714E+00,1.618608610780647E+00,1.617522290501720E+00,1.616591803810408E+00,1.615816852658952E+00,1.615197188979711E+00,1.614732614518298E+00,1.614422980700322E+00,1.614268188531580E+00,1.614268188531581E+00,1.614422980700322E+00,1.614732614518298E+00,1.615197188979711E+00,1.615816852658952E+00,1.616591803810408E+00,1.617522290501720E+00,1.618608610780647E+00,1.619851112875715E+00,1.621250195430904E+00};
      i += hgddAmp.test_2to2_amp2([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH3);
      i += hgddAmp.test_2to2_amp2_rotations([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH3);
      i += hgddAmp.test_2to2_amp2_boosts([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH3);
      i += hgddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hgddAmp.amp2(); }, mh,0,md,md,pspatial,dataCH3);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += hgddAmp.test_2to2_amp2([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH3);
      i += hgddAmp.test_2to2_amp2_rotations([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH3);
      i += hgddAmp.test_2to2_amp2_boosts([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH3);
      i += hgddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hgddAmp.amp2_gplus(); }, mh,0,md,md,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
