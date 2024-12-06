
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

//File:  SPINAS/SM/eeZh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eeZh.h"

namespace spinas {

  eeZh::eeZh(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), prope(masse,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);  
    p1=particle(me);
    p2=particle(me);
    p3=particle(MZ);
    p4=particle(mh);
    //<12>,[12],<23>,[23],<13>,[13]
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    //[312>,[213>,[343>,[143>,[341>
    s312a = sproduct(SQUARE,&p3,&p1,&p2,2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3,2);
    s343a = sproduct(SQUARE,&p3,&p4,&p3,2);
    s143a = sproduct(SQUARE,&p1,&p4,&p3,2);
    s341a = sproduct(SQUARE,&p3,&p4,&p1,2);
    s243a = sproduct(SQUARE,&p2,&p4,&p3,2);
    s342a = sproduct(SQUARE,&p3,&p4,&p2,2);
    //Couplings
    preTU = sqrt2*e*e*me/(4.0*MW*MW*SW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void eeZh::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(me);
    p2.set_mass(me);
    p3.set_mass(MZ);
    p4.set_mass(mh);
    prope.set_mass(me);
    propZ.set_mass(MZ);
    //Couplings
    preTU = sqrt2*e*e*me/(4.0*MW*MW*SW*SW);
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void eeZh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    //[312>,[213>,[343>,[143>,[341>
    s312a.update();
    s213a.update();
    s343a.update();
    s143a.update();
    s341a.update();
    s243a.update();
    s342a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propZ.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=prope.denominator(propTP);
    pDenU=prope.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eeZh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;

    //Symmetrize the Z-Boson Spin indices
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      //S-Channel Z
      //preZ = e*e/(4.0*MW*MW*SW*SW);
      //all ingoing: 
      //preZ (-Me (gLe-gRe) (<12>-[12]) (MZ [33]-[343>) + 2 MZ^2(gLe[23] <13> + gRe[13] <23>)))/(s-MZ^2)
      //34 outgoing:
      //preZ (-Me (gLe-gRe) (<12>-[12]) (MZ [33]-[343>) - 2 MZ^2(gLe[23] <13> + gRe[13] <23>)))/(s-MZ^2)
      //[33]=0 since it is symmetrized
      //preZ ( Me (gLe-gRe) (<12>-[12]) [343> - 2 MZ^2(gLe[23] <13> + gRe[13] <23>)))/(s-MZ^2)
      amplitude += normFactor*preZ*( me*(gL-gR)*(a12a.v(ds1,ds2)-s12s.v(ds1,ds2))*s343a.v(ds3a,ds3b) - 2.0*MZ*MZ*(gL*s23s.v(ds2,ds3a)*a13a.v(ds1,ds3b) + gR*a23a.v(ds2,ds3a)*s13s.v(ds1,ds3b)))/pDenS;
      
      //T-Channel e
      //preTU = e*e*me/(4.0*MW*MW*SW*SW);
      //all ingoing:
      //preh*(gLe <13> (Me [23]-[312>+MZ <23>)+gRe [13] (MZ [23]-[213>+Me <23>)))/(t-Me^2)
      //34 outgoing:
      //- preh*( gLe <13> (Me [23]-[312>-MZ <23>) + gRe [13] (Me <23>-[213>-MZ [23])) )/(t-Me^2)
      amplitude += normFactor*preTU*(
				     gR*s13s.v(ds1,ds3a)*(-2.0*me*a23a.v(ds2,ds3b)+s243a.v(ds2,ds3b))
				     +gL*a13a.v(ds1,ds3a)*(-2.0*me*s23s.v(ds2,ds3b)+s342a.v(ds3b,ds2))
				     )/pDenT;

      //U-Channel e
      //preTU = e*e*me/(4.0*MW*MW*SW*SW);
      //all ingoing:
      //preh (gLe [23] ([143>+2 Me <13>)+gRe <23> (2 Me [13]+[341>) )/(u-Me^2)
      //34 outgoing:
      //preh (gLe [23] ([143>-2 Me <13>)+gRe <23> ([341> - 2 Me [13]) )/(u-Me^2)
      amplitude += normFactor*preTU*(
				     gL*s23s.v(ds2,ds3a)*(s143a.v(ds1,ds3b)-2.0*me*a13a.v(ds1,ds3b))
				     + gR*a23a.v(ds2,ds3a)*(s341a.v(ds3b,ds1)-2.0*me*s13s.v(ds1,ds3b))
				     )/pDenU;
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eeZh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2^2=1/4
    return amp2/4.0;
  }
  



  //  Tests
  int test_eeZh(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> Z , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      eeZh eeZhAmp = eeZh(EE,me,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {9.271793604400342E-04,1.317945762319436E-03,1.665293677057919E-03,1.969223102757450E-03,2.229734039165127E-03,2.446826486208593E-03,2.620500443859534E-03,2.750755912104809E-03,2.837592890937753E-03,2.881011380354956E-03,2.881011380354956E-03,2.837592890937754E-03,2.750755912104810E-03,2.620500443859535E-03,2.446826486208595E-03,2.229734039165128E-03,1.969223102757452E-03,1.665293677057921E-03,1.317945762319438E-03,9.271793604400367E-04};
      i += eeZhAmp.test_2to2_amp2([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH);
      i += eeZhAmp.test_2to2_amp2_rotations([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH);
      i += eeZhAmp.test_2to2_amp2_boosts([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH);
      i += eeZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {3.633990947962139E-03,3.782973777474111E-03,3.915402959110596E-03,4.031278492991858E-03,4.130600379155819E-03,4.213368617617154E-03,4.279583208382356E-03,4.329244151454567E-03,4.362351446835388E-03,4.378905094525628E-03,4.378905094525628E-03,4.362351446835388E-03,4.329244151454567E-03,4.279583208382356E-03,4.213368617617155E-03,4.130600379155819E-03,4.031278492991858E-03,3.915402959110596E-03,3.782973777474111E-03,3.633990947962139E-03};
      i += eeZhAmp.test_2to2_amp2([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH2);
      i += eeZhAmp.test_2to2_amp2_rotations([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH2);
      i += eeZhAmp.test_2to2_amp2_boosts([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH2);
      i += eeZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, mh=125, MW=80.385, pspatial=95\n";
      me = 125.1;
      mh = 125;
      pspatial = 95;
      eeZhAmp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {1.568597011273475E-01,1.483365051238541E-01,1.406769388877999E-01,1.340773415251518E-01,1.285599012354351E-01,1.240833528407889E-01,1.205893102542819E-01,1.180219407445950E-01,1.163362083912132E-01,1.155011686277072E-01,1.155011686277072E-01,1.163362083912132E-01,1.180219407445950E-01,1.205893102542819E-01,1.240833528407889E-01,1.285599012354352E-01,1.340773415251518E-01,1.406769388877999E-01,1.483365051238542E-01,1.568597011273475E-01};
      i += eeZhAmp.test_2to2_amp2([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH4);
      i += eeZhAmp.test_2to2_amp2_rotations([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH4);
      i += eeZhAmp.test_2to2_amp2_boosts([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH4);
      i += eeZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH4);
      //std::cout<<"\n# me=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      me = 125;
      mh = 0.0005;
      pspatial = 125.1;
      eeZhAmp.set_masses(me,mh,MW);
      ldouble dataCH3[20] = {2.610266484110628E-01,2.275216119944078E-01,2.020741496294103E-01,1.827273604496982E-01,1.680074138067813E-01,1.568848794310616E-01,1.486534249093469E-01,1.428342582313605E-01,1.391123553751880E-01,1.372965629574767E-01,1.372965629574767E-01,1.391123553751881E-01,1.428342582313605E-01,1.486534249093469E-01,1.568848794310616E-01,1.680074138067814E-01,1.827273604496982E-01,2.020741496294103E-01,2.275216119944078E-01,2.610266484110628E-01};
      i += eeZhAmp.test_2to2_amp2([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH3);
      i += eeZhAmp.test_2to2_amp2_rotations([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH3);
      i += eeZhAmp.test_2to2_amp2_boosts([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH3);
      i += eeZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeZhAmp.amp2(); }, me,me,MZ,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
