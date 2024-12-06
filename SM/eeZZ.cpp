
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

//File:  SPINAS/SM/eeZZ.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eeZZ.h"

namespace spinas {

  eeZZ::eeZZ(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), prope(masse,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    proph = propagator(mh,wh);  
    p1=particle(me);
    p2=particle(me);
    p3=particle(MZ);
    p4=particle(MZ);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s34s = sproduct(SQUARE,&p3,&p4,2);
    a34a = sproduct(ANGLE,&p3,&p4,2);
    s24s = sproduct(SQUARE,&p2,&p4,2);
    a24a = sproduct(ANGLE,&p2,&p4,2);
    s14s = sproduct(SQUARE,&p1,&p4,2);
    a14a = sproduct(ANGLE,&p1,&p4,2);
    //[314>,[413>
    s314a = sproduct(SQUARE,&p3,&p1,&p4,2);
    s413a = sproduct(SQUARE,&p4,&p1,&p3,2);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
    preh = e*e*me/(2.0*MW*MW*SW*SW);
  }
  void eeZZ::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(me);
    p2.set_mass(me);
    p3.set_mass(MZ);
    p4.set_mass(MZ);
    prope.set_mass(me);
    proph.set_mass(mh);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
    preh = e*e*me/(2.0*MW*MW*SW*SW);
  }
  void eeZZ::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s34s.update();
    a34a.update();
    s24s.update();
    a24a.update();
    s14s.update();
    a14a.update();
    //[314>,[413>
    s314a.update();
    s413a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=proph.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=prope.denominator(propTP);
    pDenU=prope.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eeZZ::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);
      
      //S-Channel h
      //preh = e*e*me/(2.0*MW*MW*SW*SW);
      //all ingoing: 
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //34 outgoing:
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      amplitude += normFactor*preh*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))/pDenS;
      
      //T-Channel e
      //preTU = e*e/(2.0*MW*MW*SW*SW);
      //all ingoing:
      //+preTU gLe^2 [24] <13> (MZ <34>+[314>))/(Me^2-t)
      //+preTU gRe^2 [13] <24> (MZ [34]+[413>))/(Me^2-t)
      //+preTU gLe gRe Me ([13] [24] <34>+[34] <13> <24>))/(Me^2-t)
      //34 outgoing:
      //+preTU gLe^2 [24] <13> (MZ <34>-[314>))/(t-Me^2)
      //+preTU gRe^2 [13] <24> (MZ [34]-[413>))/(t-Me^2)
      //-preTU gLe gRe Me ([13] [24] <34>+[34] <13> <24>))/(t-Me^2)
      amplitude += normFactor*preTU*gL*gL*s24s.v(ds2,ds4a)*a13a.v(ds1,ds3a)*(MZ*a34a.v(ds3b,ds4b)-s314a.v(ds3b,ds4b))/pDenT
      	+          normFactor*preTU*gR*gR*a24a.v(ds2,ds4a)*s13s.v(ds1,ds3a)*(MZ*s34s.v(ds3b,ds4b)-s413a.v(ds4b,ds3b))/pDenT
      	-          normFactor*preTU*gL*gR*me*(s13s.v(ds1,ds3a)*s24s.v(ds2,ds4a)*a34a.v(ds3b,ds4b)+a13a.v(ds1,ds3a)*a24a.v(ds2,ds4a)*s34s.v(ds3b,ds4b))/pDenT;

      //U-Channel e
      //preTU = e*e/(2.0*MW*MW*SW*SW);
      //all ingoing:
      //+preTU gLe^2 [23] <14> (MZ <34>-[413>)/(u-Me^2)
      //+preTU gRe^2 [14] <23> (MZ [34]-[314>)/(u-Me^2)
      //+preTU gLe gRe Me ([14] [23] <34>+[34] <14> <23>)/(u-Me^2)
      //34 outgoing:
      //-preTU gLe^2 [23] <14> (MZ <34>+[413>)/(u-Me^2)
      //-preTU gRe^2 [14] <23> (MZ [34]+[314>)/(u-Me^2)
      //+preTU gLe gRe Me ([14] [23] <34>+[34] <14> <23>)/(u-Me^2)
      amplitude += -normFactor*preTU*gL*gL*s23s.v(ds2,ds3a)*a14a.v(ds1,ds4a)*(MZ*a34a.v(ds3b,ds4b)+s413a.v(ds4b,ds3b))/pDenU
      	-           normFactor*preTU*gR*gR*s14s.v(ds1,ds4a)*a23a.v(ds2,ds3a)*(MZ*s34s.v(ds3b,ds4b)+s314a.v(ds3b,ds4b))/pDenU
      	+           normFactor*preTU*gL*gR*me*(s14s.v(ds1,ds4a)*s23s.v(ds2,ds3a)*a34a.v(ds3b,ds4b)+a14a.v(ds1,ds4a)*a23a.v(ds2,ds3a)*s34s.v(ds3b,ds4b))/pDenU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eeZZ::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2^2=1/4
    //Symmetry Factor 1/2
    return amp2/8.0;
  }
  



  //  Tests
  int test_eeZZ(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> Z , Z       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      eeZZ eeZZAmp = eeZZ(EE,me,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {5.442211900637754E-02,1.868054814813991E-02,1.105414103484835E-02,7.790196240887772E-03,6.021116705930904E-03,4.948703454572987E-03,4.264845455102644E-03,3.827764746157513E-03,3.565896158217935E-03,3.442762308623416E-03,3.442762308623423E-03,3.565896158217951E-03,3.827764746157509E-03,4.264845455102645E-03,4.948703454572981E-03,6.021116705930900E-03,7.790196240887781E-03,1.105414103484834E-02,1.868054814813991E-02,5.442211900637747E-02};
      i += eeZZAmp.test_2to2_amp2([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH);
      i += eeZZAmp.test_2to2_amp2_rotations([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH);
      i += eeZZAmp.test_2to2_amp2_boosts([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH);
      i += eeZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {3.258950762978131E-02,2.142583969855804E-02,1.588029093209552E-02,1.273679122662298E-02,1.077465036653927E-02,9.482225317816721E-03,8.613595120711226E-03,8.039190132955377E-03,7.687590710029004E-03,7.520289758877131E-03,7.520289758877132E-03,7.687590710029001E-03,8.039190132955374E-03,8.613595120711222E-03,9.482225317816719E-03,1.077465036653927E-02,1.273679122662298E-02,1.588029093209552E-02,2.142583969855804E-02,3.258950762978130E-02};
      i += eeZZAmp.test_2to2_amp2([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH2);
      i += eeZZAmp.test_2to2_amp2_rotations([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH2);
      i += eeZZAmp.test_2to2_amp2_boosts([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH2);
      i += eeZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, mh=125, MW=80.385, pspatial=95\n";
      me = 125.1;
      mh = 125;
      pspatial = 95;
      eeZZAmp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {5.372710974321059E-02,4.114944307347489E-02,3.144917921221758E-02,2.390420680192585E-02,1.803433498128241E-02,1.351243143054706E-02,1.011282197721534E-02,7.680308488648192E-03,6.111210884096534E-03,5.341830647622328E-03,5.341830647622331E-03,6.111210884096532E-03,7.680308488648188E-03,1.011282197721534E-02,1.351243143054705E-02,1.803433498128241E-02,2.390420680192584E-02,3.144917921221758E-02,4.114944307347488E-02,5.372710974321056E-02};
      i += eeZZAmp.test_2to2_amp2([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH4);
      i += eeZZAmp.test_2to2_amp2_rotations([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH4);
      i += eeZZAmp.test_2to2_amp2_boosts([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH4);
      i += eeZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH4);
      //std::cout<<"\n# me=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      me = 125;
      mh = 0.0005;
      pspatial = 125.1;
      eeZZAmp.set_masses(me,mh,MW);
      ldouble dataCH3[20] = {7.669164781419931E-02,5.564217455921953E-02,4.086806859978095E-02,3.015468167297349E-02,2.224935111984145E-02,1.639845803401539E-02,1.212965956036441E-02,9.141056319706103E-03,7.241596845325212E-03,6.318208787414945E-03,6.318208787414933E-03,7.241596845325210E-03,9.141056319706093E-03,1.212965956036439E-02,1.639845803401538E-02,2.224935111984144E-02,3.015468167297349E-02,4.086806859978094E-02,5.564217455921951E-02,7.669164781419929E-02};
      i += eeZZAmp.test_2to2_amp2([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH3);
      i += eeZZAmp.test_2to2_amp2_rotations([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH3);
      i += eeZZAmp.test_2to2_amp2_boosts([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH3);
      i += eeZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeZZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
