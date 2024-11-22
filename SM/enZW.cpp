
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

//File:  SPINAS/SM/enZW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/enZW.h"

namespace spinas {

  enZW::enZW(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), me(masse), MW(massW), SW(sinW), WW(widthW) {
    constexpr ldouble one=1, two=2;//, sqrt2=std::sqrt(2);
    CW=std::sqrt(1-SW*SW);
    MZ=MW/CW;
    propW = propagator(MW,WW);
    prope = propagator(me,0);
    propne = propagator(0,0);
    p1=particle(me);
    p2=particle(0);
    p3=particle(MZ);
    p4=particle(MW);    
    //<34>,[34],<23>,[14],<24>,[13],<12>,[132>
    a34a=sproduct(ANGLE,&p3,&p4);
    s34s=sproduct(SQUARE,&p3,&p4);
    a23a=sproduct(ANGLE,&p2,&p3);
    s14s=sproduct(SQUARE,&p1,&p4);
    a24a=sproduct(ANGLE,&p2,&p4);
    s13s=sproduct(SQUARE,&p1,&p3);
    a12a=sproduct(ANGLE,&p1,&p2);
    s132a=sproduct(SQUARE,&p1,&p3,&p2);
    //<13>,[413>
    a13a=sproduct(ANGLE,&p1,&p3);
    s413a=sproduct(SQUARE,&p4,&p1,&p3);
    //[314>
    s314a=sproduct(SQUARE,&p3,&p1,&p4);
    //prefactor
    preW = e*e/(std::sqrt(2.0)*MW*MW*MZ*MZ*SW*SW);
    pree = e*e/(std::sqrt(2.0)*MW*MW*SW*SW);
    pren = e*e/(std::sqrt(2.0)*MW*MW*SW*SW);
    
    gL=two*SW*SW-one;
    gR=two*SW*SW;
  }
  void enZW::set_masses(const ldouble& masse, const ldouble& massW){
    constexpr ldouble two=2;//, sqrt2=std::sqrt(2);
    me=masse;
    MW=massW;
    MZ=MW/CW;
    p1.set_mass(me);
    p3.set_mass(MZ);
    p4.set_mass(MW);
    propW.set_mass(MW);
    prope.set_mass(me);
    //prefactors
    preW = e*e/(std::sqrt(2.0)*MW*MW*MZ*MZ*SW*SW);
    pree = e*e/(std::sqrt(2.0)*MW*MW*SW*SW);
    pren = e*e/(std::sqrt(2.0)*MW*MW*SW*SW);
  }
  void enZW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);   
    //<34>,[34],<23>,[14],<24>,[13],<12>,[132>
    a34a.update();
    s34s.update();
    a23a.update();
    s14s.update();
    a24a.update();
    s13s.update();
    a12a.update();
    s132a.update();
    //<13>,[413>
    a13a.update();
    s413a.update();
    //[314>
    s314a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propW.denominator(propSP);
    pDenT=prope.denominator(propTP);
    pDenU=propne.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble enZW::amp(const int& ds1, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds4a, ds4b;
    constexpr ldouble two=2;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);

            
      //S W Diagram
      //preW = e*e/(sqrt(2)*MW*MW*MZ*MZ*SW*SW);//=pree/MZ^2
      //EnZW- all ingoing:
      // - (+2(MZ<34>+MW[34])MW^2 <23>[14]+2(MW<34>+MZ[34])MW^2 <24>[13]+((2MW^2-MZ^2 )Me<12>+2MW^2 [132>)<34>[34])/(s-MW^2)
      //34 out:
      // + (+2(MZ<34>+MW[34])MW^2 <23>[14]+2(MW<34>+MZ[34])MW^2 <24>[13]+((MZ^2-2MW^2)Me<12>+2MW^2 [132>)<34>[34])/(s-MW^2)
      amplitude += normFactor*preW*(
				    +two*(MZ*a34a.v(ds3b,ds4a)+MW*s34s.v(ds3b,ds4a))*MW*MW*a23a.v(ds3a)*s14s.v(ds1,ds4b)
				    +two*(MW*a34a.v(ds3a,ds4b)+MZ*s34s.v(ds3a,ds4b))*MW*MW*a24a.v(ds4a)*s13s.v(ds1,ds3b)
				    +((MZ*MZ-two*MW*MW)*me*a12a.v(ds1)+two*MW*MW*s132a.v(ds1))*a34a.v(ds3a,ds4a)*s34s.v(ds3b,ds4b)
				    )/pDenS;
      

      //T e Diagram
      //pree = e*e/(sqrt(2)*MW*MW*SW*SW);
      //EnZW- all ingoing:
      // + <24>(+gReMe<13>[34]+(MZ[34]+[413>)gLe[13]))/(t-Me^2)
      //34 out:
      // + <24>(+gReMe<13>[34]+(-MZ[34]+[413>)gLe[13]))/(t-Me^2)
      amplitude += normFactor*pree*a24a.v(ds4a)*(
						 +gR*me*a13a.v(ds1,ds3b)*s34s.v(ds3a,ds4b)
						 +gL*(-MZ*s34s.v(ds3b,ds4b)+s413a.v(ds4b,ds3b))*s13s.v(ds1,ds3a)
						 )/pDenT;
      


      //U ne Diagram
      //pren = e*e/(sqrt(2)*MW*MW*SW*SW);//=pree
      //EnZW- all in:
      // - <23>[14](MW[34]-[314>)/u
      //34 out:
      // + <23>[14](MW[34]+[314>)/u
      amplitude += normFactor*pren*a23a.v(ds3a)*s14s.v(ds1,ds4a)*(MW*s34s.v(ds3b,ds4b)+s314a.v(ds3b,ds4b))/pDenU;

      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble enZW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-2;j4<=2;j4+=2){
	  M = amp(j1,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2
    return amp2/2.0;
  }



  



  //  Tests
  int test_enZW(){
    int n=0;//Number of fails
    std::cout<<"\t* E , ne -> Z , W+      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,me=0.0005,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;
      enZW enZWAmp = enZW(EE,me,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {6.665289742065860E-01,1.708495894679631E-01,7.575245572222450E-02,4.065577368194025E-02,2.527522950577811E-02,1.852623905446649E-02,1.622248729348434E-02,1.653989005763949E-02,1.868687004952238E-02,2.243858975621485E-02,2.797941399456004E-02,3.589880258845217E-02,4.731607663106926E-02,6.419327502696730E-02,9.003205973163220E-02,1.315224282510682E-01,2.030104363871889E-01,3.415497383358215E-01,6.807129764303653E-01,2.347588168703631E+00};
      i += enZWAmp.test_2to2_amp2([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH);
      i += enZWAmp.test_2to2_amp2_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH);
      i += enZWAmp.test_2to2_amp2_boosts([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH);
      i += enZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=90\n";
      pspatial = 90;
      ldouble dataCH2[20] = {6.486596550006811E-03,9.818951887122500E-03,1.340164496552791E-02,1.751498830505822E-02,2.231351126732690E-02,2.789816939491223E-02,3.435165608404032E-02,4.175703869619048E-02,5.020870337532984E-02,5.981986298091222E-02,7.072878099114051E-02,8.310488510598175E-02,9.715545243837975E-02,1.131332341334962E-01,1.313450341000228E-01,1.521605851297870E-01,1.760195186209685E-01,2.034304797929873E-01,2.349471163159359E-01,2.710815470691769E-01};
      i += enZWAmp.test_2to2_amp2([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH2);
      i += enZWAmp.test_2to2_amp2_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH2);
      i += enZWAmp.test_2to2_amp2_boosts([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH2);
      i += enZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH2);
      //std::cout<<"\n# me=80, MW=80.385, pspatial=250\n";
      me=80;
      enZWAmp.set_masses(me,MW);
      pspatial=250;
      ldouble dataCH3[20] = {1.040563890770294E+00,4.695503335656731E-01,2.818650893513940E-01,1.939731199748143E-01,1.464242010359667E-01,1.191225712648560E-01,1.035115231924938E-01,9.547331491414797E-02,9.295520434867405E-02,9.504954949435925E-02,1.016338022767402E-01,1.132915461218871E-01,1.314518553412818E-01,1.588203915206243E-01,2.003838923209407E-01,2.658248570343514E-01,3.761285553744222E-01,5.860872926630450E-01,1.096815486586743E+00,3.716068106672870E+00};
      i += enZWAmp.test_2to2_amp2([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH3);
      i += enZWAmp.test_2to2_amp2_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH3);
      i += enZWAmp.test_2to2_amp2_boosts([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH3);
      i += enZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH3);
      //std::cout<<"\n# me=80, MW=80.385, pspatial=70\n";
      pspatial = 70;
      ldouble dataCH4[20] = {2.653260510852911E-01,2.789733604088587E-01,2.942352618421913E-01,3.112737759133471E-01,3.302859808647262E-01,3.515113063600091E-01,3.752410998707403E-01,4.018312562795899E-01,4.317190482279561E-01,4.654458190063651E-01,5.036880090530165E-01,5.473002654657079E-01,5.973764549912751E-01,6.553378520097083E-01,7.230637080364070E-01,8.030899879254230E-01,8.989217020106980E-01,1.015542513225669E+00,1.160284036746834E+00,1.344390442874699E+00};
      i += enZWAmp.test_2to2_amp2([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH4);
      i += enZWAmp.test_2to2_amp2_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH4);
      i += enZWAmp.test_2to2_amp2_boosts([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH4);
      i += enZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH4);
      //std::cout<<"\n# me=80, MW=1, pspatial=250\n";
      MW=1;
      MZ=MW/CW;
      enZWAmp.set_masses(me,MW);
      pspatial=250;
      ldouble dataCH5[20] = {9.539988765154548E+06,4.524240922650066E+06,2.852888658843392E+06,2.017329592093593E+06,1.516039872453747E+06,1.181872411626792E+06,9.431997280104372E+05,7.642106493389856E+05,6.250114274746421E+05,5.136668436477780E+05,4.225828407924721E+05,3.466981118411907E+05,2.825107250187083E+05,2.275226070799705E+05,1.799076888075536E+05,1.383079696524892E+05,1.017113589881972E+05,6.940539966802474E+04,4.113436818104974E+04,1.989913255141225E+04};
      i += enZWAmp.test_2to2_amp2([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH5);
      i += enZWAmp.test_2to2_amp2_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH5);
      i += enZWAmp.test_2to2_amp2_boosts([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH5);
      i += enZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH5);
      //std::cout<<"\n# me=80, MW=1, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH6[20] = {1.995100953908726E+04,1.975471622952148E+04,1.961072204763000E+04,1.952841805087461E+04,1.951967006245201E+04,1.959969063993257E+04,1.978830762577570E+04,2.011185693638084E+04,2.060608991494881E+04,2.132079179856836E+04,2.232741366352662E+04,2.373229088374542E+04,2.570087962891479E+04,2.850544798540224E+04,3.262777202051880E+04,3.900836518683051E+04,4.976172237848568E+04,7.082945016090829E+04,1.276359455918856E+05,6.500626554151255E+05};
      i += enZWAmp.test_2to2_amp2([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH6);
      i += enZWAmp.test_2to2_amp2_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH6);
      i += enZWAmp.test_2to2_amp2_boosts([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH6);
      i += enZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enZWAmp.amp2(); }, me,0,MZ,MW,pspatial,dataCH6);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }



}
