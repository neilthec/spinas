
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2024 Neil Christensen

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

//File:  SPINAS/SM/ZnWe.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ZnWe.h"

namespace spinas {

  ZnWe::ZnWe(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), me(masse), MW(massW), SW(sinW), WW(widthW) {
    constexpr ldouble one=1, two=2;//, sqrt2=std::sqrt(2);
    CW=std::sqrt(1-SW*SW);
    MZ=MW/CW;
    propW = propagator(MW,WW);
    prope = propagator(me,0);
    propne = propagator(0,0);
    p1=particle(MZ);
    p2=particle(0);
    p3=particle(MW);
    p4=particle(me);    
    //<13>,[13],<21>,[43],<23>,[41],<42>,[412>
    a13a=sproduct(ANGLE,&p1,&p3);
    s13s=sproduct(SQUARE,&p1,&p3);
    a21a=sproduct(ANGLE,&p2,&p1);
    s43s=sproduct(SQUARE,&p4,&p3);
    a23a=sproduct(ANGLE,&p2,&p3);
    s41s=sproduct(SQUARE,&p4,&p1);
    a42a=sproduct(ANGLE,&p4,&p2);
    s412a=sproduct(SQUARE,&p4,&p1,&p2);
    //<41>,[341>
    a41a=sproduct(ANGLE,&p4,&p1);
    s341a=sproduct(SQUARE,&p3,&p4,&p1);
    //[143>
    s143a=sproduct(SQUARE,&p1,&p4,&p3);
    //prefactor
    preW = e*e/(std::sqrt(2.0)*MW*MW*MZ*MZ*SW*SW);
    pree = e*e/(std::sqrt(2.0)*MW*MW*SW*SW);
    pren = e*e/(std::sqrt(2.0)*MW*MW*SW*SW);
    
    gL=two*SW*SW-one;
    gR=two*SW*SW;
  }
  void ZnWe::set_masses(const ldouble& masse, const ldouble& massW){
    constexpr ldouble two=2;//, sqrt2=std::sqrt(2);
    me=masse;
    MW=massW;
    MZ=MW/CW;
    p1.set_mass(MZ);
    p3.set_mass(MW);
    p4.set_mass(me);
    propW.set_mass(MW);
    prope.set_mass(me);
    //prefactors
    preW = e*e/(std::sqrt(2.0)*MW*MW*MZ*MZ*SW*SW);
    pree = e*e/(std::sqrt(2.0)*MW*MW*SW*SW);
    pren = e*e/(std::sqrt(2.0)*MW*MW*SW*SW);
  }
  void ZnWe::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);   
    //<13>,[13],<21>,[43],<23>,[41],<42>,[412>
    a13a.update();
    s13s.update();
    a21a.update();
    s43s.update();
    a23a.update();
    s41s.update();
    a42a.update();
    s412a.update();
    //<41>,[341>
    a41a.update();
    s341a.update();
    //[143>
    s143a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT=propW.denominator(propTP);
    pDenU=prope.denominator(propUP);
    pDenS=propne.denominator(propSP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ZnWe::amp(const int& ds1, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds1a, ds1b, ds3a, ds3b;
    constexpr ldouble two=2;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds1,ds3);
    ldouble normFactor=get_spin_normalization(ds1,ds3);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds1,ds1a,ds1b, ds3,ds3a,ds3b, i);

            
      //T W Diagram
      //preW = e*e/(sqrt(2)*MW*MW*MZ*MZ*SW*SW);//=pree/MZ^2
      //EnZW- all ingoing:
      // - (+2(+MZ<34>+MW[34])MW^2 <23>[14]+2(MW<34>+MZ[34])MW^2 <24>[13]+((2MW^2-MZ^2 )Me<12>+2MW^2 [132>)<34>[34])/(s-MW^2)
      //ZnW-E: 1->4->3->1
      // - (+2(+MZ<13>+MW[13])MW^2 <21>[43]+2(MW<13>+MZ[13])MW^2 <23>[41]+((2MW^2-MZ^2 )Me<42>+2MW^2 [412>)<13>[13])/(t-MW^2)
      //34 out:
      // - (+2(-MZ<13>+MW[13])MW^2 <21>[43]+2(MW<13>-MZ[13])MW^2 <23>[41]+((2MW^2-MZ^2 )Me<42>-2MW^2 [412>)<13>[13])/(t-MW^2)
      amplitude += - normFactor*preW*(
				      +two*(-MZ*a13a.v(ds1b,ds3a)+MW*s13s.v(ds1b,ds3a))*MW*MW*a21a.v(ds1a)*s43s.v(ds4,ds3b)
				      +two*(MW*a13a.v(ds1a,ds3b)-MZ*s13s.v(ds1a,ds3b))*MW*MW*a23a.v(ds3a)*s41s.v(ds4,ds1b)
				      +((-MZ*MZ+two*MW*MW)*me*a42a.v(ds4)-two*MW*MW*s412a.v(ds4))*a13a.v(ds1a,ds3a)*s13s.v(ds1b,ds3b)
				      )/pDenT;
      

      //U e Diagram
      //pree = e*e/(sqrt(2)*MW*MW*SW*SW);
      //EnZW- all ingoing:
      // + <24>(+gReMe<13>[34]+(MZ[34]+[413>)gLe[13]))/(t-Me^2)
      //ZnW-E: 1->4->3->1
      // + <23>(+gReMe<41>[13]+(MZ[13]+[341>)gLe[41]))/(u-Me^2)
      //34 out:
      // - <23>(-gReMe<41>[13]+(MZ[13]-[341>)gLe[41]))/(u-Me^2)
      amplitude += - normFactor*pree*a23a.v(ds3a)*(
						   -gR*me*a41a.v(ds4,ds1b)*s13s.v(ds1a,ds3b)
						   +gL*(MZ*s13s.v(ds1b,ds3b)-s341a.v(ds3b,ds1b))*s41s.v(ds4,ds1a)
						   )/pDenU;
      


      //S ne Diagram
      //pren = e*e/(sqrt(2)*MW*MW*SW*SW);//=pree
      //EnZW- all in:
      // - <23>[14](MW[34]-[314>)/u
      //ZnW-E: 1->4->3->1
      // - <21>[43](MW[13]-[143>)/s
      //34 out:
      // - <21>[43](MW[13]-[143>)/s
      amplitude += - normFactor*pren*a21a.v(ds1a)*s43s.v(ds4,ds3a)*(MW*s13s.v(ds1b,ds3b)-s143a.v(ds1b,ds3b))/pDenS;

      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ZnWe::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }



  



  //  Tests
  int test_ZnWe(){
    int n=0;//Number of fails
    std::cout<<"\t* Z , ne -> W+, e       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,me=0.0005,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;
      ZnWe ZnWeAmp = ZnWe(EE,me,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {9.111264562661387E+01,2.169880593431707E+01,9.109177446054824E+00,4.831281528225641E+00,2.914912778884934E+00,1.909447773461601E+00,1.324921724132261E+00,9.598325340990970E-01,7.194957375858165E-01,5.549148543507745E-01,4.388281427625484E-01,3.552099892331410E-01,2.942657721796780E-01,2.499055614945518E-01,2.184642130297781E-01,1.981964211845729E-01,1.896129757344849E-01,1.979787629922848E-01,2.473384247380807E-01,5.771534519162178E-01};
      i += ZnWeAmp.test_2to2_amp2([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH);
      i += ZnWeAmp.test_2to2_amp2_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH);
      i += ZnWeAmp.test_2to2_amp2_boosts([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH);
      i += ZnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=90\n";
      pspatial = 90;
      ldouble dataCH2[20] = {6.671960998878340E+00,4.151577475540546E+00,2.776606252948298E+00,1.954080841244484E+00,1.428558848187979E+00,1.075871123701323E+00,8.300393946808109E-01,6.535270445022771E-01,5.238002044156221E-01,4.267453343550342E-01,3.532345324633259E-01,2.972495403392433E-01,2.548420553130143E-01,2.236090409002334E-01,2.026169074993097E-01,1.930892313326462E-01,2.016099029347257E-01,2.574732349631615E-01,6.186339468124504E-01,7.466111007204932E+00};
      i += ZnWeAmp.test_2to2_amp2([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH2);
      i += ZnWeAmp.test_2to2_amp2_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH2);
      i += ZnWeAmp.test_2to2_amp2_boosts([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH2);
      i += ZnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH2);
      //std::cout<<"\n# me=80, MW=80.385, pspatial=250\n";
      me=80;
      ZnWeAmp.set_masses(me,MW);
      pspatial=250;
      ldouble dataCH3[20] = {9.171900651466250E+01,2.226565615473914E+01,9.481565153058591E+00,5.102577399027488E+00,3.128725355195438E+00,2.087621720664467E+00,1.479683582550997E+00,1.098723187348498E+00,8.475907630671730E-01,6.759380094113949E-01,5.557763766423711E-01,4.707580123778137E-01,4.110875376351287E-01,3.710085203689428E-01,3.476450363246049E-01,3.408101680387377E-01,3.540945054719988E-01,3.992586637713471E-01,5.145007259605117E-01,8.884738256545744E-01};
      i += ZnWeAmp.test_2to2_amp2([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH3);
      i += ZnWeAmp.test_2to2_amp2_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH3);
      i += ZnWeAmp.test_2to2_amp2_boosts([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH3);
      i += ZnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH3);
      //std::cout<<"\n# me=80, MW=80.385, pspatial=70\n";
      pspatial = 70;
      ldouble dataCH4[20] = {2.153679019464761E+00,1.727496073068624E+00,1.406882845922504E+00,1.160883252249998E+00,9.689908032253364E-01,8.172121475765428E-01,6.957531817140463E-01,5.976076887188275E-01,5.176685744437661E-01,4.521524840061849E-01,3.982182413487275E-01,3.537085497339562E-01,3.169721281010859E-01,2.867396772163485E-01,2.620368561985143E-01,2.421235625502474E-01,2.264528082285404E-01,2.146453119397642E-01,2.064782272936087E-01,2.018887616654572E-01};
      i += ZnWeAmp.test_2to2_amp2([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH4);
      i += ZnWeAmp.test_2to2_amp2_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH4);
      i += ZnWeAmp.test_2to2_amp2_boosts([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH4);
      i += ZnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH4);
      //std::cout<<"\n# me=80, MW=1, pspatial=250\n";
      MW=1;
      MZ=MW/CW;
      ZnWeAmp.set_masses(me,MW);
      pspatial=250;
      ldouble dataCH5[20] = {3.430188426046772E+05,3.545760832214800E+05,3.729210205954698E+05,3.942250063333533E+05,4.184703477876615E+05,4.460797510579286E+05,4.777053277372484E+05,5.142382947423547E+05,5.568835768573765E+05,6.072921396694937E+05,6.677787756305943E+05,7.416873826044572E+05,8.340306718051287E+05,9.526769472493868E+05,1.110719099240985E+06,1.331672297510883E+06,1.662419186540622E+06,2.211828049865389E+06,3.303771235797895E+06,6.525489011887079E+06};
      i += ZnWeAmp.test_2to2_amp2([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH5);
      i += ZnWeAmp.test_2to2_amp2_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH5);
      i += ZnWeAmp.test_2to2_amp2_boosts([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH5);
      i += ZnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH5);
      //std::cout<<"\n# me=80, MW=1, pspatial=50\n";
      pspatial = 50;
      ldouble dataCH6[20] = {3.519174482859629E+05,3.427777381003698E+05,3.460143732936716E+05,3.512483206991062E+05,3.573256712517936E+05,3.639422121611105E+05,3.709914612855216E+05,3.784339267738197E+05,3.862586917651230E+05,3.944693566515283E+05,4.030781263380275E+05,4.121031044491637E+05,4.215670180343098E+05,4.314966517830988E+05,4.419226720683584E+05,4.528796902982329E+05,4.644064929049826E+05,4.765464042147471E+05,4.893477697440190E+05,5.028645605297663E+05};
      i += ZnWeAmp.test_2to2_amp2([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH6);
      i += ZnWeAmp.test_2to2_amp2_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH6);
      i += ZnWeAmp.test_2to2_amp2_boosts([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH6);
      i += ZnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZnWeAmp.amp2(); }, MZ,0,MW,me,pspatial,dataCH6);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }



}
