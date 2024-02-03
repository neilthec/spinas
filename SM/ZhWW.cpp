
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

//File:  SPINAS/SM/ZhWW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ZhWW.h"

namespace spinas {

  ZhWW::ZhWW(const ldouble& echarge, const ldouble& massh, const ldouble& massW, const ldouble& sinW, const ldouble& widthW, const ldouble& widthZ):
    e(echarge), mh(massh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WW(widthW), WZ(widthZ) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    MZ = MW/CW;
    propW = propagator(MW,WW);
    propZ = propagator(MZ,WZ);
    p1=particle(MZ);
    p2=particle(mh);
    p3=particle(MW);
    p4=particle(MW);
    //<34>,[34],<14>,[14],<13>,[13]
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //[131>,[141>,[431>,[424>,[341>,[323>,[123>,[124>,[321>,[421>
    s131a = sproduct(SQUARE,&p1,&p3,&p1);
    s141a = sproduct(SQUARE,&p1,&p4,&p1);
    s431a = sproduct(SQUARE,&p4,&p3,&p1);
    s424a = sproduct(SQUARE,&p4,&p2,&p4);
    s341a = sproduct(SQUARE,&p3,&p4,&p1);
    s323a = sproduct(SQUARE,&p3,&p2,&p3);
    s123a = sproduct(SQUARE,&p1,&p2,&p3);
    s124a = sproduct(SQUARE,&p1,&p2,&p4);
    s321a = sproduct(SQUARE,&p3,&p2,&p1);
    s421a = sproduct(SQUARE,&p4,&p2,&p1);
    //Couplings
    pre = sqrt2*e*e/(2.0*MW*MW*SW*SW);
    preS = pre;
    preTU = pre/(MZ*MZ);
  }
  void ZhWW::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    MZ=MW/CW;
    p1.set_mass(MZ);
    p2.set_mass(mh);
    p3.set_mass(MW);
    p4.set_mass(MW);
    propW.set_mass(MW);
    //Couplings
    pre = sqrt2*e*e/(2.0*MW*MW*SW*SW);
    preS = pre;
    preTU = pre/(MZ*MZ);
  }
  void ZhWW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<34>,[34],<14>,[14],<13>,[13]
    s34s.update();
    a34a.update();
    s14s.update();
    a14a.update();
    s13s.update();
    a13a.update();
    //[131>,[141>,[431>,[424>,[341>,[323>,[123>,[124>,[321>,[421>
    s131a.update();
    s141a.update();
    s431a.update();
    s424a.update();
    s341a.update();
    s323a.update();
    s123a.update();
    s124a.update();
    s321a.update();
    s421a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propZ.denominator(propSP);
    pDenT=propW.denominator(propTP);
    pDenU=propW.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ZhWW::amp(const int& ds1, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1;
    int ds1a, ds1b, ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the Z- & W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds1,ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds1,ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds1,ds1a,ds1b, ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);

      //pre = e^2 / 2 MW^2 SW^2

      //S-Channel Z
      //preS = pre;
      //All ingoing
      //- ( [34]<34>([131>-[141>) + 2MW([13]<14>+[14]<13>)([34]+<34>) )/(s-MZ^2)
      //34 outgoing:
      //  ( [34]<34>([131>-[141>) + 2MW([13]<14>+[14]<13>)([34]+<34>) )/(s-MZ^2)
      amplitude += normFactor*preS*( s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s131a.v(ds1a,ds1b)-s141a.v(ds1a,ds1b))
				     +2*MW*(s13s.v(ds1a,ds3a)*a14a.v(ds1b,ds4a)+s14s.v(ds1a,ds4a)*a13a.v(ds1b,ds3a))*(s34s.v(ds3b,ds4b)+a34a.v(ds3b,ds4b))
				     )/pDenS;
      
      //T-Channel W
      //preTU = pre / MZ^2
      //All ingoing:
      // - ( 2MW^2[13]<34>[431> - MZ^2[13]<13>[424> + 2MW^3[14]<13><34> + 2MW^2MZ[13][14]<34> + 2MW^2MZ[34]<13><14> )/(t-MW^2)
      //34 outgoing:
      //   ( 2MW^2[13]<34>[431> + MZ^2[13]<13>[424> + 2MW^3[14]<13><34> - 2MW^2MZ[13][14]<34> - 2MW^2MZ[34]<13><14> )/(t-MW^2)
      amplitude += normFactor*preTU*(
				     2.0*MW*MW*s13s.v(ds1a,ds3a)*a34a.v(ds3b,ds4a)*s431a.v(ds4b,ds1b)
				     + MZ*MZ*s13s.v(ds1a,ds3a)*a13a.v(ds1b,ds3b)*s424a.v(ds4a,ds4b)
				     + 2.0*MW*MW*MW*s14s.v(ds1a,ds4a)*a13a.v(ds1b,ds3a)*a34a.v(ds3b,ds4b)
				     - 2.0*MW*MW*MZ*s13s.v(ds1a,ds3a)*s14s.v(ds1b,ds4a)*a34a.v(ds3b,ds4b)
				     - 2.0*MW*MW*MZ*a13a.v(ds1a,ds3a)*a14a.v(ds1b,ds4a)*s34s.v(ds3b,ds4b)
				     )/pDenT;

      //U-Channel W
      //All ingoing:
      // - ( 2MW^2[14]<34>[341> + MZ^2[14]<14>[323> + 2MW^3[13]<14><34> + 2MW^2MZ[13][14]<34> + 2MW^2MZ[34]<13><14> )/(u-MW^2)
      //34 outgoing:
      //   ( 2MW^2[14]<34>[341> - MZ^2[14]<14>[323> + 2MW^3[13]<14><34> - 2MW^2MZ[13][14]<34> - 2MW^2MZ[34]<13><14> )/(u-MW^2)
      amplitude += normFactor*preTU*(
				     2.0*MW*MW*s14s.v(ds1a,ds4a)*a34a.v(ds3a,ds4b)*s341a.v(ds3b,ds1b)
				     - MZ*MZ*s14s.v(ds1a,ds4a)*a14a.v(ds1b,ds4b)*s323a.v(ds3a,ds3b)
				     + 2.0*MW*MW*MW*s13s.v(ds1a,ds3a)*a14a.v(ds1b,ds4a)*a34a.v(ds3b,ds4b)
				     - 2.0*MW*MW*MZ*s14s.v(ds1a,ds4a)*s13s.v(ds1b,ds3a)*a34a.v(ds3b,ds4b)
				     - 2.0*MW*MW*MZ*a14a.v(ds1a,ds4a)*a13a.v(ds1b,ds3a)*s34s.v(ds3b,ds4b)
				     )/pDenU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ZhWW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-2;j4<=2;j4+=2){
	  M = amp(j1,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }


  



  //  Tests
  int test_ZhWW(){
    int n=0;//Number of fails
    std::cout<<"\t* Z , h  -> W+, W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW = std::sqrt(1.0-SW*SW);
      ldouble MZ = MW/CW;
      ZhWW ZhWWAmp = ZhWW(EE,mh,MW,SW,0,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {5.131584406575832E+01,1.219745048923644E+01,5.148410656682649E+00,2.753316560440737E+00,1.682647581200057E+00,1.126779697335775E+00,8.129754279311715E-01,6.300640662813631E-01,5.273305654880298E-01,4.808398923126528E-01,4.808398923126523E-01,5.273305654880299E-01,6.300640662813662E-01,8.129754279311716E-01,1.126779697335776E+00,1.682647581200056E+00,2.753316560440739E+00,5.148410656682652E+00,1.219745048923645E+01,5.131584406575829E+01};
      i += ZhWWAmp.test_2to2_amp2([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH);
      i += ZhWWAmp.test_2to2_amp2_rotations([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH);
      i += ZhWWAmp.test_2to2_amp2_boosts([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH);
      i += ZhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {1.307275431513298E+01,6.303415538197669E+00,3.671469047979958E+00,2.396314550688745E+00,1.694146269486510E+00,1.276299463908552E+00,1.017235760888747E+00,8.561371228324864E-01,7.617051329926237E-01,7.179139073114253E-01,7.179139073114253E-01,7.617051329926229E-01,8.561371228324858E-01,1.017235760888747E+00,1.276299463908552E+00,1.694146269486510E+00,2.396314550688743E+00,3.671469047979956E+00,6.303415538197666E+00,1.307275431513297E+01};
      i += ZhWWAmp.test_2to2_amp2([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH2);
      i += ZhWWAmp.test_2to2_amp2_rotations([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH2);
      i += ZhWWAmp.test_2to2_amp2_boosts([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH2);
      i += ZhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {1.467636732984107E+00,1.467515105738395E+00,1.467407002030080E+00,1.467312418541106E+00,1.467231352368653E+00,1.467163801024945E+00,1.467109762437089E+00,1.467069234946935E+00,1.467042217310964E+00,1.467028708700213E+00,1.467028708700213E+00,1.467042217310965E+00,1.467069234946934E+00,1.467109762437089E+00,1.467163801024945E+00,1.467231352368653E+00,1.467312418541106E+00,1.467407002030080E+00,1.467515105738395E+00,1.467636732984107E+00};
      i += ZhWWAmp.test_2to2_amp2([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH4);
      i += ZhWWAmp.test_2to2_amp2_rotations([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH4);
      i += ZhWWAmp.test_2to2_amp2_boosts([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH4);
      i += ZhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH4);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=95\n";
      mh = 0.0005;
      pspatial = 95;
      ZhWWAmp.set_masses(mh,MW);
      ldouble dataCH3[20] = {3.475799289563991E+00,2.414233264039611E+00,1.791300616024454E+00,1.397765543584689E+00,1.137015918368530E+00,9.595790380132839E-01,8.382150028280996E-01,7.572207213385362E-01,7.074173884466769E-01,6.836716032381422E-01,6.836716032381422E-01,7.074173884466770E-01,7.572207213385360E-01,8.382150028280995E-01,9.595790380132839E-01,1.137015918368530E+00,1.397765543584689E+00,1.791300616024453E+00,2.414233264039610E+00,3.475799289563990E+00};
      i += ZhWWAmp.test_2to2_amp2([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH3);
      i += ZhWWAmp.test_2to2_amp2_rotations([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH3);
      i += ZhWWAmp.test_2to2_amp2_boosts([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH3);
      i += ZhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZhWWAmp.amp2(); }, MZ,mh,MW,MW,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }
  
  

}
