
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

//File:  SPINAS/SM/AhWW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AhWW.h"

namespace spinas {

  AhWW::AhWW(const ldouble& echarge, const ldouble& massh, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), mh(massh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WW(widthW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    propW = propagator(MW,WW);
    p1=particle(0);
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
    //[123>,[124>,[321>,[421>
    s123a = sproduct(SQUARE,&p1,&p2,&p3);
    s124a = sproduct(SQUARE,&p1,&p2,&p4);
    s321a = sproduct(SQUARE,&p3,&p2,&p1);
    s421a = sproduct(SQUARE,&p4,&p2,&p1);
    //Couplings
    pre = sqrt2*e*e/(MW*SW);
  }
  void AhWW::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    p2.set_mass(mh);
    p3.set_mass(MW);
    p4.set_mass(MW);
    propW.set_mass(MW);
    //Couplings
    pre = sqrt2*e*e/(MW*SW);
  }
  void AhWW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    //[123>,[124>,[321>,[421>
    s123a.update();
    s124a.update();
    s321a.update();
    s421a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT=propW.den(propTP);
    pDenU=propW.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AhWW::amp(const int& ds1, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1;
    int ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);

      //pre = sqrt2*e*e/(2.0*MW*SW);

      if(ds1>0){
      
	//T&U-Channel W
	//all ingoing:
	//+pre <34>( Mh^2[13][14] - MW([14][123>+[13][124>) )/(t-MW^2)(u-MW^2)
	//34 outgoing:
	//+pre <34>( Mh^2[13][14] + MW([14][123>+[13][124>) )/(t-MW^2)(u-MW^2)
	amplitude += normFactor*pre*a34a.v(ds3a,ds4a)*( mh*mh*s13s.v(ds3b)*s14s.v(ds4b)
							+MW*(s14s.v(ds4b)*s123a.v(ds3b)+s13s.v(ds3b)*s124a.v(ds4b))
							)/pDenT/pDenU;
      }
      else if(ds1<0){
	amplitude += normFactor*pre*s34s.v(ds3a,ds4a)*( mh*mh*a13a.v(ds3b)*a14a.v(ds4b)
							+MW*(a14a.v(ds4b)*s321a.v(ds3b)+a13a.v(ds3b)*s421a.v(ds4b))
							)/pDenT/pDenU;
      }
      else{
	std::cout<<"Photon cannot have helicity 0!\n";
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AhWW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-2;j4<=2;j4+=2){
	  M = amp(j1,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2
    return amp2/2.0;
  }

  //set_momenta(...) must be called before amp2().
  ldouble AhWW::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=2)
      for(int j4=-2;j4<=2;j4+=2){
	M = amp(2,j3,j4);
	amp2 += std::pow(std::abs(M),2);
      }
    return amp2;
  }


  



  //  Tests
  int test_AhWW(){
    int n=0;//Number of fails
    std::cout<<"\t* A , h  -> W+, W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      AhWW AhWWAmp = AhWW(EE,mh,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.917147257364860E+01,4.838048962905846E+00,2.195060642782982E+00,1.274529720739073E+00,8.531606582732509E-01,6.296229453892305E-01,5.010458815349438E-01,4.249504774531291E-01,3.817274190384819E-01,3.620319157112367E-01,3.620319157112366E-01,3.817274190384817E-01,4.249504774531288E-01,5.010458815349436E-01,6.296229453892301E-01,8.531606582732505E-01,1.274529720739072E+00,2.195060642782980E+00,4.838048962905840E+00,1.917147257364855E+01};
      i += AhWWAmp.test_2to2_amp2([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH);
      i += AhWWAmp.test_2to2_amp2_rotations([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH);
      i += AhWWAmp.test_2to2_amp2_boosts([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH);
      i += AhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH);
      i += AhWWAmp.test_2to2_amp2([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH);
      i += AhWWAmp.test_2to2_amp2_rotations([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH);
      i += AhWWAmp.test_2to2_amp2_boosts([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH);
      i += AhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {4.558315957642572E+00,2.298315239020693E+00,1.413895445830032E+00,9.807745716152709E-01,7.394684567168639E-01,5.942513728527474E-01,5.033170245778100E-01,4.463055674524968E-01,4.126838786345337E-01,3.970343099554570E-01,3.970343099554570E-01,4.126838786345335E-01,4.463055674524967E-01,5.033170245778099E-01,5.942513728527474E-01,7.394684567168638E-01,9.807745716152705E-01,1.413895445830031E+00,2.298315239020692E+00,4.558315957642567E+00};
      i += AhWWAmp.test_2to2_amp2([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH2);
      i += AhWWAmp.test_2to2_amp2_rotations([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH2);
      i += AhWWAmp.test_2to2_amp2_boosts([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH2);
      i += AhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH2);
      i += AhWWAmp.test_2to2_amp2([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH2);
      i += AhWWAmp.test_2to2_amp2_rotations([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH2);
      i += AhWWAmp.test_2to2_amp2_boosts([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH2);
      i += AhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=40\n";
      pspatial = 40;
      ldouble dataCH4[20] = {4.867219172723948E-01,4.735798868258109E-01,4.624696816325183E-01,4.531600731633016E-01,4.454687244273943E-01,4.392526646175126E-01,4.344013179385315E-01,4.308314223992898E-01,4.284833876176196E-01,4.273187873681147E-01,4.273187873681147E-01,4.284833876176197E-01,4.308314223992898E-01,4.344013179385315E-01,4.392526646175126E-01,4.454687244273944E-01,4.531600731633016E-01,4.624696816325183E-01,4.735798868258108E-01,4.867219172723948E-01};
      i += AhWWAmp.test_2to2_amp2([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH4);
      i += AhWWAmp.test_2to2_amp2_rotations([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH4);
      i += AhWWAmp.test_2to2_amp2_boosts([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH4);
      i += AhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH4);
      i += AhWWAmp.test_2to2_amp2([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH4);
      i += AhWWAmp.test_2to2_amp2_rotations([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH4);
      i += AhWWAmp.test_2to2_amp2_boosts([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH4);
      i += AhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH4);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=95\n";
      mh = 0.0005;
      pspatial = 95;
      AhWWAmp.set_masses(mh,MW);
      ldouble dataCH3[20] = {7.888370840996252E-01,6.777616340680471E-01,5.962053585293215E-01,5.353852703418480E-01,4.897237993703849E-01,4.555722771826334E-01,4.305005710721793E-01,4.128860250717157E-01,4.016697914364733E-01,3.962123791066582E-01,3.962123791066582E-01,4.016697914364732E-01,4.128860250717156E-01,4.305005710721793E-01,4.555722771826334E-01,4.897237993703849E-01,5.353852703418480E-01,5.962053585293215E-01,6.777616340680470E-01,7.888370840996251E-01};
      i += AhWWAmp.test_2to2_amp2([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH3);
      i += AhWWAmp.test_2to2_amp2_rotations([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH3);
      i += AhWWAmp.test_2to2_amp2_boosts([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH3);
      i += AhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AhWWAmp.amp2(); }, 0,mh,MW,MW,pspatial,dataCH3);
      i += AhWWAmp.test_2to2_amp2([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH3);
      i += AhWWAmp.test_2to2_amp2_rotations([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH3);
      i += AhWWAmp.test_2to2_amp2_boosts([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH3);
      i += AhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AhWWAmp.amp2_Aplus(); }, 0,mh,MW,MW,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }
  
  

}
