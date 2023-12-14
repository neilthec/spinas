
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

//File:  SPINAS/SM/AZWW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AZWW.h"

namespace spinas {

  AZWW::AZWW(const ldouble& echarge, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WW(widthW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    MZ = MW/CW;
    propW = propagator(MW,WW);
    p1=particle(0);
    p2=particle(MZ);
    p3=particle(MW);
    p4=particle(MW);
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //Couplings
    pre = 2.0*e*e/(CW*SW);
  }
  void AZWW::set_masses(const ldouble& massW){
    MW=massW;
    MZ=MW/CW;
    p2.set_mass(MZ);
    p3.set_mass(MW);
    p4.set_mass(MW);
    propW.set_mass(MW);
    //Couplings
    pre = 2.0*e*e/(CW*SW);
  }
  void AZWW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s24s.update();
    a24a.update();
    s34s.update();
    a34a.update();
    s14s.update();
    a14a.update();
    s13s.update();
    a13a.update();
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
  cdouble AZWW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds2a, ds2b, ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds2,ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds2,ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);

      if(ds1>0){
	
	//pre = 2*e*e/CW/SW;
	//all ingoing:
	//pre*( CW^2[12]^2<34>^2 + CW^2[13]^2<24>^2 + CW^2[14]^2<23>^2 + (CW^2-SW^2)[13][14]<23><24> + CW[12][13]<24><34> - CW[12][14]<23><34> )/(t-MW^2)(u-MW^2)
	//34 outgoing:
	//pre*( CW^2[12]^2<34>^2 + CW^2[13]^2<24>^2 + CW^2[14]^2<23>^2 + (CW^2-SW^2)[13][14]<23><24> - CW[12][13]<24><34> + CW[12][14]<23><34> )/(t-MW^2)(u-MW^2)
	amplitude += normFactor*pre*(
				     + CW*CW*s12s.v(ds2a)*s12s.v(ds2b)*a34a.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)
				     + CW*CW*s13s.v(ds3a)*s13s.v(ds3b)*a24a.v(ds2a,ds4a)*a24a.v(ds2b,ds4b)
				     + CW*CW*s14s.v(ds4a)*s14s.v(ds4b)*a23a.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)
				     + (CW*CW-SW*SW)*s13s.v(ds3a)*s14s.v(ds4a)*a23a.v(ds2a,ds3b)*a24a.v(ds2b,ds4b)
				     - CW*s12s.v(ds2a)*s13s.v(ds3a)*a24a.v(ds2b,ds4a)*a34a.v(ds3b,ds4b)
				     + CW*s12s.v(ds2a)*s14s.v(ds4a)*a23a.v(ds2b,ds3a)*a34a.v(ds3b,ds4b)
				     )/pDenT/pDenU;	
	
      }
      else if(ds1<0){
	amplitude += normFactor*pre*(
				     + CW*CW*a12a.v(ds2a)*a12a.v(ds2b)*s34s.v(ds3a,ds4a)*s34s.v(ds3b,ds4b)
				     + CW*CW*a13a.v(ds3a)*a13a.v(ds3b)*s24s.v(ds2a,ds4a)*s24s.v(ds2b,ds4b)
				     + CW*CW*a14a.v(ds4a)*a14a.v(ds4b)*s23s.v(ds2a,ds3a)*s23s.v(ds2b,ds3b)
				     + (CW*CW-SW*SW)*a13a.v(ds3a)*a14a.v(ds4a)*s23s.v(ds2a,ds3b)*s24s.v(ds2b,ds4b)
				     - CW*a12a.v(ds2a)*a13a.v(ds3a)*s24s.v(ds2b,ds4a)*s34s.v(ds3b,ds4b)
				     + CW*a12a.v(ds2a)*a14a.v(ds4a)*s23s.v(ds2b,ds3a)*s34s.v(ds3b,ds4b)
				     )/pDenT/pDenU;
      }
      else{
	std::cout<<"Photon cannot have helicity 0!\n";
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AZWW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2*1/3=1/6
    return amp2/6.0;
  }

  //set_momenta(...) must be called before amp2().
  ldouble AZWW::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-2;j4<=2;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }




  



  //  Tests
  int test_AZWW(){
    int n=0;//Number of fails
    std::cout<<"\t* A , Z  -> W+, W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble MZ = MW/std::sqrt(1.0-SW*SW);
      AZWW AZWWAmp = AZWW(EE,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.360252016669019E+01,1.123582767021352E+01,5.125595111507945E+00,2.986477572200235E+00,2.006859385480133E+00,1.487805874997211E+00,1.189782348111386E+00,1.013718001096369E+00,9.138575068748113E-01,8.683970563943164E-01,8.683970563943210E-01,9.138575068748105E-01,1.013718001096368E+00,1.189782348111394E+00,1.487805874997211E+00,2.006859385480138E+00,2.986477572200245E+00,5.125595111507953E+00,1.123582767021352E+01,4.360252016669016E+01};
      i += AZWWAmp.test_2to2_amp2([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH);
      i += AZWWAmp.test_2to2_amp2_rotations([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH);
      i += AZWWAmp.test_2to2_amp2_boosts([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH);
      i += AZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH);
      i += AZWWAmp.test_2to2_amp2([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH);
      i += AZWWAmp.test_2to2_amp2_rotations([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH);
      i += AZWWAmp.test_2to2_amp2_boosts([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH);
      i += AZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {8.444387738837968E+00,4.514438726337164E+00,2.829958839150267E+00,1.966152944442793E+00,1.472636959092429E+00,1.171438387147334E+00,9.813610993171955E-01,8.617007658613701E-01,7.909888606006791E-01,7.580458714882805E-01,7.580458714882806E-01,7.909888606006785E-01,8.617007658613703E-01,9.813610993171947E-01,1.171438387147334E+00,1.472636959092428E+00,1.966152944442792E+00,2.829958839150266E+00,4.514438726337163E+00,8.444387738837962E+00};
      i += AZWWAmp.test_2to2_amp2([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH2);
      i += AZWWAmp.test_2to2_amp2_rotations([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH2);
      i += AZWWAmp.test_2to2_amp2_boosts([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH2);
      i += AZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH2);
      i += AZWWAmp.test_2to2_amp2([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH2);
      i += AZWWAmp.test_2to2_amp2_rotations([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH2);
      i += AZWWAmp.test_2to2_amp2_boosts([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH2);
      i += AZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=40\n";
      pspatial = 60;
      ldouble dataCH4[20] = {1.160311729288126E+00,1.081972780736168E+00,1.016570538031790E+00,9.623850390956039E-01,9.180641309904688E-01,8.825511459822193E-01,8.550315842595916E-01,8.348940401898753E-01,8.217020841633094E-01,8.151748518694183E-01,8.151748518694184E-01,8.217020841633094E-01,8.348940401898753E-01,8.550315842595917E-01,8.825511459822193E-01,9.180641309904689E-01,9.623850390956040E-01,1.016570538031790E+00,1.081972780736167E+00,1.160311729288126E+00};
      i += AZWWAmp.test_2to2_amp2([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH4);
      i += AZWWAmp.test_2to2_amp2_rotations([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH4);
      i += AZWWAmp.test_2to2_amp2_boosts([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH4);
      i += AZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZWWAmp.amp2(); }, 0,MZ,MW,MW,pspatial,dataCH4);
      i += AZWWAmp.test_2to2_amp2([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH4);
      i += AZWWAmp.test_2to2_amp2_rotations([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH4);
      i += AZWWAmp.test_2to2_amp2_boosts([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH4);
      i += AZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZWWAmp.amp2_Aplus(); }, 0,MZ,MW,MW,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

}
