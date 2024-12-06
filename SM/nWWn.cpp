
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

//File:  SPINAS/SM/nWWn.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/nWWn.h"

namespace spinas {

  nWWn::nWWn(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), prope(0,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);
    p1=particle(0);
    p2=particle(MW);
    p3=particle(MW);
    p4=particle(0);
    //[23],<13>,<34>,[34],[24],<14>
    s34s = sproduct(SQUARE,&p3,&p4,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s24s = sproduct(SQUARE,&p2,&p4,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    //[314>,[231>
    s312a = sproduct(SQUARE,&p3,&p1,&p2,2);
    s431a = sproduct(SQUARE,&p4,&p3,&p1,2);
    //Couplings
    pree = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preZ = 2.0*e*e/(2.0*MW*MW*SW*SW);
  }
  void nWWn::set_masses(const ldouble& masse, const ldouble& massW){
    me=masse;
    prope.set_mass(me);
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p2.set_mass(MW);
    p3.set_mass(MW);
    propZ.set_mass(MZ);
    //Couplings
    pree = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preZ = 2.0*e*e/(2.0*MW*MW*SW*SW);
  }
  void nWWn::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //[23],<13>,<34>,[34],[24],<14>
    s34s.update();
    a13a.update();
    s23s.update();
    a23a.update();
    s24s.update();
    a12a.update();
    //[314>,[231>
    s312a.update();
    s431a.update();
    //Propagator Momentum
    ldouble propUP[4], propTP[4];
    for(int j=0;j<4;j++){
      propUP[j] = mom1[j]-mom4[j];
      propTP[j] = mom1[j]-mom3[j];
    }
    pDeneT=prope.denominator(propTP);
    pDenZU=propZ.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble nWWn::amp(const int& ds3, const int& ds2){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds2a, ds2b;
    constexpr ldouble two=2;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds2);
    ldouble normFactor=get_spin_normalization(ds3,ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds2,ds2a,ds2b, i);
           
      //T-Channel e
      //pree = e*e/(2.0*MW*MW*SW*SW);
      //nNW-W+ all ingoing:
      // - pree [24] <13> ([314>+MW <34>))/t
      //nW+W-N: 2<->4
      // + pree [24] <13> ([312>-MW <23>))/t
      //34 out:
      // - pree [24] <13> ([312>+MW <23>))/t
      amplitude += - normFactor*pree*s24s.v(ds2a)*a13a.v(ds3a)*(MW*a23a.v(ds2b,ds3b)+s312a.v(ds3b,ds2b))/pDeneT;

      
      //U-Channel Z
      //preZ = e*e/(2.0*MW*MW*SW*SW); //=pree
      //nNW-W+ all ingoing:
      //- preZ ( [34]<34>[231> + MW([23]<14>+[24]<13>)([34]+<34>) )/(s-MZ^2)
      //nW+W-N: 2<->4
      //- preZ ( [23]<23>[431> + MW([34]<12>+[24]<13>)([23]+<23>) )/(u-MZ^2)
      //34 out:
      //- preZ ( [23]<23>[431> + MW([34]<12>-[24]<13>)([23]-<23>) )/(u-MZ^2)
      amplitude += - normFactor*preZ*( 
				      s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*s431a.v()
				      +MW*(s23s.v(ds2a,ds3a)-a23a.v(ds2a,ds3a))*(s34s.v(ds3b)*a12a.v(ds2b)-s24s.v(ds2b)*a13a.v(ds3b))
				      )/pDenZU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble nWWn::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-2;j3<=2;j3+=2){
	M = amp(j2,j3);
	amp2 += std::pow(std::abs(M),2);
      }
    //Average over initial spins 1/3
    return amp2/3.0;
  }

  
  

  //  Tests
  int test_nWWn(){
    int n=0;//Number of fails
    std::cout<<"\t* ne, W+ -> W+, ne      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=250\n";
      ldouble me=0.0005;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      nWWn nWWnAmp = nWWn(EE,me,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.690256623190841E+00,9.684753518676835E-01,6.531784462490348E-01,5.325680131209674E-01,4.790228346756270E-01,4.590575326176777E-01,4.609523550146575E-01,4.806487925204465E-01,5.178801689751824E-01,5.751434440917650E-01,6.579095398719862E-01,7.758698453293275E-01,9.457465460355305E-01,1.197269315499389E+00,1.586595462103536E+00,2.229622629765312E+00,3.397374964441344E+00,5.851496811806930E+00,1.247569186064108E+01,4.306593086432805E+01};
      i += nWWnAmp.test_2to2_amp2([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH);
      i += nWWnAmp.test_2to2_amp2_rotations([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH);
      i += nWWnAmp.test_2to2_amp2_boosts([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH);
      i += nWWnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {6.111410474517026E+00,1.116053596357353E+00,6.732960819375647E-01,5.250985899324258E-01,4.608458297445206E-01,4.340580754284715E-01,4.293180641964791E-01,4.407893562907422E-01,4.666913876088977E-01,5.075912307263950E-01,5.659930744826498E-01,6.466067541718749E-01,7.572551956298835E-01,9.107452946446253E-01,1.128572333989559E+00,1.448569907109354E+00,1.941969594436709E+00,2.755705870220049E+00,4.233635576827009E+00,7.343502428465944E+00};
      i += nWWnAmp.test_2to2_amp2([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH2);
      i += nWWnAmp.test_2to2_amp2_rotations([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH2);
      i += nWWnAmp.test_2to2_amp2_boosts([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH2);
      i += nWWnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, MW=80.385, pspatial=95\n";
      me = 125.1;
      pspatial = 95;
      nWWnAmp.set_masses(me,MW);
      ldouble dataCH4[20] = {3.489705805196268E-01,3.347838462119878E-01,3.281639639792740E-01,3.278348990969592E-01,3.331772452918454E-01,3.440362108327760E-01,3.606384683299387E-01,3.835778559040821E-01,4.138563586107625E-01,4.529832036206782E-01,5.031499766760598E-01,5.675198695574364E-01,6.507030445722877E-01,7.595536788642447E-01,9.045523662000936E-01,1.102312623551724E+00,1.380382958622811E+00,1.787094110105507E+00,2.413547412221757E+00,3.448407490169949E+00};
      i += nWWnAmp.test_2to2_amp2([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH4);
      i += nWWnAmp.test_2to2_amp2_rotations([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH4);
      i += nWWnAmp.test_2to2_amp2_boosts([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH4);
      i += nWWnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH4);
      //std::cout<<"\n# me=125, MW=8.0, pspatial=125.1\n";
      me = 125;
      MW=8;
      pspatial = 125.1;
      nWWnAmp.set_masses(me,MW);
      ldouble dataCH3[20] = {3.658622495864547E+02,7.120513748798941E+02,8.730326076174214E+02,9.432034305316038E+02,9.662702271616200E+02,9.641125020214077E+02,9.483343786895200E+02,9.253613589333661E+02,8.988675011467785E+02,8.710050149986337E+02,8.430609615720616E+02,8.158264965494266E+02,7.898190518400831E+02,7.654349779902275E+02,7.430902361215672E+02,7.234304345479684E+02,7.078362706319579E+02,7.002182650029436E+02,7.173237607673050E+02,9.799724833776991E+02};
      i += nWWnAmp.test_2to2_amp2([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH3);
      i += nWWnAmp.test_2to2_amp2_rotations([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH3);
      i += nWWnAmp.test_2to2_amp2_boosts([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH3);
      i += nWWnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nWWnAmp.amp2(); }, 0,MW,MW,0,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}
