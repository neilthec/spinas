
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

//File:  SPINAS/SM/nnZZ.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/nnZZ.h"

namespace spinas {

  nnZZ::nnZZ(const ldouble& echarge, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propne(0,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    p1=particle(0);
    p2=particle(0);
    p3=particle(MZ);
    p4=particle(MZ);
    //[23], [24], <13>, <14>, <34>
    s23s = sproduct(SQUARE,&p2,&p3,2);
    s24s = sproduct(SQUARE,&p2,&p4,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    a14a = sproduct(ANGLE,&p1,&p4,2);
    a34a = sproduct(ANGLE,&p3,&p4,2);
    //[314>, [413>
    s314a = sproduct(SQUARE,&p3,&p1,&p4,2);
    s413a = sproduct(SQUARE,&p4,&p1,&p3,2);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
  }
  void nnZZ::set_masses(const ldouble& massW){
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p3.set_mass(MZ);
    p4.set_mass(MZ);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
  }
  void nnZZ::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //[23], [24], <13>, <14>, <34>
    s23s.update();
    s24s.update();
    a13a.update();
    a14a.update();
    a34a.update();
    //[314>, [413>
    s314a.update();
    s413a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT=propne.denominator(propTP);
    pDenU=propne.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble nnZZ::amp(const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);
      
      //T-Channel ne
      //preTU = e*e/(2.0*MW*MW*SW*SW);
      //all ingoing:
      //-preTU [24] <13> (MZ <34>+[314>))/t
      //34 outgoing:
      //+preTU [24] <13> (MZ <34>-[314>))/t
      amplitude += normFactor*preTU*s24s.v(ds4a)*a13a.v(ds3a)*(MZ*a34a.v(ds3b,ds4b)-s314a.v(ds3b,ds4b))/pDenT;

      //U-Channel ne
      //preTU = e*e/(2.0*MW*MW*SW*SW);
      //all ingoing:
      //+preTU [23] <14> (MZ <34>-[413>)/u
      //34 outgoing:
      //-preTU [23] <14> (MZ <34>+[413>)/u
      amplitude += -normFactor*preTU*s23s.v(ds3a)*a14a.v(ds4a)*(MZ*a34a.v(ds3b,ds4b)+s413a.v(ds4b,ds3b))/pDenU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble nnZZ::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=2)
      for(int j4=-2;j4<=2;j4+=2){
	M = amp(j3,j4);
	amp2 += std::pow(std::abs(M),2);
      }
    //Nothing to average over.
    //Symmetrize over final particles 1/2
    return amp2/2.0;
  }
  



  //  Tests
  int test_nnZZ(){
    int n=0;//Number of fails
    std::cout<<"\t* ne, Ne -> Z , Z       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      nnZZ nnZZAmp = nnZZ(EE,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.640346702831654E+00,5.630537017159620E-01,3.331848176861154E-01,2.348056810659994E-01,1.814835422899338E-01,1.491597450359460E-01,1.285474602684927E-01,1.153733334086609E-01,1.074803060397669E-01,1.037689069272044E-01,1.037689069272046E-01,1.074803060397662E-01,1.153733334086603E-01,1.285474602684923E-01,1.491597450359454E-01,1.814835422899336E-01,2.348056810659990E-01,3.331848176861150E-01,5.630537017159617E-01,1.640346702831651E+00};
      i += nnZZAmp.test_2to2_amp2([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH);
      i += nnZZAmp.test_2to2_amp2_rotations([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH);
      i += nnZZAmp.test_2to2_amp2_boosts([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH);
      i += nnZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH);
      //std::cout<<"\n# MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {9.822861064751632E-01,6.458000193151047E-01,4.786506543117693E-01,3.839018743801569E-01,3.247606400973517E-01,2.858054283921481E-01,2.596238921743521E-01,2.423106499675259E-01,2.317130296496787E-01,2.266703821266336E-01,2.266703821266336E-01,2.317130296496787E-01,2.423106499675258E-01,2.596238921743521E-01,2.858054283921480E-01,3.247606400973517E-01,3.839018743801568E-01,4.786506543117693E-01,6.458000193151044E-01,9.822861064751629E-01};
      i += nnZZAmp.test_2to2_amp2([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH2);
      i += nnZZAmp.test_2to2_amp2_rotations([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH2);
      i += nnZZAmp.test_2to2_amp2_boosts([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH2);
      i += nnZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH2);
      //std::cout<<"\n# MW=80.385, pspatial=95\n";
      pspatial = 95;
      ldouble dataCH4[20] = {6.690232281190653E-01,6.361857576511289E-01,6.090829278100771E-01,5.868535026792131E-01,5.688246751810152E-01,5.544781343819791E-01,5.434208737346482E-01,5.353624475959498E-01,5.300982401214256E-01,5.274978592459960E-01,5.274978592459960E-01,5.300982401214256E-01,5.353624475959498E-01,5.434208737346482E-01,5.544781343819791E-01,5.688246751810152E-01,5.868535026792130E-01,6.090829278100771E-01,6.361857576511288E-01,6.690232281190653E-01};
      i += nnZZAmp.test_2to2_amp2([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH4);
      i += nnZZAmp.test_2to2_amp2_rotations([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH4);
      i += nnZZAmp.test_2to2_amp2_boosts([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH4);
      i += nnZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH4);
      //std::cout<<"\n# MW=8.0, pspatial=125.1\n";
      pspatial = 125.1;
      MW=8;
      MZ=MW/CW;
      nnZZAmp.set_masses(MW);
      ldouble dataCH3[20] = {1.557907814488639E+00,4.959180460721435E-01,2.855195425055701E-01,1.970542256410944E-01,1.494916621374697E-01,1.207899239778720E-01,1.025384874149928E-01,9.089404047706429E-02,8.392522969946670E-02,8.065040967480019E-02,8.065040967475556E-02,8.392522970013594E-02,9.089404047738148E-02,1.025384874150277E-01,1.207899239780361E-01,1.494916621372391E-01,1.970542256411349E-01,2.855195425056500E-01,4.959180460721174E-01,1.557907814488632E+00};
      i += nnZZAmp.test_2to2_amp2([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH3);
      i += nnZZAmp.test_2to2_amp2_rotations([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH3);
      i += nnZZAmp.test_2to2_amp2_boosts([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH3);
      i += nnZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nnZZAmp.amp2(); }, 0,0,MZ,MZ,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
