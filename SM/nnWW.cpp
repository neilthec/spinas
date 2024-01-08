
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

//File:  SPINAS/SM/nnWW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/nnWW.h"

namespace spinas {

  nnWW::nnWW(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), prope(0,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);
    p1=particle(0);
    p2=particle(0);
    p3=particle(MW);
    p4=particle(MW);
    //[23],<13>,<34>,[34],[24],<14>
    s23s = sproduct(SQUARE,&p2,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s24s = sproduct(SQUARE,&p2,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    //[314>,[231>
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    s231a = sproduct(SQUARE,&p2,&p3,&p1);
    //Couplings
    pree = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preZ = 2.0*e*e/(2.0*MW*MW*SW*SW);
  }
  void nnWW::set_masses(const ldouble& masse, const ldouble& massW){
    me=masse;
    prope.set_mass(me);
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p3.set_mass(MW);
    p4.set_mass(MW);
    propZ.set_mass(MZ);
    //Couplings
    pree = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preZ = 2.0*e*e/(2.0*MW*MW*SW*SW);
  }
  void nnWW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //[23],<13>,<34>,[34],[24],<14>
    s23s.update();
    a13a.update();
    s34s.update();
    a34a.update();
    s24s.update();
    a14a.update();
    //[314>,[231>
    s314a.update();
    s231a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
    }
    pDeneT=prope.denominator(propTP);
    pDenZS=propZ.denominator(propSP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble nnWW::amp(const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds4a, ds4b;
    constexpr ldouble two=2;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);
           
      //T-Channel e
      //pree = e*e/(2.0*MW*MW*SW*SW);
      //all ingoing:
      // pree [24] <13> ([314>+MW <34>))/t
      //34 outgoing:
      // pree [24] <13> (MW <34> - [314>)/t
      amplitude += - normFactor*pree*s24s.v(ds4a)*a13a.v(ds3a)*(MW*a34a.v(ds3b,ds4b)-s314a.v(ds3b,ds4b))/pDeneT;

      
      //S-Channel Z
      //preZ = e*e/(2.0*MW*MW*SW*SW); //=pree
      //all ingoing
      //+ preZ ( [34]<34>[231> + MW([23]<14>+[24]<13>)([34]+<34>) )/(s-MZ^2)
      //34 outgoing
      //- preZ ( [34]<34>[231> + MW([34]+<34>)([23]<14>+[24]<13>) )/(s-MZ^2)
      amplitude += - normFactor*preZ*( 
				      s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*s231a.v()
				      +MW*(s34s.v(ds3a,ds4a)+a34a.v(ds3a,ds4a))*(s23s.v(ds3b)*a14a.v(ds4b)+s24s.v(ds4b)*a13a.v(ds3b))
				      )/pDenZS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble nnWW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=2)
      for(int j4=-2;j4<=2;j4+=2){
	M = amp(j3,j4);
	amp2 += std::pow(std::abs(M),2);
      }
    //Nothing to average
    return amp2;
  }

  
  

  //  Tests
  int test_nnWW(){
    int n=0;//Number of fails
    std::cout<<"\t* ne, Ne -> W+, W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=250\n";
      ldouble me=0.0005;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      nnWW nnWWAmp = nnWW(EE,me,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {7.496715821185181E+00,2.259525477726942E+00,1.191074842252847E+00,7.457231876469343E-01,5.093980060267348E-01,3.676145446833791E-01,2.760329739824261E-01,2.138720982109863E-01,1.700975526978052E-01,1.383071880477013E-01,1.145312006977507E-01,9.617719370337260E-02,8.148046112051939E-02,6.919807809296361E-02,5.842932727255536E-02,4.850551069177717E-02,3.891976134287202E-02,2.928088045947702E-02,1.928212360077718E-02,8.679577088154368E-03};
      i += nnWWAmp.test_2to2_amp2([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH);
      i += nnWWAmp.test_2to2_amp2_rotations([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH);
      i += nnWWAmp.test_2to2_amp2_boosts([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH);
      i += nnWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] ={5.032259812953072E+00,2.493806608677821E+00,1.539221209074415E+00,1.058278109926491E+00,7.746599311848482E-01,5.910078449460399E-01,4.645020185316381E-01,3.734067691971000E-01,3.055263418563311E-01,2.534917115016790E-01,2.125888751285799E-01,1.796623880315274E-01,1.525206304913928E-01,1.295940186006691E-01,1.097289916895628E-01,9.205877111085393E-02,7.591947145725975E-02,6.079402611424305E-02,4.627373300906300E-02,3.203128333255355E-02};
      i += nnWWAmp.test_2to2_amp2([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH2);
      i += nnWWAmp.test_2to2_amp2_rotations([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH2);
      i += nnWWAmp.test_2to2_amp2_boosts([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH2);
      i += nnWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, MW=80.385, pspatial=95\n";
      me = 125.1;
      pspatial = 95;
      nnWWAmp.set_masses(me,MW);
      ldouble dataCH4[20] = {6.530788845510108E-01,5.845065052967613E-01,5.236124081043461E-01,4.692849291038127E-01,4.205941850174836E-01,3.767565467050234E-01,3.371070903701099E-01,3.010780227811426E-01,2.681816287979061E-01,2.379966771414200E-01,2.101574959237704E-01,1.843451276350165E-01,1.602801173655830E-01,1.377165938776011E-01,1.164373816403554E-01,9.624994071057200E-02,7.698297571370875E-02,5.848358896811811E-02,4.061487871806973E-02,2.325390348066774E-02};
      i += nnWWAmp.test_2to2_amp2([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH4);
      i += nnWWAmp.test_2to2_amp2_rotations([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH4);
      i += nnWWAmp.test_2to2_amp2_boosts([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH4);
      i += nnWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH4);
      //std::cout<<"\n# me=125, MW=8.0, pspatial=125.1\n";
      me = 125;
      MW=8;
      pspatial = 125.1;
      nnWWAmp.set_masses(me,MW);
      ldouble dataCH3[20] = {1.068983262467667E+03,1.974312165400613E+03,2.289712777015868E+03,2.331939209846004E+03,2.243502720559766E+03,2.093246440077219E+03,1.916064595842957E+03,1.730162012232057E+03,1.545128852900057E+03,1.365944844177143E+03,1.195059277628614E+03,1.033512216628156E+03,8.815567511002222E+02,7.390120790809661E+02,6.054669767982549E+02,4.803980488255356E+02,3.632384708836291E+02,2.534175339459436E+02,1.503827858482094E+02,5.361174443595100E+01};
      i += nnWWAmp.test_2to2_amp2([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH3);
      i += nnWWAmp.test_2to2_amp2_rotations([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH3);
      i += nnWWAmp.test_2to2_amp2_boosts([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH3);
      i += nnWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nnWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}
