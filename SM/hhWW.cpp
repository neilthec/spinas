
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

//File:  SPINAS/SM/hhWW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/hhWW.h"

namespace spinas {

  hhWW::hhWW(const ldouble& echarge, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WW(widthW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    propW = propagator(MW,WW);
    proph = propagator(mh,wh);  
    p1=particle(mh);
    p2=particle(mh);
    p3=particle(MW);
    p4=particle(MW);
    //<34>,[34]
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    //[314>,[413>
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
    //Couplings
    pre = e*e/(2.0*MW*MW*SW*SW);
    preS = 3.0*pre*mh*mh;
  }
  void hhWW::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    p1.set_mass(mh);
    p2.set_mass(mh);
    p3.set_mass(MW);
    p4.set_mass(MW);
    propW.set_mass(MW);
    proph.set_mass(mh);
    //Couplings
    pre = e*e/(2.0*MW*MW*SW*SW);
    preS = 3.0*pre*mh*mh;
  }
  void hhWW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<34>,[34]
    s34s.update();
    a34a.update();
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
    pDenT=propW.denominator(propTP);
    pDenU=propW.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble hhWW::amp(const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1;
    int ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);


      //pre = e*e/(2.0*MW*MW*SW*SW);
      
      //4-Point      
      amplitude += - normFactor*pre*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b);
          
      //S-Channel h
      //preS = 3.0*pre*mh*mh;
      //all ingoing: 
      //preS [34] <34> /(s-Mh^2)
      //34 outgoing:
      //preS [34] <34> /(s-Mh^2)
      amplitude += normFactor*preS*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)/pDenS;
      
      //T-Channel Z
      //all ingoing:
      //+pre ( 2MW^2[34]<34> + MW([34][314>+<34>[413>) + [314>[413> )/(t-MW^2)
      //34 outgoing:
      //+pre ( 2MW^2[34]<34> - MW([34][314>+<34>[413>) + [314>[413> )/(t-MW^2)
      amplitude += normFactor*pre*2.0*MW*MW*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)/pDenT
      	-          normFactor*pre*MW*(s34s.v(ds3a,ds4a)*s314a.v(ds3b,ds4b)+a34a.v(ds3a,ds4a)*s413a.v(ds4b,ds3b))/pDenT
      	+          normFactor*pre*s314a.v(ds3a,ds4a)*s413a.v(ds4b,ds3b)/pDenT;
      
      //U-Channel Z
      //all ingoing:
      //+pre ( 2MW^2[34]<34> - MW([34][413>+<34>[314>) + [314>[413> )/(u-MW^2)
      //34 outgoing:
      //+pre ( 2MW^2[34]<34> + MW([34][413>+<34>[314>) + [314>[413> )/(u-MW^2)
      amplitude += normFactor*pre*2.0*MW*MW*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)/pDenU
      	+          normFactor*pre*MW*(s34s.v(ds3a,ds4a)*s413a.v(ds4b,ds3b)+a34a.v(ds3a,ds4a)*s314a.v(ds3b,ds4b))/pDenU
      	+          normFactor*pre*s314a.v(ds3a,ds4a)*s413a.v(ds4b,ds3b)/pDenU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble hhWW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=2)
      for(int j4=-2;j4<=2;j4+=2){
	M = amp(j3,j4);
	amp2 += std::pow(std::abs(M),2);
      }
    return amp2;
  }

  



  //  Tests
  int test_hhWW(){
    int n=0;//Number of fails
    std::cout<<"\t* h , h  -> W+, W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      hhWW hhWWAmp = hhWW(EE,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.150277119756263E+01,4.974236429948228E+00,2.079273190872931E+00,1.127134676724903E+00,7.157433258543224E-01,5.092771492164054E-01,3.963895476553964E-01,3.324081349869982E-01,2.972538492020933E-01,2.815674313550883E-01,2.815674313550908E-01,2.972538492020919E-01,3.324081349869970E-01,3.963895476553972E-01,5.092771492164029E-01,7.157433258543201E-01,1.127134676724904E+00,2.079273190872930E+00,4.974236429948234E+00,2.150277119756258E+01};
      i += hhWWAmp.test_2to2_amp2([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH);
      i += hhWWAmp.test_2to2_amp2_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH);
      i += hhWWAmp.test_2to2_amp2_boosts([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH);
      i += hhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {5.833223442291408E+00,2.787967955327683E+00,1.622205650744733E+00,1.069513687614552E+00,7.723295983050089E-01,5.996269470598655E-01,4.948802672960181E-01,4.309583377202177E-01,3.940260389823919E-01,3.770541855960114E-01,3.770541855960104E-01,3.940260389823912E-01,4.309583377202172E-01,4.948802672960183E-01,5.996269470598664E-01,7.723295983050097E-01,1.069513687614553E+00,1.622205650744731E+00,2.787967955327679E+00,5.833223442291399E+00};
      i += hhWWAmp.test_2to2_amp2([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH2);
      i += hhWWAmp.test_2to2_amp2_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH2);
      i += hhWWAmp.test_2to2_amp2_boosts([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH2);
      i += hhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=95\n";
      pspatial = 1;
      ldouble dataCH4[20] = {9.131466230088834E-01,9.130755122584028E-01,9.130123085138192E-01,9.129570097235838E-01,9.129096140928833E-01,9.128701200835221E-01,9.128385264138245E-01,9.128148320585502E-01,9.127990362488253E-01,9.127911384720966E-01,9.127911384720963E-01,9.127990362488245E-01,9.128148320585502E-01,9.128385264138243E-01,9.128701200835224E-01,9.129096140928830E-01,9.129570097235834E-01,9.130123085138192E-01,9.130755122584021E-01,9.131466230088842E-01};
      i += hhWWAmp.test_2to2_amp2([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH4);
      i += hhWWAmp.test_2to2_amp2_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH4);
      i += hhWWAmp.test_2to2_amp2_boosts([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH4);
      i += hhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH4);
      //std::cout<<"\n# mh=5, MW=80.385, pspatial=95\n";
      mh = 5;
      pspatial = 95;
      hhWWAmp.set_masses(mh,MW);
      ldouble dataCH6[20] = {9.944258803238153E-01,8.228482270991687E-01,6.993434068957985E-01,6.088717182061071E-01,5.420137905479475E-01,4.926853572199853E-01,4.568797410300052E-01,4.319476002803191E-01,4.161738466701093E-01,4.085287557427670E-01,4.085287557427669E-01,4.161738466701093E-01,4.319476002803190E-01,4.568797410300052E-01,4.926853572199851E-01,5.420137905479475E-01,6.088717182061069E-01,6.993434068957984E-01,8.228482270991684E-01,9.944258803238147E-01};
      i += hhWWAmp.test_2to2_amp2([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH6);
      i += hhWWAmp.test_2to2_amp2_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH6);
      i += hhWWAmp.test_2to2_amp2_boosts([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH6);
      i += hhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH6);
      //std::cout<<"\n# mh=5, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH7[20] = {1.898186157808347E+01,5.096596570375555E+00,2.386182889434564E+00,1.427602131572994E+00,9.869106029882752E-01,7.530921498680503E-01,6.188119739939866E-01,5.395038171907973E-01,4.945391624414562E-01,4.740757426734735E-01,4.740757426734779E-01,4.945391624414551E-01,5.395038171907941E-01,6.188119739939900E-01,7.530921498680496E-01,9.869106029882753E-01,1.427602131572995E+00,2.386182889434568E+00,5.096596570375555E+00,1.898186157808344E+01};
      i += hhWWAmp.test_2to2_amp2([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH7);
      i += hhWWAmp.test_2to2_amp2_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH7);
      i += hhWWAmp.test_2to2_amp2_boosts([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH7);
      i += hhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH7);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=125.1\n";
      mh = 0.0005;
      pspatial = 95;
      hhWWAmp.set_masses(mh,MW);
      ldouble dataCH3[20] = {9.901508707141952E-01,8.206042338525565E-01,6.983357833813413E-01,6.086365996548956E-01,5.422702024932342E-01,4.932576374104955E-01,4.576547125172795E-01,4.328497818690999E-01,4.171503632371886E-01,4.095395398686569E-01,4.095395398686568E-01,4.171503632371885E-01,4.328497818690998E-01,4.576547125172795E-01,4.932576374104953E-01,5.422702024932341E-01,6.086365996548955E-01,6.983357833813412E-01,8.206042338525564E-01,9.901508707141947E-01};
      i += hhWWAmp.test_2to2_amp2([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH3);
      i += hhWWAmp.test_2to2_amp2_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH3);
      i += hhWWAmp.test_2to2_amp2_boosts([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH3);
      i += hhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH3);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH5[20] = {1.897757431524353E+01,5.097314370211791E+00,2.387004797608648E+00,1.428300479953908E+00,9.874992741623664E-01,7.536004288261182E-01,6.192638157622773E-01,5.399175836363327E-01,4.949293977127790E-01,4.744547320923306E-01,4.744547320923332E-01,4.949293977127784E-01,5.399175836363310E-01,6.192638157622764E-01,7.536004288261170E-01,9.874992741623650E-01,1.428300479953909E+00,2.387004797608646E+00,5.097314370211780E+00,1.897757431524348E+01};
      i += hhWWAmp.test_2to2_amp2([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH5);
      i += hhWWAmp.test_2to2_amp2_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH5);
      i += hhWWAmp.test_2to2_amp2_boosts([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH5);
      i += hhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH5);
      //std::cout<<"\n# mh=0, MW=80.385, pspatial=95\n";
      mh = 0;
      pspatial = 95;
      hhWWAmp.set_masses(mh,MW);
      ldouble dataCH8[20] = {9.901508706714814E-01,8.206042338301200E-01,6.983357833712667E-01,6.086365996525551E-01,5.422702024958206E-01,4.932576374162513E-01,4.576547125250712E-01,4.328497818781704E-01,4.171503632470066E-01,4.095395398788197E-01,4.095395398788197E-01,4.171503632470066E-01,4.328497818781702E-01,4.576547125250712E-01,4.932576374162513E-01,5.422702024958206E-01,6.086365996525550E-01,6.983357833712667E-01,8.206042338301199E-01,9.901508706714811E-01};
      i += hhWWAmp.test_2to2_amp2([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH8);
      i += hhWWAmp.test_2to2_amp2_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH8);
      i += hhWWAmp.test_2to2_amp2_boosts([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH8);
      i += hhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH8);
      //std::cout<<"\n# mh=0, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH9[20] = {1.897757431520066E+01,5.097314370218974E+00,2.387004797616870E+00,1.428300479960896E+00,9.874992741682559E-01,7.536004288312030E-01,6.192638157667946E-01,5.399175836404710E-01,4.949293977166815E-01,4.744547320961225E-01,4.744547320961254E-01,4.949293977166837E-01,5.399175836404699E-01,6.192638157667976E-01,7.536004288312043E-01,9.874992741682556E-01,1.428300479960892E+00,2.387004797616869E+00,5.097314370218969E+00,1.897757431520061E+01};
      i += hhWWAmp.test_2to2_amp2([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH9);
      i += hhWWAmp.test_2to2_amp2_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH9);
      i += hhWWAmp.test_2to2_amp2_boosts([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH9);
      i += hhWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhWWAmp.amp2(); }, mh,mh,MW,MW,pspatial,dataCH9);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  


}
