
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

//File:  SPINAS/SM/nenenene2.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/nenenene2.h"

namespace spinas {
  //Constructors
  nenenene2::nenenene2(const ldouble& echarge, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    p1(particle(0)), p2(particle(0)),
    p3(particle(0)), p4(particle(0)),
    a12a(sproduct(ANGLE,&p1,&p2,2)),
    s34s(sproduct(SQUARE,&p3,&p4,2))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    preZ = e*e/(4.0*CW*CW*SW*SW);
  }
  void nenenene2::set_masses(const ldouble& massW){
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
  }
  void nenenene2::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a12a.update();
    s34s.update();
    //Propagator Momentum
    ldouble propPU[4], propPT[4];
    for(int j=0;j<4;j++){
      propPU[j] = mom1[j]-mom4[j];
      propPT[j] = mom1[j]-mom3[j];
    }
    pDenUZ = propZ.denominator(propPU);
    pDenTZ = propZ.denominator(propPT);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble nenenene2::amp(){
    constexpr ldouble one = 1, two = 2;
    cdouble amplitude(0,0);
    
    //Z Boson
    //Defined above:
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //neNeNene all in:
    //+ preZ [23] <14> /(s-MZ^2)
    //neneNeNe: 2<->4
    //- preZ [34] <12> /(u-MZ^2)
    //34 out:
    //- preZ [34] <12> /(u-MZ^2)
    amplitude += - two*preZ*s34s.v()*a12a.v()*(one/pDenUZ+one/pDenTZ);

    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble nenenene2::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    M = amp();
    amp2 = std::pow(std::abs(M),2);

    //Symmetry factor 1/2
    return amp2/2.0;
  }

  



  //  Tests
  int test_nenenene2(){
    int n=0;//Number of fails
    std::cout<<"\t* ne, ne -> ne, ne      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### MW=80.385, pspatial=250\n";
      ldouble MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      nenenene2 nenenene2Amp = nenenene2(0.31333,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.305574412242003E+01,4.190938086086727E+00,2.183949472905202E+00,1.412709832435260E+00,1.035960803253861E+00,8.262770622216340E-01,7.011942053930355E-01,6.251258954079334E-01,5.810897153009329E-01,5.607964149857158E-01,5.607964149857156E-01,5.810897153009328E-01,6.251258954079333E-01,7.011942053930355E-01,8.262770622216339E-01,1.035960803253861E+00,1.412709832435259E+00,2.183949472905202E+00,4.190938086086724E+00,1.305574412242001E+01};
      i += nenenene2Amp.test_2to2_amp2([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += nenenene2Amp.test_2to2_amp2_rotations([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += nenenene2Amp.test_2to2_amp2_boosts([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += nenenene2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      //std::cout<<"########### MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20,3.658343410537527E-20};
      i += nenenene2Amp.test_2to2_amp2([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += nenenene2Amp.test_2to2_amp2_rotations([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += nenenene2Amp.test_2to2_amp2_boosts([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += nenenene2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      //std::cout<<"########### MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      nenenene2Amp.set_masses(MW);
      ldouble dataCH3[20] = {2.976395526005115E-05,2.976345326221849E-05,2.976300705258446E-05,2.976261662738587E-05,2.976228198333008E-05,2.976200311759496E-05,2.976178002782879E-05,2.976161271215026E-05,2.976150116914843E-05,2.976144539788267E-05,2.976144539788267E-05,2.976150116914843E-05,2.976161271215026E-05,2.976178002782879E-05,2.976200311759496E-05,2.976228198333008E-05,2.976261662738587E-05,2.976300705258446E-05,2.976345326221849E-05,2.976395526005115E-05};
      i += nenenene2Amp.test_2to2_amp2([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += nenenene2Amp.test_2to2_amp2_rotations([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += nenenene2Amp.test_2to2_amp2_boosts([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += nenenene2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      //std::cout<<"########### MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      nenenene2Amp.set_masses(MW);
      ldouble dataCH4[20] = {2.544048796249899E+00,1.523251618966690E+00,1.055919921697098E+00,8.033286230307900E-01,6.525673705522055E-01,5.571415718202755E-01,4.951034384521930E-01,4.551314719085170E-01,4.311138500652692E-01,4.198117795866771E-01,4.198117795866771E-01,4.311138500652690E-01,4.551314719085168E-01,4.951034384521929E-01,5.571415718202755E-01,6.525673705522055E-01,8.033286230307898E-01,1.055919921697098E+00,1.523251618966690E+00,2.544048796249898E+00};
      i += nenenene2Amp.test_2to2_amp2([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += nenenene2Amp.test_2to2_amp2_rotations([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += nenenene2Amp.test_2to2_amp2_boosts([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += nenenene2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return nenenene2Amp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
