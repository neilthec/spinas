
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

//File:  SPINAS/SM/nenenene.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/nenenene.h"

namespace spinas {
  //Constructors
  nenenene::nenenene(const ldouble& echarge, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propZ(MZ,WZ),
    p1(particle(0)), p2(particle(0)),
    p3(particle(0)), p4(particle(0)),
    a14a(sproduct(ANGLE,&p1,&p4)),
    s23s(sproduct(SQUARE,&p2,&p3))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    preZ = e*e/(4.0*CW*CW*SW*SW);
  }
  void nenenene::set_masses(const ldouble& massW){
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
  }
  void nenenene::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a14a.update();
    s23s.update();
    //Propagator Momentum
    ldouble propPS[4], propPT[4];
    for(int j=0;j<4;j++){
      propPS[j] = mom1[j]+mom2[j];
      propPT[j] = mom1[j]-mom3[j];
    }
    pDenSZ = propZ.denominator(propPS);
    pDenTZ = propZ.denominator(propPT);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble nenenene::amp(){
    constexpr ldouble one = 1, two = 2;
    cdouble amplitude(0,0);
    
    //Z Boson
    //Defined above:
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //all in:
    //+ EE^2 [23] <14> /(4 CW^2 SW^2 (s-MZ^2))
    //34 out:
    //- EE^2 [23] <14>/(4 CW^2 SW^2 (s-MZ^2))
    //=-preZ [23] <14> /(s-MZ^2)
    amplitude += - two*preZ*s23s.v()*a14a.v()*(one/pDenSZ+one/pDenTZ);

    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble nenenene::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    M = amp();
    amp2 = std::pow(std::abs(M),2);

    return amp2;
  }

  



  //  Tests
  int test_nenenene(){
    int n=0;//Number of fails
    std::cout<<"\t* ne, Ne -> ne, Ne      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### MW=80.385, pspatial=250\n";
      ldouble MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      nenenene neneneneAmp = nenenene(0.31333,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.958493373158536E+01,4.564067975806970E+00,1.695673147256786E+00,7.663104846912916E-01,3.837199257474553E-01,2.036219323333062E-01,1.115873108117309E-01,6.206061689742803E-02,3.455956884848822E-02,1.904065848672037E-02,1.025477826726943E-02,5.325542238590881E-03,2.621223687515258E-03,1.193825531557359E-03,4.850042791868497E-04,1.650449607465791E-04,4.145180103746178E-05,5.474792574218010E-06,3.608047448592609E-08,9.072308172590717E-08};
      i += neneneneAmp.test_2to2_amp2([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += neneneneAmp.test_2to2_amp2_rotations([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += neneneneAmp.test_2to2_amp2_boosts([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += neneneneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      //std::cout<<"########### MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {6.955425415877487E-20,6.260340167066258E-20,5.601838352426686E-20,4.979919971956137E-20,4.394585025651977E-20,3.845833513511574E-20,3.333665435532293E-20,2.858080791711500E-20,2.419079582046562E-20,2.016661806534845E-20,1.650827465173716E-20,1.321576557960541E-20,1.028909084892685E-20,7.728250459675151E-21,5.533244411823982E-21,3.704072705347001E-21,2.240735340217871E-21,1.143232316410256E-21,4.115636338978161E-22,4.572929265421565E-23};
      i += neneneneAmp.test_2to2_amp2([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += neneneneAmp.test_2to2_amp2_rotations([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += neneneneAmp.test_2to2_amp2_boosts([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += neneneneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      //std::cout<<"########### MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      neneneneAmp.set_masses(MW);
      ldouble dataCH3[20] = {5.815121827943336E-05,5.230414784629009E-05,4.677051560811897E-05,4.154967602040299E-05,3.664098574942839E-05,3.204380366319705E-05,2.775749082238245E-05,2.378141047132900E-05,2.011492802909438E-05,1.675741108053487E-05,1.370822936743317E-05,1.096675477966878E-05,8.532361346430314E-06,6.404425227470013E-06,4.582324704399745E-06,3.065440172028601E-06,1.853154129741696E-06,9.448511729200186E-07,3.399179844010727E-07,3.774332598010512E-08};
      i += neneneneAmp.test_2to2_amp2([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += neneneneAmp.test_2to2_amp2_rotations([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += neneneneAmp.test_2to2_amp2_boosts([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += neneneneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      //std::cout<<"########### MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      neneneneAmp.set_masses(MW);
      ldouble dataCH4[20] = {2.678596472792997E+00,1.143294284758584E+00,5.532724276706362E-01,2.869444709756265E-01,1.546230954935643E-01,8.483918381147589E-02,4.667735175630105E-02,2.540542737242124E-02,1.349071654126106E-02,6.875840989338449E-03,3.290346272559204E-03,1.429397330219876E-03,5.309675731820125E-04,1.478525720383952E-04,1.993937095602714E-05,3.322935658205920E-07,1.148325345910448E-05,1.856600236206932E-05,1.302825278201604E-05,2.289945387885474E-06};
      i += neneneneAmp.test_2to2_amp2([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += neneneneAmp.test_2to2_amp2_rotations([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += neneneneAmp.test_2to2_amp2_boosts([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += neneneneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neneneneAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
