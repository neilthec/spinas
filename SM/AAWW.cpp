
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

//File:  SPINAS/SM/AAWW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AAWW.h"

namespace spinas {

  AAWW::AAWW(const ldouble& echarge, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WW(widthW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    propW = propagator(MW,WW);
    p1=particle(0);
    p2=particle(0);
    p3=particle(MW);
    p4=particle(MW);
    //<12>,[12],<34>,[34],<14>,[14],<13>,[13],[24],<24>,[23],<23>
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
  }
  void AAWW::set_masses(const ldouble& massW){
    MW=massW;
    p3.set_mass(MW);
    p4.set_mass(MW);
    propW.set_mass(MW);
  }
  void AAWW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<34>,[34],<14>,[14],<13>,[13],[24],<24>,[23],<23>
    s12s.update();
    a12a.update();
    s34s.update();
    a34a.update();
    s14s.update();
    a14a.update();
    s13s.update();
    a13a.update();
    s24s.update();
    a24a.update();
    s23s.update();
    a23a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT=propW.denominator(propTP);
    pDenU=propW.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AAWW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble two=2;
    int ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);

      if(ds1>0 && ds2>0){
      
	//TU Diagram
	//all ingoing:
	// - 2e^2 [12]^2 <34>^2 /(t-MW^2)(u-MW^2)
	//34 outgoing:
	// - 2e^2 [12]^2 <34>^2 /(t-MW^2)(u-MW^2)
	amplitude += - normFactor*2.0*e*e*s12s.v()*s12s.v()*a34a.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)/pDenT/pDenU;
	
	
      }
      else if(ds1>0 && ds2<0){

	//TU Diagram
	//all ingoing:
	//-2e^2 (<24>[13]+<23>[14])^2 /(t-MW^2)(u-MW^2)
	//34 out:
	//-2e^2 (<24>[13]+<23>[14])^2 /(t-MW^2)(u-MW^2)
	amplitude += - normFactor*2.0*e*e*(a24a.v(ds4a)*s13s.v(ds3a)+a23a.v(ds3a)*s14s.v(ds4a))*(a24a.v(ds4b)*s13s.v(ds3b)+a23a.v(ds3b)*s14s.v(ds4b))/pDenT/pDenU;

	
      }
      else if(ds1<0 && ds2>0){
	amplitude += - normFactor*2.0*e*e*(s24s.v(ds4a)*a13a.v(ds3a)+s23s.v(ds3a)*a14a.v(ds4a))*(s24s.v(ds4b)*a13a.v(ds3b)+s23s.v(ds3b)*a14a.v(ds4b))/pDenT/pDenU;
      }
      else if(ds1<0 && ds2<0){
	amplitude += - normFactor*2.0*e*e*a12a.v()*a12a.v()*s34s.v(ds3a,ds4a)*s34s.v(ds3b,ds4b)/pDenT/pDenU;
      }
      else{
	std::cout<<"Photon cannot have helicity 0!\n";
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AAWW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2^2=1/4
    return amp2/4.0;
  }

  //set_momenta(...) must be called before amp2().
  ldouble AAWW::amp2_Aplus_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=2)
      for(int j4=-2;j4<=2;j4+=2){
	M = amp(2,2,j3,j4);
	amp2 += std::pow(std::abs(M),2);
      }
    return amp2;
  }

  //set_momenta(...) must be called before amp2().
  ldouble AAWW::amp2_Aplus_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=2)
      for(int j4=-2;j4<=2;j4+=2){
	M = amp(2,-2,j3,j4);
	amp2 += std::pow(std::abs(M),2);
      }
    return amp2;
  }



  



  //  Tests
  int test_AAWW(){
    int n=0;//Number of fails
    std::cout<<"\t* A , A  -> W+, W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      AAWW AAWWAmp = AAWW(EE,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.535590404082744E+01,4.131068661007738E+00,1.927358664831104E+00,1.143902386021833E+00,7.819834389991045E-01,5.891241795646013E-01,4.779490828796939E-01,4.120851686366945E-01,3.746576984708690E-01,3.576005427899943E-01,3.576005427899948E-01,3.746576984708695E-01,4.120851686366942E-01,4.779490828796945E-01,5.891241795646016E-01,7.819834389991039E-01,1.143902386021832E+00,1.927358664831104E+00,4.131068661007737E+00,1.535590404082740E+01};
      i += AAWWAmp.test_2to2_amp2([&]() { return AAWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH);
      i += AAWWAmp.test_2to2_amp2_rotations([&]() { return AAWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH);
      i += AAWWAmp.test_2to2_amp2_boosts([&]() { return AAWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH);
      i += AAWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH);
      ldouble dataCHpp[20] = {1.522546645669118E+01,4.468778672802030E+00,2.256326152960536E+00,1.436564847661838E+00,1.043699269892456E+00,8.275469268699851E-01,6.995537547515084E-01,6.220901929580858E-01,5.773830423457349E-01,5.568154797402186E-01,5.568154797402196E-01,5.773830423457358E-01,6.220901929580852E-01,6.995537547515096E-01,8.275469268699860E-01,1.043699269892455E+00,1.436564847661836E+00,2.256326152960536E+00,4.468778672802031E+00,1.522546645669116E+01};
      i += AAWWAmp.test_2to2_amp2([&]() { return AAWWAmp.amp2_Aplus_Aplus(); }, 0,0,MW,MW,pspatial,dataCHpp);
      i += AAWWAmp.test_2to2_amp2_rotations([&]() { return AAWWAmp.amp2_Aplus_Aplus(); }, 0,0,MW,MW,pspatial,dataCHpp);
      i += AAWWAmp.test_2to2_amp2_boosts([&]() { return AAWWAmp.amp2_Aplus_Aplus(); }, 0,0,MW,MW,pspatial,dataCHpp);
      i += AAWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAWWAmp.amp2_Aplus_Aplus(); }, 0,0,MW,MW,pspatial,dataCHpp);
      ldouble dataCHpm[20] = {1.548634162496369E+01,3.793358649213447E+00,1.598391176701673E+00,8.512399243818279E-01,5.202676081057530E-01,3.507014322592174E-01,2.563444110078794E-01,2.020801443153033E-01,1.719323545960031E-01,1.583856058397699E-01,1.583856058397700E-01,1.719323545960032E-01,2.020801443153032E-01,2.563444110078794E-01,3.507014322592172E-01,5.202676081057527E-01,8.512399243818276E-01,1.598391176701672E+00,3.793358649213443E+00,1.548634162496366E+01};
      i += AAWWAmp.test_2to2_amp2([&]() { return AAWWAmp.amp2_Aplus_Aminus(); }, 0,0,MW,MW,pspatial,dataCHpm);
      i += AAWWAmp.test_2to2_amp2_rotations([&]() { return AAWWAmp.amp2_Aplus_Aminus(); }, 0,0,MW,MW,pspatial,dataCHpm);
      i += AAWWAmp.test_2to2_amp2_boosts([&]() { return AAWWAmp.amp2_Aplus_Aminus(); }, 0,0,MW,MW,pspatial,dataCHpm);
      i += AAWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAWWAmp.amp2_Aplus_Aminus(); }, 0,0,MW,MW,pspatial,dataCHpm);
      //std::cout<<"\n# MW=80.385, pspatial=81\n";
      pspatial = 81;
      ldouble dataCH4[20] = {3.790147018650531E-01,3.760730895346712E-01,3.734813507801401E-01,3.712311749487209E-01,3.693154061250488E-01,3.677279954516163E-01,3.664639612138025E-01,3.655193562361521E-01,3.648912422234288E-01,3.645776707604673E-01,3.645776707604673E-01,3.648912422234289E-01,3.655193562361521E-01,3.664639612138025E-01,3.677279954516163E-01,3.693154061250488E-01,3.712311749487208E-01,3.734813507801401E-01,3.760730895346712E-01,3.790147018650531E-01};
      i += AAWWAmp.test_2to2_amp2([&]() { return AAWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH4);
      i += AAWWAmp.test_2to2_amp2_rotations([&]() { return AAWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH4);
      i += AAWWAmp.test_2to2_amp2_boosts([&]() { return AAWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH4);
      i += AAWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAWWAmp.amp2(); }, 0,0,MW,MW,pspatial,dataCH4);
      ldouble dataCH4pp[20] = {1.249076564786190E-01,1.242208499852037E-01,1.236150996972334E-01,1.230886954805995E-01,1.226401635127734E-01,1.222682570168767E-01,1.219719484988745E-01,1.217504234026163E-01,1.216030751137465E-01,1.215295012586249E-01,1.215295012586249E-01,1.216030751137465E-01,1.217504234026163E-01,1.219719484988745E-01,1.222682570168768E-01,1.226401635127734E-01,1.230886954805995E-01,1.236150996972333E-01,1.242208499852037E-01,1.249076564786190E-01};
      i += AAWWAmp.test_2to2_amp2([&]() { return AAWWAmp.amp2_Aplus_Aplus(); }, 0,0,MW,MW,pspatial,dataCH4pp);
      i += AAWWAmp.test_2to2_amp2_rotations([&]() { return AAWWAmp.amp2_Aplus_Aplus(); }, 0,0,MW,MW,pspatial,dataCH4pp);
      i += AAWWAmp.test_2to2_amp2_boosts([&]() { return AAWWAmp.amp2_Aplus_Aplus(); }, 0,0,MW,MW,pspatial,dataCH4pp);
      i += AAWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAWWAmp.amp2_Aplus_Aplus(); }, 0,0,MW,MW,pspatial,dataCH4pp);
      //std::cout<<"A+A-:\n";
      ldouble dataCH4pm[20] = {6.331217472514873E-01,6.279253290841387E-01,6.233476018630468E-01,6.193736544168422E-01,6.159906487373242E-01,6.131877338863559E-01,6.109559739287305E-01,6.092882890696879E-01,6.081794093331112E-01,6.076258402623098E-01,6.076258402623098E-01,6.081794093331112E-01,6.092882890696879E-01,6.109559739287305E-01,6.131877338863559E-01,6.159906487373243E-01,6.193736544168421E-01,6.233476018630468E-01,6.279253290841387E-01,6.331217472514873E-01};
      i += AAWWAmp.test_2to2_amp2([&]() { return AAWWAmp.amp2_Aplus_Aminus(); }, 0,0,MW,MW,pspatial,dataCH4pm);
      i += AAWWAmp.test_2to2_amp2_rotations([&]() { return AAWWAmp.amp2_Aplus_Aminus(); }, 0,0,MW,MW,pspatial,dataCH4pm);
      i += AAWWAmp.test_2to2_amp2_boosts([&]() { return AAWWAmp.amp2_Aplus_Aminus(); }, 0,0,MW,MW,pspatial,dataCH4pm);
      i += AAWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAWWAmp.amp2_Aplus_Aminus(); }, 0,0,MW,MW,pspatial,dataCH4pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }



}
