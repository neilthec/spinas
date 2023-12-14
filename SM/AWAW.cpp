
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

//File:  SPINAS/SM/AWAW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AWAW.h"

namespace spinas {

  AWAW::AWAW(const ldouble& echarge, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WW(widthW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    propW = propagator(MW,WW);
    p1=particle(0);
    p2=particle(MW);
    p3=particle(0);
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
  void AWAW::set_masses(const ldouble& massW){
    MW=massW;
    p2.set_mass(MW);
    p4.set_mass(MW);
    propW.set_mass(MW);
  }
  void AWAW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    ldouble propSP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propW.den(propSP);
    pDenU=propW.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AWAW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble two=2;
    int ds4a, ds4b, ds2a, ds2b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds4,ds2);
    ldouble normFactor=get_spin_normalization(ds4,ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds4,ds4a,ds4b, ds2,ds2a,ds2b, i);

      if(ds1>0 && ds3>0){
      
	//TU Diagram
	//AAW-W+ all ingoing:
	// - e^2 [12]^2 <34>^2 /(t-MW^2)(u-MW^2)
	//AW+AW-: 4->2->3->4:
	// - e^2 [13]^2 <24>^2 /(u-MW^2)(s-MW^2)
	//34 out:
	// - e^2 [13]^2 <24>^2 /(u-MW^2)(s-MW^2)
	amplitude += - normFactor*2.0*e*e*s13s.v()*s13s.v()*a24a.v(ds2a,ds4a)*a24a.v(ds2b,ds4b)/pDenU/pDenS;
	
	
      }
      else if(ds1>0 && ds3<0){

	//TU Diagram
	//AAW-W+ all ingoing:
	// 2e^2*( [13]^2<24>^2 + [14]^2<23>^2 + 2[13][14]<23><24> )/(t-MW^2)(u-MW^2)
	//AW+AW-: 4->2->3->4:
	// 2e^2*( [14]^2<23>^2 + [12]^2<34>^2 - 2[14][12]<34><23> )/(u-MW^2)(s-MW^2)
	//34 out:
	// 2e^2*( [14]^2<23>^2 + [12]^2<34>^2 + 2[14][12]<34><23> )/(u-MW^2)(s-MW^2)
	amplitude += normFactor*2.0*e*e*(
				     s14s.v(ds4a)*s14s.v(ds4b)*a23a.v(ds2a)*a23a.v(ds2b)
				     +s12s.v(ds2a)*s12s.v(ds2b)*a34a.v(ds4a)*a34a.v(ds4b)
				     +two*s14s.v(ds4a)*s12s.v(ds2a)*a34a.v(ds4b)*a23a.v(ds2b)				     
				     )/pDenU/pDenS;

	
      }
      else if(ds1<0 && ds3>0){
	amplitude += normFactor*2.0*e*e*(
				     a14a.v(ds4a)*a14a.v(ds4b)*s23s.v(ds2a)*s23s.v(ds2b)
				     +a12a.v(ds2a)*a12a.v(ds2b)*s34s.v(ds4a)*s34s.v(ds4b)
				     +two*a14a.v(ds4a)*a12a.v(ds2a)*s34s.v(ds4b)*s23s.v(ds2b)				     
				     )/pDenU/pDenS;
      }
      else if(ds1<0 && ds3<0){
	amplitude += - normFactor*2.0*e*e*a13a.v()*a13a.v()*s24s.v(ds2a,ds4a)*s24s.v(ds2b,ds4b)/pDenU/pDenS;
      }
      else{
	std::cout<<"Photon cannot have helicity 0!\n";
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AWAW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/3=1/6
    return amp2/6.0;
  }

  //set_momenta(...) must be called before amp2().
  ldouble AWAW::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-2;j4<=2;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }

  //set_momenta(...) must be called before amp2().
  ldouble AWAW::amp2_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-2;j4<=2;j4+=2){
	  M = amp(-2,j2,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }



  



  //  Tests
  int test_AWAW(){
    int n=0;//Number of fails
    std::cout<<"\t* A , W+ -> A , W+      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      AWAW AWAWAmp = AWAW(EE,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.853664571705771E-02,3.869840853620973E-02,3.916253974133997E-02,3.999363676615806E-02,4.128115837745494E-02,4.315080763114322E-02,4.578221394093485E-02,4.943712939315763E-02,5.450575871929971E-02,6.158556187717782E-02,7.162080069719791E-02,8.616169831836924E-02,1.078740112960416E-01,1.416133989190397E-01,1.968975251577690E-01,2.942750504062007E-01,4.844419462304803E-01,9.197572591856272E-01,2.235617946315084E+00,9.733821628894983E+00};
      i += AWAWAmp.test_2to2_amp2([&]() { return AWAWAmp.amp2(); }, 0,MW,0,MW,pspatial,dataCH);
      i += AWAWAmp.test_2to2_amp2_rotations([&]() { return AWAWAmp.amp2(); }, 0,MW,0,MW,pspatial,dataCH);
      i += AWAWAmp.test_2to2_amp2_boosts([&]() { return AWAWAmp.amp2(); }, 0,MW,0,MW,pspatial,dataCH);
      i += AWAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWAWAmp.amp2(); }, 0,MW,0,MW,pspatial,dataCH);
      ldouble dataCHpp[20] = {3.853664571705771E-02,3.869840853620973E-02,3.916253974133997E-02,3.999363676615806E-02,4.128115837745494E-02,4.315080763114322E-02,4.578221394093485E-02,4.943712939315763E-02,5.450575871929971E-02,6.158556187717782E-02,7.162080069719791E-02,8.616169831836924E-02,1.078740112960416E-01,1.416133989190397E-01,1.968975251577690E-01,2.942750504062007E-01,4.844419462304803E-01,9.197572591856272E-01,2.235617946315084E+00,9.733821628894983E+00};
      i += AWAWAmp.test_2to2_amp2([&]() { return AWAWAmp.amp2_Aplus(); }, 0,MW,0,MW,pspatial,dataCHpp);
      i += AWAWAmp.test_2to2_amp2_rotations([&]() { return AWAWAmp.amp2_Aplus(); }, 0,MW,0,MW,pspatial,dataCHpp);
      i += AWAWAmp.test_2to2_amp2_boosts([&]() { return AWAWAmp.amp2_Aplus(); }, 0,MW,0,MW,pspatial,dataCHpp);
      i += AWAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWAWAmp.amp2_Aplus(); }, 0,MW,0,MW,pspatial,dataCHpp);
      ldouble dataCHpm[20] = {3.853664571705771E-02,3.869840853620973E-02,3.916253974133997E-02,3.999363676615806E-02,4.128115837745494E-02,4.315080763114322E-02,4.578221394093485E-02,4.943712939315763E-02,5.450575871929971E-02,6.158556187717782E-02,7.162080069719791E-02,8.616169831836924E-02,1.078740112960416E-01,1.416133989190397E-01,1.968975251577690E-01,2.942750504062007E-01,4.844419462304803E-01,9.197572591856272E-01,2.235617946315084E+00,9.733821628894983E+00};
      i += AWAWAmp.test_2to2_amp2([&]() { return AWAWAmp.amp2_Aminus(); }, 0,MW,0,MW,pspatial,dataCHpm);
      i += AWAWAmp.test_2to2_amp2_rotations([&]() { return AWAWAmp.amp2_Aminus(); }, 0,MW,0,MW,pspatial,dataCHpm);
      i += AWAWAmp.test_2to2_amp2_boosts([&]() { return AWAWAmp.amp2_Aminus(); }, 0,MW,0,MW,pspatial,dataCHpm);
      i += AWAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWAWAmp.amp2_Aminus(); }, 0,MW,0,MW,pspatial,dataCHpm);
      //std::cout<<"\n# MW=80.385, pspatial=81\n";
      pspatial = 81;
      ldouble dataCH4[20] = {3.824381063513672E-02,3.773433369366197E-02,3.739321314638147E-02,3.725272478135258E-02,3.735616776567990E-02,3.776208618853182E-02,3.855038368417938E-02,3.983133264351579E-02,4.175910417607324E-02,4.455252796361093E-02,4.852772738472456E-02,5.415085919139913E-02,6.212608727168517E-02,7.354781955893426E-02,9.017574183165374E-02,1.149577879094597E-01,1.530880048325617E-01,2.143161707153332E-01,3.185029185117990E-01,5.107944453329001E-01};
      i += AWAWAmp.test_2to2_amp2([&]() { return AWAWAmp.amp2(); }, 0,MW,0,MW,pspatial,dataCH4);
      i += AWAWAmp.test_2to2_amp2_rotations([&]() { return AWAWAmp.amp2(); }, 0,MW,0,MW,pspatial,dataCH4);
      i += AWAWAmp.test_2to2_amp2_boosts([&]() { return AWAWAmp.amp2(); }, 0,MW,0,MW,pspatial,dataCH4);
      i += AWAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWAWAmp.amp2(); }, 0,MW,0,MW,pspatial,dataCH4);
      ldouble dataCH4pp[20] = {3.824381063513672E-02,3.773433369366197E-02,3.739321314638147E-02,3.725272478135258E-02,3.735616776567990E-02,3.776208618853182E-02,3.855038368417938E-02,3.983133264351579E-02,4.175910417607324E-02,4.455252796361093E-02,4.852772738472456E-02,5.415085919139913E-02,6.212608727168517E-02,7.354781955893426E-02,9.017574183165374E-02,1.149577879094597E-01,1.530880048325617E-01,2.143161707153332E-01,3.185029185117990E-01,5.107944453329001E-01};
      i += AWAWAmp.test_2to2_amp2([&]() { return AWAWAmp.amp2_Aplus(); }, 0,MW,0,MW,pspatial,dataCH4pp);
      i += AWAWAmp.test_2to2_amp2_rotations([&]() { return AWAWAmp.amp2_Aplus(); }, 0,MW,0,MW,pspatial,dataCH4pp);
      i += AWAWAmp.test_2to2_amp2_boosts([&]() { return AWAWAmp.amp2_Aplus(); }, 0,MW,0,MW,pspatial,dataCH4pp);
      i += AWAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWAWAmp.amp2_Aplus(); }, 0,MW,0,MW,pspatial,dataCH4pp);
      //std::cout<<"A+A-:\n";
      ldouble dataCH4pm[20] = {3.824381063513672E-02,3.773433369366197E-02,3.739321314638147E-02,3.725272478135258E-02,3.735616776567990E-02,3.776208618853182E-02,3.855038368417938E-02,3.983133264351579E-02,4.175910417607324E-02,4.455252796361093E-02,4.852772738472456E-02,5.415085919139913E-02,6.212608727168517E-02,7.354781955893426E-02,9.017574183165374E-02,1.149577879094597E-01,1.530880048325617E-01,2.143161707153332E-01,3.185029185117990E-01,5.107944453329001E-01};
      i += AWAWAmp.test_2to2_amp2([&]() { return AWAWAmp.amp2_Aminus(); }, 0,MW,0,MW,pspatial,dataCH4pm);
      i += AWAWAmp.test_2to2_amp2_rotations([&]() { return AWAWAmp.amp2_Aminus(); }, 0,MW,0,MW,pspatial,dataCH4pm);
      i += AWAWAmp.test_2to2_amp2_boosts([&]() { return AWAWAmp.amp2_Aminus(); }, 0,MW,0,MW,pspatial,dataCH4pm);
      i += AWAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWAWAmp.amp2_Aminus(); }, 0,MW,0,MW,pspatial,dataCH4pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }



}
