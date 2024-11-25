
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

//File:  SPINAS/SM/eeAA.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eeAA.h"

namespace spinas {

  eeAA::eeAA(const ldouble& echarge, const ldouble& masse):
    e(echarge), me(masse), prop(masse,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(me);
    p2=particle(me);
    p3=particle(0);
    p4=particle(0);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
  }
  void eeAA::set_masses(const ldouble& masse){
    me=masse;
    p1.set_mass(me);
    p2.set_mass(me);
    prop.set_mass(me);
  }
  void eeAA::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s34s.update();
    a34a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    s24s.update();
    a24a.update();
    s23s.update();
    a23a.update();
    s14s.update();
    a14a.update();
    s314a.update();
    s413a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT = prop.denominator(propTP);
    pDenU = prop.denominator(propUP);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eeAA::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    //No sign changes due to p3 and p4 being outgoing.
    if(ds3>0&&ds4>0){
      //me[34]^2<12>
      return 2.0*e*e*me*s34s.v()*s34s.v()*a12a.v(ds1,ds2)/pDenT/pDenU;
    }
    else if(ds3<0&&ds4<0){
      //me<34>^2[12]
      return 2.0*e*e*me*a34a.v()*a34a.v()*s12s.v(ds1,ds2)/pDenT/pDenU;
    }
    else if(ds3>0&&ds4<0){
      //([13]<24>+[23]<14>)[314>
      return 2.0*e*e*(s13s.v(ds1)*a24a.v(ds2)+s23s.v(ds2)*a14a.v(ds1))*s314a.v()/pDenT/pDenU;
    }
    else if(ds3<0&&ds4>0){
      //(<13>[24]+<23>[14])*[413>
      return 2.0*e*e*(a13a.v(ds1)*s24s.v(ds2)+a23a.v(ds2)*s14s.v(ds1))*s413a.v()/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }
  //set_momenta(...) must be called before amp2().
  ldouble eeAA::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-2;j4<=2;j4+=4){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Symmetry factor for identical photons 1/2
    return amp2/8.0;
  }
  //Alternate version that is a bit more efficient but requires more care from the author.
  //However, this time we'll leave redundant calculations of the spinor products to show how the calculation of the amplitude can be made more human readable.
  //No need to call set_momenta(...) before this.
  ldouble eeAA::amp2(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    ldouble propDT=std::real(prop.denominator(propTP));
    ldouble propDU=std::real(prop.denominator(propUP));
    

    ldouble amp2 = 0;

    cdouble s34, a34, s12, a12, s13, a13, s24, a24, s23, a23, s14, a14, s314, s413, numPP, numMM, numPM, numMP;
    s34 = p3.lsquare(2)*p4.rsquare(2);
    a34 = p3.langle(2)*p4.rangle(2);
    s314 = p3.lsquare(2)*p1.umat(2)*p4.rangle(2);
    s413 = p4.lsquare(2)*p1.umat(2)*p3.rangle(2);
    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2){
	s12 = p1.lsquare(j1,2)*p2.rsquare(j2,2);
	a12 = p1.langle(j1,2)*p2.rangle(j2, 2);	
	s13 = p1.lsquare(j1,2)*p3.rsquare(2);
	a13 = p1.langle(j1,2)*p3.rangle(2);
	s24 = p2.lsquare(j2,2)*p4.rsquare(2);
	a24 = p2.langle(j2,2)*p4.rangle(2);
	s23 = p2.lsquare(j2,2)*p3.rsquare(2);
	a23 = p2.langle(j2,2)*p3.rangle(2);
	s14 = p1.lsquare(j1,2)*p4.rsquare(2);
	a14 = p1.langle(j1,2)*p4.rangle(2);

	//Sign changes due to p3 and p4 being outgoing cancel.
	numPP = me*s34*s34*a12;//std::cout<<"numPP="<<numPP<<std::endl;
	numMM = me*a34*a34*s12;
	numPM = (s13*a24+s23*a14)*s314;
	numMP = (a13*s24+a23*s14)*s413;//std::cout<<"numMP="<<numMP<<std::endl;

	amp2 += std::pow(std::abs(numPP),2) + 
	  std::pow(std::abs(numMM),2) + 
	  std::pow(std::abs(numPM),2) + 
	  std::pow(std::abs(numMP),2); 
      }
    
    //Average over initial spins 1/2*1/2=1/4
    //Symmetry factor for identical photons 1/2
    return e*e*e*e*4.0*amp2/propDT/propDT/propDU/propDU/8.0;
  }


  



  //  Tests
  int test_eeAA(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> A , A       :";
    {//amp^2
      int i=0;
      // me=0.0005, pspatial=250
      ldouble me=0.0005;
      ldouble EE=0.31333;
      eeAA eeAAAmp = eeAA(EE,me);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.761473098865852E-01,1.196559098890515E-01,6.884618493553914E-02,4.748300512537967E-02,3.599742458224401E-02,2.906647080620063E-02,2.465909507168581E-02,2.184718935306326E-02,2.016436086672569E-02,1.937355800660445E-02,1.937355800660445E-02,2.016436086672569E-02,2.184718935306325E-02,2.465909507168581E-02,2.906647080620063E-02,3.599742458224402E-02,4.748300512537965E-02,6.884618493553911E-02,1.196559098890515E-01,3.761473098865843E-01};
      i += eeAAAmp.test_2to2_amp2([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eeAAAmp.test_2to2_amp2_rotations([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eeAAAmp.test_2to2_amp2_boosts([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eeAAAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      //Close to threshold
      pspatial = 0.0001;
      ldouble dataCH2[20] = {2.081251342993773E-02,2.079750983522615E-02,2.078115831044515E-02,2.076461266895772E-02,2.074882438042710E-02,2.073456777596373E-02,2.072246020291247E-02,2.071297785405228E-02,2.070646782725441E-02,2.070315683100089E-02,2.070315683100089E-02,2.070646782725441E-02,2.071297785405228E-02,2.072246020291247E-02,2.073456777596374E-02,2.074882438042710E-02,2.076461266895772E-02,2.078115831044515E-02,2.079750983522615E-02,2.081251342993773E-02};
      i += eeAAAmp.test_2to2_amp2([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eeAAAmp.test_2to2_amp2_rotations([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eeAAAmp.test_2to2_amp2_boosts([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eeAAAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
