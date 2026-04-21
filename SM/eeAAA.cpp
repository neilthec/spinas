
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2023-2026 Neil Christensen, Nick Majestic

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

//File:  SPINAS/SM/eeAAA.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eeAAA.h"

namespace spinas {

  eeAAA::eeAAA(const ldouble& echarge, const ldouble& masse):
    e(echarge), me(masse), prop(masse,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(me);
    p2=particle(me);
    p3=particle(0);
    p4=particle(0);
    p5=particle(0);
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
    s35s = sproduct(SQUARE,&p3,&p5);
    a35a = sproduct(ANGLE,&p3,&p5);
    s45s = sproduct(SQUARE,&p4,&p5);
    a45a = sproduct(ANGLE,&p4,&p5);
    s3123s = sproduct(SQUARE,&p3,&p1,&p2,&p3);
    a3123a = sproduct(ANGLE,&p3,&p1,&p2,&p3);
    s4124s = sproduct(SQUARE,&p4,&p1,&p2,&p4);
    a4124a = sproduct(ANGLE,&p4,&p1,&p2,&p4);
    s5125s = sproduct(SQUARE,&p5,&p1,&p2,&p5);
    a5125a = sproduct(ANGLE,&p5,&p1,&p2,&p5);
  }
  void eeAAA::set_masses(const ldouble& masse){
    me=masse;
    p1.set_mass(me);
    p2.set_mass(me);
    prop.set_mass(me);
  }
  void eeAAA::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4], const ldouble mom5[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    p5.set_momentum(mom5);
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
    s35s.update();
    a35a.update();
    s45s.update();
    a45a.update();
    s3123s.update();
    a3123a.update();
    s4124s.update();
    a4124a.update();
    s5125s.update();
    a5125a.update();
    //Propagator Momentum
    ldouble propS13P[4], propS14P[4], propS15P[4], propS23P[4], propS24P[4], propS25P[4];
    for(int j=0;j<4;j++){
      propS13P[j] = mom1[j]-mom3[j];
      propS14P[j] = mom1[j]-mom4[j];
      propS15P[j] = mom1[j]-mom5[j];
      propS23P[j] = mom2[j]-mom3[j];
      propS24P[j] = mom2[j]-mom4[j];
      propS25P[j] = mom2[j]-mom5[j];
    }
    pDenS13 = prop.denominator(propS13P);
    pDenS14 = prop.denominator(propS14P);
    pDenS15 = prop.denominator(propS15P);
    pDenS23 = prop.denominator(propS23P);
    pDenS24 = prop.denominator(propS24P);
    pDenS25 = prop.denominator(propS25P);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eeAAA::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5){
    cdouble one(1,0);
    //No sign changes due to p3 and p4 being outgoing.
    if(ds3>0&&ds4>0&&ds5>0){
      //<12>([45]^2[3|p_1p_2|3](*2 denominators) + [35]^2)
      return e*e*e*me*a12a.v(ds1,ds2)*(
        s45s.v()*s45s.v()*s3123s.v()/pDenS13/pDenS23*((one/pDenS24/pDenS25)-(one/pDenS14/pDenS15)) +
        s35s.v()*s35s.v()*s4124s.v()/pDenS14/pDenS24*((one/pDenS23/pDenS25)-(one/pDenS13/pDenS15)) +
        s34s.v()*s34s.v()*s5125s.v()/pDenS15/pDenS25*((one/pDenS23/pDenS24)-(one/pDenS13/pDenS14))
      );
    }
    // else if(ds3<0&&ds4<0&&ds5<0){
    //   //me<34>^2[12]
    //   return 2.0*e*e*me*a34a.v()*a34a.v()*s12s.v(ds1,ds2)/pDenS13/pDenS14;
    // }
    // else if(ds3>0&&ds4<0){
    //   //([13]<24>+[23]<14>)[314>
    //   return 2.0*e*e*(s13s.v(ds1)*a24a.v(ds2)+s23s.v(ds2)*a14a.v(ds1))*s314a.v()/pDenT/pDenU;
    // }
    // else if(ds3<0&&ds4>0){
    //   //(<13>[24]+<23>[14])*[413>
    //   return 2.0*e*e*(a13a.v(ds1)*s24s.v(ds2)+a23a.v(ds2)*s14s.v(ds1))*s413a.v()/pDenT/pDenU;
    // }
    return cdouble(0,0);    
  }
  //set_momenta(...) must be called before amp2().
  ldouble eeAAA::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-2;j4<=2;j4+=4)
    for(int j5=-2;j5<=2;j5+=4){
	    M = amp(j1,j2,j3,j4,j5);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Symmetry factor for identical photons 1/6
    return amp2/24.0;
  }


  



  //  Tests
  int test_eeAAA(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> A , A , A   :";
    // {//amp^2
    //   int i=0;
    //   // me=0.0005, pspatial=250
    //   ldouble me=0.0005;
    //   ldouble EE=0.31333;
    //   eeAA eeAAAmp = eeAA(EE,me);
    //   ldouble pspatial=250;
    //   ldouble dataCH[20] = {3.761473098865852E-01,1.196559098890515E-01,6.884618493553914E-02,4.748300512537967E-02,3.599742458224401E-02,2.906647080620063E-02,2.465909507168581E-02,2.184718935306326E-02,2.016436086672569E-02,1.937355800660445E-02,1.937355800660445E-02,2.016436086672569E-02,2.184718935306325E-02,2.465909507168581E-02,2.906647080620063E-02,3.599742458224402E-02,4.748300512537965E-02,6.884618493553911E-02,1.196559098890515E-01,3.761473098865843E-01};
    //   i += eeAAAmp.test_2to2_amp2([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
    //   i += eeAAAmp.test_2to2_amp2_rotations([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
    //   i += eeAAAmp.test_2to2_amp2_boosts([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
    //   i += eeAAAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
    //   //Close to threshold
    //   pspatial = 0.0001;
    //   ldouble dataCH2[20] = {2.081251342993773E-02,2.079750983522615E-02,2.078115831044515E-02,2.076461266895772E-02,2.074882438042710E-02,2.073456777596373E-02,2.072246020291247E-02,2.071297785405228E-02,2.070646782725441E-02,2.070315683100089E-02,2.070315683100089E-02,2.070646782725441E-02,2.071297785405228E-02,2.072246020291247E-02,2.073456777596374E-02,2.074882438042710E-02,2.076461266895772E-02,2.078115831044515E-02,2.079750983522615E-02,2.081251342993773E-02};
    //   i += eeAAAmp.test_2to2_amp2([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
    //   i += eeAAAmp.test_2to2_amp2_rotations([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
    //   i += eeAAAmp.test_2to2_amp2_boosts([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
    //   i += eeAAAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
    //   // Done
    //   if(i==0) std::cout<<"                                         Pass"<<std::endl;
    //   else std::cout<<"                                         Fail!"<<std::endl;
    //   n+=i;
    // }
    
    return n;
  }

  
  

}
