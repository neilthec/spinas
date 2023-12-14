
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

				    /*
				      On Linux, a typical compilation of this code would look like the following.
				      Update the -I../include and -L../build according to the location of your SPINAS header files and binaries.
				      g++ -o test test.cpp eemm.cpp -I../include -L../build -lspinas -std=c++11 -DWITH_LONG_DOUBLE -O3
				     */ 
//File:  SPINAS/SM/eemm.cpp

#include <iostream>
#include "spinas.h"
#include "include/eemm.h"



//Constructor
eemm::eemm(const ldouble& echarge, const ldouble& masse, const ldouble& massmu){
  e=echarge;
  me=masse;
  mm=massmu;
  prop=spinas::propagator(0,0);
  p1=spinas::particle(me);
  p2=spinas::particle(me);
  p3=spinas::particle(mm);
  p4=spinas::particle(mm);
  a13a=spinas::sproduct(ANGLE,&p1,&p3);
  s13s=spinas::sproduct(SQUARE,&p1,&p3);
  a14a=spinas::sproduct(ANGLE,&p1,&p4);
  s14s=spinas::sproduct(SQUARE,&p1,&p4);
  a23a=spinas::sproduct(ANGLE,&p2,&p3);
  s23s=spinas::sproduct(SQUARE,&p2,&p3);
  a24a=spinas::sproduct(ANGLE,&p2,&p4);
  s24s=spinas::sproduct(SQUARE,&p2,&p4);
}

void eemm::set_masses(const ldouble& masse, const ldouble& massmu){
  me=masse;
  mm=massmu;
  p1.set_mass(me);
  p2.set_mass(me);
  p3.set_mass(mm);
  p4.set_mass(mm);
}

void eemm::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
  //Particles
  p1.set_momentum(mom1);
  p2.set_momentum(mom2);
  p3.set_momentum(mom3);
  p4.set_momentum(mom4);
  a13a.update();
  s13s.update();
  a14a.update();
  s14s.update();
  a23a.update();
  s23s.update();
  a24a.update();
  s24s.update();
  //Propagator Momentum
  ldouble propP[4];
  for(int j=0;j<4;j++)
    propP[j] = mom1[j]+mom2[j];
  pDenS = prop.den(propP);
}


//Amplitude
//set_momenta(...) must be called before amp(...).
cdouble eemm::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
  constexpr ldouble two = 2.0;
  //Sign changes due to p3 and p4 being outgoing.
  //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
  return -two*e*e*(a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3))/pDenS;
}

//Amplitude Squared
//set_momenta(...) must be called before amp2().
ldouble eemm::amp2(){
  constexpr ldouble four = 4.0;
  ldouble amp2 = 0;
  cdouble M;
  
  //Sum over spins
  for(int ds1=-1;ds1<=1;ds1+=2)
    for(int ds2=-1;ds2<=1;ds2+=2)
      for(int ds3=-1;ds3<=1;ds3+=2)
	for(int ds4=-1;ds4<=1;ds4+=2){
	  M = amp(ds1,ds2,ds3,ds4);
	  amp2 += std::pow(std::abs(M),2);
	}
  //Average over initial spins 1/2*1/2 = 1/4
  return amp2/four;
}
  



//  Tests
int test_eemm(){
  int n=0;//Number of fails
  std::cout<<"\t* e , E  -> m , M  (QED)\n";
  {
    int i=0;
    // me=0.0005, mmu=0.105, pspatial=250
    ldouble e=0.31333, me=0.0005, mmu=0.105;
    eemm eemmAmp = eemm(e,me,mmu);
    ldouble pspatial=250;
    ldouble dataCH[20] = {1.833718152338341E-02,1.660225796905280E-02,1.506010369853670E-02,1.371071871183510E-02,1.255410300894803E-02,1.159025658987546E-02,1.081917945461741E-02,1.024087160317387E-02,9.855333035544847E-03,9.662563751730333E-03,9.662563751730333E-03,9.855333035544846E-03,1.024087160317387E-02,1.081917945461741E-02,1.159025658987546E-02,1.255410300894803E-02,1.371071871183510E-02,1.506010369853669E-02,1.660225796905280E-02,1.833718152338341E-02};
    i += eemmAmp.test_2to2_amp2([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
    i += eemmAmp.test_2to2_amp2_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
    i += eemmAmp.test_2to2_amp2_boosts([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
    i += eemmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
      n+=i;
    }

  
    return n;
}

 
  
  


