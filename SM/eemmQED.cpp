
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

//File:  SPINAS/SM/eemmQED.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eemmQED.h"

namespace spinas {
  //Constructors
  eemmQED::eemmQED(const ldouble& echarge, const ldouble& masse, const ldouble& massmu):
    e(echarge), me(masse), mm(massmu), prop(0,0),
    p1(particle(me)), p2(particle(me)),
    p3(particle(mm)), p4(particle(mm)),
    a13a(sproduct(ANGLE,&p1,&p3,2)),
    s13s(sproduct(SQUARE,&p1,&p3,2)),
    a14a(sproduct(ANGLE,&p1,&p4,2)),
    s14s(sproduct(SQUARE,&p1,&p4,2)),
    a23a(sproduct(ANGLE,&p2,&p3,2)),
    s23s(sproduct(SQUARE,&p2,&p3,2)),
    a24a(sproduct(ANGLE,&p2,&p4,2)),
    s24s(sproduct(SQUARE,&p2,&p4,2))
  {}
  void eemmQED::set_masses(const ldouble& masse, const ldouble& massmu){
    me=masse;
    mm=massmu;
    p1.set_mass(me);
    p2.set_mass(me);
    p3.set_mass(mm);
    p4.set_mass(mm);
  }
  void eemmQED::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenS = prop.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eemmQED::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    //Sign changes due to p3 and p4 being outgoing.
    //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    return -2.0*e*e*(a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3))/pDenS;
  }
  //set_momenta(...) must be called before amp2().
  ldouble eemmQED::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    return amp2/4.0;
  }
  //Alternate version that is a bit more efficient but requires more care from the author.
  //No need to call set_momenta(...) before this.
  ldouble eemmQED::amp2(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]+mom2[j];
    ldouble propSD = std::real(prop.denominator(propP));
    

    ldouble amp2 = 0;

    //Calculate the scalar products
    cdouble a13[2][2], s13[2][2], a14[2][2], s14[2][2], a23[2][2], s23[2][2], a24[2][2], s24[2][2];
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
	a13[i][j] = p1.langle(2*i-1)*p3.rangle(2*j-1);
	s13[i][j] = p1.lsquare(2*i-1)*p3.rsquare(2*j-1);
	a14[i][j] = p1.langle(2*i-1)*p4.rangle(2*j-1);
	s14[i][j] = p1.lsquare(2*i-1)*p4.rsquare(2*j-1);
	a23[i][j] = p2.langle(2*i-1)*p3.rangle(2*j-1);
	s23[i][j] = p2.lsquare(2*i-1)*p3.rsquare(2*j-1);
	a24[i][j] = p2.langle(2*i-1)*p4.rangle(2*j-1);
	s24[i][j] = p2.lsquare(2*i-1)*p4.rsquare(2*j-1);
      }
    //Calculate the amp^2
    cdouble M;
    //Sum over spins
    for(int j1=0;j1<2;j1++)
      for(int j2=0;j2<2;j2++)
	for(int j3=0;j3<2;j3++)
	  for(int j4=0;j4<2;j4++){
	    //Sign changes due to p3 and p4 being outgoing.
	    //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
	    M = - (a13[j1][j3]*s24[j2][j4] + a14[j1][j4]*s23[j2][j3] + s13[j1][j3]*a24[j2][j4] + s14[j1][j4]*a23[j2][j3]);

	    amp2 += std::pow(std::abs(M),2);
	  }

    //Multiply by the overall factor including the average over the initial spins (1/4)
    return 4.0*e*e*e*e*amp2/propSD/propSD/4.0;
  }


  



  //  Tests
  int test_eemmQED(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> m , M  (QED):";
    {//amp^2
      int i=0;
      // me=0.0005, mmu=0.105, pspatial=250
      ldouble me=0.0005, mmu=0.105;
      eemmQED eemm = eemmQED(0.31333,me,mmu);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.833718152338341E-02,1.660225796905280E-02,1.506010369853670E-02,1.371071871183510E-02,1.255410300894803E-02,1.159025658987546E-02,1.081917945461741E-02,1.024087160317387E-02,9.855333035544847E-03,9.662563751730333E-03,9.662563751730333E-03,9.855333035544846E-03,1.024087160317387E-02,1.081917945461741E-02,1.159025658987546E-02,1.255410300894803E-02,1.371071871183510E-02,1.506010369853669E-02,1.660225796905280E-02,1.833718152338341E-02};
      i += eemm.test_2to2_amp2([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
      i += eemm.test_2to2_amp2_rotations([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
      i += eemm.test_2to2_amp2_boosts([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
      i += eemm.test_2to2_amp2_boosts_and_rotations([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
      // me=0.0005, mmu=0.105, pspatial=0.11
      pspatial = 0.11;
      ldouble dataCH2[20] = {1.919360703470436E-02,1.903944176218874E-02,1.890240596439707E-02,1.878249964132937E-02,1.867972279298562E-02,1.859407541936583E-02,1.852555752047000E-02,1.847416909629813E-02,1.843991014685021E-02,1.842278067212625E-02,1.842278067212625E-02,1.843991014685021E-02,1.847416909629813E-02,1.852555752047000E-02,1.859407541936583E-02,1.867972279298562E-02,1.878249964132937E-02,1.890240596439707E-02,1.903944176218874E-02,1.919360703470436E-02};
      i += eemm.test_2to2_amp2([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH2);
      i += eemm.test_2to2_amp2_rotations([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH2);
      i += eemm.test_2to2_amp2_boosts([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH2);
      i += eemm.test_2to2_amp2_boosts_and_rotations([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH2);
      // me=0.105, mmu=0.0005, pspatial=0.005
      me=0.105;
      mmu=0.0005;
      pspatial = 0.005;
      eemm.set_masses(me,mmu);
      ldouble dataCH3[20] = {1.927502326937967E-02,1.927109819107943E-02,1.926760923259032E-02,1.926455639391236E-02,1.926193967504553E-02,1.925975907598984E-02,1.925801459674529E-02,1.925670623731187E-02,1.925583399768959E-02,1.925539787787846E-02,1.925539787787846E-02,1.925583399768959E-02,1.925670623731187E-02,1.925801459674529E-02,1.925975907598984E-02,1.926193967504553E-02,1.926455639391236E-02,1.926760923259032E-02,1.927109819107943E-02,1.927502326937967E-02};
      i += eemm.test_2to2_amp2([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH3);
      i += eemm.test_2to2_amp2_rotations([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH3);
      i += eemm.test_2to2_amp2_boosts([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH3);
      i += eemm.test_2to2_amp2_boosts_and_rotations([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH3);
      // me=0.10, mmu=0.105, pspatial=0.05
      me=0.10;
      mmu=0.105;
      pspatial = 0.05;
      eemm.set_masses(me,mmu);
      ldouble dataCH4[20] = {2.605565520464624E-02,2.601471100154132E-02,2.597831615433694E-02,2.594647066303311E-02,2.591917452762983E-02,2.589642774812709E-02,2.587823032452491E-02,2.586458225682326E-02,2.585548354502217E-02,2.585093418912162E-02,2.585093418912162E-02,2.585548354502217E-02,2.586458225682326E-02,2.587823032452491E-02,2.589642774812709E-02,2.591917452762983E-02,2.594647066303311E-02,2.597831615433694E-02,2.601471100154132E-02,2.605565520464624E-02};
      i += eemm.test_2to2_amp2([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH4);
      i += eemm.test_2to2_amp2_rotations([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH4);
      i += eemm.test_2to2_amp2_boosts([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH4);
      i += eemm.test_2to2_amp2_boosts_and_rotations([&]() { return eemm.amp2(); }, me,me,mmu,mmu,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
