
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

//File:  SPINAS/SM/enenee.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/enenee.h"

namespace spinas {
  //Constructors
  enenee::enenee(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& widthW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), MW(massW), WW(widthW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    propW(MW,WW),
    p1(particle(me)), p2(particle(0)),
    p3(particle(0)), p4(particle(me)),
    a13a(sproduct(ANGLE,&p1,&p3)),
    s13s(sproduct(SQUARE,&p1,&p3)),
    a14a(sproduct(ANGLE,&p1,&p4)),
    s14s(sproduct(SQUARE,&p1,&p4)),
    a23a(sproduct(ANGLE,&p2,&p3)),
    s23s(sproduct(SQUARE,&p2,&p3)),
    a24a(sproduct(ANGLE,&p2,&p4)),
    s24s(sproduct(SQUARE,&p2,&p4)),
    s12s(sproduct(SQUARE,&p1,&p2)),
    a12a(sproduct(ANGLE,&p1,&p2)),
    s34s(sproduct(SQUARE,&p3,&p4)),
    a34a(sproduct(ANGLE,&p3,&p4))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preW = e*e/SW/SW/(4.0*MW*MW);
  }
  void enenee::set_masses(const ldouble& masse, const ldouble& massW){
    me=masse;
    MW=massW;
    propW.set_mass(MW);
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(me);
    p4.set_mass(me);
    preW = e*e/SW/SW/(4.0*MW*MW);
  }
  void enenee::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    s12s.update();
    a12a.update();
    s34s.update();
    a34a.update();
    //Propagator Momentum
    ldouble propPU[4], propPT[4];
    for(int j=0;j<4;j++){
      propPU[j] = mom1[j]-mom4[j];
      propPT[j] = mom1[j]-mom3[j];
    }
    pDenUZ = propZ.denominator(propPU);
    pDenTW = propW.denominator(propPT);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble enenee::amp(const int& ds1, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Z Boson
    //Defined above:
    //gL=2.0*SW*SW-1.0;
    //gR=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //eENene all in:
    //+ preZ (gL [23] <14> + gR [13] <24>)/(s-MZ^2)
    //eneNeE: 2<->4
    //- preZ (gL [34] <12> + gR [13] <24>)/(u-MZ^2)
    //34 out:
    //- preZ (gL [34] <12> - gR [13] <24>)/(u-MZ^2)
    amplitude = - two*preZ*(
			    + gL*s34s.v(ds4)*a12a.v(ds1)
			    - gR*s13s.v(ds1)*a24a.v(ds4)
			    )/pDenUZ;

    //W Boson
    //preW = e*e/(4.0*MW*MW*SW*SW);
    //eENene all ingoing:
    //+ preW  ( 2 MW^2 [23] <14> + Me^2 [13] <24> )/(t-MW^2)
    //eneNee: 2<->4
    //- preW  ( 2 MW^2 [34] <12> + Me^2 [13] <24> )/(t-MW^2)
    //34 out:
    //- preW  ( 2 MW^2 [34] <12> - Me^2 [13] <24> )/(t-MW^2)
    amplitude +=  - two*preW*(
			      + 2.0*MW*MW*s34s.v(ds4)*a12a.v(ds1)
			      - me*me*s13s.v(ds1)*a24a.v(ds4)
			      )/pDenTW;

    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble enenee::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2){
	    M = amp(j1,j2);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2
    return amp2/2.0;
  }

  



  //  Tests
  int test_enenee(){
    int n=0;//Number of fails
    std::cout<<"\t* e , ne -> ne , e      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### me=0.0005, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, MW=80.385, WW=0, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      enenee eneneeAmp = enenee(0.31333,me,MW,WW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.561552780270702E+01,8.698838074582023E+00,3.715467757381103E+00,1.990015639906303E+00,1.201869004506706E+00,7.795971357385616E-01,5.284327649681230E-01,3.676057256082515E-01,2.589093415224978E-01,1.825403998035229E-01,1.276284847646253E-01,8.816112449797403E-02,6.123133047919482E-02,4.658326044260976E-02,4.742175098162000E-02,7.366450494797795E-02,1.524503134345560E-01,3.673107607225687E-01,1.057550059613573E+00,4.710923706351793E+00};
      i += eneneeAmp.test_2to2_amp2([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH);
      i += eneneeAmp.test_2to2_amp2_rotations([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH);
      i += eneneeAmp.test_2to2_amp2_boosts([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH);
      i += eneneeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH);
      //std::cout<<"########### me=0.0005, MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {2.083366453828203E-20,2.089102437432628E-20,2.095761774200504E-20,2.103344464131961E-20,2.111850507227135E-20,2.121279903486158E-20,2.131632652909162E-20,2.142908755496282E-20,2.155108211247649E-20,2.168231020163397E-20,2.182277182243657E-20,2.197246697488565E-20,2.213139565898253E-20,2.229955787472853E-20,2.247695362212497E-20,2.266358290117320E-20,2.285944571187455E-20,2.306454205423034E-20,2.327887192824190E-20,2.350243533391055E-20};
      i += eneneeAmp.test_2to2_amp2([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH2);
      i += eneneeAmp.test_2to2_amp2_rotations([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH2);
      i += eneneeAmp.test_2to2_amp2_boosts([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH2);
      i += eneneeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH2);
      //std::cout<<"########### me=0.105, MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      eneneeAmp.set_masses(me,MW);
      ldouble dataCH3[20] = {1.734729776543379E-05,1.734175287058722E-05,1.734381813126238E-05,1.735352362274460E-05,1.737089954731727E-05,1.739597623481794E-05,1.742878414319763E-05,1.746935385908307E-05,1.751771609834189E-05,1.757390170665091E-05,1.763794166006751E-05,1.770986706560395E-05,1.778970916180485E-05,1.787749931932781E-05,1.797326904152699E-05,1.807704996504001E-05,1.818887386037784E-05,1.830877263251798E-05,1.843677832150073E-05,1.857292310302870E-05};
      i += eneneeAmp.test_2to2_amp2([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH3);
      i += eneneeAmp.test_2to2_amp2_rotations([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH3);
      i += eneneeAmp.test_2to2_amp2_boosts([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH3);
      i += eneneeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH3);
      //std::cout<<"########### me=0.105, MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      eneneeAmp.set_masses(me,MW);
      ldouble dataCH4[20] = {7.663985940235369E+00,3.517785565354949E+00,1.968053866034394E+00,1.230312454110716E+00,8.249186459151489E-01,5.799150840417299E-01,4.215487872388805E-01,3.140518640488217E-01,2.384799728316964E-01,1.841543946638129E-01,1.448348300199237E-01,1.168974226939512E-01,9.847500463838942E-02,8.915935618645079E-02,9.018670037035630E-02,1.053294176982803E-01,1.433283113956195E-01,2.245828232134934E-01,4.022047097809696E-01,8.454683078221448E-01};
      i += eneneeAmp.test_2to2_amp2([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH4);
      i += eneneeAmp.test_2to2_amp2_rotations([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH4);
      i += eneneeAmp.test_2to2_amp2_boosts([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH4);
      i += eneneeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eneneeAmp.amp2(); }, me,0,0,me,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
