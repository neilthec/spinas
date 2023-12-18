
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

//File:  SPINAS/SM/enmenm.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/enmenm.h"

namespace spinas {
  //Constructors
  enmenm::enmenm(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propZ(MZ,WZ),
    p1(particle(me)), p2(particle(0)),
    p3(particle(me)), p4(particle(0)),
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
  }
  void enmenm::set_masses(const ldouble& masse, const ldouble& massW){
    me=masse;
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(me);
    p3.set_mass(me);
  }
  void enmenm::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]-mom3[j];
    pDenTZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble enmenm::amp(const int& ds1, const int& ds3){
    constexpr ldouble two = 2;
    
    //Z Boson
    //Defined above:
    //gL=2.0*SW*SW-1.0;
    //gR=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //eENn all in:
    //- preZ (gL [23] <14> + gR [13] <24>)/(s-MZ^2)
    //enEN: 4->2->3->4
    //- preZ (gL [34] <12> - gR [14] <23>)/(t-MZ^2)
    //34 out:
    //- preZ (gL [34] <12> + gR [14] <23>)/(t-MZ^2)
    return - two*preZ*( gL*s34s.v(ds3)*a12a.v(ds1) + gR*s14s.v(ds1)*a23a.v(ds3) )/pDenTZ;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble enmenm::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-1;j3<=1;j3+=2){
	    M = amp(j1,j3);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2
    return amp2/2.0;
  }

  



  //  Tests
  int test_enmenm(){
    int n=0;//Number of fails
    std::cout<<"\t* e , nm -> e , nm      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### me=0.0005, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      enmenm enmenmAmp = enmenm(0.31333,me,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {5.776913136679560E+00,1.610170209830663E+00,7.250279185079357E-01,4.030808863203872E-01,2.525410707564306E-01,1.709532127422690E-01,1.222037683338240E-01,9.098357684653040E-02,6.992576154733761E-02,5.514012303190326E-02,4.442069642792457E-02,3.644394348042211E-02,3.037849989831604E-02,2.568191637374954E-02,2.198872916800686E-02,1.904588000063991E-02,1.667394723921192E-02,1.474306318308950E-02,1.315750244957108E-02,1.184555816480042E-02};
      i += enmenmAmp.test_2to2_amp2([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH);
      i += enmenmAmp.test_2to2_amp2_rotations([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH);
      i += enmenmAmp.test_2to2_amp2_boosts([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH);
      i += enmenmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH);
      //std::cout<<"########### me=0.0005, MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {5.091683423895216E-21,4.919491388629454E-21,4.756532885019918E-21,4.602807913065277E-21,4.458316472764204E-21,4.323058564115367E-21,4.197034187117438E-21,4.080243341769087E-21,3.972686028068984E-21,3.874362246015801E-21,3.785271995608207E-21,3.705415276844874E-21,3.634792089724471E-21,3.573402434245670E-21,3.521246310407140E-21,3.478323718207551E-21,3.444634657645577E-21,3.420179128719885E-21,3.404957131429147E-21,3.398968665772033E-21};
      i += enmenmAmp.test_2to2_amp2([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH2);
      i += enmenmAmp.test_2to2_amp2_rotations([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH2);
      i += enmenmAmp.test_2to2_amp2_boosts([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH2);
      i += enmenmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH2);
      //std::cout<<"########### me=0.105, MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      enmenmAmp.set_masses(me,MW);
      ldouble dataCH3[20] = {4.196596050045032E-06,4.049093258551143E-06,3.909581680582814E-06,3.778028781159763E-06,3.654402158799200E-06,3.538669544906312E-06,3.430798803167878E-06,3.330757928949001E-06,3.238515048692930E-06,3.154038419323964E-06,3.077296427653411E-06,3.008257589788597E-06,2.946890550544893E-06,2.893164082860752E-06,2.847047087215736E-06,2.808508591051519E-06,2.777517748195842E-06,2.754043838289409E-06,2.738056266215707E-06,2.729524561533731E-06};
      i += enmenmAmp.test_2to2_amp2([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH3);
      i += enmenmAmp.test_2to2_amp2_rotations([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH3);
      i += enmenmAmp.test_2to2_amp2_boosts([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH3);
      i += enmenmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH3);
      //std::cout<<"########### me=0.105, MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      enmenmAmp.set_masses(me,MW);
      ldouble dataCH4[20] = {1.110557425712026E+00,5.849262056123748E-01,3.552875194221835E-01,2.358360779255416E-01,1.663683579560873E-01,1.227078692698114E-01,9.365272332910785E-02,7.345042399052876E-02,5.890797572359789E-02,4.814160686402221E-02,3.998393936435937E-02,3.368164024545587E-02,2.873177293118895E-02,2.478872522235299E-02,2.160901146858384E-02,1.901733873168859E-02,1.688508462579726E-02,1.511626907526542E-02,1.363818724248903E-02,1.239501900972245E-02};
      i += enmenmAmp.test_2to2_amp2([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH4);
      i += enmenmAmp.test_2to2_amp2_rotations([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH4);
      i += enmenmAmp.test_2to2_amp2_boosts([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH4);
      i += enmenmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enmenmAmp.amp2(); }, me,0,me,0,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
