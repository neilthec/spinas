
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

//File:  SPINAS/SM/eeAZ.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eeAZ.h"

namespace spinas {

  eeAZ::eeAZ(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW):
    e(echarge), me(masse), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), prope(masse,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    p1=particle(me);
    p2=particle(me);
    p3=particle(0);
    p4=particle(MZ);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    //[3123]
    s3123s = sproduct(SQUARE,&p3,&p1,&p2,&p3);
    //Couplings
    preTU = sqrt2*e*e/(2.0*MW*me*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
  }
  void eeAZ::set_masses(const ldouble& masse, const ldouble& massW){
    me=masse;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(me);
    p2.set_mass(me);
    p3.set_mass(0);
    p4.set_mass(MZ);
    prope.set_mass(me);
    //Couplings
    preTU = sqrt2*e*e/(2.0*MW*me*SW);
  }
  void eeAZ::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s34s.update();
    a34a.update();
    s24s.update();
    a24a.update();
    s14s.update();
    a14a.update();
    //[3123]
    s3123s.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT=prope.den(propTP);
    pDenU=prope.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eeAZ::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds4a, ds4b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=1;
    ldouble combFactor=1;
    if(ds4<0){
      ds4a=-1;
      ds4b=-1;
    }
    else if(ds4>0){
      ds4a=1;
      ds4b=1;
    }
    else {
      nCombs*=2;
      combFactor/=sqrt2;
    }
    //Start the loop
    for(int i=0;i<nCombs;i++){
      if(ds4==0){
	if(i==0||i==1){
	  ds4a=1;
	  ds4b=-1;
	}
	else{
	  ds4a=-1;
	  ds4b=1;
	}
      }
      
      if(ds3>0){
	//-preTU [3123] (gRe <24> (me [14]+[431>)+gLe me [24] <14>))/((t-me^2) (u-me^2))
	

      }
      else if(ds3<0){
	//

      }


      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eeAZ::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2^2=1/4
    return amp2/4.0;
  }
  



  //  Tests
  int test_eeAZ(){
    int n=0;//Number of fails
    std::cout<<"\t* e, E -> A, Z      :";
    {//amp^2
      int i=0;
      std::cout<<"\n# me=0.0005, MW=80.385, pspatial=250\n";
      ldouble me=0.0005;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      eeAZ eeAZAmp = eeAZ(EE,me,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {5.442211901086163E-02,1.868054815261074E-02,1.105414103931645E-02,7.790196245354712E-03,6.021116710397218E-03,4.948703459038923E-03,4.264845459568338E-03,3.827764750623053E-03,3.565896162683383E-03,3.442762313088820E-03,3.442762313088828E-03,3.565896162683399E-03,3.827764750623049E-03,4.264845459568340E-03,4.948703459038917E-03,6.021116710397213E-03,7.790196245354721E-03,1.105414103931645E-02,1.868054815261074E-02,5.442211901086155E-02};
      i += eeAZAmp.test_2to2_amp2([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH);
      //i += eeAZAmp.test_2to2_amp2_rotations([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH);
      //i += eeAZAmp.test_2to2_amp2_boosts([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH);
      //i += eeAZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH);
      std::cout<<"\n# me=0.0005, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {3.258950763116162E-02,2.142583969979962E-02,1.588029093328496E-02,1.273679122778556E-02,1.077465036768588E-02,9.482225318953104E-03,8.613595121840859E-03,8.039190134080599E-03,7.687590711151545E-03,7.520289759998401E-03,7.520289759998404E-03,7.687590711151542E-03,8.039190134080595E-03,8.613595121840855E-03,9.482225318953102E-03,1.077465036768588E-02,1.273679122778556E-02,1.588029093328495E-02,2.142583969979962E-02,3.258950763116161E-02};
      i += eeAZAmp.test_2to2_amp2([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH2);
      //i += eeAZAmp.test_2to2_amp2_rotations([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH2);
      //i += eeAZAmp.test_2to2_amp2_boosts([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH2);
      //i += eeAZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH2);
      std::cout<<"\n# me=125.1, MW=80.385, pspatial=95\n";
      me = 125.1;
      pspatial = 95;
      eeAZAmp.set_masses(me,MW);
      ldouble dataCH4[20] = {2.810197124203936E-02,3.315269028042438E-02,3.679491077716625E-02,3.948289177375699E-02,4.148957626322836E-02,4.298641105831872E-02,4.408424142740629E-02,4.485550813579095E-02,4.534677403871230E-02,4.558588235330223E-02,4.558588235330224E-02,4.534677403871229E-02,4.485550813579094E-02,4.408424142740630E-02,4.298641105831873E-02,4.148957626322837E-02,3.948289177375699E-02,3.679491077716624E-02,3.315269028042439E-02,2.810197124203937E-02};
      i += eeAZAmp.test_2to2_amp2([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH4);
      //i += eeAZAmp.test_2to2_amp2_rotations([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH4);
      //i += eeAZAmp.test_2to2_amp2_boosts([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH4);
      //i += eeAZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH4);
      std::cout<<"\n# me=125, MW=80.385, pspatial=125.1\n";
      me = 125;
      pspatial = 125.1;
      eeAZAmp.set_masses(me,MW);
      ldouble dataCH3[20] = {4.857259033804483E-02,5.746120555747627E-02,6.320283428763762E-02,6.712514381881783E-02,6.989476120200203E-02,7.187849311291551E-02,7.329110139627615E-02,7.426291472011433E-02,7.487330179705647E-02,7.516800458681017E-02,7.516800458681017E-02,7.487330179705645E-02,7.426291472011431E-02,7.329110139627615E-02,7.187849311291551E-02,6.989476120200204E-02,6.712514381881783E-02,6.320283428763762E-02,5.746120555747631E-02,4.857259033804481E-02};
      i += eeAZAmp.test_2to2_amp2([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH3);
      //i += eeAZAmp.test_2to2_amp2_rotations([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH3);
      //i += eeAZAmp.test_2to2_amp2_boosts([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH3);
      //i += eeAZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAZAmp.amp2(); }, me,me,MZ,MZ,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                             Pass"<<std::endl;
      else std::cout<<"                                             Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
