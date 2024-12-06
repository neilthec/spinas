
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

//File:  SPINAS/SM/uuZZ.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uuZZ.h"

namespace spinas {

  uuZZ::uuZZ(const ldouble& echarge, const ldouble& massf, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), Qf(2.0/3.0), mf(massf), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propf(massf,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    proph = propagator(mh,wh);  
    p1=particle(mf);
    p2=particle(mf);
    p3=particle(MZ);
    p4=particle(MZ);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s34s = sproduct(SQUARE,&p3,&p4,2);
    a34a = sproduct(ANGLE,&p3,&p4,2);
    s24s = sproduct(SQUARE,&p2,&p4,2);
    a24a = sproduct(ANGLE,&p2,&p4,2);
    s14s = sproduct(SQUARE,&p1,&p4,2);
    a14a = sproduct(ANGLE,&p1,&p4,2);
    //[314>,[413>
    s314a = sproduct(SQUARE,&p3,&p1,&p4,2);
    s413a = sproduct(SQUARE,&p4,&p1,&p3,2);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
    gL=-2.0*Qf*SW*SW+1.0;
    gR=-2.0*Qf*SW*SW;
    preh = e*e*mf/(2.0*MW*MW*SW*SW);
  }
  void uuZZ::set_masses(const ldouble& massf, const ldouble& massh, const ldouble& massW){
    mf=massf;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(mf);
    p2.set_mass(mf);
    p3.set_mass(MZ);
    p4.set_mass(MZ);
    propf.set_mass(mf);
    proph.set_mass(mh);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
    preh = e*e*mf/(2.0*MW*MW*SW*SW);
  }
  void uuZZ::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    //[314>,[413>
    s314a.update();
    s413a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=proph.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=propf.denominator(propTP);
    pDenU=propf.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uuZZ::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);
      
      //S-Channel h
      //preh = e*e*mf/(2.0*MW*MW*SW*SW);
      //all ingoing: 
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //34 outgoing:
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      amplitude += normFactor*preh*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))/pDenS;
      
      //T-Channel e
      //preTU = e*e/(2.0*MW*MW*SW*SW);
      //all ingoing:
      //+preTU gLe^2 [24] <13> (MZ <34>+[314>))/(Mf^2-t)
      //+preTU gRe^2 [13] <24> (MZ [34]+[413>))/(Mf^2-t)
      //+preTU gLe gRe Mf ([13] [24] <34>+[34] <13> <24>))/(Mf^2-t)
      //34 outgoing:
      //+preTU gLe^2 [24] <13> (MZ <34>-[314>))/(t-Mf^2)
      //+preTU gRe^2 [13] <24> (MZ [34]-[413>))/(t-Mf^2)
      //-preTU gLe gRe Mf ([13] [24] <34>+[34] <13> <24>))/(t-Mf^2)
      amplitude += normFactor*preTU*gL*gL*s24s.v(ds2,ds4a)*a13a.v(ds1,ds3a)*(MZ*a34a.v(ds3b,ds4b)-s314a.v(ds3b,ds4b))/pDenT
      	+          normFactor*preTU*gR*gR*a24a.v(ds2,ds4a)*s13s.v(ds1,ds3a)*(MZ*s34s.v(ds3b,ds4b)-s413a.v(ds4b,ds3b))/pDenT
      	-          normFactor*preTU*gL*gR*mf*(s13s.v(ds1,ds3a)*s24s.v(ds2,ds4a)*a34a.v(ds3b,ds4b)+a13a.v(ds1,ds3a)*a24a.v(ds2,ds4a)*s34s.v(ds3b,ds4b))/pDenT;

      //U-Channel u
      //preTU = e*e/(2.0*MW*MW*SW*SW);
      //all ingoing:
      //+preTU gLe^2 [23] <14> (MZ <34>-[413>)/(u-Mf^2)
      //+preTU gRe^2 [14] <23> (MZ [34]-[314>)/(u-Mf^2)
      //+preTU gLe gRe Mf ([14] [23] <34>+[34] <14> <23>)/(u-Mf^2)
      //34 outgoing:
      //-preTU gLe^2 [23] <14> (MZ <34>+[413>)/(u-Mf^2)
      //-preTU gRe^2 [14] <23> (MZ [34]+[314>)/(u-Mf^2)
      //+preTU gLe gRe Mf ([14] [23] <34>+[34] <14> <23>)/(u-Mf^2)
      amplitude += -normFactor*preTU*gL*gL*s23s.v(ds2,ds3a)*a14a.v(ds1,ds4a)*(MZ*a34a.v(ds3b,ds4b)+s413a.v(ds4b,ds3b))/pDenU
      	-           normFactor*preTU*gR*gR*s14s.v(ds1,ds4a)*a23a.v(ds2,ds3a)*(MZ*s34s.v(ds3b,ds4b)+s314a.v(ds3b,ds4b))/pDenU
      	+           normFactor*preTU*gL*gR*mf*(s14s.v(ds1,ds4a)*s23s.v(ds2,ds3a)*a34a.v(ds3b,ds4b)+a14a.v(ds1,ds4a)*a23a.v(ds2,ds3a)*s34s.v(ds3b,ds4b))/pDenU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uuZZ::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3^2=1/9
    //Symmetry factor 1/2
    return amp2/72.0;
  }
  



  //  Tests
  int test_uuZZ(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> Z , Z       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mf=0.0042, mh=125, MW=80.385, pspatial=250\n";
      ldouble mf=0.0042, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      uuZZ uuZZAmp = uuZZ(EE,mf,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.400256878537088E-02,1.167147910081938E-02,6.906551941380427E-03,4.867261489726477E-03,3.761952658617381E-03,3.091916171214088E-03,2.664646355852423E-03,2.391561309004831E-03,2.227947600723886E-03,2.151014411042451E-03,2.151014411042445E-03,2.227947600723874E-03,2.391561309004842E-03,2.664646355852452E-03,3.091916171214094E-03,3.761952658617385E-03,4.867261489726467E-03,6.906551941380418E-03,1.167147910081937E-02,3.400256878537077E-02};
      i += uuZZAmp.test_2to2_amp2([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH);
      i += uuZZAmp.test_2to2_amp2_rotations([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH);
      i += uuZZAmp.test_2to2_amp2_boosts([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH);
      i += uuZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH);
      //std::cout<<"\n# mf=0.0042, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {2.036170210852353E-02,1.338671850852420E-02,9.921897459893671E-03,7.957860282903757E-03,6.731928055516502E-03,5.924429699515233E-03,5.381715476473840E-03,5.022831157592048E-03,4.803154236068024E-03,4.698625742785657E-03,4.698625742785656E-03,4.803154236068021E-03,5.022831157592047E-03,5.381715476473843E-03,5.924429699515232E-03,6.731928055516502E-03,7.957860282903753E-03,9.921897459893669E-03,1.338671850852420E-02,2.036170210852353E-02};
      i += uuZZAmp.test_2to2_amp2([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH2);
      i += uuZZAmp.test_2to2_amp2_rotations([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH2);
      i += uuZZAmp.test_2to2_amp2_boosts([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH2);
      i += uuZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH2);
      //std::cout<<"\n# mf=125.1, mh=125, MW=80.385, pspatial=95\n";
      mf = 125.1;
      mh = 125;
      pspatial = 95;
      uuZZAmp.set_masses(mf,mh,MW);
      ldouble dataCH4[20] = {2.033095893871820E-02,1.599669183481019E-02,1.263641091710642E-02,1.001268209007929E-02,7.965593583996088E-03,6.385198285662170E-03,5.195132810350378E-03,4.342618360809686E-03,3.792268941442023E-03,3.522291226690466E-03,3.522291226690467E-03,3.792268941442022E-03,4.342618360809685E-03,5.195132810350380E-03,6.385198285662169E-03,7.965593583996086E-03,1.001268209007928E-02,1.263641091710642E-02,1.599669183481018E-02,2.033095893871819E-02};
      i += uuZZAmp.test_2to2_amp2([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH4);
      i += uuZZAmp.test_2to2_amp2_rotations([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH4);
      i += uuZZAmp.test_2to2_amp2_boosts([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH4);
      i += uuZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH4);
      //std::cout<<"\n# mf=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      mf = 125;
      mh = 0.0005;
      pspatial = 125.1;
      uuZZAmp.set_masses(mf,mh,MW);
      ldouble dataCH3[20] = {2.866841703907730E-02,2.140393199897433E-02,1.625257765559415E-02,1.249157180005816E-02,9.703167017937702E-03,7.632417536056093E-03,6.117930088097471E-03,5.055815718252405E-03,4.380001442244305E-03,4.051252140957807E-03,4.051252140957801E-03,4.380001442244305E-03,5.055815718252402E-03,6.117930088097467E-03,7.632417536056094E-03,9.703167017937700E-03,1.249157180005816E-02,1.625257765559415E-02,2.140393199897433E-02,2.866841703907729E-02};
      i += uuZZAmp.test_2to2_amp2([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH3);
      i += uuZZAmp.test_2to2_amp2_rotations([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH3);
      i += uuZZAmp.test_2to2_amp2_boosts([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH3);
      i += uuZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
