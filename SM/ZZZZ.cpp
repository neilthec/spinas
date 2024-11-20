
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

//File:  SPINAS/SM/ZZZZ.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ZZZZ.h"

namespace spinas {

  ZZZZ::ZZZZ(const ldouble& echarge, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW):
    e(echarge), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    proph = propagator(mh,wh);  
    p1=particle(MZ);
    p2=particle(MZ);
    p3=particle(MZ);
    p4=particle(MZ);
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //Couplings
    pre = e*e/(MW*MW*SW*SW);
  }
  void ZZZZ::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(MZ);
    p2.set_mass(MZ);
    p3.set_mass(MZ);
    p4.set_mass(MZ);
    proph.set_mass(mh);
    //Couplings
    pre = e*e/(MW*MW*SW*SW);
  }
  void ZZZZ::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s24s.update();
    a24a.update();
    s34s.update();
    a34a.update();
    s14s.update();
    a14a.update();
    s13s.update();
    a13a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=proph.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=proph.denominator(propTP);
    pDenU=proph.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ZZZZ::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1;
    int ds1a, ds1b, ds2a, ds2b, ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds1,ds2,ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds1,ds2,ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds1,ds1a,ds1b, ds2,ds2a,ds2b, ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);

      //pre = e*e/(MW*MW*SW*SW);

      //all ingoing = 34 outgoing:
      //S-Channel h
      //pre [12]<12>[34]<34>/(s-Mh^2)
      //T-Channel h
      //pre [13]<13>[24]<24>/(t-Mh^2)
      //U-Channel h
      //pre [14]<14>[23]<23>/(u-Mh^2)
      amplitude += normFactor*pre*s12s.v(ds1a,ds2a)*a12a.v(ds1b,ds2b)*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)/pDenS;
      amplitude += normFactor*pre*s13s.v(ds1a,ds3a)*a13a.v(ds1b,ds3b)*s24s.v(ds2a,ds4a)*a24a.v(ds2b,ds4b)/pDenT;
      amplitude += normFactor*pre*s14s.v(ds1a,ds4a)*a14a.v(ds1b,ds4b)*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)/pDenU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ZZZZ::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over spins 1/3^2=1/9
    //Symmetry factor 1/2
    return amp2/18.0;
  }



  //  Tests
  int test_ZZZZ(){
    int n=0;//Number of fails
    std::cout<<"\t* Z , Z  -> Z , Z       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      ZZZZ ZZZZAmp = ZZZZ(EE,mh,0,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {7.941279098231038E-02,8.005258243639021E-02,8.123604544834155E-02,8.214984739214962E-02,8.281214332188951E-02,8.328807260180823E-02,8.362602636291921E-02,8.385761261540871E-02,8.400258525292180E-02,8.407243035102346E-02,8.407243035102357E-02,8.400258525292199E-02,8.385761261540860E-02,8.362602636291930E-02,8.328807260180818E-02,8.281214332188945E-02,8.214984739214974E-02,8.123604544834173E-02,8.005258243639042E-02,7.941279098231062E-02};
      i += ZZZZAmp.test_2to2_amp2([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH);
      i += ZZZZAmp.test_2to2_amp2_rotations([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH);
      i += ZZZZAmp.test_2to2_amp2_boosts([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH);
      i += ZZZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=2500\n";
      pspatial = 2500;
      ldouble dataCH1[20] = {8.719882999897610E-02,8.761082604677561E-02,8.769532622017007E-02,8.773114837491658E-02,8.775048001847185E-02,8.776217025829292E-02,8.776961369880576E-02,8.777436653498683E-02,8.777721241114511E-02,8.777855013592317E-02,8.777855013475033E-02,8.777721240963410E-02,8.777436653598603E-02,8.776961369702407E-02,8.776217025684874E-02,8.775048002004127E-02,8.773114837624885E-02,8.769532622172527E-02,8.761082604851556E-02,8.719882999892725E-02};
      i += ZZZZAmp.test_2to2_amp2([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH1);
      i += ZZZZAmp.test_2to2_amp2_rotations([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH1);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {8.148411866052552E-02,8.032767043695974E-02,7.992784113844699E-02,7.985673181227587E-02,7.992491962546096E-02,8.004252277818570E-02,8.016483744542846E-02,8.026898257665077E-02,8.034315347769726E-02,8.038146049084957E-02,8.038146049084956E-02,8.034315347769726E-02,8.026898257665077E-02,8.016483744542847E-02,8.004252277818569E-02,7.992491962546096E-02,7.985673181227589E-02,7.992784113844699E-02,8.032767043695972E-02,8.148411866052555E-02};
      i += ZZZZAmp.test_2to2_amp2([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH2);
      i += ZZZZAmp.test_2to2_amp2_rotations([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH2);
      i += ZZZZAmp.test_2to2_amp2_boosts([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH2);
      i += ZZZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {1.024987963092666E-01,1.024987961428384E-01,1.024987959949021E-01,1.024987958654579E-01,1.024987957545058E-01,1.024987956620456E-01,1.024987955880775E-01,1.024987955326014E-01,1.024987954956174E-01,1.024987954771253E-01,1.024987954771253E-01,1.024987954956174E-01,1.024987955326014E-01,1.024987955880775E-01,1.024987956620456E-01,1.024987957545058E-01,1.024987958654579E-01,1.024987959949021E-01,1.024987961428384E-01,1.024987963092666E-01};
      i += ZZZZAmp.test_2to2_amp2([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH4);
      i += ZZZZAmp.test_2to2_amp2_rotations([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH4);
      i += ZZZZAmp.test_2to2_amp2_boosts([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH4);
      i += ZZZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZZZAmp.amp2(); }, MZ,MZ,MZ,MZ,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
  
  

}
