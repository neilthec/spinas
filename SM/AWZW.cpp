
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

//File:  SPINAS/SM/AWZW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AWZW.h"

namespace spinas {

  AWZW::AWZW(const ldouble& echarge, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WW(widthW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    MZ = MW/CW;
    propW = propagator(MW,WW);
    p1=particle(0);
    p2=particle(MW);
    p3=particle(MZ);
    p4=particle(MW);
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
    pre = 2.0*e*e/(CW*SW);
  }
  void AWZW::set_masses(const ldouble& massW){
    MW=massW;
    MZ=MW/CW;
    p2.set_mass(MW);
    p3.set_mass(MZ);
    p4.set_mass(MW);
    propW.set_mass(MW);
    //Couplings
    pre = 2.0*e*e/(CW*SW);
  }
  void AWZW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
  cdouble AWZW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds4a, ds4b, ds2a, ds2b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4,ds2);
    ldouble normFactor=get_spin_normalization(ds3,ds4,ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, ds2,ds2a,ds2b, i);

      if(ds1>0){
	
	//pre = 2*e*e/CW/SW;
	//AZW-W+ all ingoing:
	//+ pre*( CW^2[12]^2<34>^2 + CW^2[13]^2<24>^2 + CW^2[14]^2<23>^2 + (CW^2-SW^2)[13][14]<23><24> + CW[12][13]<24><34> - CW[12][14]<23><34> )/(t-MW^2)(u-MW^2)
	//AW+ZW-: 4->2->3->4:
	//+ pre*( CW^2[13]^2<24>^2 + CW^2[14]^2<23>^2 + CW^2[12]^2<34>^2 - (CW^2-SW^2)[14][12]<34><23> + CW[13][14]<23><24> + CW[13][12]<34><24> )/(u-MW^2)(s-MW^2)
	//34 out:
	//+ pre*( CW^2[13]^2<24>^2 + CW^2[14]^2<23>^2 + CW^2[12]^2<34>^2 + (CW^2-SW^2)[14][12]<34><23> + CW[13][14]<23><24> - CW[13][12]<34><24> )/(u-MW^2)(s-MW^2)
	amplitude += normFactor*pre*(
				     + CW*CW*s13s.v(ds3a)*s13s.v(ds3b)*a24a.v(ds2a,ds4a)*a24a.v(ds2b,ds4b)
				     + CW*CW*s14s.v(ds4a)*s14s.v(ds4b)*a23a.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)
				     + CW*CW*s12s.v(ds2a)*s12s.v(ds2b)*a34a.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)
				     + (CW*CW-SW*SW)*s14s.v(ds4a)*s12s.v(ds2a)*a34a.v(ds3a,ds4b)*a23a.v(ds2b,ds3b)
				     + CW*s13s.v(ds3a)*s14s.v(ds4a)*a23a.v(ds2a,ds3b)*a24a.v(ds2b,ds4b)
				     - CW*s13s.v(ds3a)*s12s.v(ds2a)*a34a.v(ds3b,ds4a)*a24a.v(ds2b,ds4b)
				     )/pDenU/pDenS;	
	
      }
      else if(ds1<0){
	amplitude += normFactor*pre*(
				     + CW*CW*a13a.v(ds3a)*a13a.v(ds3b)*s24s.v(ds2a,ds4a)*s24s.v(ds2b,ds4b)
				     + CW*CW*a14a.v(ds4a)*a14a.v(ds4b)*s23s.v(ds2a,ds3a)*s23s.v(ds2b,ds3b)
				     + CW*CW*a12a.v(ds2a)*a12a.v(ds2b)*s34s.v(ds3a,ds4a)*s34s.v(ds3b,ds4b)
				     + (CW*CW-SW*SW)*a14a.v(ds4a)*a12a.v(ds2a)*s34s.v(ds3a,ds4b)*s23s.v(ds2b,ds3b)
				     + CW*a13a.v(ds3a)*a14a.v(ds4a)*s23s.v(ds2a,ds3b)*s24s.v(ds2b,ds4b)
				     - CW*a13a.v(ds3a)*a12a.v(ds2a)*s34s.v(ds3b,ds4a)*s24s.v(ds2b,ds4b)
				     )/pDenU/pDenS;
      }
      else{
	std::cout<<"Photon cannot have helicity 0!\n";
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AWZW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2*1/3=1/6
    return amp2/6.0;
  }

  //set_momenta(...) must be called before amp2().
  ldouble AWZW::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-2;j4<=2;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }




  



  //  Tests
  int test_AWZW(){
    int n=0;//Number of fails
    std::cout<<"\t* A , W+ -> Z , W+      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble MZ = MW/std::sqrt(1.0-SW*SW);
      AWZW AWZWAmp = AWZW(EE,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {9.460402761480981E-02,9.557592325821845E-02,9.777740828319822E-02,1.014843915949987E-01,1.070786493996535E-01,1.150961614562064E-01,1.263021199456995E-01,1.418104064265187E-01,1.632796591557312E-01,1.932462844134336E-01,2.357131876413291E-01,2.972410106652790E-01,3.890887910740202E-01,5.317144603611019E-01,7.650934221076464E-01,1.175178859528767E+00,1.972812186704114E+00,3.786114322456728E+00,9.194221377754333E+00,3.889602896192213E+01};
      i += AWZWAmp.test_2to2_amp2([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH);
      i += AWZWAmp.test_2to2_amp2_rotations([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH);
      i += AWZWAmp.test_2to2_amp2_boosts([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH);
      i += AWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH);
      i += AWZWAmp.test_2to2_amp2([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH);
      i += AWZWAmp.test_2to2_amp2_rotations([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH);
      i += AWZWAmp.test_2to2_amp2_boosts([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH);
      i += AWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {9.802570324269670E-02,9.893565462326039E-02,1.008879986760644E-01,1.041256220072035E-01,1.089751947935632E-01,1.158810619337128E-01,1.254555141118484E-01,1.385549023277239E-01,1.563975762483043E-01,1.807514780612333E-01,2.142415937962524E-01,2.608715237985808E-01,3.269445743770044E-01,4.227694104859003E-01,5.660032835718314E-01,7.886762753998552E-01,1.153289340988025E+00,1.794089635048234E+00,3.040200909272057E+00,5.872354550595943E+00};
      i += AWZWAmp.test_2to2_amp2([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH2);
      i += AWZWAmp.test_2to2_amp2_rotations([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH2);
      i += AWZWAmp.test_2to2_amp2_boosts([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH2);
      i += AWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH2);
      i += AWZWAmp.test_2to2_amp2([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH2);
      i += AWZWAmp.test_2to2_amp2_rotations([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH2);
      i += AWZWAmp.test_2to2_amp2_boosts([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH2);
      i += AWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=70\n";
      pspatial = 70;
      ldouble dataCH4[20] = {1.651422729805644E-01,1.700275663118974E-01,1.754842472613868E-01,1.815722326300015E-01,1.883589214974382E-01,1.959202881036823E-01,2.043421607288232E-01,2.137217230936769E-01,2.241692829953306E-01,2.358103631631344E-01,2.487881822651545E-01,2.632666103947099E-01,2.794337042589451E-01,2.975059540623076E-01,3.177334085914486E-01,3.404058897904081E-01,3.658605666155240E-01,3.944912349277522E-01,4.267597522119178E-01,4.632102122520777E-01};
      i += AWZWAmp.test_2to2_amp2([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH4);
      i += AWZWAmp.test_2to2_amp2_rotations([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH4);
      i += AWZWAmp.test_2to2_amp2_boosts([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH4);
      i += AWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWZWAmp.amp2(); }, 0,MW,MZ,MW,pspatial,dataCH4);
      i += AWZWAmp.test_2to2_amp2([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH4);
      i += AWZWAmp.test_2to2_amp2_rotations([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH4);
      i += AWZWAmp.test_2to2_amp2_boosts([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH4);
      i += AWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWZWAmp.amp2_Aplus(); }, 0,MW,MZ,MW,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

}
