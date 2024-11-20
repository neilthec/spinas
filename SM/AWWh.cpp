
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

//File:  SPINAS/SM/AWWh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AWWh.h"

namespace spinas {

  AWWh::AWWh(const ldouble& echarge, const ldouble& massh, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), mh(massh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WW(widthW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    propW = propagator(MW,WW);
    p1=particle(0);
    p2=particle(MW);
    p3=particle(MW);
    p4=particle(mh);
    //<23>,[23],<12>,[12],<13>,[13]
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //[143>,[142>,[341>,[241>
    s143a = sproduct(SQUARE,&p1,&p4,&p3);
    s142a = sproduct(SQUARE,&p1,&p4,&p2);
    s341a = sproduct(SQUARE,&p3,&p4,&p1);
    s241a = sproduct(SQUARE,&p2,&p4,&p1);
    //Couplings
    pre = sqrt2*e*e/(MW*SW);
  }
  void AWWh::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    p2.set_mass(MW);
    p3.set_mass(MW);
    p4.set_mass(mh);
    propW.set_mass(MW);
    //Couplings
    pre = sqrt2*e*e/(MW*SW);
  }
  void AWWh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<23>,[23],<12>,[12],<13>,[13]
    s23s.update();
    a23a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    //[143>,[142>,[341>,[241>
    s143a.update();
    s142a.update();
    s341a.update();
    s241a.update();
    //Propagator Momentum
    ldouble propTP[4], propSP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propSP[j] = mom1[j]+mom2[j];
    }
    pDenT=propW.denominator(propTP);
    pDenS=propW.denominator(propSP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AWWh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1;
    int ds3a, ds3b, ds2a, ds2b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds2);
    ldouble normFactor=get_spin_normalization(ds3,ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds2,ds2a,ds2b, i);

      //pre = sqrt2*e*e/(2.0*MW*SW);

      if(ds1>0){
      
	//T&U-Channel W
	//AhW-W+ all ingoing:
	//+pre <34>( Mh^2[13][14] - MW([14][123>+[13][124>) )/(t-MW^2)(u-MW^2)
	//AW+W-h: 2<->4:
	//-pre <23>( Mh^2[13][12] - MW([12][143>+[13][142>) )/(t-MW^2)(s-MW^2)
	//34 out:
	//+pre <23>( Mh^2[13][12] - MW([12][143>-[13][142>) )/(t-MW^2)(s-MW^2)
	amplitude += normFactor*pre*a23a.v(ds2a,ds3a)*( mh*mh*s13s.v(ds3b)*s12s.v(ds2b)
							-MW*(s12s.v(ds2b)*s143a.v(ds3b)-s13s.v(ds3b)*s142a.v(ds2b))
							)/pDenT/pDenS;
      }
      else if(ds1<0){
	amplitude += normFactor*pre*s23s.v(ds2a,ds3a)*( mh*mh*a13a.v(ds3b)*a12a.v(ds2b)
							-MW*(a12a.v(ds2b)*s341a.v(ds3b)-a13a.v(ds3b)*s241a.v(ds2b))
							)/pDenT/pDenS;
      }
      else{
	std::cout<<"Photon cannot have helicity 0!\n";
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AWWh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2*1/3
    return amp2/6.0;
  }

  //set_momenta(...) must be called before amp2().
  ldouble AWWh::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-2;j3<=2;j3+=2){
	M = amp(2,j2,j3);
	amp2 += std::pow(std::abs(M),2);
      }
    //Average over initial spins 1/3
    return amp2/3.0;
  }


  



  //  Tests
  int test_AWWh(){
    int n=0;//Number of fails
    std::cout<<"\t* A , W+ -> W+, h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      AWWh AWWhAmp = AWWh(EE,mh,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {5.409950341836966E+00,1.329914120081493E+00,5.546257118435115E-01,2.888939917465141E-01,1.700187300456206E-01,1.079660452497632E-01,7.212301962553640E-02,4.988189533516871E-02,3.532679473574713E-02,2.540532426450842E-02,1.842435400729740E-02,1.338962125525697E-02,9.689734056125617E-03,6.934185585332866E-03,4.865494240076958E-03,3.309247882421935E-03,2.144540593111295E-03,1.285865826629810E-03,6.716904331315391E-04,2.570534586040455E-04};
      i += AWWhAmp.test_2to2_amp2([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH);
      i += AWWhAmp.test_2to2_amp2_rotations([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH);
      i += AWWhAmp.test_2to2_amp2_boosts([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH);
      i += AWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH);
      i += AWWhAmp.test_2to2_amp2([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH);
      i += AWWhAmp.test_2to2_amp2_rotations([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH);
      i += AWWhAmp.test_2to2_amp2_boosts([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH);
      i += AWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {6.646507979456209E-01,4.008845448080308E-01,2.624516177347133E-01,1.816744810264004E-01,1.308856506457957E-01,9.712731255258195E-02,7.370259204440020E-02,5.688129244465427E-02,4.445953908264935E-02,3.507171239720982E-02,2.783768159929141E-02,2.217117555313992E-02,1.767053958690804E-02,1.405385792828710E-02,1.111904604496788E-02,8.718541716428493E-03,6.742819042791249E-03,5.109392471970163E-03,3.755326400885328E-03,2.632035460161731E-03};
      i += AWWhAmp.test_2to2_amp2([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH2);
      i += AWWhAmp.test_2to2_amp2_rotations([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH2);
      i += AWWhAmp.test_2to2_amp2_boosts([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH2);
      i += AWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH2);
      i += AWWhAmp.test_2to2_amp2([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH2);
      i += AWWhAmp.test_2to2_amp2_rotations([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH2);
      i += AWWhAmp.test_2to2_amp2_boosts([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH2);
      i += AWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=90\n";
      pspatial = 90;
      ldouble dataCH4[20] = {9.021736702755855E-02,8.207758443429436E-02,7.482435990476147E-02,6.833675210999440E-02,6.251350625024814E-02,5.726935918310927E-02,5.253212564235061E-02,4.824038285815122E-02,4.434161751221099E-02,4.079073274406904E-02,3.754883764982127E-02,3.458225998075762E-02,3.186173636402718E-02,2.936174459975845E-02,2.705995034006038E-02,2.493674637054253E-02,2.297486726111448E-02,2.115906567015370E-02,1.947583932480323E-02,1.791319984546122E-02};
      i += AWWhAmp.test_2to2_amp2([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH4);
      i += AWWhAmp.test_2to2_amp2_rotations([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH4);
      i += AWWhAmp.test_2to2_amp2_boosts([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH4);
      i += AWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH4);
      i += AWWhAmp.test_2to2_amp2([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH4);
      i += AWWhAmp.test_2to2_amp2_rotations([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH4);
      i += AWWhAmp.test_2to2_amp2_boosts([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH4);
      i += AWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH4);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=95\n";
      mh = 0.0005;
      pspatial = 95;
      AWWhAmp.set_masses(mh,MW);
      ldouble dataCH3[20] = {4.414111142628276E-01,2.606518348282668E-01,1.693046719790523E-01,1.170451666617239E-01,8.447263273817678E-02,6.285118042043714E-02,4.778930285950597E-02,3.688968278998235E-02,2.875687475478902E-02,2.253572244784899E-02,1.767958226773389E-02,1.382648487729594E-02,1.072947261513810E-02,8.215673432081875E-03,6.161347547401282E-03,4.476175976011579E-03,3.093086108520585E-03,1.961497629236441E-03,1.042739104534081E-03,3.068757151444159E-04};
      i += AWWhAmp.test_2to2_amp2([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH3);
      i += AWWhAmp.test_2to2_amp2_rotations([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH3);
      i += AWWhAmp.test_2to2_amp2_boosts([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH3);
      i += AWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWWhAmp.amp2(); }, 0,MW,MW,mh,pspatial,dataCH3);
      i += AWWhAmp.test_2to2_amp2([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH3);
      i += AWWhAmp.test_2to2_amp2_rotations([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH3);
      i += AWWhAmp.test_2to2_amp2_boosts([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH3);
      i += AWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWWhAmp.amp2_Aplus(); }, 0,MW,MW,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }
  
  

}
