
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

//File:  SPINAS/SM/udAW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/udAW.h"

namespace spinas {

  udAW::udAW(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), mu(massu), md(massd), MW(massW), SW(sinW), WW(widthW) {
    //constexpr ldouble sqrt2 = std::sqrt(2);
    propW = propagator(MW,WW);
    propu = propagator(mu,0);
    propd = propagator(md,0);
    p1=particle(mu);
    p2=particle(md);
    p3=particle(0);
    p4=particle(MW);
    //<42>,[41],[34],[31],[3123]
    a42a = sproduct(ANGLE,&p4,&p2,2);
    s41s = sproduct(SQUARE,&p4,&p1,2);
    s34s = sproduct(SQUARE,&p3,&p4,2);
    s31s = sproduct(SQUARE,&p3,&p1,2);
    s3123s = sproduct(SQUARE,&p3,&p1,&p2,&p3,2);
    a3123a = sproduct(ANGLE,&p3,&p1,&p2,&p3,2);
    a34a = sproduct(ANGLE,&p3,&p4,2);
    a32a = sproduct(ANGLE,&p3,&p2,2);
    //prefactor
    pre = std::sqrt(2.0)*e*e/(MW*SW);
  }
  void udAW::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massW){
    //constexpr ldouble sqrt2 = std::sqrt(2);
    mu=massu;
    md=massd;
    MW=massW;
    p1.set_mass(mu);
    p2.set_mass(md);
    p4.set_mass(MW);
    propW.set_mass(MW);
    propu.set_mass(mu);
    propd.set_mass(md);
    pre = std::sqrt(2.0)*e*e/(MW*SW);
  }
  void udAW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    constexpr ldouble one=1, two=2, three=3;
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<42>,[41],[34],[31],[3123]
    a42a.update();
    s41s.update();
    s34s.update();
    s31s.update();
    s3123s.update();
    a3123a.update();
    a34a.update();
    a32a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propW.denominator(propSP);
    pDenT=propu.denominator(propTP);
    pDenU=propd.denominator(propUP);
    prop = (one/three/pDenS/pDenU-two/three/pDenS/pDenT);
    propAp = (one/pDenS+two/three/pDenT);
    propAm = (one/pDenS+one/three/pDenU);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble udAW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds4a, ds4b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds4);
    ldouble normFactor=get_spin_normalization(ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds4,ds4a,ds4b, i);

      if(ds3>0){
      
	//STU Diagram
	//pre = sqrt(2)*e*e/(MW*SW);
	//AW+Ud all ingoing:
	//S W: +    pre<24>[23][1341]/(s-MW^2)/(t-Mu^2)  +     pre[12][13]<24>/(s-MW^2)
	//S W: -    pre<24>[23][1341]/(s-MW^2)/(u-Md^2)  +     pre[12][13]<24>/(s-MW^2)
	//T u: +2/3*pre<24>[23][1341]/(t-Mu^2)/(s-MW^2)  + 2/3*pre[12][13]<24>/(t-Mu^2)
	//U d: -1/3*pre<24>[23][1341]/(u-Md^2)/(s-MW^2)
	//UdAW+: 1<->3 and 2<->4
	//S W: +    pre<42>[41][3123]/(s-MW^2)/(t-Mu^2)  +     pre[34][31]<42>/(s-MW^2)
	//S W: -    pre<42>[41][3123]/(s-MW^2)/(u-Md^2)  +     pre[34][31]<42>/(s-MW^2)
	//T u: +2/3*pre<42>[41][3123]/(t-Mu^2)/(s-MW^2)  + 2/3*pre[34][31]<42>/(t-Mu^2)
	//U d: -1/3*pre<42>[41][3123]/(u-Md^2)/(s-MW^2)
	//34 out:
	//S W: -    pre<42>[41][3123]/(s-MW^2)/(t-Mu^2)  -     pre[34][31]<42>/(s-MW^2)
	//S W: +    pre<42>[41][3123]/(s-MW^2)/(u-Md^2)  -     pre[34][31]<42>/(s-MW^2)
	//T u: -2/3*pre<42>[41][3123]/(t-Mu^2)/(s-MW^2)  - 2/3*pre[34][31]<42>/(t-Mu^2)
	//U d: +1/3*pre<42>[41][3123]/(u-Md^2)/(s-MW^2)
	//  prop = (one/three/pDenS/pDenU-two/three/pDenS/pDenT);
	//  propAp = (one/pDenS+two/three/pDenT);
	amplitude += normFactor*pre*a42a.v(ds4a,ds2)*s41s.v(ds4b,ds1)*s3123s.v()*prop
	  -normFactor*pre*s34s.v(ds4a)*s31s.v(ds1)*a42a.v(ds4b,ds2)*propAp;
	
	
      }
      else if(ds3<0){

	//STU Diagram
	//pre = sqrt(2)*e*e/(MW*SW);
	//AW+Ud all in:
	//S W: -    pre<24>[23]<1341>/(s-MW^2)(u-Md^2)  +     pre<12><14>[23]/(s-MW^2)
	//S W: +    pre<24>[23]<1341>/(s-MW^2)(t-Mu^2)  +     pre<12><14>[23]/(s-MW^2)
	//U d: -1/3*pre<24>[23]<1341>/(u-Md^2)(s-MW^2)  + 1/3*pre<12><14>[23]/(u-Md^2)
	//T u: +2/3*pre<24>[23]<1341>/(t-Mu^2)(s-MW^2)
	//UdAW+: 1<->3 and 2<->4
	//S W: -    pre<42>[41]<3123>/(s-MW^2)(u-Md^2)  +     pre<34><32>[41]/(s-MW^2)
	//S W: +    pre<42>[41]<3123>/(s-MW^2)(t-Mu^2)  +     pre<34><32>[41]/(s-MW^2)
	//U d: -1/3*pre<42>[41]<3123>/(u-Md^2)(s-MW^2)  + 1/3*pre<34><32>[41]/(u-Md^2)
	//T u: +2/3*pre<42>[41]<3123>/(t-Mu^2)(s-MW^2)
	//34 out:
	//S W: +    pre<42>[41]<3123>/(s-MW^2)(u-Md^2)  -     pre<34><32>[41]/(s-MW^2)
	//S W: -    pre<42>[41]<3123>/(s-MW^2)(t-Mu^2)  -     pre<34><32>[41]/(s-MW^2)
	//U d: +1/3*pre<42>[41]<3123>/(u-Md^2)(s-MW^2)  - 1/3*pre<34><32>[41]/(u-Md^2)
	//T u: -2/3*pre<42>[41]<3123>/(t-Mu^2)(s-MW^2)
	//  prop = (one/three/pDenS/pDenU-two/three/pDenS/pDenT);
	//  propAm = (one/pDenS+one/three/pDenU);
	amplitude += normFactor*pre*a42a.v(ds4a,ds2)*s41s.v(ds4b,ds1)*a3123a.v()*prop
	  -normFactor*pre*a34a.v(ds4a)*a32a.v(ds2)*s41s.v(ds4b,ds1)*propAm;

	
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble udAW::amp2(){
    constexpr ldouble three=3;
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    //Color factor 3
	    amp2 += three*std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3^2=1/9
    return amp2/36.0;
  }



  



  //  Tests
  int test_udAW(){
    int n=0;//Number of fails
    std::cout<<"\t* U , d  -> A , W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,mu=0.0042,md=0.0075,MW=80.385, SW=0.474;
      udAW udAWAmp = udAW(EE,mu,md,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.214654980030611E-01,3.303813892486947E-02,1.602908672060694E-02,9.166879475001317E-03,5.644683815138391E-03,3.607308528597372E-03,2.342630731371288E-03,1.520039877699170E-03,9.666760778313188E-04,5.853276672231173E-04,3.197725820336874E-04,1.390818140263917E-04,3.102122214560715E-05,1.393593490470620E-06,8.001725528131994E-05,3.396053973849570E-04,9.506588566945047E-04,2.371166667327498E-03,6.298284346237135E-03,2.804625829997848E-02};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH);
      i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH);
      i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH);
      i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=81\n";
      pspatial = 81;
      ldouble dataCH2[20] = {2.195590649454362E-01,6.234149939547436E-02,3.160711108242041E-02,1.889022591504900E-02,1.214106131876322E-02,8.075468892034155E-03,5.432814568674048E-03,3.627835710525935E-03,2.354527452948571E-03,1.440801046252570E-03,7.871294971062667E-04,3.387607918489249E-04,7.403746692304418E-05,3.231893477082901E-06,1.791299024366945E-04,7.304518650379627E-04,1.959026526826431E-03,4.675607995169821E-03,1.188458210830164E-02,5.069595724728110E-02};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH2);
      i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH2);
      i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH2);
      i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH2);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=250\n";
      mu=80;
      udAWAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH3[20] = {8.256818546187714E-02,3.472084943584747E-02,1.931092877135298E-02,1.196511353034598E-02,7.811120820899880E-03,5.226644501488082E-03,3.518395397884914E-03,2.342133124669676E-03,1.509670464272961E-03,9.121536502018401E-04,4.859420217252450E-04,1.968041592758696E-04,3.364968127686558E-05,9.270232079551694E-06,1.699535814176212E-04,6.234609844991037E-04,1.619158981270846E-03,3.824055043236203E-03,9.706752124639347E-03,4.153622145759785E-02};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH3);
      i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH3);
      i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH3);
      i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH3);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {4.437373833746103E-02,5.185938535997558E-02,6.025638587082559E-02,6.973045569826253E-02,8.049007345173763E-02,9.280122584919707E-02,1.070087062955788E-01,1.235676363956481E-01,1.430914328034084E-01,1.664271824766789E-01,1.947786210080023E-01,2.299159144048951E-01,2.745532558121124E-01,3.330749879356295E-01,4.130537915554323E-01,5.287929146779445E-01,7.109399818420717E-01,1.039169622811190E+00,1.805555797014838E+00,5.636576504334291E+00};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH4);
      i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH4);
      i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH4);
      i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH4);
      //std::cout<<"\n# mu=80, md=40, MW=1, pspatial=25\n";
      mu=80;
      md=40;
      MW=30;
      udAWAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH5[20] = {3.136396109698285E-01,1.373163853761221E-01,7.966453097183732E-02,5.139956504153692E-02,3.483345601025350E-02,2.409868381700405E-02,1.669292147863683E-02,1.137509091089483E-02,7.466926326010541E-03,4.575691565631318E-03,2.470824171295090E-03,1.026505073892502E-03,1.981727052517375E-04,2.326994518801740E-05,6.513528600429410E-04,2.432795863164998E-03,6.168995664035803E-03,1.394696233930737E-02,3.306922513711885E-02,1.177677572534412E-01};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH5);
      i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH5);
      i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH5);
      i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH5);
      //std::cout<<"\n# mu=80, md=0.0075, MW=1, pspatial=1\n";
      mu=80;
      md=0.0075;
      MW=1;
      udAWAmp.set_masses(mu,md,MW);
      pspatial = 1;
      ldouble dataCH6[20] = {9.305785701150859E+01,9.850421223466601E+01,1.045717424204576E+02,1.113733797148777E+02,1.190511996427738E+02,1.277864703626617E+02,1.378141681687003E+02,1.494444602847809E+02,1.630953958127214E+02,1.793442767472553E+02,1.990114729307221E+02,2.233033986467391E+02,2.540698664167933E+02,2.942989883893278E+02,3.491519439439225E+02,4.283773555734648E+02,5.528642791164643E+02,7.769206909771283E+02,1.299650063935370E+03,3.912133494177679E+03};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH6);
      i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH6);
      i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH6);
      i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH6);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }



}
