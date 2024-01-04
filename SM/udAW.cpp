
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
    constexpr ldouble sqrt2 = std::sqrt(2);
    propW = propagator(MW,WW);
    propu = propagator(mu,0);
    propd = propagator(md,0);
    p1=particle(mu);
    p2=particle(md);
    p3=particle(0);
    p4=particle(MW);
    //<14>,[34],[23],[312>,[314>
    a14a = sproduct(ANGLE,&p1,&p4);
    s34s = sproduct(SQUARE,&p3,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    //[24],<34>,<13>,[123>,[413>
    s24s = sproduct(SQUARE,&p2,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    a13a = sproduct(ANGLE,&p1,&p3);
    s123a = sproduct(SQUARE,&p1,&p2,&p3);
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
    //prefactor
    pre = sqrt2*e*e/(MW*SW);
  }
  void udAW::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massW){
    constexpr ldouble sqrt2 = std::sqrt(2);
    mu=massu;
    md=massd;
    MW=massW;
    p1.set_mass(mu);
    p2.set_mass(md);
    p4.set_mass(MW);
    propW.set_mass(MW);
    propu.set_mass(mu);
    propd.set_mass(md);
    pre = sqrt2*e*e/(MW*SW);
  }
  void udAW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    constexpr ldouble one=1;
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<14>,[34],[23],[312>,[314>
    a14a.update();
    s34s.update();
    s23s.update();
    s312a.update();
    s314a.update();
    //[24],<34>,<13>,[123>,[413>
    s24s.update();
    a34a.update();
    a13a.update();
    s123a.update();
    s413a.update();
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
    prop=(one/pDenT/pDenU);
    std::cout<<"pDenS="<<pDenS<<std::endl;
    std::cout<<"pDenT="<<pDenT<<std::endl;
    std::cout<<"pDenU="<<pDenU<<std::endl;
    std::cout<<"prop="<<prop<<std::endl;
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
	//uDAW- all ingoing:
	// -<14>([34](mu^2[23]+md[312>)-MW[23][314>)/(s-MW^2)(t-mu^2)
	//34 out:
	// +<14>([34](mu^2[23]+md[312>)+MW[23][314>)/(s-MW^2)(t-mu^2)
	amplitude += normFactor*pre*a14a.v(ds1,ds4a)*(s34s.v(ds4b)*(mu*mu*s23s.v(ds2)+md*s312a.v(ds2))+MW*s23s.v(ds2)*s314a.v(ds4b))*prop;
	
	
      }
      else if(ds3<0){

	//ST Diagram
	//pre = sqrt(2)*e*e/(MW*SW);
	//uDAW- all ingoing:
	// [24](<34>(<13>(md^2-MW^2)+mu[123>)+MW<13>[413>))/(s-MW^2)(t-mu^2)
	//34 out:
	// -[24](<34>(<13>(md^2-MW^2)+mu[123>)-MW<13>[413>))/(s-MW^2)(t-mu^2)
	amplitude += normFactor*pre*s24s.v(ds2,ds4a)*(-a34a.v(ds4b)*((md*md-MW*MW)*a13a.v(ds1)+mu*s123a.v(ds1))+MW*a13a.v(ds1)*s413a.v(ds4b))*prop;

	
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
    std::cout<<"\t* u , D  -> A , W+      :";
    {//amp^2
      int i=0;
      std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,mu=0.0042,md=0.0075,MW=80.385, SW=0.474;
      udAW udAWAmp = udAW(EE,mu,md,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.214654980030611E-01,3.303813892486947E-02,1.602908672060694E-02,9.166879475001317E-03,5.644683815138391E-03,3.607308528597372E-03,2.342630731371288E-03,1.520039877699170E-03,9.666760778313188E-04,5.853276672231173E-04,3.197725820336874E-04,1.390818140263917E-04,3.102122214560715E-05,1.393593490470620E-06,8.001725528131994E-05,3.396053973849570E-04,9.506588566945047E-04,2.371166667327498E-03,6.298284346237135E-03,2.804625829997848E-02};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH);
      //i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH);
      //i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH);
      //i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH);
      std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=81\n";
      pspatial = 81;
      ldouble dataCH2[20] = {2.195590649454362E-01,6.234149939547436E-02,3.160711108242041E-02,1.889022591504900E-02,1.214106131876322E-02,8.075468892034155E-03,5.432814568674048E-03,3.627835710525935E-03,2.354527452948571E-03,1.440801046252570E-03,7.871294971062667E-04,3.387607918489249E-04,7.403746692304418E-05,3.231893477082901E-06,1.791299024366945E-04,7.304518650379627E-04,1.959026526826431E-03,4.675607995169821E-03,1.188458210830164E-02,5.069595724728110E-02};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH2);
      //i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH2);
      //i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH2);
      //i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH2);
      std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=250\n";
      mu=80;
      udAWAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH3[20] = {8.256818546187714E-02,3.472084943584747E-02,1.931092877135298E-02,1.196511353034598E-02,7.811120820899880E-03,5.226644501488082E-03,3.518395397884914E-03,2.342133124669676E-03,1.509670464272961E-03,9.121536502018401E-04,4.859420217252450E-04,1.968041592758696E-04,3.364968127686558E-05,9.270232079551694E-06,1.699535814176212E-04,6.234609844991037E-04,1.619158981270846E-03,3.824055043236203E-03,9.706752124639347E-03,4.153622145759785E-02};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH3);
      //i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH3);
      //i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH3);
      //i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH3);
      std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {4.437373833746103E-02,5.185938535997558E-02,6.025638587082559E-02,6.973045569826253E-02,8.049007345173763E-02,9.280122584919707E-02,1.070087062955788E-01,1.235676363956481E-01,1.430914328034084E-01,1.664271824766789E-01,1.947786210080023E-01,2.299159144048951E-01,2.745532558121124E-01,3.330749879356295E-01,4.130537915554323E-01,5.287929146779445E-01,7.109399818420717E-01,1.039169622811190E+00,1.805555797014838E+00,5.636576504334291E+00};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH4);
      //i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH4);
      //i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH4);
      //i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH4);
      std::cout<<"\n# mu=80, md=0.0075, MW=1, pspatial=250\n";
      MW=1;
      udAWAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH5[20] = {1.838075074003835E+02,8.113253514570728E+01,4.759092543407630E+01,3.102351737893039E+01,2.120874076396370E+01,1.477155686710830E+01,1.027654476756248E+01,7.013183420784458E+00,4.593874539626313E+00,2.794737683574226E+00,1.484979365592721E+00,5.942151212721001E-01,9.947499468705701E-02,2.660351071631339E-02,4.697952607650677E-01,1.648581881427617E+00,4.071126665300199E+00,9.096790693768298E+00,2.175670939799966E+01,8.743638154857922E+01};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH5);
      //i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH5);
      //i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH5);
      //i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH5);
      std::cout<<"\n# mu=80, md=0.0075, MW=1, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH6[20] = {9.305785701150859E+01,9.850421223466601E+01,1.045717424204576E+02,1.113733797148777E+02,1.190511996427738E+02,1.277864703626617E+02,1.378141681687003E+02,1.494444602847809E+02,1.630953958127214E+02,1.793442767472553E+02,1.990114729307221E+02,2.233033986467391E+02,2.540698664167933E+02,2.942989883893278E+02,3.491519439439225E+02,4.283773555734648E+02,5.528642791164643E+02,7.769206909771283E+02,1.299650063935370E+03,3.912133494177679E+03};
      i += udAWAmp.test_2to2_amp2([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH6);
      //i += udAWAmp.test_2to2_amp2_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH6);
      //i += udAWAmp.test_2to2_amp2_boosts([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH6);
      //i += udAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udAWAmp.amp2(); }, mu,md,0,MW,pspatial,dataCH6);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }



}
