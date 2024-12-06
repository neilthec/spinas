
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

//File:  SPINAS/SM/hnWe.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/hnWe.h"

namespace spinas {
  
  hnWe::hnWe(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), me(masse), mh(massh), MW(massW), WW(widthW), SW(sinW), prope(masse,0), propW(massW,widthW) {
    p1=particle(mh);
    p2=particle(0);
    p3=particle(MW);
    p4=particle(me);
    //<23>, [43], <42>, [313>, [314>
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s43s = sproduct(SQUARE,&p4,&p3,2);
    a42a = sproduct(ANGLE,&p4,&p2,2);
    s313a = sproduct(SQUARE,&p3,&p1,&p3,2);
    s314a = sproduct(SQUARE,&p3,&p1,&p4,2);
    //Couplings
    preW = e*e/(2.0*MW*MW*SW*SW);
    pree = me*preW;
  }
  void hnWe::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    mh=massh;
    MW=massW;
    p1.set_mass(mh);
    p3.set_mass(MW);
    p4.set_mass(me);
    prope.set_mass(me);
    propW.set_mass(MW);
    //Couplings
    preW = e*e/(2.0*MW*MW*SW*SW);
    pree = me*preW;
  }
  void hnWe::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);    
    //<23>, [13], <12>, [343>, [341>
    a23a.update();
    s43s.update();
    a42a.update();
    s313a.update();
    s314a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      //propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenWT=propW.denominator(propTP);
    pDeneU=prope.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble hnWe::amp(const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;
    constexpr ldouble two=2;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      //S-Channel W
      //preW = e*e/(2.0*MW*MW*SW*SW);
      //EnW-h all in:
      // - preW (2MW^2 <23>[13]-Me<12>[343>)/(s-MW^2)
      //hnW-E: 1<->4
      // - preW (2MW^2 <23>[43]-Me<42>[313>)/(t-MW^2)
      //34 out:
      // + preW (2MW^2 <23>[43]+Me<42>[313>)/(t-MW^2)
      amplitude += normFactor*preW*(two*MW*MW*a23a.v(ds3a)*s43s.v(ds4,ds3b)+me*a42a.v(ds4)*s313a.v(ds3a,ds3b))/pDenWT;
      
      //U-Channel e
      //pree = e*e*me/(2.0*MW*MW*SW*SW);
      //EnW-h all ingoing:
      // - pree <23>(2Me[13]+[341>)/(u-Me^2)
      //hnW-e: 1<->4
      // - pree <23>(2Me[43]+[314>)/(u-Me^2)
      //34 out:
      // + pree <23>(2Me[43]-[314>)/(u-Me^2)
      amplitude += normFactor*pree*a23a.v(ds3a)*(two*me*s43s.v(ds4,ds3b)-s314a.v(ds3b,ds4))/pDeneU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble hnWe::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(j3,j4);
	amp2 += std::pow(std::abs(M),2);
      }
    //
    return amp2;
  }

  
  

  //  Tests
  int test_hnWe(){
    int n=0;//Number of fails
    std::cout<<"\t* h , ne -> W+, e       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=2500\n";
      ldouble me=0.0005, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      hnWe hnWeAmp = hnWe(EE,me,mh,MW,0,SW);
      ldouble pspatial=2500;
      ldouble dataCH[20] = {1.460317872601314E+02,1.560623329585625E+01,5.329351855164090E+00,2.566798665545301E+00,1.459663535964208E+00,9.145128733202343E-01,6.098213166580460E-01,4.242302196653767E-01,3.039314565761930E-01,2.222032599280088E-01,1.646050868573055E-01,1.228044612264239E-01,9.173469828417953E-02,6.817995453110717E-02,5.002484839767410E-02,3.583513827795396E-02,2.461284394169340E-02,1.564775580774758E-02,8.424414654210058E-03,2.562160676237956E-03};
      i += hnWeAmp.test_2to2_amp2([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH);
      i += hnWeAmp.test_2to2_amp2_rotations([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH1[20] = {4.043690666667543E+01,9.611888461189922E+00,4.052020889945705E+00,2.155707370566055E+00,1.300617839675368E+00,8.483353633971600E-01,5.829588839827347E-01,4.154699336713075E-01,3.038924073570483E-01,2.264057739580365E-01,1.707928121803789E-01,1.298008924291866E-01,9.891574932611168E-02,7.521541926536067E-02,5.674535398509362E-02,4.216095678116444E-02,3.051460172890541E-02,2.112421942106937E-02,1.348970799250565E-02,7.238327129418925E-03};
      i += hnWeAmp.test_2to2_amp2([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH1);
      i += hnWeAmp.test_2to2_amp2_rotations([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH1);
      i += hnWeAmp.test_2to2_amp2_boosts([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH1);
      i += hnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH1);
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {8.855629303259658E+00,4.157280056751604E+00,2.349637888580353E+00,1.477562326453237E+00,9.953490859685301E-01,7.031522586062190E-01,5.140666965024892E-01,3.854992526326875E-01,2.946435958766370E-01,2.284254852547406E-01,1.789312222321683E-01,1.411536103837635E-01,1.118054864076977E-01,8.866058175582683E-02,7.017036919549603E-02,5.523256800170041E-02,4.304653844890716E-02,3.302028829052352E-02,2.470912048049240E-02,1.777422769163612E-02};
      i += hnWeAmp.test_2to2_amp2([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH2);
      i += hnWeAmp.test_2to2_amp2_rotations([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH2);
      i += hnWeAmp.test_2to2_amp2_boosts([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH2);
      i += hnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, mh=125, MW=80.385, pspatial=95\n";
      me = 125.1;
      mh = 125;
      pspatial = 95;
      hnWeAmp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {6.706283535772274E+00,4.520821667595869E+00,3.335992096316690E+00,2.621389209745765E+00,2.158125630836939E+00,1.842246112947276E+00,1.619155526685041E+00,1.457947514272191E+00,1.340101010923711E+00,1.254034905765683E+00,1.192293800373036E+00,1.150012135707604E+00,1.124042773020609E+00,1.112453582481667E+00,1.114242700583175E+00,1.129196681948998E+00,1.157856489996697E+00,1.201583146653225E+00,1.262738702861306E+00,1.345028146867413E+00};
      i += hnWeAmp.test_2to2_amp2([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH4);
      i += hnWeAmp.test_2to2_amp2_rotations([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH4);
      i += hnWeAmp.test_2to2_amp2_boosts([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH4);
      i += hnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH4);
      //std::cout<<"\n# me=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      me = 125;
      mh = 0.0005;
      pspatial = 125.1;
      hnWeAmp.set_masses(me,mh,MW);
      ldouble dataCH3[20] = {4.852427353600977E+00,3.662239154969347E+00,2.899388121833132E+00,2.381740794033985E+00,2.015363427894203E+00,1.747790921336854E+00,1.547861346532991E+00,1.396153832059081E+00,1.280104996070471E+00,1.191362318517207E+00,1.124280802818376E+00,1.075038949483560E+00,1.041109278355061E+00,1.020944651709623E+00,1.013807274245253E+00,1.019704801901706E+00,1.039423382843637E+00,1.074670270405519E+00,1.128367545303186E+00,1.205186126286275E+00};
      i += hnWeAmp.test_2to2_amp2([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH3);
      i += hnWeAmp.test_2to2_amp2_rotations([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH3);
      i += hnWeAmp.test_2to2_amp2_boosts([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH3);
      i += hnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH3);
      //std::cout<<"\n# me=125, mh=0.0005, MW=0.0004, pspatial=125.1\n";
      me = 125;
      mh = 0.0005;
      MW = 0.0004;
      pspatial = 125.1;
      hnWeAmp.set_masses(me,mh,MW);
      ldouble dataCH5[20] = {5.928286884026059E+19,6.407513106844272E+19,6.947423960387744E+19,7.558721467297795E+19,8.254575827875514E+19,9.051341449404218E+19,9.969527621382698E+19,1.103513282812634E+20,1.228150715220422E+20,1.375199588383285E+20,1.550376258938821E+20,1.761343375743486E+20,2.018562896506845E+20,2.336619431236146E+20,2.736335443043545E+20,3.248270043921423E+20,3.918741391369385E+20,4.820689525285011E+20,6.074396048554935E+20,7.889784559936135E+20};
      i += hnWeAmp.test_2to2_amp2([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH5);
      i += hnWeAmp.test_2to2_amp2_rotations([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH5);
      i += hnWeAmp.test_2to2_amp2_boosts([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH5);
      i += hnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hnWeAmp.amp2(); }, mh,0,MW,me,pspatial,dataCH5);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}
