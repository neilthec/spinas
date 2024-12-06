
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

//File:  SPINAS/SM/enWh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/enWh.h"

namespace spinas {
  
  enWh::enWh(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), me(masse), mh(massh), MW(massW), WW(widthW), SW(sinW), prope(masse,0), propW(massW,widthW) {
    p1=particle(me);
    p2=particle(0);
    p3=particle(MW);
    p4=particle(mh);
    //<23>, [13], <12>, [343>, [341>
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s343a = sproduct(SQUARE,&p3,&p4,&p3,2);
    s341a = sproduct(SQUARE,&p3,&p4,&p1,2);
    //Couplings
    preW = e*e/(2.0*MW*MW*SW*SW);
    pree = me*preW;
  }
  void enWh::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    mh=massh;
    MW=massW;
    p1.set_mass(me);
    p3.set_mass(MW);
    p4.set_mass(mh);
    prope.set_mass(me);
    propW.set_mass(MW);
    //Couplings
    preW = e*e/(2.0*MW*MW*SW*SW);
    pree = me*preW;
  }
  void enWh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);    
    //<23>, [13], <12>, [343>, [341>
    a23a.update();
    s13s.update();
    a12a.update();
    s343a.update();
    s341a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      //propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenWS=propW.denominator(propSP);
    pDeneU=prope.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble enWh::amp(const int& ds1, const int& ds3){//Double Spin
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
      //all ingoing:
      // - preW (2MW^2 <23>[13]-Me<12>[343>)/(s-MW^2)
      //34 outgoing:
      // + preW (2MW^2 <23>[13]+Me<12>[343>)/(s-MW^2)
      amplitude += normFactor*preW*(two*MW*MW*a23a.v(ds3a)*s13s.v(ds1,ds3b)+me*a12a.v(ds1)*s343a.v(ds3a,ds3b))/pDenWS;
      
      //U-Channel e
      //pree = e*e*me/(2.0*MW*MW*SW*SW);
      //all ingoing:
      // - pree <23>(2Me[13]+[341>)/(u-Me^2)
      //34 outgoing:
      // + pree <23>(2Me[13]-[341>)/(u-Me^2)
      amplitude += normFactor*pree*a23a.v(ds3a)*(two*me*s13s.v(ds1,ds3b)-s341a.v(ds3b,ds1))/pDeneU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble enWh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-2;j3<=2;j3+=2){
	M = amp(j1,j3);
	amp2 += std::pow(std::abs(M),2);
      }
    //Average over initial spins 1/2
    return amp2/2.0;
  }

  
  

  //  Tests
  int test_enWh(){
    int n=0;//Number of fails
    std::cout<<"\t* E , ne -> W+, h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=2500\n";
      ldouble me=0.0005, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      enWh enWhAmp = enWh(EE,me,mh,MW,0,SW);
      ldouble pspatial=2500;
      ldouble dataCH[20] = {1.186767934131281E-03,3.332144888703923E-03,5.239146626103111E-03,6.907773146329908E-03,8.338024449385655E-03,9.529900535272058E-03,1.048340140399133E-02,1.119852705554640E-02,1.167527748994120E-02,1.191365270718116E-02,1.191365270727403E-02,1.167527749023117E-02,1.119852705607000E-02,1.048340140481868E-02,9.529900536525839E-03,8.338024451283342E-03,6.907773149288023E-03,5.239146631051662E-03,3.332144898450782E-03,1.186767968373887E-03};
      i += enWhAmp.test_2to2_amp2([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH);
      i += enWhAmp.test_2to2_amp2_rotations([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH1[20] = {3.611380881222951E-03,5.478027452858238E-03,7.137268849869218E-03,8.589105072256966E-03,9.833536120022836E-03,1.087056199316855E-02,1.170018269169634E-02,1.232239821560914E-02,1.273720856491092E-02,1.294461373960712E-02,1.294461373970547E-02,1.273720856521731E-02,1.232239821615995E-02,1.170018269256123E-02,1.087056199446898E-02,9.833536121972766E-03,8.589105075262312E-03,7.137268854821578E-03,5.478027462359839E-03,3.611380911317014E-03};
      i += enWhAmp.test_2to2_amp2([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH1);
      i += enWhAmp.test_2to2_amp2_rotations([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH1);
      i += enWhAmp.test_2to2_amp2_boosts([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH1);
      i += enWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH1);
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {1.271067266058783E-02,1.355412879612740E-02,1.430386758327713E-02,1.495988902203780E-02,1.552219311241035E-02,1.599077985439592E-02,1.636564924799593E-02,1.664680129321215E-02,1.683423599004683E-02,1.692795333850284E-02,1.692795333858392E-02,1.683423599029501E-02,1.664680129364270E-02,1.636564924863597E-02,1.599077985528703E-02,1.552219311361232E-02,1.495988902363193E-02,1.430386758535959E-02,1.355412879873379E-02,1.271067266308490E-02};
      i += enWhAmp.test_2to2_amp2([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH2);
      i += enWhAmp.test_2to2_amp2_rotations([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH2);
      i += enWhAmp.test_2to2_amp2_boosts([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH2);
      i += enWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, mh=125, MW=80.385, pspatial=95\n";
      me = 125.1;
      mh = 125;
      pspatial = 95;
      enWhAmp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {3.075279318865571E-02,3.441668909660767E-02,3.838636509296271E-02,4.269207235844051E-02,4.736727771150183E-02,5.244903335666629E-02,5.797837438335404E-02,6.400073391425362E-02,7.056635580999675E-02,7.773066804111671E-02,8.555455167256468E-02,9.410439304348364E-02,1.034517265320720E-01,1.136721380862298E-01,1.248428616712936E-01,1.370380807268284E-01,1.503201900589582E-01,1.647238772145887E-01,1.802272315497096E-01,1.966988883992521E-01};
      i += enWhAmp.test_2to2_amp2([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH4);
      i += enWhAmp.test_2to2_amp2_rotations([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH4);
      i += enWhAmp.test_2to2_amp2_boosts([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH4);
      i += enWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH4);
      //std::cout<<"\n# me=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      me = 125;
      mh = 0.0005;
      pspatial = 125.1;
      enWhAmp.set_masses(me,mh,MW);
      ldouble dataCH3[20] = {3.103507322622199E-02,3.572786745415579E-02,4.091762852003847E-02,4.669321148801709E-02,5.316002704969990E-02,6.044404667671945E-02,6.869701339665645E-02,7.810329546693059E-02,8.888900583416909E-02,1.013342849800070E-01,1.157900542043697E-01,1.327011574368732E-01,1.526387115172555E-01,1.763457639800118E-01,2.048019632333628E-01,2.393141357713542E-01,2.816367992158181E-01,3.341041444144451E-01,3.996451819758785E-01,4.810333584133815E-01};
      i += enWhAmp.test_2to2_amp2([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH3);
      i += enWhAmp.test_2to2_amp2_rotations([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH3);
      i += enWhAmp.test_2to2_amp2_boosts([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH3);
      i += enWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH3);
      //std::cout<<"\n# me=125, mh=0.0005, MW=0.0004, pspatial=125.1\n";
      me = 125;
      mh = 0.0005;
      MW = 0.0004;
      pspatial = 125.1;
      enWhAmp.set_masses(me,mh,MW);
      ldouble dataCH5[20] = {6.813312401116914E+19,7.552753558876905E+19,8.354326442788661E+19,9.226115386510693E+19,1.017764858262178E+20,1.122022808850074E+20,1.236735147790681E+20,1.363525450007007E+20,1.504361412276982E+20,1.661646421906618E+20,1.838339142342936E+20,2.038109299137046E+20,2.265537998379429E+20,2.526366113938465E+20,2.827773600906687E+20,3.178602906359190E+20,3.589216114634847E+20,4.069941910970704E+20,4.624487872872736E+20,5.224510587696056E+20};
      i += enWhAmp.test_2to2_amp2([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH5);
      i += enWhAmp.test_2to2_amp2_rotations([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH5);
      i += enWhAmp.test_2to2_amp2_boosts([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH5);
      i += enWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enWhAmp.amp2(); }, me,0,MW,mh,pspatial,dataCH5);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}
