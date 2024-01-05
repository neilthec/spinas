
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

//File:  SPINAS/SM/dWgu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dWgu.h"

namespace spinas {

  dWgu::dWgu(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& massd, const ldouble& massW, const ldouble& sinW):
    e(echarge), gs(gscharge), mu(massu), md(massd), MW(massW), SW(sinW) {
    constexpr ldouble sqrt2 = std::sqrt(2);
    propu = propagator(mu,0);
    propd = propagator(md,0);
    p1=particle(md);
    p2=particle(MW);
    p3=particle(0);
    p4=particle(mu);
    a21a = sproduct(ANGLE,&p2,&p1);
    s24s = sproduct(SQUARE,&p2,&p4);
    s32s = sproduct(SQUARE,&p3,&p2);
    s34s = sproduct(SQUARE,&p3,&p4);
    s3413s = sproduct(SQUARE,&p3,&p4,&p1,&p3);
    a32a = sproduct(ANGLE,&p3,&p2);
    a31a = sproduct(ANGLE,&p3,&p1);
    a3413a = sproduct(ANGLE,&p3,&p4,&p1,&p3);
    //prefactor
    pre = sqrt2*e*gs/(MW*SW);
  }
  void dWgu::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massW){
    constexpr ldouble sqrt2 = std::sqrt(2);
    mu=massu;
    md=massd;
    MW=massW;
    p1.set_mass(md);
    p2.set_mass(MW);
    p4.set_mass(mu);
    propu.set_mass(mu);
    propd.set_mass(md);
    pre = sqrt2*e*gs/(MW*SW);
  }
  void dWgu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    constexpr ldouble one=1;
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a21a.update();
    s24s.update();
    s32s.update();
    s34s.update();
    s3413s.update();
    a32a.update();
    a31a.update();
    a3413a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propu.denominator(propSP);
    pDenT=propd.denominator(propTP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble dWgu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    constexpr ldouble one=1, two=2, three=3;
    cdouble amplitude(0,0);
    int ds2a, ds2b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds2);
    ldouble normFactor=get_spin_normalization(ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, i);

      if(ds3>0){
      
	//ST Diagram
	//pre = sqrt(2)*e*gs/(MW*SW);
	//gW+Ud all ingoing:
	//T u: +pre<24>[23][1341]/(t-Mu^2)/(u-Md^2)  - pre[12][13]<24>/(t-Mu^2)
	//U d: +pre<24>[23][1341]/(u-Md^2)/(t-Mu^2)
	//dW+gU all in: 4->1->3->4
	//T u: +pre<21>[24][3413]/(s-Mu^2)/(t-Md^2)  - pre[32][34]<21>/(s-Mu^2)
	//U d: +pre<21>[24][3413]/(t-Md^2)/(s-Mu^2)
	//34 out:
	//T u: -pre<21>[24][3413]/(s-Mu^2)/(t-Md^2)  - pre[32][34]<21>/(s-Mu^2)
	//U d: -pre<21>[24][3413]/(t-Md^2)/(s-Mu^2)
	amplitude += -normFactor*pre*a21a.v(ds2a,ds1)*s24s.v(ds2b,ds4)*s3413s.v()/pDenS/pDenT
	  -normFactor*pre*s32s.v(ds2a)*s34s.v(ds4)*a21a.v(ds2b,ds1)/pDenS;
	
	
      }
      else if(ds3<0){

	//TU Diagram
	amplitude += -normFactor*pre*a21a.v(ds2a,ds1)*s24s.v(ds2b,ds4)*a3413a.v()/pDenS/pDenT
	  +normFactor*pre*a32a.v(ds2a)*a31a.v(ds1)*s24s.v(ds2b,ds4)/pDenT;
	
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble dWgu::amp2(){
    constexpr ldouble three=3, four=4;
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    //Color factor Tr(Ta,Ta) = 4
	    amp2 += four*std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/3=1/6
    //Average over initial colors 1/3
    return amp2/18.0;
  }


  



  //  Tests
  int test_dWgu(){
    int n=0;//Number of fails
    std::cout<<"\t* d , W+ -> g , u       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,gs=1.238,mu=0.0042,md=0.0075,MW=80.385, SW=0.474;
      dWgu dWguAmp = dWgu(EE,gs,mu,md,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.164249809165811E+01,3.909947591246862E+00,2.375050699707527E+00,1.725532901635787E+00,1.371141468613759E+00,1.150898383466671E+00,1.002889024139098E+00,8.982199003197011E-01,8.215944583074006E-01,7.641567985881479E-01,7.204246976518766E-01,6.868227923391377E-01,6.609198357011132E-01,6.410049502550751E-01,6.258392247202893E-01,6.145034668561200E-01,6.063013189044545E-01,6.006955905947083E-01,5.972652409039657E-01,5.956755961986067E-01};
      i += dWguAmp.test_2to2_amp2([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH);
      i += dWguAmp.test_2to2_amp2_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH);
      i += dWguAmp.test_2to2_amp2_boosts([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH);
      i += dWguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=81\n";
      pspatial = 81;
      ldouble dataCH2[20] = {1.040552451907005E+01,3.552343832851400E+00,2.191592836378876E+00,1.615474687758129E+00,1.300900824781999E+00,1.105210728683495E+00,9.735349688158675E-01,8.802678107756812E-01,8.118532774889030E-01,7.604431716864986E-01,7.211790855278402E-01,6.908924921565178E-01,6.674287933080597E-01,6.492717902025321E-01,6.353235457916456E-01,6.247694616197416E-01,6.169924176443427E-01,6.115163498324709E-01,6.079681269148903E-01,6.060511573697984E-01};
      i += dWguAmp.test_2to2_amp2([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH2);
      i += dWguAmp.test_2to2_amp2_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH2);
      i += dWguAmp.test_2to2_amp2_boosts([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH2);
      i += dWguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH2);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=250\n";
      mu=80;
      dWguAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH3[20] = {1.795841265414397E+01,5.818687466897673E+00,3.408106752856043E+00,2.387403831417540E+00,1.829993508439551E+00,1.483170730423001E+00,1.249741245914482E+00,1.084347737537981E+00,9.629775078478467E-01,8.717284688444406E-01,8.019945660495255E-01,7.481631552369947E-01,7.064176383324505E-01,6.740722611350112E-01,6.491821671541569E-01,6.303043983483849E-01,6.163458046861608E-01,6.064630990002365E-01,5.999953264838043E-01,5.964171128032177E-01};
      i += dWguAmp.test_2to2_amp2([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH3);
      i += dWguAmp.test_2to2_amp2_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH3);
      i += dWguAmp.test_2to2_amp2_boosts([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH3);
      i += dWguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH3);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {4.424282508268378E+02,1.521132402414763E+02,9.404186619828083E+01,6.915354842960294E+01,5.532671533459799E+01,4.652793376943807E+01,4.043659839441104E+01,3.596974441440589E+01,3.255403184527842E+01,2.985752460440131E+01,2.767473730915981E+01,2.587165716468395E+01,2.435715514544003E+01,2.306710322916243E+01,2.195506376690550E+01,2.098658294588407E+01,2.013555934482097E+01,1.938185754929196E+01,1.870969572569816E+01,1.810652932391350E+01};
      i += dWguAmp.test_2to2_amp2([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH4);
      i += dWguAmp.test_2to2_amp2_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH4);
      i += dWguAmp.test_2to2_amp2_boosts([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH4);
      i += dWguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH4);
      //std::cout<<"\n# mu=80, md=0.0075, MW=1, pspatial=250\n";
      MW=1;
      dWguAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH5[20] = {3.823358879675073E+04,1.150382701277157E+04,6.195985583228986E+03,3.948418514850055E+03,2.720942902336895E+03,1.957145224786719E+03,1.443020281603106E+03,1.078698979102433E+03,8.113094441318962E+02,6.102416448551188E+02,4.565465139395658E+02,3.378659639917373E+02,2.457964955825594E+02,1.744245350388548E+02,1.194678393553596E+02,7.774926049559087E+01,4.686186812864959E+01,2.494888742769942E+01,1.055500839321855E+01,2.523446153722787E+00};
      i += dWguAmp.test_2to2_amp2([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH5);
      i += dWguAmp.test_2to2_amp2_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH5);
      i += dWguAmp.test_2to2_amp2_boosts([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH5);
      i += dWguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH5);
      //std::cout<<"\n# mu=80, md=0.0075, MW=1, pspatial=50\n";
      pspatial = 50;
      ldouble dataCH6[20] = {4.028381679452666E+05,1.268062226879477E+05,7.163792198318374E+04,4.802158211019632E+04,3.492256351779988E+04,2.660414573137429E+04,2.085989768256391E+04,1.666015150579117E+04,1.345978890704364E+04,1.094321507933873E+04,8.915061806372356E+03,7.247915032764791E+03,5.855133200157719E+03,4.675746322631529E+03,3.665603720478234E+03,2.791948437627673E+03,2.029964598165805E+03,1.360508526317948E+03,7.685757194018622E+02,2.422395214054864E+02};
      i += dWguAmp.test_2to2_amp2([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH6);
      i += dWguAmp.test_2to2_amp2_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH6);
      i += dWguAmp.test_2to2_amp2_boosts([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH6);
      i += dWguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH6);
      //std::cout<<"\n# mu=0.0042, md=0.0042, MW=80.385, pspatial=250\n";
      mu = 0.0042;
      md = 0.0042;
      MW = 80.385;
      dWguAmp.set_masses(mu,md,MW);
      pspatial = 250;
      ldouble dataCH7[20] = {1.164249812153867E+01,3.909947587191118E+00,2.375050696103047E+00,1.725532898943162E+00,1.371141466601783E+00,1.150898381944400E+00,1.002889022975779E+00,8.982198994255801E-01,8.215944576191240E-01,7.641567980596362E-01,7.204246972486816E-01,6.868227920348855E-01,6.609198354751922E-01,6.410049500910808E-01,6.258392246049592E-01,6.145034667785859E-01,6.063013188557005E-01,6.006955905671744E-01,5.972652408912499E-01,5.956755961952424E-01};
      i += dWguAmp.test_2to2_amp2([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH7);
      i += dWguAmp.test_2to2_amp2_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH7);
      i += dWguAmp.test_2to2_amp2_boosts([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH7);
      i += dWguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dWguAmp.amp2(); }, md,MW,0,mu,pspatial,dataCH7);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  



}
