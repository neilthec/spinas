
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

//File:  SPINAS/SM/gguu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/gguu.h"

namespace spinas {

  gguu::gguu(const ldouble& gcharge, const ldouble& massu):
    gs(gcharge), mu(massu), propu(massu,0), propg(0,0){
    constexpr ldouble two=2;
    sqrt2 = std::sqrt(two);
    p1=particle(0);
    p2=particle(0);
    p3=particle(mu);
    p4=particle(mu);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s132a = sproduct(SQUARE,&p1,&p3,&p2);
    s231a = sproduct(SQUARE,&p2,&p3,&p1);
  }
  void gguu::set_masses(const ldouble& massu){
    mu=massu;
    p3.set_mass(mu);
    p4.set_mass(mu);
    propu.set_mass(mu);
  }
  void gguu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s34s.update();
    a34a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    s24s.update();
    a24a.update();
    s23s.update();
    a23a.update();
    s14s.update();
    a14a.update();
    s132a.update();
    s231a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS = std::real(propg.denominator(propSP));
    pDenT = std::real(propu.denominator(propTP));
    pDenU = std::real(propu.denominator(propUP));
    
  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble gguu::amp_Num(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds1>0&&ds2>0){
      //eeAA:   me[34]^2<12>
      //AAee:   me[12]^2<34>
      //34 out: me[12]^2<34>
      return 2.0*gs*gs*mu*s12s.v()*s12s.v()*a34a.v(ds3,ds4);
    }
    else if(ds1<0&&ds2<0){
      //eeAA:   me<34>^2[12]
      //AAee:   me<12>^2[34]
      //34 out: me<12>^2[34]
      return 2.0*gs*gs*mu*a12a.v()*a12a.v()*s34s.v(ds3,ds4);
    }
    else if(ds1>0&&ds2<0){
      //eeAA:   ([13]<24>+[23]<14>)[314>
      //AAee:   ([31]<42>+[41]<32>)[132>
      //34 out: -([13]<24>+[14]<23>)[132>
      return -2.0*gs*gs*(s13s.v(ds3)*a24a.v(ds4)+s14s.v(ds4)*a23a.v(ds3))*s132a.v();
    }
    else if(ds1<0&&ds2>0){
      //eeAA:   (<13>[24]+<23>[14])[413>
      //AAee:   (<31>[42]+<41>[32])[231>
      //34 out: -(<13>[24]+<14>[23])[231>
      return -2.0*gs*gs*(a13a.v(ds3)*s24s.v(ds4)+a14a.v(ds4)*s23s.v(ds3))*s231a.v();
    }
    return cdouble(0,0);    
  }

  //set_momenta(...) must be called before amp2_Feynman().
  //This expression comes from Ellis, Stirling and Weber: QCD and Collider Physics, Table 7.1
  //This agrees with CalcHEP.
  //This assumes massless quarks (up quark here)
  //Ours are not massless, so we do not expect agreement.
  //But, with very small quark masses, we might get close.
  ldouble gguu::amp2_Feynman(){
    ldouble pre = gs*gs*gs*gs;
    ldouble oneSixth=1.0/6.0, threeEights = 3.0/8.0;
    cdouble M2 =  pre*(pDenT*pDenT+pDenU*pDenU)*(oneSixth/pDenT/pDenU - threeEights/pDenS/pDenS);
    return std::abs(M2);
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble gguu::amp2(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M_Num = amp_Num(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenS*pDenT,2)+6.0/std::pow(pDenS*pDenU,2)-2.0/3.0/std::pow(pDenT*pDenU,2));
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/8*1/8=1/64
    return amp2/256.0;
  }

  //g+, g+ -> e, E
  ldouble gguu::amp2_gplus_gplus(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M_Num = amp_Num(2,2,j3,j4);
	amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenS*pDenT,2)+6.0/std::pow(pDenS*pDenU,2)-2.0/3.0/std::pow(pDenT*pDenU,2));
      }
    //Average over initial colors 1/8*1/8=1/64
    return amp2/64.0;
  }

  //g+, g- -> e, E
  ldouble gguu::amp2_gplus_gminus(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M_Num = amp_Num(2,-2,j3,j4);
	amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenS*pDenT,2)+6.0/std::pow(pDenS*pDenU,2)-2.0/3.0/std::pow(pDenT*pDenU,2));
      }
    //Average over initial colors 1/8*1/8=1/64
    return amp2/64.0;
  }
  




  //  Tests
  int test_gguu(){
    int n=0;//Number of fails
    std::cout<<"\t* g , g  -> u , U       :";
    {//amp^2
      int i=0;
      // mu=0.0042, pspatial=250
      ldouble mu=0.0042;
      ldouble gg=1.238;
      gguu gguuAmp = gguu(gg,mu);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.444059151504384E+01,4.101585509807672E+00,2.108242850779164E+00,1.302164654548335E+00,8.884905505642698E-01,6.510096134501014E-01,5.072237794654312E-01,4.194347333956981E-01,3.686983203722325E-01,3.453858701080085E-01,3.453858701080085E-01,3.686983203722323E-01,4.194347333956979E-01,5.072237794654310E-01,6.510096134501013E-01,8.884905505642697E-01,1.302164654548334E+00,2.108242850779164E+00,4.101585509807675E+00,1.444059151504386E+01};
      i += gguuAmp.test_2to2_amp2([&]() { return gguuAmp.amp2_Feynman(); }, 0,0,mu,mu,pspatial,dataCH);
      i += gguuAmp.test_2to2_amp2([&]() { return gguuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      i += gguuAmp.test_2to2_amp2_rotations([&]() { return gguuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      i += gguuAmp.test_2to2_amp2_boosts([&]() { return gguuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      i += gguuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gguuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      //Helicities g+, g+ -> u, U
      ldouble dataCHpp[20] = {8.788894583813790E-08,9.687419494936052E-09,3.481778227155410E-09,1.789535021890268E-09,1.104102586701322E-09,7.663906419208555E-10,5.813594342800148E-10,4.753820332715521E-10,4.164566128474038E-10,3.899293132846905E-10,3.899292838757066E-10,4.164564885436248E-10,4.753818531584662E-10,5.813592700504907E-10,7.663906206433879E-10,1.104102580711105E-09,1.789534967910553E-09,3.481778209794622E-09,9.687419496304858E-09,8.788894604819868E-08};
      i += gguuAmp.test_2to2_amp2([&]() { return gguuAmp.amp2_gplus_gplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      i += gguuAmp.test_2to2_amp2_rotations([&]() { return gguuAmp.amp2_gplus_gplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      i += gguuAmp.test_2to2_amp2_boosts([&]() { return gguuAmp.amp2_gplus_gplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      i += gguuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gguuAmp.amp2_gplus_gplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      //Helicities g+, g- -> u, U
      ldouble dataCHpm[20] = {2.888118294219873E+01,8.203171009927924E+00,4.216485698076549E+00,2.604329307307134E+00,1.776981100024437E+00,1.302019226133812E+00,1.014447558349503E+00,8.388694663160141E-01,7.373966403280083E-01,6.907717398260877E-01,6.907717398260876E-01,7.373966403280080E-01,8.388694663160139E-01,1.014447558349503E+00,1.302019226133812E+00,1.776981100024437E+00,2.604329307307134E+00,4.216485698076549E+00,8.203171009927930E+00,2.888118294219878E+01};
      i += gguuAmp.test_2to2_amp2([&]() { return gguuAmp.amp2_gplus_gminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      i += gguuAmp.test_2to2_amp2_rotations([&]() { return gguuAmp.amp2_gplus_gminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      i += gguuAmp.test_2to2_amp2_boosts([&]() { return gguuAmp.amp2_gplus_gminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      i += gguuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gguuAmp.amp2_gplus_gminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      //Close to threshold
      pspatial = 0.005;
      ldouble dataCH2[20] = {8.266968024506537E-01,7.531073093494767E-01,6.888970829412314E-01,6.345762242887008E-01,5.897577931447362E-01,5.537668931419482E-01,5.259047902792747E-01,5.055624056815474E-01,4.922669524550842E-01,4.856992963966429E-01,4.856992963966428E-01,4.922669524550842E-01,5.055624056815473E-01,5.259047902792747E-01,5.537668931419482E-01,5.897577931447363E-01,6.345762242887008E-01,6.888970829412312E-01,7.531073093494767E-01,8.266968024506536E-01};
      //i += gguuAmp.test_2to2_amp2([&]() { return gguuAmp.amp2_Feynman(); }, 0,0,mu,mu,pspatial,dataCH); //This test fails because the Ellis, Stirling, Weber result assumes massless quarks.
      i += gguuAmp.test_2to2_amp2([&]() { return gguuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      i += gguuAmp.test_2to2_amp2_rotations([&]() { return gguuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      i += gguuAmp.test_2to2_amp2_boosts([&]() { return gguuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      i += gguuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gguuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      //Helicities g+, g+ -> u, U
      ldouble dataCH2pp[20] = {1.556935692060359E+00,1.285616160791518E+00,1.090129202686146E+00,9.466685919208927E-01,8.404142149918430E-01,7.618337173834990E-01,7.046679234105000E-01,6.647860157631863E-01,6.395171061469740E-01,6.272587707445658E-01,6.272587707445656E-01,6.395171061469740E-01,6.647860157631861E-01,7.046679234104999E-01,7.618337173834990E-01,8.404142149918431E-01,9.466685919208927E-01,1.090129202686146E+00,1.285616160791518E+00,1.556935692060358E+00};
      i += gguuAmp.test_2to2_amp2([&]() { return gguuAmp.amp2_gplus_gplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      i += gguuAmp.test_2to2_amp2_rotations([&]() { return gguuAmp.amp2_gplus_gplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      i += gguuAmp.test_2to2_amp2_boosts([&]() { return gguuAmp.amp2_gplus_gplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      i += gguuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gguuAmp.amp2_gplus_gplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      //Helicities g+, g- -> u, U
      ldouble dataCH2pm[20] = {9.645791284094850E-02,2.205984579074350E-01,2.876649631963165E-01,3.224838566565088E-01,3.391013712976295E-01,3.457000689003972E-01,3.471416571480495E-01,3.463387955999085E-01,3.450167987631943E-01,3.441398220487201E-01,3.441398220487201E-01,3.450167987631943E-01,3.463387955999085E-01,3.471416571480495E-01,3.457000689003972E-01,3.391013712976295E-01,3.224838566565089E-01,2.876649631963164E-01,2.205984579074352E-01,9.645791284094876E-02};
      i += gguuAmp.test_2to2_amp2([&]() { return gguuAmp.amp2_gplus_gminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      i += gguuAmp.test_2to2_amp2_rotations([&]() { return gguuAmp.amp2_gplus_gminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      i += gguuAmp.test_2to2_amp2_boosts([&]() { return gguuAmp.amp2_gplus_gminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      i += gguuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gguuAmp.amp2_gplus_gminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }




}
