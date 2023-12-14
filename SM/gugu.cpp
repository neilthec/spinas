
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

//File:  SPINAS/SM/gugu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/gugu.h"

namespace spinas {

  gugu::gugu(const ldouble& gcharge, const ldouble& massu):
    gs(gcharge), mu(massu), propu(massu,0), propg(0,0){
    constexpr ldouble two=2;
    sqrt2 = std::sqrt(two);
    p1=particle(0);
    p2=particle(mu);
    p3=particle(0);
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
    s143a = sproduct(SQUARE,&p1,&p4,&p3);
    s341a = sproduct(SQUARE,&p3,&p4,&p1);
  }
  void gugu::set_masses(const ldouble& massu){
    mu=massu;
    p2.set_mass(mu);
    p4.set_mass(mu);
    propu.set_mass(mu);
  }
  void gugu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    s143a.update();
    s341a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT = std::real(propg.den(propTP));
    pDenU = std::real(propu.den(propUP));
    pDenS = std::real(propu.den(propSP));
    
  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble gugu::amp_Num(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds1>0&&ds3>0){
      //ggUu all in:      + mu[12]^2<34>
      //gugU: 4->2->3->4: - mu[13]^2<24>
      //34 out:           + mu[13]^2<24>
      return + 2.0*gs*gs*mu*s13s.v()*s13s.v()*a24a.v(ds2,ds4);
    }
    else if(ds1<0&&ds3<0){
      //ggUu all in:      + mu<12>^2[34]
      //gugU: 4->2->3->4: - mu<13>^2[24]
      //34 out:           - mu<13>^2[24]
      return - 2.0*gs*gs*mu*a13a.v()*a13a.v()*s24s.v(ds2,ds4);
    }
    else if(ds1>0&&ds3<0){
      //ggUu all in:      + ([31]<42>+[41]<32>)[132>
      //gugU: 4->2->3->4: - ([14]<23>-[12]<34>)[143>
      //34 out:           + ([14]<23>+[12]<34>)[143>
      return + 2.0*gs*gs*(s14s.v(ds4)*a23a.v(ds2)+s12s.v(ds2)*a34a.v(ds4))*s143a.v();
    }
    else if(ds1<0&&ds3>0){
      //ggUu all in:      + (<31>[42]+<41>[32])[231>
      //gugU: 4->2->3->4: - (<14>[23]-<12>[34])[341>
      //34 out:           - (<14>[23]+<12>[34])[341>
      return - 2.0*gs*gs*(a14a.v(ds4)*s23s.v(ds2)+a12a.v(ds2)*s34s.v(ds4))*s341a.v();
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble gugu::amp2(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-1;j4<=1;j4+=2){
	    M_Num = amp_Num(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenT*pDenU,2)+6.0/std::pow(pDenT*pDenS,2)-2.0/3.0/std::pow(pDenU*pDenS,2));
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/8*1/3=1/24
    return amp2/96.0;
  }

  //g+, u -> g, u
  ldouble gugu::amp2_gplus(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-1;j4<=1;j4+=2){
	  M_Num = amp_Num(2,j2,j3,j4);
	  amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenT*pDenU,2)+6.0/std::pow(pDenT*pDenS,2)-2.0/3.0/std::pow(pDenU*pDenS,2));
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/8*1/3=1/24
    return amp2/48.0;
  }

  //g-, u -> g, u
  ldouble gugu::amp2_gminus(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-1;j4<=1;j4+=2){
	  M_Num = amp_Num(-2,j2,j3,j4);
	  amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenT*pDenU,2)+6.0/std::pow(pDenT*pDenS,2)-2.0/3.0/std::pow(pDenU*pDenS,2));
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/8*1/3=1/24
    return amp2/48.0;
  }
  



  //  Tests
  int test_gugu(){
    int n=0;//Number of fails
    std::cout<<"\t* g , u  -> g , u       :";
    {//amp^2
      int i=0;
      // mu=0.0042, pspatial=250
      ldouble mu=0.0042;
      ldouble gg=1.238;
      gugu guguAmp = gugu(gg,mu);
      ldouble pspatial=250;
      ldouble dataCH[20] = {7.333310273672108E+03,7.770025650469427E+02,2.675433730399247E+02,1.310339897437831E+02,7.642511969828286E+01,4.958452767864433E+01,3.462305364143432E+01,2.555187422951544E+01,1.972050321857333E+01,1.581728982433201E+01,1.313910801503922E+01,1.128817621558588E+01,1.003456987818882E+01,9.251709470837460E+00,8.890385702217413E+00,8.983807689799983E+00,9.705338377394890E+00,1.159850870319251E+01,1.675908731700997E+01,4.425860552874457E+01};
      i += guguAmp.test_2to2_amp2([&]() { return guguAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guguAmp.test_2to2_amp2_rotations([&]() { return guguAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guguAmp.test_2to2_amp2_boosts([&]() { return guguAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guguAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH);
      //Helicities g+, u -> g, u
      i += guguAmp.test_2to2_amp2([&]() { return guguAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guguAmp.test_2to2_amp2_rotations([&]() { return guguAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guguAmp.test_2to2_amp2_boosts([&]() { return guguAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guguAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH);
      //Helicities g-, u -> g, u
      i += guguAmp.test_2to2_amp2([&]() { return guguAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guguAmp.test_2to2_amp2_rotations([&]() { return guguAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guguAmp.test_2to2_amp2_boosts([&]() { return guguAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guguAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH);
      //Close to threshold
      pspatial = 0.005;
      ldouble dataCH2[20] = {9.714451875342194E+03,1.020618194519480E+03,3.477384000410041E+02,1.681235134498774E+02,9.653711415283021E+01,6.147377759338010E+01,4.198507534930147E+01,3.018911975258480E+01,2.260217609265208E+01,1.750043946676486E+01,1.395705731722020E+01,1.144177027285896E+01,9.637540739244280E+00,8.350497813958661E+00,7.464200369661594E+00,6.917888854649910E+00,6.701772867899306E+00,6.872913255024232E+00,7.615691138492007E+00,9.445088107397801E+00};
      i += guguAmp.test_2to2_amp2([&]() { return guguAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guguAmp.test_2to2_amp2_rotations([&]() { return guguAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guguAmp.test_2to2_amp2_boosts([&]() { return guguAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guguAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH2);
      //Helicities g+, u -> g, u
      i += guguAmp.test_2to2_amp2([&]() { return guguAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guguAmp.test_2to2_amp2_rotations([&]() { return guguAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guguAmp.test_2to2_amp2_boosts([&]() { return guguAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guguAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH2);
      //Helicities g-, u -> g, u
      i += guguAmp.test_2to2_amp2([&]() { return guguAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guguAmp.test_2to2_amp2_rotations([&]() { return guguAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guguAmp.test_2to2_amp2_boosts([&]() { return guguAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guguAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guguAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH2);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


}
