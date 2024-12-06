
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

//File:  SPINAS/SM/uAuh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uAuh.h"

namespace spinas {

  uAuh::uAuh(const ldouble& echarge, const ldouble& massu, const ldouble& massh, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qu(2.0/3.0), mu(massu), mh(massh), prop(massu,0), MW(massW), SW(sinW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(mu);
    p2=particle(0);
    p3=particle(mu);
    p4=particle(mh);
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s243a = sproduct(SQUARE,&p2,&p4,&p3,2);
    a243s = sproduct(ANGLE,&p2,&p4,&p3,2);
    s241a = sproduct(SQUARE,&p2,&p4,&p1,2);
    a241s = sproduct(ANGLE,&p2,&p4,&p1,2);
    s2342s = sproduct(SQUARE,&p2,&p3,&p4,&p2,2);
    a2342a = sproduct(ANGLE,&p2,&p3,&p4,&p2,2);
    pre = sqrt2*e*e*Qu*mu/(2.0*MW*SW);
  }
  void uAuh::set_masses(const ldouble& massu, const ldouble& massh, const ldouble& massW){
    mu=massu;
    mh=massh;
    p1.set_mass(mu);
    p3.set_mass(mu);
    p4.set_mass(mh);
    prop.set_mass(mu);
    pre = sqrt2*e*e*Qu*mu/(2.0*MW*SW);
  }
  void uAuh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s243a.update();
    a243s.update();
    s241a.update();
    a241s.update();
    s2342s.update();
    a2342a.update();
    //Propagator Momentum
    ldouble propSP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=prop.denominator(propSP);
    pDenU=prop.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uAuh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    /*mh^2[12][23] - mu[12][243> + mu[23][241> - <13>[2342]
      mh^2<12><23> - mu<12><243] + mu<23><241] - [13]<2342>
      Becomes after a sign change due to p3 and p4 being outgoing:*/
    //pre = sqrt2*e*e*Qu*mu/(2.0*MW*SW);
    if(ds2>0){
      //mh^2[12][23] - mu[12][243> - mu[23][241> + <13>[2342]
      return pre*(mh*mh*s12s.v(ds1)*s23s.v(ds3) - mu*s12s.v(ds1)*s243a.v(ds3) - mu*s23s.v(ds3)*s241a.v(ds1) + a13a.v(ds1,ds3)*s2342s.v())/pDenS/pDenU;
    }
    else if(ds2<0){
      //-mh^2<12><23> + mu<12><243] + mu<23><241] - [13]<2342>
      return pre*(- mh*mh*a12a.v(ds1)*a23a.v(ds3) + mu*a12a.v(ds1)*a243s.v(ds3) + mu*a23a.v(ds3)*a241s.v(ds1) - s13s.v(ds1,ds3)*a2342a.v())/pDenS/pDenU;
    }
    return cdouble(0,0);    
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uAuh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-1;j3<=1;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3
    return amp2/12.0;
  }
  
  //set_momenta(...) must be called before amp2_Aplus().
  //Positive Helicity Photon Only
  ldouble uAuh::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-1;j3<=1;j3+=2){
	M = amp(j1,2,j3);
	amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
      }
    //Average over initial spins 1/2^1
    //Average over initial colors 1/3
    return amp2/6.0;
  }
  



  //  Tests
  int test_uAuh(){
    int n=0;//Number of fails
    std::cout<<"\t* u , A  -> u , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, mh=125, pspatial=250\n";
      ldouble mu=0.0042, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      uAuh uAuhAmp = uAuh(EE,mu,mh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {6.342847948664793E-14,1.327982421471252E-13,2.798052418175084E-13,5.185653397711458E-13,8.668372001902914E-13,1.347278429100123E-12,1.989259997273308E-12,2.831551586128737E-12,3.926407983664147E-12,5.345990563778535E-12,7.192844436744968E-12,9.617771235310047E-12,1.285200012967031E-11,1.726905516395889E-11,2.351411423251879E-11,3.280584669115379E-11,4.775568687317498E-11,7.515336404286398E-11,1.398945521524449E-10,4.660403162172609E-10};
      i += uAuhAmp.test_2to2_amp2([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uAuhAmp.test_2to2_amp2_rotations([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uAuhAmp.test_2to2_amp2_boosts([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uAuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      //(CH produces the same M^2s since it is averaging over spins and the to helicity cases give the same result -- averaging them gaves the same as each individually.)    
      i += uAuhAmp.test_2to2_amp2([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uAuhAmp.test_2to2_amp2_rotations([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uAuhAmp.test_2to2_amp2_boosts([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uAuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, mh=125, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {1.114287808927230E-12,1.227300641848175E-12,1.409026345345722E-12,1.671958168612582E-12,2.031813425148426E-12,2.508645239072314E-12,3.128448400907201E-12,3.925536034140715E-12,4.946156261730248E-12,6.254175539181699E-12,7.940351461149953E-12,1.013815107857792E-11,1.305222386953996E-11,1.701315745238615E-11,2.259196680409112E-11,3.086723576959871E-11,4.415021544380359E-11,6.845015866749497E-11,1.258009924806883E-10,4.145080588285423E-10};
      i += uAuhAmp.test_2to2_amp2([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uAuhAmp.test_2to2_amp2_rotations([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uAuhAmp.test_2to2_amp2_boosts([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uAuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH2);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += uAuhAmp.test_2to2_amp2([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uAuhAmp.test_2to2_amp2_rotations([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uAuhAmp.test_2to2_amp2_boosts([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uAuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, mh=125, pspatial=95\n";
      mu = 125.1;
      mh = 125;
      pspatial = 95;
      uAuhAmp.set_masses(mu,mh,MW);
      ldouble dataCH4[20] = {4.036408890809026E-03,4.267462260633547E-03,4.499202280169918E-03,4.731433212637942E-03,4.963939752394659E-03,5.196485242436028E-03,5.428809718307906E-03,5.660627759781720E-03,5.891626129429998E-03,6.121461174724579E-03,6.349755967433104E-03,6.576097150857806E-03,6.800031461788515E-03,7.021061889861501E-03,7.238643432251182E-03,7.452178396181511E-03,7.661011195522977E-03,7.864422580615755E-03,8.061623232283809E-03,8.251746641608666E-03};
      i += uAuhAmp.test_2to2_amp2([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uAuhAmp.test_2to2_amp2_rotations([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uAuhAmp.test_2to2_amp2_boosts([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uAuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH4);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += uAuhAmp.test_2to2_amp2([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uAuhAmp.test_2to2_amp2_rotations([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uAuhAmp.test_2to2_amp2_boosts([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uAuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, mh=0.0005, pspatial=125.1\n";
      mu = 125;
      mh = 0.0005;
      pspatial = 125.1;
      uAuhAmp.set_masses(mu,mh,MW);
      ldouble dataCH3[20] = {4.067056268666889E-04,1.293870417939113E-03,2.288835362223577E-03,3.404475217492258E-03,4.655688689126602E-03,6.059772200141682E-03,7.636859979665393E-03,9.410431236463145E-03,1.140787080406892E-02,1.366103671766918E-02,1.620671547141002E-02,1.908668630629939E-02,2.234676272800088E-02,2.603337791148048E-02,3.018439408509305E-02,3.480615343833617E-02,3.981651479037487E-02,4.489859977833998E-02,4.909909929703737E-02,4.960239948660988E-02};
      i += uAuhAmp.test_2to2_amp2([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uAuhAmp.test_2to2_amp2_rotations([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uAuhAmp.test_2to2_amp2_boosts([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uAuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAuhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH3);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += uAuhAmp.test_2to2_amp2([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uAuhAmp.test_2to2_amp2_rotations([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uAuhAmp.test_2to2_amp2_boosts([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uAuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAuhAmp.amp2_Aplus(); }, mu,0,mu,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
