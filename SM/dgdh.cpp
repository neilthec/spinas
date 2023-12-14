
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

//File:  SPINAS/SM/dgdh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dgdh.h"

namespace spinas {

  dgdh::dgdh(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd, const ldouble& massh, const ldouble& massW, const ldouble& sinW):
    e(echarge), gs(gscharge), md(massd), mh(massh), prop(massd,0), MW(massW), SW(sinW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(md);
    p2=particle(0);
    p3=particle(md);
    p4=particle(mh);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s243a = sproduct(SQUARE,&p2,&p4,&p3);
    a243s = sproduct(ANGLE,&p2,&p4,&p3);
    s241a = sproduct(SQUARE,&p2,&p4,&p1);
    a241s = sproduct(ANGLE,&p2,&p4,&p1);
    s2342s = sproduct(SQUARE,&p2,&p3,&p4,&p2);
    a2342a = sproduct(ANGLE,&p2,&p3,&p4,&p2);
    pre = sqrt2*e*gs*md/(2.0*MW*SW);
  }
  void dgdh::set_masses(const ldouble& massd, const ldouble& massh, const ldouble& massW){
    md=massd;
    mh=massh;
    p1.set_mass(md);
    p3.set_mass(md);
    p4.set_mass(mh);
    prop.set_mass(md);
    pre = sqrt2*e*gs*md/(2.0*MW*SW);
  }
  void dgdh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenS=prop.den(propSP);
    pDenU=prop.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble dgdh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    /*mh^2[12][23] - md[12][243> + md[23][241> - <13>[2342]
      mh^2<12><23> - md<12><243] + md<23><241] - [13]<2342>
      Becomes after a sign change due to p3 and p4 being outgoing:*/
    //pre = sqrt2*e*gs*md/(2.0*MW*SW);
    if(ds2>0){
      //mh^2[12][23] - md[12][243> - md[23][241> + <13>[2342]
      return pre*(mh*mh*s12s.v(ds1)*s23s.v(ds3) - md*s12s.v(ds1)*s243a.v(ds3) - md*s23s.v(ds3)*s241a.v(ds1) + a13a.v(ds1,ds3)*s2342s.v())/pDenS/pDenU;
    }
    else if(ds2<0){
      //-mh^2<12><23> + md<12><243] + md<23><241] - [13]<2342>
      return pre*(- mh*mh*a12a.v(ds1)*a23a.v(ds3) + md*a12a.v(ds1)*a243s.v(ds3) + md*a23a.v(ds3)*a241s.v(ds1) - s13s.v(ds1,ds3)*a2342a.v())/pDenS/pDenU;
    }
    return cdouble(0,0);    
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble dgdh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-1;j3<=1;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(TaTa) = 4 //Both channels have the same color factor
	}
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3*1/8=1/24
    return amp2/96.0;
  }
  
  //set_momenta(...) must be called before amp2_Aplus().
  //Positive Helicity Photon Only
  ldouble dgdh::amp2_gplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-1;j3<=1;j3+=2){
	M = amp(j1,2,j3);
	amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(TaTa) = 4
      }
    //Average over initial spins 1/2^1
    //Average over initial colors 1/3*1/8=1/24
    return amp2/48.0;
  }
  



  //  Tests
  int test_dgdh(){
    int n=0;//Number of fails
    std::cout<<"\t* d , g  -> d , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# md=0.0075, mh=125, pspatial=250\n";
      ldouble md=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble gs=1.238;
      dgdh dgdhAmp = dgdh(EE,gs,md,mh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.184068558085889E-12,2.479047659482750E-12,5.223341194094240E-12,9.680460876979496E-12,1.618192144641906E-11,2.515068998226030E-11,3.713505714352313E-11,5.285876660005170E-11,7.329729895023429E-11,9.979774644883093E-11,1.342743980891829E-10,1.795423847433803E-10,2.399182404144058E-10,3.223748277206474E-10,4.389561820897720E-10,6.124121482079199E-10,8.914923932266923E-10,1.402946051601714E-09,2.611519948669723E-09,8.699935479443086E-09};
      i += dgdhAmp.test_2to2_amp2([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH);
      i += dgdhAmp.test_2to2_amp2_rotations([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH);
      i += dgdhAmp.test_2to2_amp2_boosts([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH);
      i += dgdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      //(CH produces the same M^2s since it is averaging over spins and the to helicity cases give the same result -- averaging them gaves the same as each individually.)    
      i += dgdhAmp.test_2to2_amp2([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH);
      i += dgdhAmp.test_2to2_amp2_rotations([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH);
      i += dgdhAmp.test_2to2_amp2_boosts([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH);
      i += dgdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH);
      //std::cout<<"\n# md=0.0042, mh=125, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {2.080127353619655E-11,2.291097168188958E-11,2.630338620549200E-11,3.121173826644254E-11,3.792943512979227E-11,4.683082399570715E-11,5.840116975318592E-11,7.328102205554832E-11,9.233373043528112E-11,1.167515394844657E-10,1.482286915174884E-10,1.892567193484229E-10,2.436559732899148E-10,3.175977884577619E-10,4.217417436778984E-10,5.762225989279183E-10,8.241862683613764E-10,1.277812128424104E-09,2.348424556865058E-09,7.737942869194889E-09};
      i += dgdhAmp.test_2to2_amp2([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH2);
      i += dgdhAmp.test_2to2_amp2_rotations([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH2);
      i += dgdhAmp.test_2to2_amp2_boosts([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH2);
      i += dgdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH2);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += dgdhAmp.test_2to2_amp2([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH2);
      i += dgdhAmp.test_2to2_amp2_rotations([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH2);
      i += dgdhAmp.test_2to2_amp2_boosts([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH2);
      i += dgdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH2);
      //std::cout<<"\n# md=125.1, mh=125, pspatial=95\n";
      md = 125.1;
      mh = 125;
      pspatial = 95;
      dgdhAmp.set_masses(md,mh,MW);
      ldouble dataCH4[20] = {2.363000175002506E-02,2.498263764024288E-02,2.633929332486995E-02,2.769882291889219E-02,2.905996597698369E-02,3.042133705837990E-02,3.178141427551381E-02,3.313852671726539E-02,3.449084062467721E-02,3.583634418228123E-02,3.717283077151347E-02,3.849788051377494E-02,3.980883990920000E-02,4.110279935272074E-02,4.237656828112354E-02,4.362664767294498E-02,4.484919958663686E-02,4.604001338071557E-02,4.719446821174936E-02,4.830749135102964E-02};
      i += dgdhAmp.test_2to2_amp2([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH4);
      i += dgdhAmp.test_2to2_amp2_rotations([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH4);
      i += dgdhAmp.test_2to2_amp2_boosts([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH4);
      i += dgdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH4);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += dgdhAmp.test_2to2_amp2([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH4);
      i += dgdhAmp.test_2to2_amp2_rotations([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH4);
      i += dgdhAmp.test_2to2_amp2_boosts([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH4);
      i += dgdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH4);
      //std::cout<<"\n# md=125, mh=0.0005, pspatial=125.1\n";
      md = 125;
      mh = 0.0005;
      pspatial = 125.1;
      dgdhAmp.set_masses(md,mh,MW);
      ldouble dataCH3[20] = {2.380941806091071E-03,7.574594414808863E-03,1.339933219798798E-02,1.993052674381923E-02,2.725539826307912E-02,3.547520371886823E-02,4.470781319878121E-02,5.509067901729470E-02,6.678411785238592E-02,7.997463346211629E-02,9.487758185836782E-02,1.117375476619851E-01,1.308227329427511E-01,1.524049674472023E-01,1.767059047655688E-01,2.037626734328878E-01,2.330944013919951E-01,2.628460148702128E-01,2.874366382839632E-01,2.903830653388305E-01};
      i += dgdhAmp.test_2to2_amp2([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH3);
      i += dgdhAmp.test_2to2_amp2_rotations([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH3);
      i += dgdhAmp.test_2to2_amp2_boosts([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH3);
      i += dgdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dgdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH3);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += dgdhAmp.test_2to2_amp2([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH3);
      i += dgdhAmp.test_2to2_amp2_rotations([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH3);
      i += dgdhAmp.test_2to2_amp2_boosts([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH3);
      i += dgdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dgdhAmp.amp2_gplus(); }, md,0,md,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
