
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

//File:  SPINAS/SM/uguh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uguh.h"

namespace spinas {

  uguh::uguh(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& massh, const ldouble& massW, const ldouble& sinW):
    e(echarge), gs(gscharge), mu(massu), mh(massh), prop(massu,0), MW(massW), SW(sinW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(mu);
    p2=particle(0);
    p3=particle(mu);
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
    pre = sqrt2*e*gs*mu/(2.0*MW*SW);
  }
  void uguh::set_masses(const ldouble& massu, const ldouble& massh, const ldouble& massW){
    mu=massu;
    mh=massh;
    p1.set_mass(mu);
    p3.set_mass(mu);
    p4.set_mass(mh);
    prop.set_mass(mu);
    pre = sqrt2*e*gs*mu/(2.0*MW*SW);
  }
  void uguh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
  cdouble uguh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    /*mh^2[12][23] - mu[12][243> + mu[23][241> - <13>[2342]
      mh^2<12><23> - mu<12><243] + mu<23><241] - [13]<2342>
      Becomes after a sign change due to p3 and p4 being outgoing:*/
    //pre = sqrt2*e*gs*mu/(2.0*MW*SW);
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
  ldouble uguh::amp2(){
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
  ldouble uguh::amp2_gplus(){
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
  int test_uguh(){
    int n=0;//Number of fails
    std::cout<<"\t* u , g  -> u , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, mh=125, pspatial=250\n";
      ldouble mu=0.0042, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble gs=1.238;
      uguh uguhAmp = uguh(EE,gs,mu,mh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.713238975079426E-13,7.774293386089188E-13,1.638039785530422E-12,3.035792512051096E-12,5.074650540019468E-12,7.887256345584432E-12,1.164555387930521E-11,1.657650915591395E-11,2.298603289090959E-11,3.129657321520380E-11,4.210845115690813E-11,5.630449175731681E-11,7.523836007966027E-11,1.010967458412861E-10,1.376566585535739E-10,1.920524495147572E-10,2.795720143606642E-10,4.399638817308233E-10,8.189726565384397E-10,2.728299779755178E-09};
      i += uguhAmp.test_2to2_amp2([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uguhAmp.test_2to2_amp2_rotations([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uguhAmp.test_2to2_amp2_boosts([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uguhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      //(CH produces the same M^2s since it is averaging over spins and the to helicity cases give the same result -- averaging them gaves the same as each individually.)    
      i += uguhAmp.test_2to2_amp2([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uguhAmp.test_2to2_amp2_rotations([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uguhAmp.test_2to2_amp2_boosts([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH);
      i += uguhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, mh=125, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {6.523279377106051E-12,7.184880694499340E-12,8.248741866108954E-12,9.788001047229202E-12,1.189467075580607E-11,1.468614627334976E-11,1.831460666847438E-11,2.298092831173487E-11,2.895585761559382E-11,3.661328248293655E-11,4.648451730200710E-11,5.935090676128161E-11,7.641051271620056E-11,9.959866585586640E-11,1.322582100967924E-10,1.807034061675203E-10,2.584648127637925E-10,4.007218824597834E-10,7.364659411097801E-10,2.426618913118772E-09};
      i += uguhAmp.test_2to2_amp2([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uguhAmp.test_2to2_amp2_rotations([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uguhAmp.test_2to2_amp2_boosts([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uguhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH2);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += uguhAmp.test_2to2_amp2([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uguhAmp.test_2to2_amp2_rotations([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uguhAmp.test_2to2_amp2_boosts([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH2);
      i += uguhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, mh=125, pspatial=95\n";
      mu = 125.1;
      mh = 125;
      pspatial = 95;
      uguhAmp.set_masses(mu,mh,MW);
      ldouble dataCH4[20] = {2.363000175002506E-02,2.498263764024288E-02,2.633929332486995E-02,2.769882291889219E-02,2.905996597698369E-02,3.042133705837990E-02,3.178141427551381E-02,3.313852671726539E-02,3.449084062467721E-02,3.583634418228123E-02,3.717283077151347E-02,3.849788051377494E-02,3.980883990920000E-02,4.110279935272074E-02,4.237656828112354E-02,4.362664767294498E-02,4.484919958663686E-02,4.604001338071557E-02,4.719446821174936E-02,4.830749135102964E-02};
      i += uguhAmp.test_2to2_amp2([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uguhAmp.test_2to2_amp2_rotations([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uguhAmp.test_2to2_amp2_boosts([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uguhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH4);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += uguhAmp.test_2to2_amp2([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uguhAmp.test_2to2_amp2_rotations([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uguhAmp.test_2to2_amp2_boosts([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH4);
      i += uguhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, mh=0.0005, pspatial=125.1\n";
      mu = 125;
      mh = 0.0005;
      pspatial = 125.1;
      uguhAmp.set_masses(mu,mh,MW);
      ldouble dataCH3[20] = {2.380941806091071E-03,7.574594414808863E-03,1.339933219798798E-02,1.993052674381923E-02,2.725539826307912E-02,3.547520371886823E-02,4.470781319878121E-02,5.509067901729470E-02,6.678411785238592E-02,7.997463346211629E-02,9.487758185836782E-02,1.117375476619851E-01,1.308227329427511E-01,1.524049674472023E-01,1.767059047655688E-01,2.037626734328878E-01,2.330944013919951E-01,2.628460148702128E-01,2.874366382839632E-01,2.903830653388305E-01};
      i += uguhAmp.test_2to2_amp2([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uguhAmp.test_2to2_amp2_rotations([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uguhAmp.test_2to2_amp2_boosts([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uguhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uguhAmp.amp2(); }, mu,0,mu,mh,pspatial,dataCH3);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += uguhAmp.test_2to2_amp2([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uguhAmp.test_2to2_amp2_rotations([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uguhAmp.test_2to2_amp2_boosts([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH3);
      i += uguhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uguhAmp.amp2_gplus(); }, mu,0,mu,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
