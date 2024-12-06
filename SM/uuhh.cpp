
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

//File:  SPINAS/SM/uuhh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uuhh.h"

namespace spinas {
  //Constructors
  uuhh::uuhh(const ldouble& echarge, const ldouble& massu, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW):
    e(echarge), mu(massu), mh(massh), wh(widthh), MW(massW), SW(sinW),
    propu(mu,0), proph(mh,wh),
    p1(particle(mu)), p2(particle(mu)),
    p3(particle(mh)), p4(particle(mh)),
    s12s(sproduct(SQUARE,&p1,&p2,2)),
    a12a(sproduct(ANGLE,&p1,&p2,2)),
    s132a(sproduct(SQUARE,&p1,&p3,&p2,2)),
    s231a(sproduct(SQUARE,&p2,&p3,&p1,2))
  {
    ldouble sqrt2 = std::sqrt(2);
    prehS = 3.0*e*e*mu*mh*mh/(4.0*MW*MW*SW*SW);
    prehTU = e*e*mu*mu/(4.0*MW*MW*SW*SW);
  }
  void uuhh::set_masses(const ldouble& massu, const ldouble& massh, const ldouble& massW){
    ldouble sqrt2 = std::sqrt(2);
    mu=massu;
    propu.set_mass(mu);
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    p1.set_mass(mu);
    p2.set_mass(mu);
    p3.set_mass(mh);
    p4.set_mass(mh);
    prehS = 3.0*e*e*mu*mh*mh/(4.0*MW*MW*SW*SW);
    prehTU = e*e*mu*mu/(4.0*MW*MW*SW*SW);
  }
  void uuhh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s12s.update();
    a12a.update();
    s132a.update();
    s231a.update();
    //Propagator Momentum
    ldouble propPS[4], propPT[4], propPU[4];
    for(int j=0;j<4;j++){
      propPS[j] = mom1[j]+mom2[j];
      propPT[j] = mom1[j]-mom3[j];
      propPU[j] = mom1[j]-mom4[j];
    }
    pDenSh = proph.denominator(propPS);
    pDenTu = propu.denominator(propPT);
    pDenUu = propu.denominator(propPU);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uuhh::amp(const int& ds1, const int& ds2){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Higgs
    //prehS = 3*e*e*mu*mh*mh/(4*MW*MW*SW*SW);
    //S Channel
    //all ingoing: prehS*([12]+<12>)/(s-mh^2)
    amplitude += prehS*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))/pDenSh;

    //T Channel
    //prehTU = e*e*mu*mu/(4*MW*MW*SW*SW);
    //all ingoing:  prehTU * ( 2mu([12]+<12>) + [132> + [231> )/(t-mu^2)
    //34 outgoing:  prehTU * ( 2mu([12]+<12>) - [132> - [231> )/(t-mu^2)
    amplitude += prehTU*( two*mu*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2)) - s132a.v(ds1,ds2) - s231a.v(ds2,ds1) )/pDenTu;

    //U Channel
    //all ingoing:  prehTU * ( 2mu([12]+<12>) + [142> + [241> )/(u-mu^2)
    //all ingoing:  prehTU * ( 2mu([12]+<12>) - [132> - [231> )/(u-mu^2)
    //34 outgoing:  prehTU * ( 2mu([12]+<12>) + [132> + [231> )/(u-mu^2)
    amplitude += prehTU*( two*mu*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2)) + s132a.v(ds1,ds2) + s231a.v(ds2,ds1) )/pDenUu;

    //amplitude += prehTU*two*mu*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))*(one/pDenTu+one/pDenUu);
    //amplitude += prehTU*(s132a.v(ds1,ds2) + s231a.v(ds2,ds1))*(-one/pDenTu+one/pDenUu);
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uuhh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2){
	M = amp(j1,j2);
	amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over colors 1/3^2 = 1/9
    //Symmetry factor 1/2
    return amp2/72.0;
  }

  



  //  Tests
  int test_uuhh(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> h , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mu=0.0042, mh=125, MW=80.385, pspatial=250\n";
      ldouble mu=0.0042, mh=125, wh=0, MW=80.385, SW=0.474;//Set width to 0 for comparison with Feynman diagrams.
      uuhh uuhhAmp = uuhh(0.31333,mu,mh,wh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.201378802934037E-12,4.201378762779703E-12,4.201378748939706E-12,4.201378742360639E-12,4.201378738623534E-12,4.201378736298322E-12,4.201378734791739E-12,4.201378733819027E-12,4.201378733232569E-12,4.201378732955854E-12,4.201378732955854E-12,4.201378733232569E-12,4.201378733819027E-12,4.201378734791739E-12,4.201378736298322E-12,4.201378738623534E-12,4.201378742360639E-12,4.201378748939706E-12,4.201378762779703E-12,4.201378802934037E-12};
      i += uuhhAmp.test_2to2_amp2([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH);
      i += uuhhAmp.test_2to2_amp2_rotations([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH);
      i += uuhhAmp.test_2to2_amp2_boosts([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH);
      i += uuhhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH);
      //std::cout<<"\n########### mu=0.0042, mh=125, MW=80.385, pspatial=126\n";
      pspatial = 126;
      ldouble dataCH2[20] = {2.557321507774824E-11,2.557321507984529E-11,2.557321508166234E-11,2.557321508321725E-11,2.557321508452478E-11,2.557321508559705E-11,2.557321508644376E-11,2.557321508707241E-11,2.557321508748852E-11,2.557321508769569E-11,2.557321508769569E-11,2.557321508748852E-11,2.557321508707241E-11,2.557321508644376E-11,2.557321508559705E-11,2.557321508452478E-11,2.557321508321725E-11,2.557321508166234E-11,2.557321507984529E-11,2.557321507774824E-11};
      i += uuhhAmp.test_2to2_amp2([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH2);
      i += uuhhAmp.test_2to2_amp2_rotations([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH2);
      i += uuhhAmp.test_2to2_amp2_boosts([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH2);
      i += uuhhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH2);
      //std::cout<<"\n########### mu=0.1, mh=0.15, MW=80.385, pspatial=0.2\n";
      mu=0.10;
      mh=0.15;
      pspatial = 0.2;
      uuhhAmp.set_masses(mu,mh,MW);
      ldouble dataCH4[20] = {2.125152568394375E-14,1.299664724857766E-14,8.510336016507414E-15,5.751261332755809E-15,3.931952383290586E-15,2.688052005359802E-15,1.830253155688700E-15,1.252949479523918E-15,8.955061328505184E-16,7.243198866008038E-16,7.243198866008024E-16,8.955061328505171E-16,1.252949479523917E-15,1.830253155688699E-15,2.688052005359801E-15,3.931952383290583E-15,5.751261332755804E-15,8.510336016507411E-15,1.299664724857766E-14,2.125152568394374E-14};
      i += uuhhAmp.test_2to2_amp2([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH4);
      i += uuhhAmp.test_2to2_amp2_rotations([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH4);
      i += uuhhAmp.test_2to2_amp2_boosts([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH4);
      i += uuhhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH4);
      //std::cout<<"\n########### mu=0.1, mh=0.15, MW=0.11, pspatial=0.2\n";
      mu=0.10;
      mh=0.15;
      MW=0.11;
      pspatial = 0.2;
      uuhhAmp.set_masses(mu,mh,MW);
      ldouble dataCH5[20] = {6.060653260141943E-03,3.706471417133078E-03,2.427034956945145E-03,1.640183451517083E-03,1.121340669132170E-03,7.665967795443336E-04,5.219639248437857E-04,3.573245733414366E-04,2.553864717410923E-04,2.065664248016549E-04,2.065664248016545E-04,2.553864717410920E-04,3.573245733414363E-04,5.219639248437855E-04,7.665967795443335E-04,1.121340669132169E-03,1.640183451517081E-03,2.427034956945144E-03,3.706471417133079E-03,6.060653260141940E-03};
      i += uuhhAmp.test_2to2_amp2([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH5);
      i += uuhhAmp.test_2to2_amp2_rotations([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH5);
      i += uuhhAmp.test_2to2_amp2_boosts([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH5);
      i += uuhhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH5);
      //std::cout<<"\n########### mu=0.1, mh=0.15, MW=0.006, pspatial=0.2\n";
      mu=0.10;
      mh=0.15;
      MW=0.006;
      pspatial = 0.2;
      uuhhAmp.set_masses(mu,mh,MW);
      ldouble dataCH6[20] = {6.846761140566218E+02,4.187225927333750E+02,2.741837870727921E+02,1.852926382226976E+02,1.266786167960192E+02,8.660295871380083E+01,5.896661900955144E+01,4.036719967817881E+01,2.885118312315844E+01,2.333594927099560E+01,2.333594927099555E+01,2.885118312315839E+01,4.036719967817877E+01,5.896661900955141E+01,8.660295871380082E+01,1.266786167960192E+02,1.852926382226974E+02,2.741837870727920E+02,4.187225927333751E+02,6.846761140566214E+02};
      i += uuhhAmp.test_2to2_amp2([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH6);
      i += uuhhAmp.test_2to2_amp2_rotations([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH6);
      i += uuhhAmp.test_2to2_amp2_boosts([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH6);
      i += uuhhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuhhAmp.amp2(); }, mu,mu,mh,mh,pspatial,dataCH6);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
