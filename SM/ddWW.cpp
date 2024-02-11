
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

//File:  SPINAS/SM/ddWW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ddWW.h"

namespace spinas {

  ddWW::ddWW(const ldouble& echarge, const ldouble& massd, const ldouble& massu, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), Qf(-1.0/3.0), md(massd), mu(massu), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propu(mu,0), propA(0,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);
    proph = propagator(mh,wh);  
    p1=particle(md);
    p2=particle(md);
    p3=particle(MW);
    p4=particle(MW);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    //[413>,[314>,[231>,[132>
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    s231a = sproduct(SQUARE,&p2,&p3,&p1);
    s132a = sproduct(SQUARE,&p1,&p3,&p2);
    //Couplings
    pred = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preh = 2.0*e*e*md/(4.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
    gL=-2.0*Qf*SW*SW-1.0;
    gR=-2.0*Qf*SW*SW;
  }
  void ddWW::set_masses(const ldouble& massd, const ldouble& massu, const ldouble& massh, const ldouble& massW){
    md=massd;
    mu=massu;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(md);
    p2.set_mass(md);
    p3.set_mass(MW);
    p4.set_mass(MW);
    propu.set_mass(mu);
    proph.set_mass(mh);
    propZ.set_mass(MZ);
    //Couplings
    pred = e*e/(MW*MW*SW*SW);
    preh = e*e*md/(2.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
  }
  void ddWW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s34s.update();
    a34a.update();
    s24s.update();
    a24a.update();
    s14s.update();
    a14a.update();
    //[413>,[314>,[231>,[132>
    s413a.update();
    s314a.update();
    s231a.update();
    s132a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenhS=proph.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenuT=propu.denominator(propTP);
    pDenuU=propu.denominator(propUP);
    pDenZS=propZ.denominator(propSP);
    pDenAS=propA.denominator(propSP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ddWW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds4a, ds4b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);
      
      //S-Channel h
      //preh = e*e*md/(4.0*MW*MW*SW*SW);
      //all ingoing: 
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //34 outgoing:
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      amplitude += normFactor*preh*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))/pDenhS;
      
      //T-Channel u
      //pred = e*e/(MW*MW*SW*SW);
      //all in:
      // - pred <24>[13]( MW[34]+[413>)/(u-mu^2)
      // 34 out:
      // - pred <24>[13](-MW[34]+[413>)/(u-mu^2)
      amplitude += - normFactor*pred*a24a.v(ds2,ds4a)*s13s.v(ds1,ds3a)*(-MW*s34s.v(ds3b,ds4b)+s413a.v(ds4b,ds3b))/pDenuT;

      //S-Channel A
      //all in:
      //((2(<24>[13]+<23>[14]+<14>[23]+<13>[24])(<34>+[34])EE^2 )/(3MWs12) +(2([132>+[231>)EE^2 <34>[34])/(3MW^2 s12) )
      amplitude +=
	- normFactor*e*e*Qf*two/MW*(
				    a13a.v(ds1,ds3a)*s24s.v(ds2,ds4a) + s13s.v(ds1,ds3a)*a24a.v(ds2,ds4a)
				    + a14a.v(ds1,ds4a)*s23s.v(ds2,ds3a) + s14s.v(ds1,ds4a)*a23a.v(ds2,ds3a)
				    )*(s34s.v(ds3b,ds4b)+a34a.v(ds3b,ds4b))/pDenAS
	- normFactor*e*e*Qf*two/MW/MW*(
				       s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s231a.v(ds2,ds1)+s132a.v(ds1,ds2))
				       )/pDenAS;
      
      //S-Channel Z
      //preZ = e*e/(2.0*MW*MW*SW*SW); //=pred
      //all ingoing
      // = preZ ( (gR-gL)md[34]<34>([12]-<12>)
      //          + 2[34]<34>(gR[231>+gL[132>)
      //          + 2MW([34]+<34>)( gR([23]<14>+[24]<13>) + gL([13]<24>+[14]<23>) )
      //        )/(s-MZ^2)
      //34 outgoing
      //preZ ( (gR-gL)md[34]<34>([12]-<12>)
      //        - 2[34]<34>(gR[231>+gL[132>)
      //        - 2MW([34]+<34>)( gR([23]<14>+[24]<13>) + gL([13]<24>+[14]<23>) )
      //      )/(s-MZ^2)
      amplitude += normFactor*preZ*( (gR-gL)*md*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s12s.v(ds1,ds2)-a12a.v(ds1,ds2))
				       -two*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(gR*s231a.v(ds2,ds1)+gL*s132a.v(ds1,ds2))
				       -two*MW*(s34s.v(ds3a,ds4a)+a34a.v(ds3a,ds4a))*( gR*(s23s.v(ds2,ds3b)*a14a.v(ds1,ds4b)+s24s.v(ds2,ds4b)*a13a.v(ds1,ds3b))
										       + gL*(s13s.v(ds1,ds3b)*a24a.v(ds2,ds4b)+s14s.v(ds1,ds4b)*a23a.v(ds2,ds3b)) )
										       )/pDenZS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ddWW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3^2=1/9
    return amp2/36.0;
  }

  
  

  //  Tests
  int test_ddWW(){
    int n=0;//Number of fails
    std::cout<<"\t* D , d  -> W+, W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# md=0.0075, mu=0.0042, mh=125, MW=80.385, pspatial=250\n";
      ldouble md=0.0075, mu=0.0042, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      ddWW ddWWAmp = ddWW(EE,md,mu,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {6.247097351524106E-01,1.885450765693643E-01,9.958530248122271E-02,6.251856021462568E-02,4.285684978731518E-02,3.106473119831614E-02,2.344913178986098E-02,1.827918007480336E-02,1.463582012949063E-02,1.198581450749802E-02,9.998516633628227E-03,8.458004628637380E-03,7.217295321824076E-03,6.172864907320081E-03,5.249692192739634E-03,4.392080985800483E-03,3.557814071160806E-03,2.714308297361198E-03,1.836014755615100E-03,9.026177431286505E-04};
      i += ddWWAmp.test_2to2_amp2([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH);
      i += ddWWAmp.test_2to2_amp2_rotations([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH);
      i += ddWWAmp.test_2to2_amp2_boosts([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH);
      i += ddWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH);
      //std::cout<<"\n# md=0.0075, mu=0.0042, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {4.159858443977861E-01,2.070217674308176E-01,1.282211664402882E-01,8.848586648543240E-02,6.503998661921200E-02,4.984972364514384E-02,3.937985139345815E-02,3.183530757663921E-02,2.620868715814048E-02,2.189119761562672E-02,1.849334814969863E-02,1.575440008172255E-02,1.349324432738830E-02,1.158016763915240E-02,9.919833135872193E-03,8.440602351394036E-03,7.087603852369650E-03,5.818099775258145E-03,4.598308173874532E-03,3.401174162584671E-03};
      i += ddWWAmp.test_2to2_amp2([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH2);
      i += ddWWAmp.test_2to2_amp2_rotations([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH2);
      i += ddWWAmp.test_2to2_amp2_boosts([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH2);
      i += ddWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH2);
      //std::cout<<"\n# md=125.1, mh=125, MW=80.385, pspatial=95\n";
      md = 125.1;
      mh = 125;
      pspatial = 95;
      ddWWAmp.set_masses(md,mu,mh,MW);
      ldouble dataCH4[20] = {1.008318920694584E+00,6.871652702789750E-01,4.869695018300897E-01,3.670549181898975E-01,2.887306192845187E-01,2.340150127385204E-01,1.938546034187554E-01,1.632598231335534E-01,1.392701516306541E-01,1.200225404520403E-01,1.042876782564961E-01,9.122232882860383E-02,8.022910966204171E-02,7.087312707269798E-02,6.283032424395209E-02,5.585432824333812E-02,4.975451172383762E-02,4.438108446922223E-02,3.961472243345002E-02,3.535920219254238E-02};
      i += ddWWAmp.test_2to2_amp2([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH4);
      i += ddWWAmp.test_2to2_amp2_rotations([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH4);
      i += ddWWAmp.test_2to2_amp2_boosts([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH4);
      i += ddWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH4);
      //std::cout<<"\n# md=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      md = 125;
      mh = 0.0005;
      pspatial = 125.1;
      ddWWAmp.set_masses(md,mu,mh,MW);
      ldouble dataCH3[20] = {1.661242989082962E+00,8.453598922542708E-01,5.398263675443259E-01,3.854729742574838E-01,2.931228148801160E-01,2.320038001062697E-01,1.887661336742280E-01,1.566970578014326E-01,1.320559247700869E-01,1.125943383478283E-01,9.687990465295451E-02,8.395699891426170E-02,7.316370999645008E-02,6.402708417342735E-02,5.620020642583147E-02,4.942284679360014E-02,4.349598214024042E-02,3.826480358533749E-02,3.360708632142176E-02,2.942504649798380E-02};
      i += ddWWAmp.test_2to2_amp2([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH3);
      i += ddWWAmp.test_2to2_amp2_rotations([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH3);
      i += ddWWAmp.test_2to2_amp2_boosts([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH3);
      i += ddWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddWWAmp.amp2(); }, md,md,MW,MW,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}
