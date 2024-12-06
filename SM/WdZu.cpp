
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

//File:  SPINAS/SM/WdZu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/WdZu.h"

namespace spinas {

  WdZu::WdZu(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), Qu(2.0/3.0), Qd(-1.0/3.0), mu(massu), md(massd), MW(massW), SW(sinW), WW(widthW) {
    //constexpr ldouble sqrt2 = std::sqrt(2);
    CW=std::sqrt(1.0-SW*SW);
    MZ=MW/CW;
    propW = propagator(MW,WW);
    propu = propagator(mu,0);
    propd = propagator(md,0);
    p1=particle(MW);
    p2=particle(md);
    p3=particle(MZ);
    p4=particle(mu);
    //Spinor Products
    a42a= sproduct(ANGLE,&p4,&p2,2);
    a43a= sproduct(ANGLE,&p4,&p3,2);
    a23a= sproduct(ANGLE,&p2,&p3,2);
    a21a= sproduct(ANGLE,&p2,&p1,2);
    a31a= sproduct(ANGLE,&p3,&p1,2);
    s42s= sproduct(SQUARE,&p4,&p2,2);
    s43s= sproduct(SQUARE,&p4,&p3,2);
    s41s= sproduct(SQUARE,&p4,&p1,2);
    s23s= sproduct(SQUARE,&p2,&p3,2);
    s31s= sproduct(SQUARE,&p3,&p1,2);
    s432a= sproduct(SQUARE,&p4,&p3,&p2,2);
    s341a= sproduct(SQUARE,&p3,&p4,&p1,2);
    s143a= sproduct(SQUARE,&p1,&p4,&p3,2);
    //prefactor
    preud = e*e/(std::sqrt(2.0)*MW*MW*SW*SW);
    preW = preud/(MZ*MZ);
    gLu=-2.0*Qu*SW*SW+1.0;
    gRu=-2.0*Qu*SW*SW;
    gLd=-2.0*Qd*SW*SW-1.0;
    gRd=-2.0*Qd*SW*SW;
  }
  void WdZu::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massW){
    //constexpr ldouble sqrt2 = std::sqrt(2);
    mu=massu;
    md=massd;
    MW=massW;
    MZ=MW/CW;
    p1.set_mass(MW);
    p2.set_mass(md);
    p3.set_mass(MZ);
    p4.set_mass(mu);
    propW.set_mass(MW);
    propu.set_mass(mu);
    propd.set_mass(md);
    preud = e*e/(std::sqrt(2.0)*MW*MW*SW*SW);
    preW = preud/(MZ*MZ);
  }
  void WdZu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    constexpr ldouble one=1, two=2, three=3;
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //Spinor Products
    a42a.update();
    a43a.update();
    a23a.update();
    a21a.update();
    a31a.update();
    s42s.update();
    s43s.update();
    s41s.update();
    s23s.update();
    s31s.update();
    s432a.update();
    s341a.update();
    s143a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propu.denominator(propSP);
    pDenT=propW.denominator(propTP);
    pDenU=propd.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble WdZu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    constexpr ldouble two=2;
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds1a, ds1b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds1);
    ldouble normFactor=get_spin_normalization(ds3,ds1);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds1,ds1a,ds1b, i);
      
      
      //W Diagram
      //preW = e*e/(sqrt2*MW*MW*MZ*MZ*SW*SW);
      //UdZW all ingoing:
      // - preW (2(MZ<34>+MW[34])MW^2 <23>[14]+2(MW<34>+MZ[34])MW^2 <24>[13]+<34>[34]((2MW^2 -MZ^2)(Mu<12>-Md[12])+2MW^2 [132>)))/(s-MW^2)
      //WdZU all ingoing: 1<->4
      // - preW (2(MZ<31>+MW[31])MW^2 <23>[41]+2(MW<31>+MZ[31])MW^2 <21>[43]+<31>[31]((2MW^2 -MZ^2)(Mu<42>-Md[42])+2MW^2 [432>)))/(t-MW^2)
      //34 out:
      // - preW (
      //         +2(+MZ<31>-MW[31])MW^2 <23>[41]
      //         +2(-MW<31>+MZ[31])MW^2 <21>[43]
      //         +<31>[31]((2MW^2 -MZ^2)(Mu<42>+Md[42])+2MW^2 [432>))
      //        )/(t-MW^2)
      amplitude += - normFactor*preW*(
				      +two*(+MZ*a31a.v(ds3a,ds1b)-MW*s31s.v(ds3a,ds1b))*MW*MW*a23a.v(ds2,ds3b)*s41s.v(ds4,ds1a)
				      +two*(-MW*a31a.v(ds3a,ds1b)+MZ*s31s.v(ds3a,ds1b))*MW*MW*a21a.v(ds2,ds1a)*s43s.v(ds4,ds3b)
				      +a31a.v(ds3a,ds1a)*s31s.v(ds3b,ds1b)*(
									    +(two*MW*MW-MZ*MZ)*(mu*a42a.v(ds4,ds2)+md*s42s.v(ds4,ds2))
									    +two*MW*MW*s432a.v(ds4,ds2)
									    )
				      )/pDenT;
      

      //u Diagram
      //preu = e*e/(sqrt2*MW*MW*SW*SW);
      //UdZW all in:
      // - preu<24>(gRu*mu<13>[34]+(MZ[34]+[413>)gLu[13]))/(t-mu^2)
      //WdZU all in: 1<->4
      // - preu<21>(gRu*mu<43>[31]+(MZ[31]+[143>)gLu[43]))/(s-mu^2)
      //34 out:
      // - preu<21>(gRu*mu<43>[31]+(MZ[31]+[143>)gLu[43]))/(s-mu^2)
      amplitude += - normFactor*preud*a21a.v(ds2,ds1a)*(
						       +gRu*mu*a43a.v(ds4,ds3b)*s31s.v(ds3a,ds1b)
						       +gLu*s43s.v(ds4,ds3a)*(MZ*s31s.v(ds3b,ds1b)+s143a.v(ds1b,ds3b))
						       )/pDenS;
      
      //d Diagram
      //pred = e*e/(sqrt2*MW*MW*SW*SW);
      //UdZW all in:
      // + pred[14](gRd*md<34>[23]+(MW[34]-[314>)gLd<23>))/(u-md^2)
      //WdZU all in: 1<->4 
      // + pred[41](gRd*md<31>[23]+(MW[31]-[341>)gLd<23>))/(u-md^2)
      //34 out:
      // - pred[41](gRd*md<31>[23]+(MW[31]+[341>)gLd<23>))/(u-md^2)
      amplitude += - normFactor*preud*s41s.v(ds4,ds1a)*(
						       +gRd*md*a31a.v(ds3a,ds1b)*s23s.v(ds2,ds3b)
						       +gLd*a23a.v(ds2,ds3a)*(MW*s31s.v(ds3b,ds1b)+s341a.v(ds3b,ds1b))
						       )/pDenU;
      
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble WdZu::amp2(){
    constexpr ldouble three=3;
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    //Color factor 3
	    amp2 += three*std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/3=1/6
    //Average over initial colors 1/3
    return amp2/18.0;
  }



  



  //  Tests
  int test_WdZu(){
    int n=0;//Number of fails
    std::cout<<"\t* W+, d  -> Z , u       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,mu=0.0042,md=0.0075,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;
      WdZu WdZuAmp = WdZu(EE,mu,md,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.567839010866377E+01,1.115027151393418E+01,4.779407678341179E+00,2.588001565836544E+00,1.595450630435823E+00,1.069192860370316E+00,7.601593534717510E-01,5.653120707543174E-01,4.359736845576661E-01,3.468397070176106E-01,2.837902759428861E-01,2.385355392226418E-01,2.060686561432385E-01,1.834055969150515E-01,1.690093933725599E-01,1.627291752274924E-01,1.665231821311190E-01,1.876410807574028E-01,2.555182037681686E-01,6.573849031087551E-01};
      i += WdZuAmp.test_2to2_amp2([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH);
      i += WdZuAmp.test_2to2_amp2_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH);
      i += WdZuAmp.test_2to2_amp2_boosts([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH);
      i += WdZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=90\n";
      pspatial = 90;
      ldouble dataCH2[20] = {3.096558839904828E+00,2.023405668681739E+00,1.408780331641348E+00,1.027270457782631E+00,7.761875699200338E-01,6.034882625355010E-01,4.806050820388220E-01,3.908615075119663E-01,3.240466791526461E-01,2.736842226133623E-01,2.355832910118321E-01,2.070532758030837E-01,1.865081986851646E-01,1.733757113521537E-01,1.684241591760230E-01,1.751405913971914E-01,2.052026399095969E-01,3.093387564084727E-01,1.022313151760880E+00,3.390691125261484E+00};    
      i += WdZuAmp.test_2to2_amp2([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH2);
      i += WdZuAmp.test_2to2_amp2_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH2);
      i += WdZuAmp.test_2to2_amp2_boosts([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH2);
      i += WdZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH2);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=250\n";
      mu=80;
      WdZuAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH3[20] = {4.534710016843481E+01,1.121814436509207E+01,4.820053179037382E+00,2.611132832804994E+00,1.610478580363463E+00,1.080760761301974E+00,7.705828797882556E-01,5.758298241802273E-01,4.473198642468194E-01,3.595215328118330E-01,2.982386386913134E-01,2.552049410712101E-01,2.255279181363698E-01,2.064542853401903E-01,1.968724338234822E-01,1.974476307873697E-01,2.118815579902903E-01,2.518934529149709E-01,3.627991072660263E-01,9.571074309198674E-01};
      i += WdZuAmp.test_2to2_amp2([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH3);
      i += WdZuAmp.test_2to2_amp2_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH3);
      i += WdZuAmp.test_2to2_amp2_boosts([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH3);
      i += WdZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH3);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=70\n";
      pspatial = 70;
      ldouble dataCH4[20] = {6.639803371742666E-01,6.218731704907989E-01,5.847457948744046E-01,5.520452501537766E-01,5.233159534551159E-01,4.981858597722506E-01,4.763562930270449E-01,4.575949880922612E-01,4.417321462317235E-01,4.286595813280469E-01,4.183333692298557E-01,4.107808783403112E-01,4.061137759031871E-01,4.045497902082069E-01,4.064480885319410E-01,4.123669721161715E-01,4.231600705517166E-01,4.401426643778438E-01,4.653939058749514E-01,5.023426302860564E-01};
      i += WdZuAmp.test_2to2_amp2([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH4);
      i += WdZuAmp.test_2to2_amp2_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH4);
      i += WdZuAmp.test_2to2_amp2_boosts([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH4);
      i += WdZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH4);
      //std::cout<<"\n# mu=80, md=40, MW=30, pspatial=250\n";
      mu=80;
      md=40;
      MW=30;
      MZ=MW/CW;
      WdZuAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH5[20] = {1.398652233467756E+02,1.736618029277250E+01,6.141633251744914E+00,3.056918936906525E+00,1.820098374109405E+00,1.216501571899369E+00,8.848254192576800E-01,6.885094542212998E-01,5.674537624590448E-01,4.924465809964916E-01,4.485133632127896E-01,4.281138805744102E-01,4.282728939825162E-01,4.496480620591736E-01,4.970587785056045E-01,5.821540924128490E-01,7.313976727121867E-01,1.012654054365424E+00,1.657324714090705E+00,4.253148179824843E+00};
      i += WdZuAmp.test_2to2_amp2([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH5);
      i += WdZuAmp.test_2to2_amp2_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH5);
      i += WdZuAmp.test_2to2_amp2_boosts([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH5);
      i += WdZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH5);
      //std::cout<<"\n# mu=80, md=40, MW=1, pspatial=70\n";
      mu=80;
      md=40;
      MW=1;
      MZ=MW/CW;
      WdZuAmp.set_masses(mu,md,MW);
      pspatial = 70;
      ldouble dataCH6[20] ={8.637852896768165E+04,8.076492075331710E+04,7.446222535071924E+04,6.832285646439328E+04,6.250605296812027E+04,5.711134584940790E+04,5.223937128636795E+04,4.801023661047519E+04,4.457680456501944E+04,4.214066592423459E+04,4.097535628054380E+04,4.146284488485294E+04,4.415398659572930E+04,4.987403438169479E+04,5.991824172557577E+04,7.644237501386610E+04,1.033199919190544E+05,1.482794209716351E+05,2.292873163232224E+05,3.998995679083028E+05};
      i += WdZuAmp.test_2to2_amp2([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH6);
      i += WdZuAmp.test_2to2_amp2_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH6);
      i += WdZuAmp.test_2to2_amp2_boosts([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH6);
      i += WdZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WdZuAmp.amp2(); }, MW,md,MZ,mu,pspatial,dataCH6);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }



}
