
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2024 Neil Christensen

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

//File:  SPINAS/SM/uhWd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uhWd.h"

namespace spinas {
  
  uhWd::uhWd(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), mu(massu), md(massd), mh(massh), MW(massW), WW(widthW), SW(sinW), propu(massu,0), propd(massd,0), propW(massW,widthW) {
    p1=particle(mu);
    p2=particle(mh);
    p3=particle(MW);
    p4=particle(md);
    //<13>,[43],<14>,[14],[323>
    a13a = sproduct(ANGLE,&p1,&p3);
    s43s = sproduct(SQUARE,&p4,&p3);
    a14a = sproduct(ANGLE,&p1,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    s323a = sproduct(SQUARE,&p3,&p2,&p3);
    //<13>,[43],[324>
    s324a = sproduct(SQUARE,&p3,&p2,&p4);
    //[43],<13>,[123>
    s123a = sproduct(SQUARE,&p1,&p2,&p3);
    //Couplings
    preW = e*e/(2.0*MW*MW*SW*SW);
    preu = mu*preW;
    pred = md*preW;
  }
  void uhWd::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& massW){
    mu=massu;
    md=massd;
    mh=massh;
    MW=massW;
    p1.set_mass(mu);
    p2.set_mass(mh);
    p3.set_mass(MW);
    p4.set_mass(md);
    propu.set_mass(mu);
    propd.set_mass(md);
    propW.set_mass(MW);
    //Couplings
    preW = e*e/(2.0*MW*MW*SW*SW);
    preu = mu*preW;
    pred = md*preW;
  }
  void uhWd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4); 
    //<13>,[43],<14>,[14],[323>
    a13a.update();
    s43s.update();
    a14a.update();
    s14s.update();
    s323a.update();
    //<13>,[43],[324>
    s324a.update();
    //[43],<13>,[123>
    s123a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenU=propW.denominator(propUP);
    pDenS=propu.denominator(propSP);
    pDenT=propd.denominator(propTP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uhWd::amp(const int& ds1, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;
    constexpr ldouble two=2;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      //U-Channel W
      //preW = e*e/(2.0*MW*MW*SW*SW);
      //uDW-h: all ingoing:
      // - preW (2MW^2 <13>[23]+(Md<12>-Mu[12])[343>))/(s-MW^2)
      //uhW-D: 2<->4
      // - preW (2MW^2 <13>[43]+(Md<14>-Mu[14])[323>))/(u-MW^2)
      //34 outgoing:
      // + preW (2MW^2 <13>[43]-(Md<14>+Mu[14])[323>))/(u-MW^2)
      amplitude += normFactor*preW*(two*MW*MW*a13a.v(ds1,ds3a)*s43s.v(ds4,ds3b)-(md*a14a.v(ds1,ds4)+mu*s14s.v(ds1,ds4))*s323a.v(ds3a,ds3b))/pDenU;
      
      //T-Channel d
      //pred = e*e*md/(2.0*MW*MW*SW*SW);
      //uDW-h: all ingoing:
      // - pred <13>(2Md[23]+[342>)/(t-Md^2)
      //uhW-D: 2<->4
      // - pred <13>(2Md[43]+[324>)/(t-Md^2)
      //34 outgoing:
      // + pred <13>(2Md[43]-[324>)/(t-Md^2)
      amplitude += normFactor*pred*a13a.v(ds1,ds3a)*(two*md*s43s.v(ds4,ds3b)-s324a.v(ds3b,ds4))/pDenT;
      
      //S-Channel u
      //preu = e*e*mu/(2.0*MW*MW*SW*SW);
      //uDW-h: all ingoing:
      // - preu [23](2Mu<13>+[143>)/(u-Mu^2)
      //uhW-D: 2<->4
      // - preu [43](2Mu<13>+[123>)/(s-Mu^2)
      //34 outgoing:
      // + preu [43](2Mu<13>+[123>)/(s-Mu^2)
      amplitude += normFactor*preu*s43s.v(ds4,ds3a)*(two*mu*a13a.v(ds1,ds3b)+s123a.v(ds1,ds3b))/pDenS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uhWd::amp2(){
    constexpr ldouble three=3;
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,j3,j4);
	  //Color Factor 3
	  amp2 += three*std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2
    //Average over colors 1/3
    return amp2/6.0;
  }

  
  

  //  Tests
  int test_uhWd(){
    int n=0;//Number of fails
    std::cout<<"\t* u , h  -> W+, d       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, md=0.0075, mh=125, MW=80.385, pspatial=2500\n";
      ldouble mu=0.0042, md=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      uhWd uhWdAmp = uhWd(EE,mu,md,mh,MW,0,SW);
      ldouble pspatial=2500;
      ldouble dataCH[20] = {1.281089040268323E-03,4.212210528617322E-03,7.823880022244186E-03,1.230642363900492E-02,1.791757056972479E-02,2.501242549080821E-02,3.408997847394889E-02,4.586735030242093E-02,6.140223175130074E-02,8.230254456568008E-02,1.111016311196162E-01,1.519657294828779E-01,2.121151110908376E-01,3.049106596821528E-01,4.572564381538664E-01,7.298317696896108E-01,1.283399334826865E+00,2.664675930272234E+00,7.803116652117112E+00,7.301589364175759E+01};
      i += uhWdAmp.test_2to2_amp2([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH);
      i += uhWdAmp.test_2to2_amp2_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH);
      //std::cout<<"\n#mu=0.0042, md=0.0075, mh=125, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH1[20] = {3.619173889304420E-03,6.744857504577909E-03,1.056211199344819E-02,1.525730264800346E-02,2.108047991270789E-02,2.837267836206452E-02,3.760771090991696E-02,4.945787588649442E-02,6.490044741196635E-02,8.539640728401239E-02,1.132028881899329E-01,1.519462049278911E-01,2.077349681482667E-01,2.914794433992084E-01,4.241676832472087E-01,6.503089215989875E-01,1.077853687385026E+00,2.026010447685343E+00,4.805944234672729E+00,2.021845334284789E+01};
      i += uhWdAmp.test_2to2_amp2([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH1);
      i += uhWdAmp.test_2to2_amp2_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH1);
      i += uhWdAmp.test_2to2_amp2_boosts([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH1);
      i += uhWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH1);
      //std::cout<<"\n# mu=0.0042, md=0.0075, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {8.887191860397957E-03,1.235456627569828E-02,1.651014737654348E-02,2.152327155935699E-02,2.761628591117684E-02,3.508518627163063E-02,4.433029241090252E-02,5.590274465297786E-02,7.057680659596022E-02,8.946561250497702E-02,1.142127440270964E-01,1.473217993745582E-01,1.927496278176588E-01,2.570333498530543E-01,3.515761310538839E-01,4.976745449569000E-01,7.387811655424102E-01,1.174818947177576E+00,2.078640032346669E+00,4.427814658204641E+00};
      i += uhWdAmp.test_2to2_amp2([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH2);
      i += uhWdAmp.test_2to2_amp2_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH2);
      i += uhWdAmp.test_2to2_amp2_boosts([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH2);
      i += uhWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, md=130, mh=125, MW=80.385, pspatial=1\n";
      mu = 125.1;
      md = 130;
      mh = 125;
      pspatial = 1;
      uhWdAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH4[20] = {9.829054515825820E-01,9.839722031459947E-01,9.850449563576438E-01,9.861237376844983E-01,9.872085737727985E-01,9.882994914493873E-01,9.893965177230490E-01,9.904996797858635E-01,9.916090050145698E-01,9.927245209719439E-01,9.938462554081896E-01,9.949742362623392E-01,9.961084916636689E-01,9.972490499331275E-01,9.983959395847765E-01,9.995491893272438E-01,1.000708828065193E+00,1.001874884900800E+00,1.003047389135252E+00,1.004226370270253E+00};
      i += uhWdAmp.test_2to2_amp2([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH4);
      i += uhWdAmp.test_2to2_amp2_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH4);
      i += uhWdAmp.test_2to2_amp2_boosts([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH4);
      i += uhWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, md=130, mh=0.0005, MW=80.385, pspatial=250\n";
      mu = 125;
      mh = 0.0005;
      pspatial = 250;
      uhWdAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH3[20] = {2.410166052470733E+00,1.507167249632516E+00,1.087778265294518E+00,8.501336684717467E-01,7.009190747811114E-01,6.020190538379218E-01,5.352188351889656E-01,4.910429556347564E-01,4.645187708794406E-01,4.534363250201062E-01,4.577344191516176E-01,4.796096588891750E-01,5.244182692915614E-01,6.029688178774822E-01,7.369896900953757E-01,9.731719595434225E-01,1.424686543489094E+00,2.423284702185716E+00,5.307638330309302E+00,2.063364838107693E+01};
      i += uhWdAmp.test_2to2_amp2([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH3);
      i += uhWdAmp.test_2to2_amp2_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH3);
      i += uhWdAmp.test_2to2_amp2_boosts([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH3);
      i += uhWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH3);
      //std::cout<<"\n# mu=125, md=130, mh=0.0005, MW=80.385, pspatial=70\n";
      pspatial = 70;
      ldouble dataCH9[20] = {4.711790968820574E-01,4.723603582884213E-01,4.742098401760201E-01,4.767448168111300E-01,4.799865286231537E-01,4.839604623619190E-01,4.886966872064016E-01,4.942302539786073E-01,5.006016661484326E-01,5.078574331378771E-01,5.160507186177905E-01,5.252420991273210E-01,5.355004515485616E-01,5.469039918801274E-01,5.595414925564576E-01,5.735137114892190E-01,5.889350733659310E-01,6.059356529215718E-01,6.246635214140601E-01,6.452875320573878E-01};
      i += uhWdAmp.test_2to2_amp2([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH9);
      i += uhWdAmp.test_2to2_amp2_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH9);
      i += uhWdAmp.test_2to2_amp2_boosts([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH9);
      i += uhWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH9);
      //std::cout<<"\n# mu=0.004, md=130, mh=0.0005, MW=0.0004, pspatial=250\n";
      mu = 0.004;
      mh = 0.0005;
      MW = 0.0004;
      pspatial = 250;
      uhWdAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH5[20] = {1.362984343411446E+21,8.349768622829917E+20,5.992642264708001E+20,4.666320733261463E+20,3.818172115899212E+20,3.229835892534492E+20,2.798076897848471E+20,2.467860055858777E+20,2.207197431881196E+20,1.996241666226298E+20,1.822031605027509E+20,1.675747049674899E+20,1.551178789948168E+20,1.443829804971412E+20,1.350363582161110E+20,1.268252699506375E+20,1.195547740933162E+20,1.130721091189968E+20,1.072558765896478E+20,1.020083891411462E+20};
      i += uhWdAmp.test_2to2_amp2([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH5);
      i += uhWdAmp.test_2to2_amp2_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH5);
      i += uhWdAmp.test_2to2_amp2_boosts([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH5);
      i += uhWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH5);
      //std::cout<<"\n# mu=0.004, md=130, mh=0.0005, MW=0.0004, pspatial=130\n";
      pspatial = 130;
      ldouble dataCH6[20] = {4.608409597839037E+20,3.548916845815366E+20,2.816842759433706E+20,2.289962200427394E+20,1.898184305697152E+20,1.598977171836805E+20,1.365317953711206E+20,1.179370152173278E+20,1.028975890561853E+20,9.056168933991540E+19,8.031811999449506E+19,7.171912591897810E+19,6.443061890871122E+19,5.819923059212668E+19,5.282999163903616E+19,4.817089350567368E+19,4.410201000756703E+19,4.052770205108339E+19,3.737094584877986E+19,3.456914843707228E+19};
      i += uhWdAmp.test_2to2_amp2([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH6);
      i += uhWdAmp.test_2to2_amp2_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH6);
      i += uhWdAmp.test_2to2_amp2_boosts([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH6);
      i += uhWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH6);
      //std::cout<<"\n# mu=0.004, md=0.007, mh=0.0005, MW=0.0004, pspatial=250\n";
      mu = 0.004;
      md = 0.007;
      mh = 0.0005;
      MW = 0.0004;
      pspatial = 2.5;
      uhWdAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH7[20] = {5.824156182170629E+04,1.848325401943568E+04,1.055108315780976E+04,7.165874397816486E+03,5.296410331615451E+03,4.115994188897035E+03,3.306679564801314E+03,2.720111206497459E+03,2.277773948367656E+03,1.934246174836015E+03,1.661457718187594E+03,1.441171864092095E+03,1.261095418341858E+03,1.112755283049519E+03,9.903340989475400E+02,8.901511909750152E+02,8.108899430020247E+02,7.557969171487584E+02,7.456701423728682E+02,1.016469481318160E+03};
      i += uhWdAmp.test_2to2_amp2([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH7);
      i += uhWdAmp.test_2to2_amp2_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH7);
      i += uhWdAmp.test_2to2_amp2_boosts([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH7);
      i += uhWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH7);
      //std::cout<<"\n# mu=0.004, md=0.007, mh=0.0005, MW=0.0004, pspatial=1\n";
      pspatial = 0.003;
      ldouble dataCH8[20] = {9.081260053998093E+02,8.535642367351599E+02,8.024835588110154E+02,7.547673319657013E+02,7.103284075356523E+02,6.691136397675330E+02,6.311108924048177E+02,5.963598252614419E+02,5.649686149747388E+02,5.371403460371773E+02,5.132158197251651E+02,4.937455663259024E+02,4.796167135651181E+02,4.722898749744272E+02,4.742752890035154E+02,4.901857786384407E+02,5.293848329300056E+02,6.139878482967816E+02,8.111392885737976E+02,1.460119294990857E+03};
      i += uhWdAmp.test_2to2_amp2([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH8);
      i += uhWdAmp.test_2to2_amp2_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH8);
      i += uhWdAmp.test_2to2_amp2_boosts([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH8);
      i += uhWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uhWdAmp.amp2(); }, mu,mh,MW,md,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}
