
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

//File:  SPINAS/SM/AnWe.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AnWe.h"

namespace spinas {

  AnWe::AnWe(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), me(masse), MW(massW), SW(sinW), WW(widthW) {
    //constexpr ldouble sqrt2 = std::sqrt(2);
    propW = propagator(MW,WW);
    prope = propagator(me,0);
    p1=particle(0);
    p2=particle(0);
    p3=particle(MW);
    p4=particle(me);
    //[13],<42>,<23>,[41],[143>
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a42a = sproduct(ANGLE,&p4,&p2,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s41s = sproduct(SQUARE,&p4,&p1,2);
    s143a = sproduct(SQUARE,&p1,&p4,&p3,2);
    //<21>,[43],<13>,[341>
    a21a = sproduct(ANGLE,&p2,&p1,2);
    s43s = sproduct(SQUARE,&p4,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s341a = sproduct(SQUARE,&p3,&p4,&p1,2);
    //prefactor
    pre = std::sqrt(2.0)*e*e/(MW*SW);
  }
  void AnWe::set_masses(const ldouble& masse, const ldouble& massW){
    //constexpr ldouble sqrt2 = std::sqrt(2);
    me=masse;
    MW=massW;
    p3.set_mass(MW);
    p4.set_mass(me);
    propW.set_mass(MW);
    prope.set_mass(me);
    pre = std::sqrt(2.0)*e*e/(MW*SW);
  }
  void AnWe::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //[13],<42>,<23>,[41],[143>
    s13s.update();
    a42a.update();
    a23a.update();
    s41s.update();
    s143a.update();
    //<21>,[43],<13>,[341>
    a21a.update();
    s43s.update();
    a13a.update();
    s341a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT=propW.denominator(propTP);
    pDenU=prope.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AnWe::amp(const int& ds1, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);

      if(ds1>0){
      
	//ST Diagram
	//pre = e*e/(2.0*MW*SW);
	//EnAW- all ingoing:
	// -(me[34]<12>+MW<24>[13])(MW[34]-[314>))/((s-MW^2)(t-me^2))
	//AnW-E: 1->4->3->1
	// -(me[13]<42>+MW<23>[41])(MW[13]-[143>))/((t-MW^2)(u-me^2))
	//34 out:
	// +(me[13]<42>+MW<23>[41])(MW[13]-[143>))/((t-MW^2)(u-me^2))
	amplitude += normFactor*pre*(me*s13s.v(ds3a)*a42a.v(ds4)+MW*a23a.v(ds3a)*s41s.v(ds4))*(MW*s13s.v(ds3b)-s143a.v(ds3b))/pDenT/pDenU;
	
	
      }
      else if(ds1<0){

	//ST Diagram
	//pre = e*e/(2.0*MW*SW);
	//EnAW- all ingoing:
	// (<23>[14](-me^2 <34>+MW[413>))/((s-MW^2)(t-Me^2))
	//AnW-E: 1->4->3->1
	// (<21>[43](-me^2 <13>+MW[341>))/((t-MW^2)(u-Me^2))
	//34 out:
	// (<21>[43](me^2 <13>-MW[341>))/((t-MW^2)(u-Me^2))
	amplitude += normFactor*pre*a21a.v()*s43s.v(ds4,ds3a)*(me*me*a13a.v(ds3b)-MW*s341a.v(ds3b))/pDenT/pDenU;

	
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AnWe::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2
    return amp2/2.0;
  }


  
  //set_momenta(...) must be called before amp2().
  ldouble AnWe::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(1,j3,j4);
	amp2 += std::pow(std::abs(M),2);
      }
    //
    return amp2;
  }


  



  //  Tests
  int test_AnWe(){
    int n=0;//Number of fails
    std::cout<<"\t* A , ne -> W+, e       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,me=0.0005,MW=80.385, SW=0.474;
      AnWe AnWeAmp = AnWe(EE,me,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.406809722358972E+01,8.799524775934520E+00,3.973227432705334E+00,2.267226986762877E+00,1.475117383345656E+00,1.045384920629694E+00,7.877720818529990E-01,6.225337991148903E-01,5.115690451031906E-01,4.349170496402392E-01,3.814440322577813E-01,3.447469157598848E-01,3.212461668509977E-01,3.093477300493950E-01,3.093066145799037E-01,3.239231611662637E-01,3.611382880869027E-01,4.434151961031652E-01,6.551933196167097E-01,1.760781642897382E+00};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH);
      ldouble dataCHp[20] = {3.412027629288433E+01,8.362713710369944E+00,3.562836025564448E+00,1.905949700095471E+00,1.153999792587975E+00,7.546452247983108E-01,5.196543320898771E-01,3.710261843953806E-01,2.718480569531195E-01,2.028797771942403E-01,1.533262207101591E-01,1.167683370732885E-01,8.920562864910626E-02,6.804659798862765E-02,5.155782340790181E-02,3.854944696746540E-02,2.819062284294757E-02,1.990646766041888E-02,1.337081497620645E-02,9.327192714610702E-03};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCHp);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCHp);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCHp);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCHp);
      //std::cout<<"\n#  me=0.0005, MW=80.385, pspatial=81\n";
      pspatial = 81;
      ldouble dataCH2[20] = {1.271854233511375E+00,9.789084319451568E-01,7.792870021798937E-01,6.377102950526774E-01,5.341821016650772E-01,4.567070800939126E-01,3.977637413898730E-01,3.524762455967830E-01,3.176158628611309E-01,2.910364063907009E-01,2.713541599439645E-01,2.577817831120852E-01,2.500821963431110E-01,2.486569628406702E-01,2.548618122473640E-01,2.718452703713859E-01,3.069112891311704E-01,3.796331199703679E-01,5.625146809951643E-01,1.509779789972140E+00};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH2);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH2);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH2);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH2);
      ldouble dataCH2p[20] = {1.599592821626530E+00,1.195014823429346E+00,9.208459320003854E-01,7.271552010076066E-01,5.857317977643410E-01,4.796848663763820E-01,3.984170145747019E-01,3.350202068206851E-01,2.848482431116003E-01,2.446973890092003E-01,2.123189775914168E-01,1.861232994118847E-01,1.650018353981934E-01,1.482333938088856E-01,1.354703787938563E-01,1.268492900069245E-01,1.234104542240376E-01,1.286370355257731E-01,1.562768953823917E-01,3.322897885286755E-01};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH2p);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH2p);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH2p);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH2p);
      //std::cout<<"\n#  me=80, MW=80.385, pspatial=81\n";
      me=80;
      AnWeAmp.set_masses(me,MW);
      pspatial=250;
      ldouble dataCH3[20] = {3.408130582284264E+01,9.081208253437987E+00,4.159011044673917E+00,2.399490399808116E+00,1.578062360412619E+00,1.131373660919391E+00,8.636440097713755E-01,6.924325288572009E-01,5.782433003605887E-01,5.003717189869644E-01,4.472999381659096E-01,4.124659647734681E-01,3.922976482082296E-01,3.853570696864339E-01,3.922075030449231E-01,4.160908194554485E-01,4.652611951334960E-01,5.603933406397523E-01,7.640491927907039E-01,1.378228802216603E+00};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH3);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH3);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH3);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH3);
      ldouble dataCH3p[20] = {3.422529379644325E+01,8.707703860941635E+00,3.798846925070202E+00,2.083635265275877E+00,1.300964782671038E+00,8.850247580952758E-01,6.414533240186562E-01,4.893384357968512E-01,3.903497828046646E-01,3.246113575972830E-01,2.811799073190585E-01,2.538812978966983E-01,2.393903106851939E-01,2.364062595010777E-01,2.455508044855353E-01,2.701175597270303E-01,3.186923925452804E-01,4.139821678832165E-01,6.317318584791487E-01,1.419533150659489E+00};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH3p);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH3p);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH3p);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH3p);
      //std::cout<<"\n#  me=80, MW=80.385, pspatial=81\n";
      pspatial = 81;
      ldouble dataCH4[20] = {5.683234210351511E-01,5.553441689932734E-01,5.430753583285957E-01,5.314770539690670E-01,5.205125208538031E-01,5.101479497522606E-01,5.003522121492351E-01,4.910966409192917E-01,4.823548339436391E-01,4.741024781926566E-01,4.663171921174914E-01,4.589783844720268E-01,4.520671279285847E-01,4.455660460624429E-01,4.394592124662713E-01,4.337320609199060E-01,4.283713056869043E-01,4.233648711400227E-01,4.187018300357254E-01,4.143723498653714E-01};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH4);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH4);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH4);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH4);
      ldouble dataCH4p[20] = {8.032415699001588E-01,7.859171954477987E-01,7.696462672858951E-01,7.543726408862493E-01,7.400449474675302E-01,7.266162018425570E-01,7.140434537464926E-01,7.022874779227727E-01,6.913124988873259E-01,6.810859468497543E-01,6.715782417564544E-01,6.627626028468194E-01,6.546148814895875E-01,6.471134154006049E-01,6.402389026431010E-01,6.339742940834949E-01,6.283047032254913E-01,6.232173325780026E-01,6.187014159330438E-01,6.147481761427832E-01};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH4p);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH4p);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH4p);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH4p);
      //std::cout<<"\n#  me=80, MW=1, pspatial=81\n";
      MW=1;
      AnWeAmp.set_masses(me,MW);
      pspatial=250;
      ldouble dataCH5[20] = {5.732153756503400E+02,2.651113430109291E+02,2.259507085603407E+02,2.166663274881091E+02,2.173869450336996E+02,2.233655537691529E+02,2.330860227513259E+02,2.461258019932797E+02,2.625927409534980E+02,2.829714384635357E+02,3.081301267165755E+02,3.394386760085387E+02,3.790261076745048E+02,4.302850800306324E+02,4.988931108296292E+02,5.950574697808853E+02,7.390919849851654E+02,9.779635462873715E+02,1.450219905094143E+03,2.820928365332030E+03};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH5);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH5);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH5);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH5);
      ldouble dataCH5p[20] = {7.031141330677245E+02,4.038636474039782E+02,3.724357992542840E+02,3.714220347036365E+02,3.812601897069321E+02,3.974204322086327E+02,4.186186164541273E+02,4.447129733551058E+02,4.761675978701174E+02,5.139313797243310E+02,5.594908917017602E+02,6.150567503441340E+02,6.839249658880804E+02,7.711392243123749E+02,8.847400729030192E+02,1.038260829874827E+03,1.256078049857171E+03,1.585528857224976E+03,2.119303047698969E+03,2.752729774967880E+03};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH5p);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH5p);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH5p);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH5p);
      //std::cout<<"\n#  me=80, MW=1, pspatial=81\n";
      pspatial = 81;
      ldouble dataCH6[20] = {3.077954176660309E+03,1.098856104335329E+03,7.324259157335989E+02,5.822097206868574E+02,5.031997425056386E+02,4.566564632129541E+02,4.278817060874249E+02,4.101155023997267E+02,3.998549842532088E+02,3.951445095449307E+02,3.948490721501614E+02,3.983107192904380E+02,4.051728933283950E+02,4.152848368031354E+02,4.286452939142676E+02,4.453622568894837E+02,4.656083356772928E+02,4.895400866812836E+02,5.171108007572594E+02,5.475916642591168E+02};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH6);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH6);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH6);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2(); }, 0,0,MW,me,pspatial,dataCH6);
      ldouble dataCH6p[20] = {3.145474913101438E+03,1.168194424969924E+03,8.018089725805911E+02,6.510925764249997E+02,5.710954340926689E+02,5.230053507547566E+02,4.919933742803370E+02,4.711149271951255E+02,4.566163912066535E+02,4.462022118780772E+02,4.382724523219340E+02,4.315243426106810E+02,4.246915125946772E+02,4.163132096103714E+02,4.044568412223710E+02,3.863031254304595E+02,3.574405857552935E+02,3.105613726296397E+02,2.328740378741121E+02,1.005717255059575E+02};
      i += AnWeAmp.test_2to2_amp2([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH6p);
      i += AnWeAmp.test_2to2_amp2_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH6p);
      i += AnWeAmp.test_2to2_amp2_boosts([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH6p);
      i += AnWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AnWeAmp.amp2_Aplus(); }, 0,0,MW,me,pspatial,dataCH6p);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }



}
