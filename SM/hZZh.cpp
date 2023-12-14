
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

//File:  SPINAS/SM/hZZh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/hZZh.h"

namespace spinas {

  hZZh::hZZh(const ldouble& echarge, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);
    proph = propagator(mh,wh);  
    p1=particle(mh);
    p2=particle(MZ);
    p3=particle(MZ);
    p4=particle(mh);
    //<23>,[23]
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    //[312>,[213>
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
    //Couplings
    pre = 2.0*e*e/(4.0*MW*MW*SW*SW);
    preS = 3.0*pre*mh*mh;
  }
  void hZZh::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(mh);
    p2.set_mass(MZ);
    p3.set_mass(MZ);
    p4.set_mass(mh);
    propZ.set_mass(MZ);
    proph.set_mass(mh);
    //Couplings
    pre = 2.0*e*e/(4.0*MW*MW*SW*SW);
    preS = 3.0*pre*mh*mh;
  }
  void hZZh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<23>,[23]
    s23s.update();
    a23a.update();
    //[312>,[213>
    s312a.update();
    s213a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propZ.den(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=propZ.den(propTP);
    pDenU=proph.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble hZZh::amp(const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1;
    int ds3a, ds3b, ds2a, ds2b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds2);
    ldouble normFactor=get_spin_normalization(ds3,ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds2,ds2a,ds2b, i);

      //pre = sqrt2*e*e/(4.0*MW*MW*SW*SW);
      
      //4-Point
      //hhZZ all in:  - [34]<34>
      //hZZh: 2<->4:  - [23]<23>
      //34 out:       + [23]<23>
      amplitude += + normFactor*pre*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b);
          
      //U-Channel h
      //preS = 3.0*pre*mh*mh;
      //hhZZ all ingoing: 
      //+ preS [34] <34> /(s-Mh^2)
      //hZZh: 2<->4:
      //+ preS [23] <23> /(u-Mh^2)
      //34 out:
      //- preS [23] <23> /(u-Mh^2)
      amplitude += - normFactor*preS*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)/pDenU;
      
      //T-Channel Z
      //hhZZ all ingoing:
      //+ pre ( 2MZ^2[34]<34> + MZ([34][314>+<34>[413>) + [314>[413> )/(t-MZ^2)
      //hZZh: 2<->4:
      //+ pre ( 2MZ^2[23]<23> - MZ([23][312>+<23>[213>) + [312>[213> )/(t-MZ^2)
      //34 out:
      //- pre ( 2MZ^2[23]<23> + MZ([23][312>+<23>[213>) + [312>[213> )/(t-MZ^2)
      amplitude +=
	- normFactor*pre*2.0*MZ*MZ*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)/pDenT
      	- normFactor*pre*MZ*(s23s.v(ds2a,ds3a)*s312a.v(ds3b,ds2b)+a23a.v(ds2a,ds3a)*s213a.v(ds2b,ds3b))/pDenT
      	- normFactor*pre*s312a.v(ds3a,ds2a)*s213a.v(ds2b,ds3b)/pDenT;
      
      //S-Channel Z
      //hhZZ all ingoing:
      //+ pre ( 2MZ^2[34]<34> - MZ([34][413>+<34>[314>) + [314>[413> )/(u-MZ^2)
      //hZZh: 2<->4:
      //+ pre ( 2MZ^2[23]<23> + MZ([23][213>+<23>[312>) + [312>[213> )/(s-MZ^2)
      //34 out:
      //- pre ( 2MZ^2[23]<23> + MZ([23][213>+<23>[312>) + [312>[213> )/(s-MZ^2)
      amplitude +=
	- normFactor*pre*2.0*MZ*MZ*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)/pDenS
      	- normFactor*pre*MZ*(s23s.v(ds2a,ds3a)*s213a.v(ds2b,ds3b)+a23a.v(ds2a,ds3a)*s312a.v(ds3b,ds2b))/pDenS
      	- normFactor*pre*s312a.v(ds3a,ds2a)*s213a.v(ds2b,ds3b)/pDenS;
      


    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble hZZh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-2;j3<=2;j3+=2){
	M = amp(j2,j3);
	amp2 += std::pow(std::abs(M),2);
      }
    //Average over initial spins 1/3
    return amp2/3.0;
  }

   


  //  Tests
  int test_hZZh(){
    int n=0;//Number of fails
    std::cout<<"\t* h , Z  -> Z , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      hZZh hZZhAmp = hZZh(EE,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {9.702185971437629E+00,2.840658500540382E+00,1.373923916406302E+00,8.333585443605196E-01,5.768249306874748E-01,4.357360328181883E-01,3.503079955259413E-01,2.950252995444098E-01,2.575298775563592E-01,2.312662175897053E-01,2.125214126555473E-01,1.991042923571574E-01,1.897170400885594E-01,1.836612238163996E-01,1.807499259542944E-01,1.814339480316751E-01,1.873559787244303E-01,2.032357546030402E-01,2.444882884831904E-01,3.829641650374800E-01};
      i += hZZhAmp.test_2to2_amp2([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH);
      i += hZZhAmp.test_2to2_amp2_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH);
      i += hZZhAmp.test_2to2_amp2_boosts([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH);
      i += hZZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {2.451115627934268E+00,1.433589042249805E+00,9.638141295310676E-01,7.100449054207139E-01,5.585044119311606E-01,4.616850313622220E-01,3.969413485620730E-01,3.524138950080675E-01,3.214601066256358E-01,3.001940140275007E-01,2.863095085912475E-01,2.784910744517876E-01,2.761237403235718E-01,2.791809636236844E-01,2.882509655072914E-01,3.047203789733652E-01,3.312088019106163E-01,3.724986391118508E-01,4.375867698806776E-01,5.446079750152306E-01};
      i += hZZhAmp.test_2to2_amp2([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH2);
      i += hZZhAmp.test_2to2_amp2_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH2);
      i += hZZhAmp.test_2to2_amp2_boosts([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH2);
      i += hZZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {1.196558899334687E+00,1.196543735850995E+00,1.196528575062413E+00,1.196513416968725E+00,1.196498261569714E+00,1.196483108865165E+00,1.196467958854862E+00,1.196452811538589E+00,1.196437666916130E+00,1.196422524987269E+00,1.196407385751791E+00,1.196392249209480E+00,1.196377115360120E+00,1.196361984203496E+00,1.196346855739391E+00,1.196331729967591E+00,1.196316606887879E+00,1.196301486500041E+00,1.196286368803860E+00,1.196271253799122E+00};
      i += hZZhAmp.test_2to2_amp2([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH4);
      i += hZZhAmp.test_2to2_amp2_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH4);
      i += hZZhAmp.test_2to2_amp2_boosts([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH4);
      i += hZZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH4);
      //std::cout<<"\n# mh=5, MW=80.385, pspatial=95\n";
      mh = 5;
      pspatial = 95;
      hZZhAmp.set_masses(mh,MW);
      ldouble dataCH6[20] = {5.874597034818223E-01,3.396984083902609E-01,2.150231703654108E-01,1.472853852834462E-01,1.087089599241926E-01,8.619434708257544E-02,7.301467433680130E-02,6.548407112691888E-02,6.147766923040730E-02,5.971965225964169E-02,5.941804099143193E-02,6.006703487737657E-02,6.133484887774605E-02,6.299687742552730E-02,6.489287910491245E-02,6.689461591681263E-02,6.886907150260431E-02,7.059698819278432E-02,7.141610622898256E-02,6.595135939862370E-02};
      i += hZZhAmp.test_2to2_amp2([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH6);
      i += hZZhAmp.test_2to2_amp2_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH6);
      i += hZZhAmp.test_2to2_amp2_boosts([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH6);
      i += hZZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH6);
      //std::cout<<"\n# mh=5, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH7[20] = {7.692661729495359E+00,2.005554963502649E+00,8.715101754992688E-01,4.805133132361269E-01,3.069529984138050E-01,2.177095882915593E-01,1.671972798044100E-01,1.366255516256181E-01,1.172027290871176E-01,1.044173012710976E-01,9.578530504493739E-02,8.986115619222992E-02,8.576496999674822E-02,8.294076496672272E-02,8.102545050754731E-02,7.977380654777712E-02,7.901185763491177E-02,7.859968371433251E-02,7.836131034872335E-02,7.738454117428552E-02};
      i += hZZhAmp.test_2to2_amp2([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH7);
      i += hZZhAmp.test_2to2_amp2_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH7);
      i += hZZhAmp.test_2to2_amp2_boosts([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH7);
      i += hZZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH7);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=95\n";
      mh = 0.0005;
      pspatial = 95;
      hZZhAmp.set_masses(mh,MW);
      ldouble dataCH3[20] = {5.869345830057631E-01,3.389960907446721E-01,2.143898425258087E-01,1.467801381803489E-01,1.083392526594626E-01,8.595376674802124E-02,7.289462898271902E-02,6.547761298971541E-02,6.158082389322628E-02,5.993204711134646E-02,5.974353724178325E-02,6.051476852667850E-02,6.192112786884610E-02,6.374873284270761E-02,6.585503204317286E-02,6.814433053919124E-02,7.055222235970626E-02,7.303548216742455E-02,7.556538246408513E-02,7.812320541453786E-02};
      i += hZZhAmp.test_2to2_amp2([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH3);
      i += hZZhAmp.test_2to2_amp2_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH3);
      i += hZZhAmp.test_2to2_amp2_boosts([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH3);
      i += hZZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH3);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH5[20] = {7.689815629042528E+00,2.004185281730475E+00,8.707173727556421E-01,4.799951489967848E-01,3.065902480150300E-01,2.174455615166891E-01,1.670014272179612E-01,1.364800024900087E-01,1.170965137211280E-01,1.043436068784935E-01,9.574010365676766E-02,8.984254885123399E-02,8.577299415219168E-02,8.297762954711485E-02,8.109638122533321E-02,7.988912442215551E-02,7.919226194373406E-02,7.889229786224733E-02,7.890927869926707E-02,7.918614977758390E-02};
      i += hZZhAmp.test_2to2_amp2([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH5);
      i += hZZhAmp.test_2to2_amp2_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH5);
      i += hZZhAmp.test_2to2_amp2_boosts([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH5);
      i += hZZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH5);
      //std::cout<<"\n# mh=0, MW=80.385, pspatial=95\n";
      mh = 0;
      pspatial = 95;
      hZZhAmp.set_masses(mh,MW);
      ldouble dataCH8[20] = {5.869345830005516E-01,3.389960907376668E-01,2.143898425194881E-01,1.467801381753089E-01,1.083392526557797E-01,8.595376674563213E-02,7.289462898153841E-02,6.547761298967426E-02,6.158082389428558E-02,5.993204711350338E-02,5.974353724507769E-02,6.051476853120384E-02,6.192112787476856E-02,6.374873285030259E-02,6.585503205289650E-02,6.814433055183375E-02,7.055222237676574E-02,7.303548219223119E-02,7.556538250667626E-02,7.812320554528798E-02};
      i += hZZhAmp.test_2to2_amp2([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH8);
      i += hZZhAmp.test_2to2_amp2_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH8);
      i += hZZhAmp.test_2to2_amp2_boosts([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH8);
      i += hZZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH8);
      //std::cout<<"\n# mh=0, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH9[20] = {7.689815629014073E+00,2.004185281716780E+00,8.707173727477151E-01,4.799951489916052E-01,3.065902480114051E-01,2.174455615140514E-01,1.670014272160054E-01,1.364800024885565E-01,1.170965137200690E-01,1.043436068777596E-01,9.574010365631952E-02,8.984254885105109E-02,8.577299415227572E-02,8.297762954748787E-02,8.109638122604831E-02,7.988912442331354E-02,7.919226194554621E-02,7.889229786519014E-02,7.890927870478022E-02,7.918614979582940E-02};
      i += hZZhAmp.test_2to2_amp2([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH9);
      i += hZZhAmp.test_2to2_amp2_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH9);
      i += hZZhAmp.test_2to2_amp2_boosts([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH9);
      i += hZZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hZZhAmp.amp2(); }, mh,MZ,MZ,mh,pspatial,dataCH9);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  

}
