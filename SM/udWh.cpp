
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

//File:  SPINAS/SM/udWh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/udWh.h"

namespace spinas {
  
  udWh::udWh(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), mu(massu), md(massd), mh(massh), MW(massW), WW(widthW), SW(sinW), propu(massu,0), propd(massd,0), propW(massW,widthW) {
    p1=particle(mu);
    p2=particle(md);
    p3=particle(MW);
    p4=particle(mh);
    //<13>,[23],<12>,[12],[343>
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s12s = sproduct(SQUARE,&p1,&p2,2);
    s343a = sproduct(SQUARE,&p3,&p4,&p3,2);
    //<13>,[23],[342>
    s342a = sproduct(SQUARE,&p3,&p4,&p2,2);
    //[23],<13>,[143>
    s143a = sproduct(SQUARE,&p1,&p4,&p3,2);
    //Couplings
    preW = e*e/(2.0*MW*MW*SW*SW);
    preu = mu*preW;
    pred = md*preW;
  }
  void udWh::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& massW){
    mu=massu;
    md=massd;
    mh=massh;
    MW=massW;
    p1.set_mass(mu);
    p2.set_mass(md);
    p3.set_mass(MW);
    p4.set_mass(mh);
    propu.set_mass(mu);
    propd.set_mass(md);
    propW.set_mass(MW);
    //Couplings
    preW = e*e/(2.0*MW*MW*SW*SW);
    preu = mu*preW;
    pred = md*preW;
  }
  void udWh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4); 
    //<13>,[23],<12>,[12],[343>
    a13a.update();
    s23s.update();
    a12a.update();
    s12s.update();
    s343a.update();
    //<13>,[23],[342>
    s342a.update();
    //[23],<13>,[143>
    s143a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propW.denominator(propSP);
    pDenU=propu.denominator(propUP);
    pDenT=propd.denominator(propTP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble udWh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;
    constexpr ldouble two=2;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      //S-Channel W
      //preW = e*e/(2.0*MW*MW*SW*SW);
      //all ingoing:
      // - preW (2MW^2 <13>[23]+(Md<12>-Mu[12])[343>))/(s-MW^2)
      //34 outgoing:
      // + preW (2MW^2 <13>[23]-(Md<12>-Mu[12])[343>))/(s-MW^2)
      amplitude += normFactor*preW*(two*MW*MW*a13a.v(ds1,ds3a)*s23s.v(ds2,ds3b)-(md*a12a.v(ds1,ds2)-mu*s12s.v(ds1,ds2))*s343a.v(ds3a,ds3b))/pDenS;
      
      //T-Channel d
      //pred = e*e*md/(2.0*MW*MW*SW*SW);
      //all ingoing:
      // - pred <13>(2Md[23]+[342>)/(t-Md^2)
      //34 outgoing:
      // + pred <13>(2Md[23]-[342>)/(t-Md^2)
      amplitude += normFactor*pred*a13a.v(ds1,ds3a)*(two*md*s23s.v(ds2,ds3b)-s342a.v(ds3b,ds2))/pDenT;
      
      //U-Channel u
      //preu = e*e*mu/(2.0*MW*MW*SW*SW);
      //all ingoing:
      // - preu [23](2Mu<13>+[143>)/(u-Mu^2)
      //34 outgoing:
      // + preu [23](2Mu<13>-[143>)/(u-Mu^2)
      amplitude += normFactor*preu*s23s.v(ds2,ds3a)*(two*mu*a13a.v(ds1,ds3b)-s143a.v(ds1,ds3b))/pDenU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble udWh::amp2(){
    constexpr ldouble three=3;
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2){
	  M = amp(j1,j2,j3);
	  //Color Factor 3
	  amp2 += three*std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2^2=1/4
    //Average over colors 1/3^2=1/9
    return amp2/36.0;
  }

  
  

  //  Tests
  int test_udWh(){
    int n=0;//Number of fails
    std::cout<<"\t* u , D  -> W+, h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, md=0.0075, mh=125, MW=80.385, pspatial=2500\n";
      ldouble mu=0.0042, md=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      udWh udWhAmp = udWh(EE,mu,md,mh,MW,0,SW);
      ldouble pspatial=2500;
      ldouble dataCH[20] = {1.977959387023590E-04,5.553578438534919E-04,8.731912850838541E-04,1.151295629087078E-03,1.389670805509331E-03,1.588316795178934E-03,1.747233590740464E-03,1.866421188783604E-03,1.945879587531015E-03,1.985608785990860E-03,1.985608783600447E-03,1.945879580067191E-03,1.866421175306207E-03,1.747233569444622E-03,1.588316762906632E-03,1.389670756662903E-03,1.151295552945240E-03,8.731911577081375E-04,5.553575929692955E-04,1.977950572976166E-04};
      i += udWhAmp.test_2to2_amp2([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH);
      i += udWhAmp.test_2to2_amp2_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH);
      //std::cout<<"\n#mu=0.0042, md=0.0075, mh=125, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH1[20] = {6.018979424438169E-04,9.130049301957049E-04,1.189544990773324E-03,1.431517620192453E-03,1.638922754371022E-03,1.811760375180865E-03,1.950030475533859E-03,2.053733052105178E-03,2.122868103148396E-03,2.157435627683442E-03,2.157435625151872E-03,2.122868095261964E-03,2.053733037927507E-03,1.950030453271502E-03,1.811760341707863E-03,1.638922704179838E-03,1.431517542834832E-03,1.189544863299574E-03,9.130046856244819E-04,6.018971678226306E-04};
      i += udWhAmp.test_2to2_amp2([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH1);
      i += udWhAmp.test_2to2_amp2_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH1);
      i += udWhAmp.test_2to2_amp2_boosts([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH1);
      i += udWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH1);
      //std::cout<<"\n# mu=0.0042, md=0.0075, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {2.118445532675537E-03,2.259021558942225E-03,2.383978003613005E-03,2.493314891808605E-03,2.587032225846304E-03,2.665130005225040E-03,2.727608229207615E-03,2.774466897198961E-03,2.805706008766915E-03,2.821325563609254E-03,2.821325561522218E-03,2.805706002378829E-03,2.774466886116577E-03,2.727608212733013E-03,2.665129982287821E-03,2.587032194907675E-03,2.493314850775818E-03,2.383977950010606E-03,2.259021491853656E-03,2.118445468400985E-03};
      i += udWhAmp.test_2to2_amp2([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH2);
      i += udWhAmp.test_2to2_amp2_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH2);
      i += udWhAmp.test_2to2_amp2_boosts([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH2);
      i += udWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, md=130, mh=125, MW=80.385, pspatial=1\n";
      mu = 125.1;
      md = 130;
      mh = 125;
      pspatial = 1;
      udWhAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH4[20] = {6.195948655716483E-02,6.196195317128930E-02,6.196447158208101E-02,6.196704179973250E-02,6.196966383450624E-02,6.197233769673457E-02,6.197506339681978E-02,6.197784094523418E-02,6.198067035252006E-02,6.198355162928985E-02,6.198648478622602E-02,6.198946983408119E-02,6.199250678367815E-02,6.199559564590989E-02,6.199873643173958E-02,6.200192915220065E-02,6.200517381839678E-02,6.200847044150193E-02,6.201181903276031E-02,6.201521960348649E-02};
      i += udWhAmp.test_2to2_amp2([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH4);
      i += udWhAmp.test_2to2_amp2_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH4);
      i += udWhAmp.test_2to2_amp2_boosts([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH4);
      i += udWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, md=130, mh=0.0005, MW=80.385, pspatial=250\n";
      mu = 125;
      mh = 0.0005;
      pspatial = 250;
      udWhAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH3[20] = {4.651311599840420E-01,3.012586631856314E-01,2.246147063029743E-01,1.810226298549274E-01,1.535863141178511E-01,1.353599869331981E-01,1.229946196170931E-01,1.147118826093797E-01,1.095221864470204E-01,1.068848757117513E-01,1.065536734777962E-01,1.085146335553705E-01,1.129847344206419E-01,1.204709486088184E-01,1.319214619246607E-01,1.490624928061836E-01,1.751829577853955E-01,2.171997082104798E-01,2.922873962234312E-01,4.575980225757061E-01};
      i += udWhAmp.test_2to2_amp2([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH3);
      i += udWhAmp.test_2to2_amp2_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH3);
      i += udWhAmp.test_2to2_amp2_boosts([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH3);
      i += udWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH3);
      //std::cout<<"\n# mu=125, md=130, mh=0.0005, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH9[20] = {1.026577502806115E-01,1.026570682392098E-01,1.026565235239810E-01,1.026561161417563E-01,1.026558461001255E-01,1.026557134074379E-01,1.026557180728019E-01,1.026558601060855E-01,1.026561395179166E-01,1.026565563196827E-01,1.026571105235316E-01,1.026578021423716E-01,1.026586311898714E-01,1.026595976804607E-01,1.026607016293300E-01,1.026619430524317E-01,1.026633219664792E-01,1.026648383889484E-01,1.026664923380769E-01,1.026682838328651E-01};
      i += udWhAmp.test_2to2_amp2([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH9);
      i += udWhAmp.test_2to2_amp2_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH9);
      i += udWhAmp.test_2to2_amp2_boosts([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH9);
      i += udWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH9);
      //std::cout<<"\n# mu=0.004, md=130, mh=0.0005, MW=0.0004, pspatial=250\n";
      mu = 0.004;
      mh = 0.0005;
      MW = 0.0004;
      pspatial = 250;
      udWhAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH5[20] = {2.622829819215927E+20,1.701745938260428E+20,1.228050792994900E+20,9.408619991385385E+19,7.484112817603979E+19,6.105322284997108E+19,5.069222533451361E+19,4.262287077730448E+19,3.616109765671067E+19,3.087029882020166E+19,2.645876859272316E+19,2.272422264068348E+19,1.952198515658280E+19,1.674585061777926E+19,1.431609750058111E+19,1.217171948433830E+19,1.026523916829466E+19,8.559155528433439E+18,7.023454777742112E+18,5.633831331831932E+18};
      i += udWhAmp.test_2to2_amp2([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH5);
      i += udWhAmp.test_2to2_amp2_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH5);
      i += udWhAmp.test_2to2_amp2_boosts([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH5);
      i += udWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH5);
      //std::cout<<"\n# mu=0.004, md=130, mh=0.0005, MW=0.0004, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH6[20] = {3.696679483958449E+17,4.208182030176835E+17,4.717831927484247E+17,5.225635885956637E+17,5.731600575676567E+17,6.235732622896219E+17,6.738038604432486E+17,7.238525039304539E+17,7.737198375935857E+17,8.234064971964276E+17,8.729131061213902E+17,9.222402697257053E+17,9.713885651717865E+17,1.020358521857921E+18,1.069150580486342E+18,1.117764997538638E+18,1.166201585975759E+18,1.214458833606051E+18,1.262529495057358E+18,1.310343295539311E+18};
      i += udWhAmp.test_2to2_amp2([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH6);
      i += udWhAmp.test_2to2_amp2_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH6);
      i += udWhAmp.test_2to2_amp2_boosts([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH6);
      i += udWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH6);
      //std::cout<<"\n# mu=0.004, md=0.007, mh=0.0005, MW=0.0004, pspatial=250\n";
      mu = 0.004;
      md = 0.007;
      mh = 0.0005;
      MW = 0.0004;
      pspatial = 2.5;
      udWhAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH7[20] = {9.940785282944042E+03,3.314793763887899E+03,1.993693992802162E+03,1.431021309530957E+03,1.121685056882137E+03,9.280486686439108E+02,7.973113765156540E+02,7.049987601831344E+02,6.383639819953642E+02,5.903049855402769E+02,5.568087592414322E+02,5.358142547470852E+02,5.267727014493285E+02,5.306948720511245E+02,5.507956751236761E+02,5.943640022612253E+02,6.781088032206662E+02,8.465631364511974E+02,1.264339451756419E+03,3.415873105222789E+03};
      i += udWhAmp.test_2to2_amp2([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH7);
      i += udWhAmp.test_2to2_amp2_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH7);
      i += udWhAmp.test_2to2_amp2_boosts([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH7);
      i += udWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH7);
      //std::cout<<"\n# mu=0.004, md=0.007, mh=0.0005, MW=0.0004, pspatial=1\n";
      pspatial = 0.001;
      ldouble dataCH8[20] = {5.446150612272197E+02,5.429214992722684E+02,5.413843233974013E+02,5.400092487253778E+02,5.388019389033711E+02,5.377680607758939E+02,5.369133336165869E+02,5.362435732464145E+02,5.357647311352385E+02,5.354829283252889E+02,5.354044837122315E+02,5.355359358514395E+02,5.358840569962923E+02,5.364558574854678E+02,5.372585778279731E+02,5.382996648207916E+02,5.395867266814922E+02,5.411274603573125E+02,5.429295417007579E+02,5.450004658208965E+02};
      i += udWhAmp.test_2to2_amp2([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH8);
      i += udWhAmp.test_2to2_amp2_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH8);
      i += udWhAmp.test_2to2_amp2_boosts([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH8);
      i += udWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udWhAmp.amp2(); }, mu,md,MW,mh,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}
