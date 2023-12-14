
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

//File:  SPINAS/SM/AAuu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AAuu.h"

namespace spinas {

  AAuu::AAuu(const ldouble& echarge, const ldouble& massu):
    e(echarge), Qu(2.0/3.0), mu(massu), prop(massu,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(0);
    p2=particle(0);
    p3=particle(mu);
    p4=particle(mu);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s132a = sproduct(SQUARE,&p1,&p3,&p2);
    s231a = sproduct(SQUARE,&p2,&p3,&p1);
  }
  void AAuu::set_masses(const ldouble& massu){
    mu=massu;
    p3.set_mass(mu);
    p4.set_mass(mu);
    prop.set_mass(mu);
  }
  void AAuu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s34s.update();
    a34a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    s24s.update();
    a24a.update();
    s23s.update();
    a23a.update();
    s14s.update();
    a14a.update();
    s132a.update();
    s231a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT = prop.den(propTP);
    pDenU = prop.den(propUP);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AAuu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds1>0&&ds2>0){
      //eeAA:   mu[34]^2<12>
      //AAuu:   mu[12]^2<34>
      //34 out: mu[12]^2<34>
      return 2.0*e*e*Qu*Qu*mu*s12s.v()*s12s.v()*a34a.v(ds3,ds4)/pDenT/pDenU;
    }
    else if(ds1<0&&ds2<0){
      //eeAA:   mu<34>^2[12]
      //AAuu:   mu<12>^2[34]
      //34 out: mu<12>^2[34]
      return 2.0*e*e*Qu*Qu*mu*a12a.v()*a12a.v()*s34s.v(ds3,ds4)/pDenT/pDenU;
    }
    else if(ds1>0&&ds2<0){
      //eeAA:   ([13]<24>+[23]<14>)[314>
      //AAuu:   ([31]<42>+[41]<32>)[132>
      //34 out: -([13]<24>+[14]<23>)[132>
      return -2.0*e*e*Qu*Qu*(s13s.v(ds3)*a24a.v(ds4)+s14s.v(ds4)*a23a.v(ds3))*s132a.v()/pDenT/pDenU;
    }
    else if(ds1<0&&ds2>0){
      //eeAA:   (<13>[24]+<23>[14])[413>
      //AAuu:   (<31>[42]+<41>[32])[231>
      //34 out: -(<13>[24]+<14>[23])[231>
      return -2.0*e*e*Qu*Qu*(a13a.v(ds3)*s24s.v(ds4)+a14a.v(ds4)*s23s.v(ds3))*s231a.v()/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble AAuu::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2*1/2=1/4
    return amp2/4.0;
  }

  //A+, A+ -> e, E
  ldouble AAuu::amp2_Aplus_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(2,2,j3,j4);
	amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
      }
    return amp2;
  }

  //A+, A- -> u, U
  ldouble AAuu::amp2_Aplus_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(2,-2,j3,j4);
	amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
      }
    return amp2;
  }



  //  Tests
  int test_AAuu(){
    int n=0;//Number of fails
    std::cout<<"\t* A , A  -> u , U       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n#  mu=0.0042, pspatial=250\n";
      ldouble mu=0.0042;
      ldouble EE=0.31333;
      AAuu AAuuAmp = AAuu(EE,mu);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.458042180482265E-01,1.418144116468877E-01,8.159547843381763E-02,5.627615422858103E-02,4.266361433001905E-02,3.444915059702072E-02,2.922559417150102E-02,2.589296517184332E-02,2.389850178067821E-02,2.296125394646710E-02,2.296125394646710E-02,2.389850178067820E-02,2.589296517184332E-02,2.922559417150101E-02,3.444915059702072E-02,4.266361433001905E-02,5.627615422858103E-02,8.159547843381765E-02,1.418144116468878E-01,4.458042180482273E-01};
      i += AAuuAmp.test_2to2_amp2([&]() { return AAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      i += AAuuAmp.test_2to2_amp2_rotations([&]() { return AAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      i += AAuuAmp.test_2to2_amp2_boosts([&]() { return AAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      i += AAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      //std::cout<<"\n#  A+ , A+ -> e , E\n";
      ldouble dataCHpp[20] = {2.713272703418789E-09,3.349474717400204E-10,1.347555206277213E-10,7.733902743356429E-11,5.301689032060475E-11,4.055471091172542E-11,3.349719552216433E-11,2.934675663032938E-11,2.699412198066772E-11,2.592250121884052E-11,2.592250085292229E-11,2.699411693912762E-11,2.934675683785245E-11,3.349719391026063E-11,4.055470769045912E-11,5.301688966076609E-11,7.733902582335466E-11,1.347555195922235E-10,3.349474720396583E-10,2.713272709022865E-09};
      i += AAuuAmp.test_2to2_amp2([&]() { return AAuuAmp.amp2_Aplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      i += AAuuAmp.test_2to2_amp2_rotations([&]() { return AAuuAmp.amp2_Aplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      i += AAuuAmp.test_2to2_amp2_boosts([&]() { return AAuuAmp.amp2_Aplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      i += AAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAuuAmp.amp2_Aplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      //std::cout<<"\n#  A+ , A- -> e , E\n";
      ldouble dataCHpm[20] = {8.916084333831802E-01,2.836288229588279E-01,1.631909567328797E-01,1.125523083798230E-01,8.532722860702120E-02,6.889830115348672E-02,5.845118830950484E-02,5.178593031433990E-02,4.779700353436230E-02,4.592250786701170E-02,4.592250786701170E-02,4.779700353436229E-02,5.178593031433989E-02,5.845118830950483E-02,6.889830115348673E-02,8.532722860702120E-02,1.125523083798230E-01,1.631909567328798E-01,2.836288229588281E-01,8.916084333831819E-01};
      i += AAuuAmp.test_2to2_amp2([&]() { return AAuuAmp.amp2_Aplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      i += AAuuAmp.test_2to2_amp2_rotations([&]() { return AAuuAmp.amp2_Aplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      i += AAuuAmp.test_2to2_amp2_boosts([&]() { return AAuuAmp.amp2_Aplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      i += AAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAuuAmp.amp2_Aplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      //Close to threshold
      //std::cout<<"\n#  mu=0.0042, pspatial=0.006\n";
      pspatial = 0.006;
      ldouble dataCH2[20] = {6.563581875682613E-02,5.745259024673335E-02,5.105176244257331E-02,4.611092464392209E-02,4.231743097242840E-02,3.943406038640045E-02,3.729158502344099E-02,3.577281534710020E-02,3.479967046295057E-02,3.432441893428717E-02,3.432441893428717E-02,3.479967046295056E-02,3.577281534710019E-02,3.729158502344099E-02,3.943406038640044E-02,4.231743097242839E-02,4.611092464392207E-02,5.105176244257328E-02,5.745259024673331E-02,6.563581875682606E-02};
      i += AAuuAmp.test_2to2_amp2([&]() { return AAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      i += AAuuAmp.test_2to2_amp2_rotations([&]() { return AAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      i += AAuuAmp.test_2to2_amp2_boosts([&]() { return AAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      i += AAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      //std::cout<<"\n#  A+ , A+ -> e , E\n";
      ldouble dataCH2pp[20] = {1.160598694233970E-01,8.477075527443655E-02,6.648073680648145E-02,5.493050257336342E-02,4.726817232757115E-02,4.204444191586187E-02,3.846461879935979E-02,3.607150622774553E-02,3.459806683137773E-02,3.389497400286161E-02,3.389497400286161E-02,3.459806683137773E-02,3.607150622774552E-02,3.846461879935979E-02,4.204444191586186E-02,4.726817232757115E-02,5.493050257336340E-02,6.648073680648142E-02,8.477075527443652E-02,1.160598694233969E-01};
      i += AAuuAmp.test_2to2_amp2([&]() { return AAuuAmp.amp2_Aplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      i += AAuuAmp.test_2to2_amp2_rotations([&]() { return AAuuAmp.amp2_Aplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      i += AAuuAmp.test_2to2_amp2_boosts([&]() { return AAuuAmp.amp2_Aplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      i += AAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAuuAmp.amp2_Aplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      //std::cout<<"\n#  A+ , A- -> e , E\n";
      ldouble dataCH2pm[20] = {1.521176809025528E-02,3.013442521903016E-02,3.562278807866518E-02,3.729134671448077E-02,3.736668961728565E-02,3.682367885693903E-02,3.611855124752220E-02,3.547412446645487E-02,3.500127409452342E-02,3.475386386571273E-02,3.475386386571273E-02,3.500127409452340E-02,3.547412446645486E-02,3.611855124752218E-02,3.682367885693902E-02,3.736668961728564E-02,3.729134671448075E-02,3.562278807866515E-02,3.013442521903011E-02,1.521176809025521E-02};
      i += AAuuAmp.test_2to2_amp2([&]() { return AAuuAmp.amp2_Aplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      i += AAuuAmp.test_2to2_amp2_rotations([&]() { return AAuuAmp.amp2_Aplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      i += AAuuAmp.test_2to2_amp2_boosts([&]() { return AAuuAmp.amp2_Aplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      i += AAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAuuAmp.amp2_Aplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
