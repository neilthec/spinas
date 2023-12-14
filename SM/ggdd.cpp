
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

//File:  SPINAS/SM/ggdd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ggdd.h"

namespace spinas {

  ggdd::ggdd(const ldouble& gcharge, const ldouble& massd):
    gs(gcharge), md(massd), propd(massd,0), propg(0,0){
    constexpr ldouble two=2;
    sqrt2 = std::sqrt(two);
    p1=particle(0);
    p2=particle(0);
    p3=particle(md);
    p4=particle(md);
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
  void ggdd::set_masses(const ldouble& massd){
    md=massd;
    p3.set_mass(md);
    p4.set_mass(md);
    propd.set_mass(md);
  }
  void ggdd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS = std::real(propg.den(propSP));
    pDenT = std::real(propd.den(propTP));
    pDenU = std::real(propd.den(propUP));
  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ggdd::amp_Num(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    ////////////////////////////////////////////////////
    //Missing color factor includig averaging over color
    //Missing symmetry factor when identical
    ////////////////////////////////////////////////////
    if(ds1>0&&ds2>0){
      //eeAA:   me[34]^2<12>
      //AAee:   me[12]^2<34>
      //34 out: me[12]^2<34>
      return 2.0*gs*gs*md*s12s.v()*s12s.v()*a34a.v(ds3,ds4);
    }
    else if(ds1<0&&ds2<0){
      //eeAA:   me<34>^2[12]
      //AAee:   me<12>^2[34]
      //34 out: me<12>^2[34]
      return 2.0*gs*gs*md*a12a.v()*a12a.v()*s34s.v(ds3,ds4);
    }
    else if(ds1>0&&ds2<0){
      //eeAA:   ([13]<24>+[23]<14>)[314>
      //AAee:   ([31]<42>+[41]<32>)[132>
      //34 out: -([13]<24>+[14]<23>)[132>
      return -2.0*gs*gs*(s13s.v(ds3)*a24a.v(ds4)+s14s.v(ds4)*a23a.v(ds3))*s132a.v();
    }
    else if(ds1<0&&ds2>0){
      //eeAA:   (<13>[24]+<23>[14])[413>
      //AAee:   (<31>[42]+<41>[32])[231>
      //34 out: -(<13>[24]+<14>[23])[231>
      return -2.0*gs*gs*(a13a.v(ds3)*s24s.v(ds4)+a14a.v(ds4)*s23s.v(ds3))*s231a.v();
    }
    return cdouble(0,0);    
  }

  //set_momenta(...) must be called before amp2_Feynman().
  //This expression comes from Ellis, Stirling and Weber: QCD and Collider Physics, Table 7.1
  //This agrees with CalcHEP.
  //This assumes massless quarks (up quark here)
  ldouble ggdd::amp2_Feynman(){
    ldouble pre = gs*gs*gs*gs;
    constexpr ldouble oneSixth=1.0/6.0, threeEights = 3.0/8.0;
    cdouble M2 =  pre*(pDenT*pDenT+pDenU*pDenU)*(oneSixth/pDenT/pDenU - threeEights/pDenS/pDenS);
    return std::abs(M2);
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble ggdd::amp2(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M_Num = amp_Num(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenS*pDenT,2)+6.0/std::pow(pDenS*pDenU,2)-2.0/3.0/std::pow(pDenT*pDenU,2));
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/8*1/8=1/64
    return amp2/256.0;
  }

  //g+, g+ -> e, E
  ldouble ggdd::amp2_gplus_gplus(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M_Num = amp_Num(2,2,j3,j4);
	amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenS*pDenT,2)+6.0/std::pow(pDenS*pDenU,2)-2.0/3.0/std::pow(pDenT*pDenU,2));
      }
    //Average over initial colors 1/8*1/8=1/64
    return amp2/64.0;
  }

  //g+, g- -> e, E
  ldouble ggdd::amp2_gplus_gminus(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M_Num = amp_Num(2,-2,j3,j4);
	amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenS*pDenT,2)+6.0/std::pow(pDenS*pDenU,2)-2.0/3.0/std::pow(pDenT*pDenU,2));
      }
    //Average over initial colors 1/8*1/8=1/64
    return amp2/64.0;
  }
  



  //  Tests
  int test_ggdd(){
    int n=0;//Number of fails
    std::cout<<"\t* g , g  -> d , D       :";
    {//amp^2
      int i=0;
      // md=0.0075, pspatial=250
      ldouble md=0.0075;
      ldouble gg=1.238;
      ggdd ggddAmp = ggdd(gg,md);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.444059143282378E+01,4.101585503869656E+00,2.108242849756263E+00,1.302164654568734E+00,8.884905508878588E-01,6.510096138660674E-01,5.072237799031505E-01,4.194347338316461E-01,3.686983208010841E-01,3.453858705319250E-01,3.453858705319250E-01,3.686983208010840E-01,4.194347338316460E-01,5.072237799031504E-01,6.510096138660673E-01,8.884905508878587E-01,1.302164654568734E+00,2.108242849756262E+00,4.101585503869654E+00,1.444059143282374E+01};
      i += ggddAmp.test_2to2_amp2([&]() { return ggddAmp.amp2_Feynman(); }, 0,0,md,md,pspatial,dataCH);
      i += ggddAmp.test_2to2_amp2([&]() { return ggddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      i += ggddAmp.test_2to2_amp2_rotations([&]() { return ggddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      i += ggddAmp.test_2to2_amp2_boosts([&]() { return ggddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      i += ggddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ggddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      //Helicities g+, g+ -> d, D
      ldouble dataCHpp[20] = {2.802581147451239E-07,3.089100592415728E-08,1.110260921101478E-08,5.706425354272925E-09,3.520735247856651E-09,2.443847655790371E-09,1.853824598007942E-09,1.515886524514282E-09,1.327986217601043E-09,1.243396828071652E-09,1.243396984820181E-09,1.327986222046272E-09,1.515886565551334E-09,1.853824603510268E-09,2.443847556856922E-09,3.520735182587680E-09,5.706425301363859E-09,1.110260924162994E-08,3.089100595050340E-08,2.802581149358418E-07};
      i += ggddAmp.test_2to2_amp2([&]() { return ggddAmp.amp2_gplus_gplus(); }, 0,0,md,md,pspatial,dataCHpp);
      i += ggddAmp.test_2to2_amp2_rotations([&]() { return ggddAmp.amp2_gplus_gplus(); }, 0,0,md,md,pspatial,dataCHpp);
      i += ggddAmp.test_2to2_amp2_boosts([&]() { return ggddAmp.amp2_gplus_gplus(); }, 0,0,md,md,pspatial,dataCHpp);
      i += ggddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ggddAmp.amp2_gplus_gplus(); }, 0,0,md,md,pspatial,dataCHpp);
      //Helicities g+, g- -> d, D
      ldouble dataCHpm[20] = {2.888118258538944E+01,8.203170976848307E+00,4.216485688409916E+00,2.604329303431043E+00,1.776981098254982E+00,1.302019225288287E+00,1.014447557952477E+00,8.388694661474058E-01,7.373966402741821E-01,6.907717398204531E-01,6.907717398204531E-01,7.373966402741818E-01,8.388694661474053E-01,1.014447557952476E+00,1.302019225288287E+00,1.776981098254982E+00,2.604329303431042E+00,4.216485688409914E+00,8.203170976848304E+00,2.888118258538937E+01};
      i += ggddAmp.test_2to2_amp2([&]() { return ggddAmp.amp2_gplus_gminus(); }, 0,0,md,md,pspatial,dataCHpm);
      i += ggddAmp.test_2to2_amp2_rotations([&]() { return ggddAmp.amp2_gplus_gminus(); }, 0,0,md,md,pspatial,dataCHpm);
      i += ggddAmp.test_2to2_amp2_boosts([&]() { return ggddAmp.amp2_gplus_gminus(); }, 0,0,md,md,pspatial,dataCHpm);
      i += ggddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ggddAmp.amp2_gplus_gminus(); }, 0,0,md,md,pspatial,dataCHpm);
      //Close to threshold
      pspatial = 0.008;
      ldouble dataCH2[20] = {4.967904333632854E-01,4.812084782526009E-01,4.669796505298057E-01,4.543188255478394E-01,4.433560925505287E-01,4.341669702215529E-01,4.267926609257171E-01,4.212536478062357E-01,4.175587456609391E-01,4.157109537043129E-01,4.157109537043129E-01,4.175587456609391E-01,4.212536478062357E-01,4.267926609257170E-01,4.341669702215528E-01,4.433560925505286E-01,4.543188255478393E-01,4.669796505298057E-01,4.812084782526009E-01,4.967904333632853E-01};
      i += ggddAmp.test_2to2_amp2([&]() { return ggddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      i += ggddAmp.test_2to2_amp2_rotations([&]() { return ggddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      i += ggddAmp.test_2to2_amp2_boosts([&]() { return ggddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      i += ggddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ggddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      //Helicities g+, g+ -> d, D
      ldouble dataCH2pp[20] = {9.704613515124714E-01,9.019322944671070E-01,8.454527556974414E-01,7.991666086663294E-01,7.616496388352150E-01,7.318116282473249E-01,7.088265881673369E-01,6.920830910376005E-01,6.811493740446898E-01,6.757496980644847E-01,6.757496980644846E-01,6.811493740446899E-01,6.920830910376005E-01,7.088265881673368E-01,7.318116282473247E-01,7.616496388352150E-01,7.991666086663293E-01,8.454527556974414E-01,9.019322944671071E-01,9.704613515124710E-01};
      i += ggddAmp.test_2to2_amp2([&]() { return ggddAmp.amp2_gplus_gplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      i += ggddAmp.test_2to2_amp2_rotations([&]() { return ggddAmp.amp2_gplus_gplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      i += ggddAmp.test_2to2_amp2_boosts([&]() { return ggddAmp.amp2_gplus_gplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      i += ggddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ggddAmp.amp2_gplus_gplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      //Helicities g+, g- -> d, D
      ldouble dataCH2pm[20] = {2.311951521409944E-02,6.048466203809481E-02,8.850654536216990E-02,1.094710424293493E-01,1.250625462658423E-01,1.365223121957808E-01,1.447587336840973E-01,1.504242045748710E-01,1.539681172771885E-01,1.556722093441412E-01,1.556722093441412E-01,1.539681172771885E-01,1.504242045748710E-01,1.447587336840973E-01,1.365223121957808E-01,1.250625462658423E-01,1.094710424293494E-01,8.850654536216990E-02,6.048466203809486E-02,2.311951521409951E-02};
      i += ggddAmp.test_2to2_amp2([&]() { return ggddAmp.amp2_gplus_gminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      i += ggddAmp.test_2to2_amp2_rotations([&]() { return ggddAmp.amp2_gplus_gminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      i += ggddAmp.test_2to2_amp2_boosts([&]() { return ggddAmp.amp2_gplus_gminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      i += ggddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ggddAmp.amp2_gplus_gminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }




  
  

}
