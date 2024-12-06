
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

//File:  SPINAS/SM/udtb.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/udtb.h"

namespace spinas {
  //Constructors
  udtb::udtb(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& masst, const ldouble& massb, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), mu(massu), md(massd), mt(masst), mb(massb), MW(massW), WW(widthW), SW(sinW), prop(massW,widthW),
    p1(particle(mu)), p2(particle(md)),
    p3(particle(mt)), p4(particle(mb)),
    //<14>, [23], [12], <12>, [34], <34>
    a14a(sproduct(ANGLE,&p1,&p4,2)),
    s23s(sproduct(SQUARE,&p2,&p3,2)),
    s12s(sproduct(SQUARE,&p1,&p2,2)),
    a12a(sproduct(ANGLE,&p1,&p2,2)),
    s34s(sproduct(SQUARE,&p3,&p4,2)),
    a34a(sproduct(ANGLE,&p3,&p4,2))
  {}
  void udtb::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& masst, const ldouble& massb, const ldouble& massW){
    mu=massu;
    md=massd;
    mt=masst;
    mb=massb;
    MW=massW;
    p1.set_mass(mu);
    p2.set_mass(md);
    p3.set_mass(mt);
    p4.set_mass(mb);
    prop.set_mass(massW);
  }
  void udtb::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a14a.update();
    s23s.update();
    s12s.update();
    a12a.update();
    s34s.update();
    a34a.update();
    //Propagator Momentum
    ldouble propPS[4];
    for(int j=0;j<4;j++){
      propPS[j] = mom1[j]+mom2[j];
    }
    pDenS = prop.denominator(propPS);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble udtb::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    //All ingoing tBUd: 2MW^2 [23]⟨14⟩ + mbmu⟨12⟩⟨34⟩ − mtmu[12]⟨34⟩ − mbmd⟨12⟩[34] + mtmd[12][34]
    //All ingoing uDTb: 2MW^2 [23]⟨14⟩ + mdmt⟨12⟩⟨34⟩ − mtmu[12]⟨34⟩ − mbmd⟨12⟩[34] + mumb[12][34]
    //Sign changes due to p3 and p4 being outgoing.
    // - ( 2MW^2 [23]⟨14⟩ - mdmt⟨12⟩⟨34⟩ + mtmu[12]⟨34⟩ + mbmd⟨12⟩[34] - mumb[12][34])
    return e*e/2.0/SW/SW*(  2.0*MW*MW*a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  - md*mt*a34a.v(ds3,ds4)*a12a.v(ds1,ds2) + mt*mu*s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
			  + mb*md*s34s.v(ds3,ds4)*a12a.v(ds1,ds2) - mu*mb*s12s.v(ds1,ds2)*s34s.v(ds3,ds4)
			  )/(MW*MW*pDenS);
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble udtb::amp2(){
    constexpr ldouble nine=9;
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += nine*std::pow(std::abs(M),2);//Color factor 9
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over initial colors 1/3^2=1/9
    return amp2/36.0;
  }
  



  //  Tests
  int test_udtb(){
    int n=0;//Number of fails
    std::cout<<"\t* u , D  -> t , B       :";
    {//amp^2
      int i=0;
      // mu=0.0042, md=0.0075, mt=172.5, mb=4.25, pspatial=250
      //When testing against CH, it is better to set the width to 0 (both here and in CH) because
      //CH sets a cutoff for using the width and interpolates between the regions.  This is necessary to
      //keep gauge invariance in Feynman diagrams.  See the CH manual.
      //Our spinor "diagrams" are trivially gauge invariant so we don't have such a constraint.
      ldouble mu=0.0042, md=0.0075, mt=172.5, mb=4.25, MW=80.385, WW=0;
      ldouble sinW=0.474;
      udtb udtbAmp = udtb(0.31333,mu,md,mt,mb,MW,WW,sinW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.225102641152419E-02,3.827911861371615E-02,3.450236905443277E-02,3.092077773367405E-02,2.753434465143998E-02,2.434306980773057E-02,2.134695320254582E-02,1.854599483588572E-02,1.594019470775028E-02,1.352955281813949E-02,1.131406916705337E-02,9.293743754491894E-03,7.468576580455081E-03,5.838567644942925E-03,4.403716947955424E-03,3.164024489492582E-03,2.119490269554397E-03,1.270114288140868E-03,6.158965452519952E-04,1.568370408877803E-04};
      i += udtbAmp.test_2to2_amp2([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH);
      i += udtbAmp.test_2to2_amp2_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH);
      i += udtbAmp.test_2to2_amp2_boosts([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH);
      i += udtbAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH);
      // mu=0.0042, md=0.0075, mt=172.5, mb=4.25, pspatial=173
      pspatial = 173;
      ldouble dataCH2[20] = {3.833383790526305E-02,3.497612170678214E-02,3.176889360832107E-02,2.871215360987981E-02,2.580590171145839E-02,2.305013791305680E-02,2.044486221467503E-02,1.799007461631309E-02,1.568577511797098E-02,1.353196371964870E-02,1.152864042134625E-02,9.675805223063622E-03,7.973458124800825E-03,6.421599126557855E-03,5.020228228334714E-03,3.769345430131402E-03,2.668950731947917E-03,1.719044133784262E-03,9.196256356404343E-04,2.706952375164343E-04};
      i += udtbAmp.test_2to2_amp2([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH2);
      i += udtbAmp.test_2to2_amp2_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH2);
      i += udtbAmp.test_2to2_amp2_boosts([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH2);
      i += udtbAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH2);
      // mu=0.0042, md=0.0075, mt=0.005, pspatial=250
      mu=0.0042;
      md=0.0075;
      mt=0.005;
      mb=0.0047;
      pspatial = 250;
      udtbAmp.set_masses(mu,md,mt,mb,MW);
      ldouble dataCH3[20] = {4.781778301785994E-02,4.303914855569133E-02,3.851202117041491E-02,3.423640086203068E-02,3.021228763053863E-02,2.643968147593877E-02,2.291858239823109E-02,1.964899039741560E-02,1.663090547349228E-02,1.386432762646116E-02,1.134925685632222E-02,9.085693163075465E-03,7.073636546720894E-03,5.313087007258507E-03,3.804044544688305E-03,2.546509159010288E-03,1.540480850224457E-03,7.859596183308101E-04,2.829454633293479E-04,3.143838522007047E-05};
      i += udtbAmp.test_2to2_amp2([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH3);
      i += udtbAmp.test_2to2_amp2_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH3);
      i += udtbAmp.test_2to2_amp2_boosts([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH3);
      i += udtbAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH3);
      // mu=0.0042, md=0.0075, mt=0.005, pspatial=0.001
      pspatial = 0.001;
      ldouble dataCH4[20] = {6.379146164736003E-18,6.258020160288300E-18,6.137971176298778E-18,6.018999212767435E-18,5.901104269694272E-18,5.784286347079290E-18,5.668545444922489E-18,5.553881563223867E-18,5.440294701983426E-18,5.327784861201165E-18,5.216352040877084E-18,5.105996241011184E-18,4.996717461603464E-18,4.888515702653925E-18,4.781390964162566E-18,4.675343246129386E-18,4.570372548554387E-18,4.466478871437569E-18,4.363662214778931E-18,4.261922578578473E-18};
      i += udtbAmp.test_2to2_amp2([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH4);
      i += udtbAmp.test_2to2_amp2_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH4);
      i += udtbAmp.test_2to2_amp2_boosts([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH4);
      i += udtbAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH4);
      // mu=0.0042, md=0.0075, mt=0.005, MW=0.006, pspatial=250
      MW=0.006;
      pspatial=250;
      udtbAmp.set_masses(mu,md,mt,mb,MW);
      ldouble dataCH5[20] = {7.741720325731080E-02,7.288240300756464E-02,6.858627645511432E-02,6.452882359995986E-02,6.071004444210122E-02,5.712993898153843E-02,5.378850721827148E-02,5.068574915230036E-02,4.782166478362509E-02,4.519625411224566E-02,4.280951713816206E-02,4.066145386137431E-02,3.875206428188240E-02,3.708134839968632E-02,3.564930621478608E-02,3.445593772718169E-02,3.350124293687314E-02,3.278522184386042E-02,3.230787444814355E-02,3.206920074972252E-02};
      i += udtbAmp.test_2to2_amp2([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH5);
      i += udtbAmp.test_2to2_amp2_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH5);
      i += udtbAmp.test_2to2_amp2_boosts([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH5);
      i += udtbAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH5);
      // mu=0.0042, md=0.0075, mt=0.005, MW=0.006, pspatial=0.001
      pspatial=0.001;
      ldouble dataCH6[20] = {4.818563315038857E-02,4.773808896709074E-02,4.729460634111719E-02,4.685518527246794E-02,4.641982576114297E-02,4.598852780714228E-02,4.556129141046589E-02,4.513811657111377E-02,4.471900328908594E-02,4.430395156438240E-02,4.389296139700315E-02,4.348603278694817E-02,4.308316573421749E-02,4.268436023881109E-02,4.228961630072898E-02,4.189893391997115E-02,4.151231309653761E-02,4.112975383042836E-02,4.075125612164339E-02,4.037681997018271E-02};
      i += udtbAmp.test_2to2_amp2([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH6);
      i += udtbAmp.test_2to2_amp2_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH6);
      i += udtbAmp.test_2to2_amp2_boosts([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH6);
      i += udtbAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH6);
      // mu=0.0042, md=0.0075, mt=0.005, MW=0.0006, pspatial=0.001
      MW=0.0006;
      udtbAmp.set_masses(mu,md,mt,mb,MW);
      ldouble dataCH7[20] = {2.742798922212869E+02,2.742801526328648E+02,2.742804153107700E+02,2.742806802550024E+02,2.742809474655621E+02,2.742812169424490E+02,2.742814886856632E+02,2.742817626952046E+02,2.742820389710733E+02,2.742823175132693E+02,2.742825983217926E+02,2.742828813966431E+02,2.742831667378209E+02,2.742834543453259E+02,2.742837442191582E+02,2.742840363593178E+02,2.742843307658046E+02,2.742846274386187E+02,2.742849263777600E+02,2.742852275832286E+02};
      i += udtbAmp.test_2to2_amp2([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH7);
      i += udtbAmp.test_2to2_amp2_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH7);
      i += udtbAmp.test_2to2_amp2_boosts([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH7);
      i += udtbAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udtbAmp.amp2(); }, mu,md,mt,mb,pspatial,dataCH7);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
