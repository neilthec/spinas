
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

//File:  SPINAS/SM/ubdt.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ubdt.h"

namespace spinas {
  //Constructors
  ubdt::ubdt(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& masst, const ldouble& massb, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), mu(massu), md(massd), mt(masst), mb(massb), MW(massW), WW(widthW), SW(sinW), prop(massW,widthW),
    p1(particle(mu)), p2(particle(mb)),
    p3(particle(md)), p4(particle(mt)),
    //<14>, [23], [12], <12>, [34], <34>
    s34s(sproduct(SQUARE,&p3,&p4,2)),
    a12a(sproduct(ANGLE,&p1,&p2,2)),
    a13a(sproduct(ANGLE,&p1,&p3,2)),
    s13s(sproduct(SQUARE,&p1,&p3,2)),
    s24s(sproduct(SQUARE,&p2,&p4,2)),
    a24a(sproduct(ANGLE,&p2,&p4,2))
  {
    preW = e*e/(4.0*MW*MW*SW*SW);
  }
  void ubdt::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& masst, const ldouble& massb, const ldouble& massW){
    mu=massu;
    md=massd;
    mt=masst;
    mb=massb;
    MW=massW;
    p1.set_mass(mu);
    p2.set_mass(mb);
    p3.set_mass(md);
    p4.set_mass(mt);
    prop.set_mass(massW);
    preW = e*e/(4.0*MW*MW*SW*SW);
  }
  void ubdt::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s34s.update();
    a12a.update();
    a13a.update();
    s13s.update();
    s24s.update();
    a24a.update();
    //Propagator Momentum
    ldouble propPS[4];
    for(int j=0;j<4;j++){
      propPS[j] = mom1[j]-mom3[j];
    }
    pDenT = prop.denominator(propPS);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ubdt::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    //preW = e*e/(4.0*MW*MW*SW*SW);
    //udtb all ingoing:
    //+ 2 preW ( 2MW^2 [23]⟨14⟩ + mdmt⟨12⟩⟨34⟩ − mtmu[12]⟨34⟩ − mbmd⟨12⟩[34] + mumb[12][34] )/(s-MW^2)
    //ubdt: 4->2->3->4
    //+ 2 preW ( 2MW^2 [34]⟨12⟩ - mdmt⟨13⟩⟨24⟩ + mtmu[13]⟨24⟩ + mbmd⟨13⟩[24] - mumb[13][24] )/(t-MW^2)
    //34 out:
    //+ 2 preW ( 2MW^2 [34]⟨12⟩ - mdmt⟨13⟩⟨24⟩ - mtmu[13]⟨24⟩ - mbmd⟨13⟩[24] - mumb[13][24] )/(t-MW^2)
    return preW*2.0*(
		     + 2.0*MW*MW*s34s.v(ds3,ds4)*a12a.v(ds1,ds2)
		     - md*mt*a13a.v(ds1,ds3)*a24a.v(ds2,ds4)
		     - mt*mu*s13s.v(ds1,ds3)*a24a.v(ds2,ds4)
		     - mb*md*a13a.v(ds1,ds3)*s24s.v(ds2,ds4)
		     - mu*mb*s13s.v(ds1,ds3)*s24s.v(ds2,ds4)
		     )/pDenT;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ubdt::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over initial colors (which must be the same) 1/3
    //Sum over final colors (which must be the same) 3
    return amp2/4.0;
  }
  



  //  Tests
  int test_ubdt(){
    int n=0;//Number of fails
    std::cout<<"\t* u , b  -> d , t       :";
    {//amp^2
      int i=0;
      // mu=0.0042, md=0.0075, mt=172.5, mb=4.25, pspatial=250
      ldouble mu=0.0042, md=0.0075, mt=172.5, mb=4.25, MW=80.385, WW=0;
      ldouble sinW=0.474;
      ubdt ubdtAmp = ubdt(0.31333,mu,md,mt,mb,MW,WW,sinW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.835322532096487E+01,4.977626481586548E+00,2.274873368042293E+00,1.297784136858694E+00,8.376775653910948E-01,5.850407501709444E-01,4.315810549746624E-01,3.314425218376886E-01,2.625057594132746E-01,2.130361331260763E-01,1.763384838307548E-01,1.483666833812242E-01,1.265586206702295E-01,1.092279105245741E-01,9.522780235226312E-02,8.375648533389614E-02,7.423967628509097E-02,6.625732014303620E-02,5.949643062635731E-02,5.371993559444489E-02};
      i += ubdtAmp.test_2to2_amp2([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH);
      i += ubdtAmp.test_2to2_amp2_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH);
      i += ubdtAmp.test_2to2_amp2_boosts([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH);
      i += ubdtAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH);
      // mu=0.0042, md=0.0075, mt=172.5, mb=4.25, pspatial=173
      pspatial = 173;
      ldouble dataCH2[20] = {6.778105887493330E+00,2.947416888291588E+00,1.640062948328296E+00,1.042870836450705E+00,7.211070739294910E-01,5.281705122139446E-01,4.034526537773075E-01,3.182089198962393E-01,2.573817328064005E-01,2.124641666035482E-01,1.783555321865865E-01,1.518467376951761E-01,1.308369152531970E-01,1.139039821010124E-01,1.000575125784040E-01,8.859054711050293E-02,7.898752966935225E-02,7.086532969322340E-02,6.393442645795351E-02,5.797274654454656E-02};
      i += ubdtAmp.test_2to2_amp2([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH2);
      i += ubdtAmp.test_2to2_amp2_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH2);
      i += ubdtAmp.test_2to2_amp2_boosts([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH2);
      i += ubdtAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH2);
      // mu=0.0042, md=0.0075, mt=0.005, pspatial=250
      mu=0.0042;
      md=0.0075;
      mt=0.005;
      mb=0.0047;
      pspatial = 250;
      ubdtAmp.set_masses(mu,md,mt,mb,MW);
      ldouble dataCH3[20] = {1.846307395909712E+01,4.693627842468835E+00,2.097786290591745E+00,1.183324589721133E+00,7.586068477418508E-01,5.274037424481904E-01,3.877921735530439E-01,2.970826554225933E-01,2.348422336033058E-01,1.902937023319428E-01,1.573158898185977E-01,1.322229288870205E-01,1.126877614631766E-01,9.718251374874448E-02,8.467041051920535E-02,7.442784741652674E-02,6.593736994857251E-02,5.882101179432095E-02,5.279749069030364E-02,4.765398016178167E-02};
      i += ubdtAmp.test_2to2_amp2([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH3);
      i += ubdtAmp.test_2to2_amp2_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH3);
      i += ubdtAmp.test_2to2_amp2_boosts([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH3);
      i += ubdtAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH3);
      // mu=0.0042, md=0.0075, mt=0.005, pspatial=0.005
      pspatial = 0.005;
      ldouble dataCH4[20] = {1.565842299687222E-17,1.565842298613354E-17,1.565842297539485E-17,1.565842296465617E-17,1.565842295391748E-17,1.565842294317879E-17,1.565842293244011E-17,1.565842292170142E-17,1.565842291096273E-17,1.565842290022405E-17,1.565842288948536E-17,1.565842287874667E-17,1.565842286800799E-17,1.565842285726930E-17,1.565842284653061E-17,1.565842283579193E-17,1.565842282505324E-17,1.565842281431455E-17,1.565842280357587E-17,1.565842279283718E-17};
      i += ubdtAmp.test_2to2_amp2([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH4);
      i += ubdtAmp.test_2to2_amp2_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH4);
      i += ubdtAmp.test_2to2_amp2_boosts([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH4);
      i += ubdtAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH4);
      // mu=0.0042, md=0.0075, mt=0.005, MW=0.006, pspatial=250
      MW=0.006;
      pspatial=250;
      ubdtAmp.set_masses(mu,md,mt,mb,MW);
      ldouble dataCH5[20] = {7.640762171129001E+01,8.518215247817096E+00,3.087062688499223E+00,1.590724735916037E+00,9.749478002130890E-01,6.632425330313395E-01,4.839658941015936E-01,3.714864028381312E-01,2.963147408449165E-01,2.436060789491521E-01,2.052266299963350E-01,1.764166504492835E-01,1.542402996459037E-01,1.368069703725236E-01,1.228545665355582E-01,1.115144791293067E-01,1.021730517559231E-01,9.438678139857519E-02,8.782868955183254E-02,8.225342518858512E-02};
      i += ubdtAmp.test_2to2_amp2([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH5);
      i += ubdtAmp.test_2to2_amp2_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH5);
      i += ubdtAmp.test_2to2_amp2_boosts([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH5);
      i += ubdtAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH5);
      // mu=0.0042, md=0.0075, mt=0.005, MW=0.006, pspatial=0.005
      pspatial=0.005;
      ldouble dataCH6[20] = {4.112546003207651E-01,3.729623639685265E-01,3.402503797228241E-01,3.120820183832028E-01,2.876508883267742E-01,2.663220600752031E-01,2.475900840837934E-01,2.310485329364645E-01,2.163675886577873E-01,2.032773346194554E-01,1.915551506928970E-01,1.810160988867555E-01,1.715055150883264E-01,1.628932466510192E-01,1.550691307256575E-01,1.479394170700120E-01,1.414239163582433E-01,1.354537105266921E-01,1.299693019968880E-01,1.249191081728294E-01};
      i += ubdtAmp.test_2to2_amp2([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH6);
      i += ubdtAmp.test_2to2_amp2_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH6);
      i += ubdtAmp.test_2to2_amp2_boosts([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH6);
      i += ubdtAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH6);
      // mu=0.0042, md=0.0075, mt=0.005, MW=0.0006, pspatial=0.005
      MW=0.0006;
      ubdtAmp.set_masses(mu,md,mt,mb,MW);
      ldouble dataCH7[20] = {1.418679042734380E+03,1.108390218517031E+03,9.342449898920364E+02,8.229445702338417E+02,7.457285096516597E+02,6.890414358048592E+02,6.456669690634877E+02,6.114131515368401E+02,5.836794170721334E+02,5.607673590156892E+02,5.415209650270815E+02,5.251262157683636E+02,5.109933838405875E+02,4.986847933763456E+02,4.878688008940016E+02,4.782895303375593E+02,4.697464115240600E+02,4.620800079497204E+02,4.551619890168397E+02,4.488878987453491E+02};
      i += ubdtAmp.test_2to2_amp2([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH7);
      i += ubdtAmp.test_2to2_amp2_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH7);
      i += ubdtAmp.test_2to2_amp2_boosts([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH7);
      i += ubdtAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ubdtAmp.amp2(); }, mu,mb,md,mt,pspatial,dataCH7);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
