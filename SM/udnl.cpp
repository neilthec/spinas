
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

//File:  SPINAS/SM/udnl.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/udnl.h"

namespace spinas {
  //Constructors
  udnl::udnl(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& masstau, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), mu(massu), md(massd), ml(masstau), MW(massW), WW(widthW), SW(sinW), prop(massW,widthW),
    p1(particle(mu)), p2(particle(md)),
    p3(particle(0)), p4(particle(ml)),
    //[14], <23>, <34>, [12], <12>
    s14s(sproduct(SQUARE,&p1,&p4,2)),
    a23a(sproduct(ANGLE,&p2,&p3,2)),
    s12s(sproduct(SQUARE,&p1,&p2,2)),
    a12a(sproduct(ANGLE,&p1,&p2,2)),
    a34a(sproduct(ANGLE,&p3,&p4,2))
  {}
  void udnl::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& masstau, const ldouble& massW){
    mu=massu;
    md=massd;
    ml=masstau;
    MW=massW;
    p1.set_mass(mu);
    p2.set_mass(md);
    p3.set_mass(0);
    p4.set_mass(ml);
    prop.set_mass(massW);
  }
  void udnl::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s14s.update();
    a23a.update();
    s12s.update();
    a12a.update();
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
  cdouble udnl::amp(const int& ds1, const int& ds2, const int& ds4){
    //All ingoing: (2MW^2 [14]<23> + mlmu<34><12> - mlmd<34>[12])/ (2MW^2 (s−MW^2))
    //Sign changes due to p3 and p4 being outgoing.
    // - (2MW^2 <14>[23] - mlmu<34>[12] + mlmd<34><12>)/ (2MW^2 (s−MW^2))
    return e*e/2.0/SW/SW*(  2.0*MW*MW*s14s.v(ds1,ds4)*a23a.v(ds2)
			  + ml*md*a34a.v(ds4)*s12s.v(ds1,ds2) - mu*ml*a12a.v(ds1,ds2)*a34a.v(ds4)
			  )/(MW*MW*pDenS);
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble udnl::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,j2,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over initial colors (which must be the same) 1/3
    return amp2/12.0;
  }
  



  //  Tests
  int test_udnl(){
    int n=0;//Number of fails
    std::cout<<"\t* u , D  -> nl, L       :";
    {//amp^2
      int i=0;
      // mu=0.0042, md=0.0075, mtau=1.777, pspatial=250
      //When testing against CH, it is better to set the width to 0 (both here and in CH) because
      //CH sets a cutoff for using the width and interpolates between the regions.  This is necessary to
      //keep gauge invariance in Feynman diagrams.  See the CH manual.
      //Our spinor "diagrams" are trivially gauge invariant so we don't have such a constraint.
      ldouble mu=0.0042, md=0.0075, mtau=1.777, MW=80.385, WW=0;
      ldouble sinW=0.474;
      udnl udnlAmp = udnl(0.31333,mu,md,mtau,MW,WW,sinW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.593906484361222E-02,1.434621633879787E-02,1.283720140848269E-02,1.141202005266669E-02,1.007067227134987E-02,8.813158064532230E-03,7.639477432213766E-03,6.549630374394482E-03,5.543616891074376E-03,4.621436982253448E-03,3.783090647931700E-03,3.028577888109130E-03,2.357898702785739E-03,1.771053091961527E-03,1.268041055636493E-03,8.488625938106384E-04,5.135177064839621E-04,2.620063936564647E-04,9.432865532814595E-05,1.048449149900575E-05};
      i += udnlAmp.test_2to2_amp2([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH);
      i += udnlAmp.test_2to2_amp2_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH);
      i += udnlAmp.test_2to2_amp2_boosts([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH);
      i += udnlAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH);
      // mu=0.0042, md=0.0075, mtau=1.777, pspatial=1.8
      pspatial = 1.8;
      ldouble dataCH2[20] = {4.649561930386768E-08,4.241091230335185E-08,3.851002356204619E-08,3.479295307995068E-08,3.125970085706533E-08,2.791026689339014E-08,2.474465118892511E-08,2.176285374367024E-08,1.896487455762553E-08,1.635071363079098E-08,1.392037096316659E-08,1.167384655475236E-08,9.611140405548282E-09,7.732252515554369E-09,6.037182884770615E-09,4.525931513197021E-09,3.198498400833586E-09,2.054883547680311E-09,1.095086953737196E-09,3.191086190042396E-10};
      i += udnlAmp.test_2to2_amp2([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH2);
      i += udnlAmp.test_2to2_amp2_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH2);
      i += udnlAmp.test_2to2_amp2_boosts([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH2);
      i += udnlAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH2);
      // mu=0.0042, md=0.0075, mtau=0.005, pspatial=250
      mu=0.0042;
      md=0.0075;
      mtau=0.005;
      pspatial = 250;
      udnlAmp.set_masses(mu,md,mtau,MW);
      ldouble dataCH3[20] = {1.593926100732559E-02,1.434638285306197E-02,1.283734039111057E-02,1.141213362147137E-02,1.007076254414439E-02,8.813227159129610E-03,7.639527466427045E-03,6.549663466036689E-03,5.543635157958545E-03,4.621442542192610E-03,3.783085618738886E-03,3.028564387597372E-03,2.357878848768069E-03,1.771029002250976E-03,1.268014848046093E-03,8.488363861534209E-04,5.134936165729594E-04,2.619865393047080E-04,9.431515434866716E-05,1.047946170483660E-05};
      i += udnlAmp.test_2to2_amp2([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH3);
      i += udnlAmp.test_2to2_amp2_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH3);
      i += udnlAmp.test_2to2_amp2_boosts([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH3);
      i += udnlAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH3);
      // mu=0.0042, md=0.0075, mtau=0.005, pspatial=0.001
      pspatial = 0.001;
      ldouble dataCH4[20] = {2.211682187158320E-18,2.155034754373015E-18,2.099116321267990E-18,2.043926887843247E-18,1.989466454098787E-18,1.935735020034607E-18,1.882732585650709E-18,1.830459150947093E-18,1.778914715923759E-18,1.728099280580706E-18,1.678012844917934E-18,1.628655408935445E-18,1.580026972633237E-18,1.532127536011311E-18,1.484957099069666E-18,1.438515661808303E-18,1.392803224227221E-18,1.347819786326422E-18,1.303565348105904E-18,1.260039909565667E-18};
      i += udnlAmp.test_2to2_amp2([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH4);
      i += udnlAmp.test_2to2_amp2_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH4);
      i += udnlAmp.test_2to2_amp2_boosts([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH4);
      i += udnlAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH4);
      // mu=0.0042, md=0.0075, mtau=0.005, MW=0.006, pspatial=250
      MW=0.006;
      pspatial=250;
      udnlAmp.set_masses(mu,md,mtau,MW);
      ldouble dataCH5[20] = {2.079582710287643E-02,1.928422701933562E-02,1.785218483490748E-02,1.649970054959202E-02,1.522677416338923E-02,1.403340567629911E-02,1.291959508832167E-02,1.188534239945690E-02,1.093064760970480E-02,1.005551071906537E-02,9.259931727538619E-03,8.543910635124540E-03,7.907447441823132E-03,7.350542147634397E-03,6.873194752558335E-03,6.475405256594946E-03,6.157173659744230E-03,5.918499962006187E-03,5.759384163380816E-03,5.679826263868119E-03};
      i += udnlAmp.test_2to2_amp2([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH5);
      i += udnlAmp.test_2to2_amp2_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH5);
      i += udnlAmp.test_2to2_amp2_boosts([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH5);
      i += udnlAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH5);
      // mu=0.0042, md=0.0075, mtau=0.005, MW=0.006, pspatial=0.001
      pspatial=0.001;
      ldouble dataCH6[20] = {1.222991260421230E-02,1.197860481409388E-02,1.173004615837524E-02,1.148423663705638E-02,1.124117625013730E-02,1.100086499761800E-02,1.076330287949848E-02,1.052848989577874E-02,1.029642604645878E-02,1.006711133153860E-02,9.840545751018200E-03,9.616729304897577E-03,9.395661993176735E-03,9.177343815855672E-03,8.961774772934389E-03,8.748954864412885E-03,8.538884090291162E-03,8.331562450569218E-03,8.126989945247052E-03,7.925166574324669E-03};
      i += udnlAmp.test_2to2_amp2([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH6);
      i += udnlAmp.test_2to2_amp2_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH6);
      i += udnlAmp.test_2to2_amp2_boosts([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH6);
      i += udnlAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH6);
      // mu=0.0042, md=0.0075, mtau=0.005, MW=0.0006, pspatial=0.001
      MW=0.0006;
      udnlAmp.set_masses(mu,md,mtau,MW);
      ldouble dataCH7[20] = {4.002015006185405E+01,4.001792809351610E+01,4.001570765918045E+01,4.001348875884709E+01,4.001127139251602E+01,4.000905556018724E+01,4.000684126186074E+01,4.000462849753654E+01,4.000241726721463E+01,4.000020757089501E+01,3.999799940857768E+01,3.999579278026264E+01,3.999358768594989E+01,3.999138412563943E+01,3.998918209933125E+01,3.998698160702537E+01,3.998478264872178E+01,3.998258522442048E+01,3.998038933412147E+01,3.997819497782475E+01};
      i += udnlAmp.test_2to2_amp2([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH7);
      i += udnlAmp.test_2to2_amp2_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH7);
      i += udnlAmp.test_2to2_amp2_boosts([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH7);
      i += udnlAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udnlAmp.amp2(); }, mu,md,0,mtau,pspatial,dataCH7);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
