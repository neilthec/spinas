
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

//File:  SPINAS/SM/dhdh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dhdh.h"

namespace spinas {
  //Constructors
  dhdh::dhdh(const ldouble& echarge, const ldouble& massd, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW):
    e(echarge), md(massd), mh(massh), wh(widthh), MW(massW), SW(sinW),
    propd(md,0), proph(mh,wh),
    p1(particle(md)), p2(particle(mh)),
    p3(particle(md)), p4(particle(mh)),
    s13s(sproduct(SQUARE,&p1,&p3,2)),
    a13a(sproduct(ANGLE,&p1,&p3,2)),
    s123a(sproduct(SQUARE,&p1,&p2,&p3,2)),
    s321a(sproduct(SQUARE,&p3,&p2,&p1,2))
  {
    prehS = 3*e*e*md*mh*mh/(4*MW*MW*SW*SW);
    prehTU = e*e*md*md/(4.0*MW*MW*SW*SW);
  }
  void dhdh::set_masses(const ldouble& massd, const ldouble& massh, const ldouble& massW){
    //constexpr ldouble sqrt2 = std::sqrt(2);
    md=massd;
    propd.set_mass(md);
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    p1.set_mass(md);
    p2.set_mass(mh);
    p3.set_mass(md);
    p4.set_mass(mh);
    prehS = 3.0*e*e*md*mh*mh/(4.0*MW*MW*SW*SW);
    prehTU = e*e*md*md/(4.0*MW*MW*SW*SW);
  }
  void dhdh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s13s.update();
    a13a.update();
    s123a.update();
    s321a.update();
    //Propagator Momentum
    ldouble propPS[4], propPT[4], propPU[4];
    for(int j=0;j<4;j++){
      propPS[j] = mom1[j]+mom2[j];
      propPT[j] = mom1[j]-mom3[j];
      propPU[j] = mom1[j]-mom4[j];
    }
    pDenTh = proph.denominator(propPT);
    pDenSd = propd.denominator(propPS);
    pDenUd = propd.denominator(propPU);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble dhdh::amp(const int& ds1, const int& ds3){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Higgs
    //prehS = 3*e*e*md*mh*mh/(4*MW*MW*SW*SW);
    //S Channel
    //dDhh all in:     + prehS*([12]+<12>)/(s-mh^2)
    //dhDh: 2<->3:     + prehS*([13]+<13>)/(t-mh^2)
    //34 out:          + prehS*([13]-<13>)/(t-mh^2)
    amplitude += + prehS*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3))/pDenTh;

    //T Channel
    //prehTU = e*e*md*md/(4*MW*MW*SW*SW);
    //dDhh all in:     + prehTU * ( 2md([12]+<12>) + [132> + [231> )/(t-md^2)
    //dhDh: 2<->3:     + prehTU * ( 2md([13]+<13>) + [123> + [321> )/(s-md^2)
    //34 out:          + prehTU * ( 2md([13]-<13>) - [123> + [321> )/(s-md^2)
    amplitude += + prehTU*( two*md*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3)) - s123a.v(ds1,ds3) + s321a.v(ds3,ds1) )/pDenSd;

    //U Channel
    //dDhh all in:     + prehTU * ( 2md([12]+<12>) - [132> - [231> )/(u-md^2)
    //dhDh: 2<->3:     + prehTU * ( 2md([13]+<13>) - [123> - [321> )/(u-md^2)
    //34 out:          + prehTU * ( 2md([13]-<13>) + [123> - [321> )/(u-md^2)
    amplitude += + prehTU*( two*md*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3)) + s123a.v(ds1,ds3) - s321a.v(ds3,ds1) )/pDenUd;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble dhdh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-1;j3<=1;j3+=2){
	M = amp(j1,j3);
	amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2
    //Average over colors 1/3
    return amp2/6.0;
  }

  



  //  Tests
  int test_dhdh(){
    int n=0;//Number of fails
    std::cout<<"\t* d , h  -> d , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### md=0.0042, mh=125, MW=80.385, pspatial=250\n";
      ldouble md=0.0075, mh=125, wh=0, MW=80.385, SW=0.474;//Set width to 0 for comparison with Feynman diagrams.
      dhdh dhdhAmp = dhdh(0.31333,md,mh,wh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.613850175759994E-10,5.605255946548275E-10,5.023970140898474E-10,4.383796662584196E-10,3.846328752560225E-10,3.411337759361986E-10,3.058296092399333E-10,2.768310093904151E-10,2.526848907072952E-10,2.323144679936743E-10,2.149231255066567E-10,1.999157699553016E-10,1.868418683405860E-10,1.753554687651109E-10,1.651872318716851E-10,1.561247187156038E-10,1.479983348563931E-10,1.406711751599823E-10,1.340315891774055E-10,1.279876949547878E-10};
      i += dhdhAmp.test_2to2_amp2([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH);
      i += dhdhAmp.test_2to2_amp2_rotations([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH);
      i += dhdhAmp.test_2to2_amp2_boosts([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH);
      i += dhdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH);
      //std::cout<<"\n########### md=0.0042, mh=125, MW=80.385, pspatial=126\n";
      pspatial = 126;
      ldouble dataCH2[20] = {1.892900344160835E-10,4.047633775552472E-10,5.050443538381120E-10,5.491037204062651E-10,5.640682216903785E-10,5.634515538907888E-10,5.543925822955673E-10,5.408292231032814E-10,5.250019570383787E-10,5.082113689914925E-10,4.912194494838085E-10,4.744711343594867E-10,4.582208622446343E-10,4.426069567259529E-10,4.276963838414671E-10,4.135122060593660E-10,4.000506814295029E-10,3.872920381559901E-10,3.752073557874912E-10,3.637626504508571E-10};
      i += dhdhAmp.test_2to2_amp2([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH2);
      i += dhdhAmp.test_2to2_amp2_rotations([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH2);
      i += dhdhAmp.test_2to2_amp2_boosts([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH2);
      i += dhdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH2);
      //std::cout<<"\n########### md=0.1, mh=0.15, MW=80.385, pspatial=0.2\n";
      md=0.10;
      mh=0.15;
      pspatial = 0.2;
      dhdhAmp.set_masses(md,mh,MW);
      ldouble dataCH4[20] = {3.529499456952278E-13,2.415383046226727E-13,1.902207104895969E-13,1.635010706027575E-13,1.486629047580430E-13,1.403097060040832E-13,1.358839101205679E-13,1.340913589125564E-13,1.342692504853799E-13,1.361108047071021E-13,1.395443696229677E-13,1.446919502533162E-13,1.518828630816734E-13,1.617274427557318E-13,1.752893438926153E-13,1.944660273991694E-13,2.228904216766091E-13,2.683853012883027E-13,3.512814823054526E-13,5.454659082625271E-13};
      i += dhdhAmp.test_2to2_amp2([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH4);
      i += dhdhAmp.test_2to2_amp2_rotations([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH4);
      i += dhdhAmp.test_2to2_amp2_boosts([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH4);
      i += dhdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH4);
      //std::cout<<"\n########### md=0.1, mh=0.15, MW=0.11, pspatial=0.2\n";
      md=0.10;
      mh=0.15;
      MW=0.11;
      pspatial = 0.2;
      dhdhAmp.set_masses(md,mh,MW);
      ldouble dataCH5[20] = {1.006566432386016E-01,6.888352088841179E-02,5.424842368123804E-02,4.662833677555567E-02,4.239668867937969E-02,4.001446718555772E-02,3.875228890014571E-02,3.824107707072800E-02,3.829180938787165E-02,3.881699622685900E-02,3.979620340053371E-02,4.126422512250993E-02,4.331497808607276E-02,4.612252163770529E-02,4.999019596666304E-02,5.545913175702459E-02,6.356539200427269E-02,7.653992825823670E-02,1.001808196091335E-01,1.555596423698312E-01};
      i += dhdhAmp.test_2to2_amp2([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH5);
      i += dhdhAmp.test_2to2_amp2_rotations([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH5);
      i += dhdhAmp.test_2to2_amp2_boosts([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH5);
      i += dhdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH5);
      //std::cout<<"\n########### md=0.1, mh=0.15, MW=0.006, pspatial=0.2\n";
      md=0.10;
      mh=0.15;
      MW=0.006;
      pspatial = 0.2;
      dhdhAmp.set_masses(md,mh,MW);
      ldouble dataCH6[20] = {1.137124933376825E+04,7.781818127524977E+03,6.128481258618874E+03,5.267634866750854E+03,4.789582707984553E+03,4.520461528269680E+03,4.377872390332047E+03,4.320120442843585E+03,4.325851707159173E+03,4.385182420967922E+03,4.495804120271712E+03,4.661647531008240E+03,4.893322485788513E+03,5.210492587173172E+03,5.647426382314148E+03,6.265255772026211E+03,7.181025496408613E+03,8.646767666889225E+03,1.131749521525712E+04,1.757367842543748E+04};
      i += dhdhAmp.test_2to2_amp2([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH6);
      i += dhdhAmp.test_2to2_amp2_rotations([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH6);
      i += dhdhAmp.test_2to2_amp2_boosts([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH6);
      i += dhdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dhdhAmp.amp2(); }, md,mh,md,mh,pspatial,dataCH6);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
