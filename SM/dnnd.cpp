
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

//File:  SPINAS/SM/dnnd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dnnd.h"

namespace spinas {
  //Constructors
  dnnd::dnnd(const ldouble& echarge, const ldouble& massd, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), md(massd), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propZ(MZ,WZ),
    p1(particle(md)), p2(particle(0)),
    p3(particle(0)), p4(particle(md)),
    s34s(sproduct(SQUARE,&p3,&p4)),
    a12a(sproduct(ANGLE,&p1,&p2)),
    s13s(sproduct(SQUARE,&p1,&p3)),
    a24a(sproduct(ANGLE,&p2,&p4))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    gLd=-1.0+2.0/3.0*SW*SW;
    gRd=2.0/3.0*SW*SW;
    gLn=1.0;
    gRn=0;
    preZ = e*e/(4.0*CW*CW*SW*SW);
  }
  void dnnd::set_masses(const ldouble& massd, const ldouble& massW){
    md=massd;
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(md);
    p4.set_mass(md);
  }
  void dnnd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s34s.update();
    a12a.update();
    s13s.update();
    a24a.update();
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]-mom4[j];
    pDenUZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble dnnd::amp(const int& ds1, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Z Boson
    //Defined above:
    //gLd=-1.0+2.0/3.0*SW*SW;
    //gRd=2.0/3.0*SW*SW;
    //gLn=1.0;
    //gRn=0;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //dDNn all in:
    // - preZ ( gLd gLn [23] <14> + gLn gRd [13] <24>)/(s-MZ^2)
    //dnND: 2<->4
    // + preZ ( gLd gLn [34] <12> + gLn gRd [13] <24>)/(u-MZ^2)
    //34 out:
    // + preZ ( gLd gLn [34] <12> - gLn gRd [13] <24>)/(u-MZ^2)
    amplitude += 
      + two*preZ*gLn*(
		      + gLd*s34s.v(ds4)*a12a.v(ds1)
		      - gRd*s13s.v(ds1)*a24a.v(ds4)
		      )/pDenUZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble dnnd::amp2(){
    ldouble amp2 = 0, two=2, three = 3;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2){
	M = amp(j1,j2);
	amp2 += three*std::pow(std::abs(M),2);// Color factor 3
      }
    //Average over initial spins 1/2
    //Average over colors 1/3
    return amp2/6.0;
  }

  



  //  Tests
  int test_dnnd(){
    int n=0;//Number of fails
    std::cout<<"\t* d , nm -> nm, d       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### md=0.0042, MW=80.385, pspatial=250\n";
      ldouble md=0.0042, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      dnnd dnndAmp = dnnd(0.31333,md,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.822888026752336E-02,3.125618250998174E-02,3.480271941659601E-02,3.899362050557115E-02,4.399401014437376E-02,5.002536758366299E-02,5.739025552125183E-02,6.651068974952962E-02,7.798947527434627E-02,9.271166813491261E-02,1.120192119072925E-01,1.380258881439003E-01,1.742178701259082E-01,2.266791789225574E-01,3.068124794502121E-01,4.380711249934854E-01,6.752153677595792E-01,1.172008716283988E+00,2.510298329877469E+00,8.682328335650210E+00};
      i += dnndAmp.test_2to2_amp2([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dnndAmp.test_2to2_amp2_rotations([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dnndAmp.test_2to2_amp2_boosts([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dnndAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      //std::cout<<"########### md=0.0005, MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {5.733054360670961E-20,5.685196568965621E-20,5.637441372051610E-20,5.589788769928942E-20,5.542238762597634E-20,5.494791350057698E-20,5.447446532309151E-20,5.400204309352005E-20,5.353064681186278E-20,5.306027647811983E-20,5.259093209229135E-20,5.212261365437750E-20,5.165532116437840E-20,5.118905462229422E-20,5.072381402812510E-20,5.025959938187120E-20,4.979641068353265E-20,4.933424793310960E-20,4.887311113060220E-20,4.841300027601061E-20};
      i += dnndAmp.test_2to2_amp2([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dnndAmp.test_2to2_amp2_rotations([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dnndAmp.test_2to2_amp2_boosts([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dnndAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      //std::cout<<"########### md=0.105, MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      dnndAmp.set_masses(md,MW);
      ldouble dataCH3[20] = {4.603900249991003E-05,4.571684015237733E-05,4.539458341244120E-05,4.507223279160570E-05,4.474978880469286E-05,4.442725196986059E-05,4.410462280862059E-05,4.378190184585648E-05,4.345908960984203E-05,4.313618663225937E-05,4.281319344821740E-05,4.249011059627035E-05,4.216693861843630E-05,4.184367806021592E-05,4.152032947061130E-05,4.119689340214478E-05,4.087337041087808E-05,4.054976105643138E-05,4.022606590200250E-05,3.990228551438636E-05};
      i += dnndAmp.test_2to2_amp2([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH3);
      i += dnndAmp.test_2to2_amp2_rotations([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH3);
      i += dnndAmp.test_2to2_amp2_boosts([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH3);
      i += dnndAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH3);
      //std::cout<<"########### md=0.105, MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      dnndAmp.set_masses(md,MW);
      ldouble dataCH4[20] = {2.090672929700221E-01,2.277143949984637E-01,2.491597003220729E-01,2.740030969934361E-01,3.030155914715371E-01,3.372016576913164E-01,3.778898619387558E-01,4.268674525757845E-01,4.865851337340404E-01,5.604772602388374E-01,6.534783834431496E-01,7.728870688769213E-01,9.298722951277453E-01,1.142233814710207E+00,1.439770257172785E+00,1.875507732969076E+00,2.551440427709048E+00,3.685003971257117E+00,5.810999794850599E+00,1.055945793982426E+01};
      i += dnndAmp.test_2to2_amp2([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH4);
      i += dnndAmp.test_2to2_amp2_rotations([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH4);
      i += dnndAmp.test_2to2_amp2_boosts([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH4);
      i += dnndAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dnndAmp.amp2(); }, md,0,0,md,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
