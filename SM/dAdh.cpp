
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

//File:  SPINAS/SM/dAdh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dAdh.h"

namespace spinas {

  dAdh::dAdh(const ldouble& echarge, const ldouble& massd, const ldouble& massh, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qd(-1.0/3.0), md(massd), mh(massh), prop(massd,0), MW(massW), SW(sinW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(md);
    p2=particle(0);
    p3=particle(md);
    p4=particle(mh);
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s243a = sproduct(SQUARE,&p2,&p4,&p3,2);
    a243s = sproduct(ANGLE,&p2,&p4,&p3,2);
    s241a = sproduct(SQUARE,&p2,&p4,&p1,2);
    a241s = sproduct(ANGLE,&p2,&p4,&p1,2);
    s2342s = sproduct(SQUARE,&p2,&p3,&p4,&p2,2);
    a2342a = sproduct(ANGLE,&p2,&p3,&p4,&p2,2);
    pre = sqrt2*e*e*Qd*md/(2.0*MW*SW);
  }
  void dAdh::set_masses(const ldouble& massd, const ldouble& massh, const ldouble& massW){
    md=massd;
    mh=massh;
    p1.set_mass(md);
    p3.set_mass(md);
    p4.set_mass(mh);
    prop.set_mass(md);
    pre = sqrt2*e*e*Qd*md/(2.0*MW*SW);
  }
  void dAdh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s243a.update();
    a243s.update();
    s241a.update();
    a241s.update();
    s2342s.update();
    a2342a.update();
    //Propagator Momentum
    ldouble propSP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=prop.denominator(propSP);
    pDenU=prop.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble dAdh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    /*mh^2[12][23] - md[12][243> + md[23][241> - <13>[2342]
      mh^2<12><23> - md<12><243] + md<23><241] - [13]<2342>
      Becomes after a sign change due to p3 and p4 being outgoing:*/
    //pre = sqrt2*e*e*Qd*md/(2.0*MW*SW);
    if(ds2>0){
      //mh^2[12][23] - md[12][243> - md[23][241> + <13>[2342]
      return pre*(mh*mh*s12s.v(ds1)*s23s.v(ds3) - md*s12s.v(ds1)*s243a.v(ds3) - md*s23s.v(ds3)*s241a.v(ds1) + a13a.v(ds1,ds3)*s2342s.v())/pDenS/pDenU;
    }
    else if(ds2<0){
      //-mh^2<12><23> + md<12><243] + md<23><241] - [13]<2342>
      return pre*(- mh*mh*a12a.v(ds1)*a23a.v(ds3) + md*a12a.v(ds1)*a243s.v(ds3) + md*a23a.v(ds3)*a241s.v(ds1) - s13s.v(ds1,ds3)*a2342a.v())/pDenS/pDenU;
    }
    return cdouble(0,0);    
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble dAdh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-1;j3<=1;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3
    return amp2/12.0;
  }
  
  //set_momenta(...) must be called before amp2_Aplus().
  //Positive Helicity Photon Only
  ldouble dAdh::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-1;j3<=1;j3+=2){
	M = amp(j1,2,j3);
	amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
      }
    //Average over initial spins 1/2^1
    //Average over initial colors 1/3
    return amp2/6.0;
  }
  



  //  Tests
  int test_dAdh(){
    int n=0;//Number of fails
    std::cout<<"\t* d , A  -> d , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# md=0.0042, mh=125, pspatial=250\n";
      ldouble md=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      dAdh dAdhAmp = dAdh(EE,md,mh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {5.056479582338834E-14,1.058659465975798E-13,2.230590274453605E-13,4.133971165589053E-13,6.910373123081751E-13,1.074042119508829E-12,1.585825896253468E-12,2.257295595221672E-12,3.130108413499155E-12,4.261790956578524E-12,5.734091558569538E-12,7.667228357848796E-12,1.024553583322246E-11,1.376678506576146E-11,1.874530792259529E-11,2.615262015229518E-11,3.807054121467661E-11,5.991180169934364E-11,1.115230803921166E-10,3.715244849553568E-10};
      i += dAdhAmp.test_2to2_amp2([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH);
      i += dAdhAmp.test_2to2_amp2_rotations([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH);
      i += dAdhAmp.test_2to2_amp2_boosts([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH);
      i += dAdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      //(CH produces the same M^2s since it is averaging over spins and the to helicity cases give the same result -- averaging them gaves the same as each individually.)    
      i += dAdhAmp.test_2to2_amp2([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH);
      i += dAdhAmp.test_2to2_amp2_rotations([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH);
      i += dAdhAmp.test_2to2_amp2_boosts([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH);
      i += dAdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH);
      //std::cout<<"\n# md=0.0075, mh=125, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {8.883034196301450E-13,9.783965609921596E-13,1.123267182344991E-12,1.332874825497414E-12,1.619749236593850E-12,1.999876643470674E-12,2.493979934059430E-12,3.129413320422986E-12,3.943045523157202E-12,4.985790489734504E-12,6.329999619165751E-12,8.082072027592136E-12,1.040515302633331E-11,1.356278504116948E-11,1.801017771618410E-11,2.460717148905063E-11,3.519628157976169E-11,5.456804754520583E-11,1.002877809849751E-10,3.304432826984178E-10};
      i += dAdhAmp.test_2to2_amp2([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH2);
      i += dAdhAmp.test_2to2_amp2_rotations([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH2);
      i += dAdhAmp.test_2to2_amp2_boosts([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH2);
      i += dAdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH2);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += dAdhAmp.test_2to2_amp2([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH2);
      i += dAdhAmp.test_2to2_amp2_rotations([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH2);
      i += dAdhAmp.test_2to2_amp2_boosts([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH2);
      i += dAdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH2);
      //std::cout<<"\n# md=125.1, mh=125, pspatial=95\n";
      md = 125.1;
      mh = 125;
      pspatial = 95;
      dAdhAmp.set_masses(md,mh,MW);
      ldouble dataCH4[20] = {1.009102222702257E-03,1.066865565158387E-03,1.124800570042480E-03,1.182858303159486E-03,1.240984938098665E-03,1.299121310609007E-03,1.357202429576976E-03,1.415156939945430E-03,1.472906532357500E-03,1.530365293681145E-03,1.587438991858276E-03,1.644024287714451E-03,1.700007865447129E-03,1.755265472465375E-03,1.809660858062796E-03,1.863044599045378E-03,1.915252798880744E-03,1.966105645153939E-03,2.015405808070952E-03,2.062936660402167E-03};
      i += dAdhAmp.test_2to2_amp2([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH4);
      i += dAdhAmp.test_2to2_amp2_rotations([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH4);
      i += dAdhAmp.test_2to2_amp2_boosts([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH4);
      i += dAdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH4);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += dAdhAmp.test_2to2_amp2([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH4);
      i += dAdhAmp.test_2to2_amp2_rotations([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH4);
      i += dAdhAmp.test_2to2_amp2_boosts([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH4);
      i += dAdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH4);
      //std::cout<<"\n# md=125, mh=0.0005, pspatial=125.1\n";
      md = 125;
      mh = 0.0005;
      pspatial = 125.1;
      dAdhAmp.set_masses(md,mh,MW);
      ldouble dataCH3[20] = {1.016764067166722E-04,3.234676044847782E-04,5.722088405558941E-04,8.511188043730646E-04,1.163922172281650E-03,1.514943050035420E-03,1.909214994916348E-03,2.352607809115786E-03,2.851967701017229E-03,3.415259179417294E-03,4.051678867852506E-03,4.771671576574847E-03,5.586690682000219E-03,6.508344477870120E-03,7.546098521273263E-03,8.701538359584043E-03,9.954128697593717E-03,1.122464994458500E-02,1.227477482425934E-02,1.240059987165247E-02};
      i += dAdhAmp.test_2to2_amp2([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH3);
      i += dAdhAmp.test_2to2_amp2_rotations([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH3);
      i += dAdhAmp.test_2to2_amp2_boosts([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH3);
      i += dAdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAdhAmp.amp2(); }, md,0,md,mh,pspatial,dataCH3);
      //std::cout<<"\nPositive Helicity Photon Only\n";
      i += dAdhAmp.test_2to2_amp2([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH3);
      i += dAdhAmp.test_2to2_amp2_rotations([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH3);
      i += dAdhAmp.test_2to2_amp2_boosts([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH3);
      i += dAdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAdhAmp.amp2_Aplus(); }, md,0,md,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
