
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

//File:  SPINAS/SM/WWWW2.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/WWWW2.h"

namespace spinas {

  WWWW2::WWWW2(const ldouble& echarge, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& widthZ, const ldouble& sinW):
    e(echarge), mh(massh), wh(widthh), MW(massW), WZ(widthZ), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW) {
    CW = std::sqrt(1.0-sinW*sinW);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    proph = propagator(mh,wh);
    propZ = propagator(MZ,WZ);
    propA = propagator(0,0);
    p1=particle(MW);
    p2=particle(MW);
    p3=particle(MW);
    p4=particle(MW);
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s32s = sproduct(SQUARE,&p3,&p2);
    a32a = sproduct(ANGLE,&p3,&p2);
    s43s = sproduct(SQUARE,&p4,&p3);
    a43a = sproduct(ANGLE,&p4,&p3);
    s42s = sproduct(SQUARE,&p4,&p2);
    a42a = sproduct(ANGLE,&p4,&p2);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //
    s231a = sproduct(SQUARE,&p2,&p3,&p1);
    s143a = sproduct(SQUARE,&p1,&p4,&p3);
    s341a = sproduct(SQUARE,&p3,&p4,&p1);
    s123a = sproduct(SQUARE,&p1,&p2,&p3);
    s142a = sproduct(SQUARE,&p1,&p4,&p2);
    s241a = sproduct(SQUARE,&p2,&p4,&p1);
    s423a = sproduct(SQUARE,&p4,&p2,&p3);
    s234a = sproduct(SQUARE,&p2,&p3,&p4);
    //
    s214a = sproduct(SQUARE,&p2,&p1,&p4);
    s412a = sproduct(SQUARE,&p4,&p1,&p2);
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
    //Couplings
    preh = e*e/(MW*MW*SW*SW);
    preA = e*e/(2.0*MW*MW*MW*MW);
    preZ = e*e/(2.0*MW*MW*MZ*MZ*SW*SW);
    pre4 = e*e/(MW*MW*MW*MW*SW*SW);
  }
  void WWWW2::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(MW);
    p2.set_mass(MW);
    p3.set_mass(MW);
    p4.set_mass(MW);
    proph.set_mass(mh);
    propZ.set_mass(MZ);
    //Couplings
    preh = e*e/(MW*MW*SW*SW);
    preA = e*e/(2.0*MW*MW*MW*MW);
    preZ = e*e/(2.0*MW*MW*MZ*MZ*SW*SW);
    pre4 = e*e/(MW*MW*MW*MW*SW*SW);
  }
  void WWWW2::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    s14s.update();
    a14a.update();
    s32s.update();
    a32a.update();
    s43s.update();
    a43a.update();
    s42s.update();
    a42a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    //
    s231a.update();
    s143a.update();
    s341a.update();
    s123a.update();
    s142a.update();
    s241a.update();
    s423a.update();
    s234a.update();
    //
    s214a.update();
    s412a.update();
    s314a.update();
    s413a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenTh=proph.denominator(propTP);
    pDenSh=proph.denominator(propSP);
    pDenTZ=propZ.denominator(propTP);
    pDenSZ=propZ.denominator(propSP);
    pDenTA=propA.denominator(propTP);
    pDenSA=propA.denominator(propSP);
    //Mandelstahms
    m14 = -2.*p1.dot(p4)+2.*MW*MW;//(p1-p4)^2
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble WWWW2::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1, two=2, three=3, four=4, six=6;
    int ds1a, ds1b, ds2a, ds2b, ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds1,ds2,ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds1,ds2,ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds1,ds1a,ds1b, ds2,ds2a,ds2b, ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);

      
      //T-Channel h
      //preh = e*e/(MW*MW*SW*SW);
      //W+W+W-W- all in:
      //preh [13]<13>[24]<24>/(t-Mh^2)
      //W+W-W-W+ all in: 2<->4
      //preh [13]<13>[24]<24>/(t-Mh^2)
      //34 out:
      //preh [13]<13>[24]<24>/(t-Mh^2)
      amplitude += normFactor*preh*s13s.v(ds1a,ds3a)*a13a.v(ds1b,ds3b)*s42s.v(ds4a,ds2a)*a42a.v(ds4b,ds2b)/pDenTh;
      
      //S-Channel h
      //preh = e*e/(MW*MW*SW*SW);
      //W+W+W-W- all in:
      //preh [14]<14>[23]<23>/(u-Mh^2)
      //W+W-W-W+ all in: 2<->4
      //preh [12]<12>[43]<43>/(s-Mh^2)
      //34 out:
      //preh [12]<12>[43]<43>/(s-Mh^2)
      amplitude += normFactor*preh*s12s.v(ds1a,ds2a)*a12a.v(ds1b,ds2b)*s43s.v(ds4a,ds3a)*a43a.v(ds4b,ds3b)/pDenSh;

      //4-Point
      //pre4 = e*e/(MW*MW*MW*MW*SW*SW);
      //W+W+W-W- all in:
      //- pre4 [12]<12>[34]<34>
      //W+W-W-W+ all in: 2<->4
      //- pre4 [14]<14>[32]<32>
      //34 out:
      //- pre4 [14]<14>[32]<32>
      amplitude += - normFactor*pre4*s14s.v(ds1a,ds4a)*a14a.v(ds1b,ds4b)*s32s.v(ds3a,ds2a)*a32a.v(ds3b,ds2b);

      //T-Channel A
      //preA = e*e/(2*MW*MW*MW*MW);
      //W+W+W-W- all in:
      //((3(<34>[12]-<23>[14]-<14>[23]+<12>[34])MW[13]+((<34>[12]-4<23>[14]-(<14>+3[14])[23]+4<12>[34])MW+3[23][124>-3[13][432>)<13>)MW[24]+(4MW^2 <34>[12][13]-((4<23>[14]+(<14>+3[14])[23]-<12>[34])MW-3[24][123>+3[23][431>)MW[13]+(-3MW^2 <23>[14]+3MW^2 <14>[23]-6MW^2 [14][23]+4MW^2 [13][24]-2s12[13][24]+3MW[24][123>+3MW[23][124>-3MW[23][431>-3MW[13][432>)<13>)<24>)
      //W+W-W-W+: 2<->4:
      //(-2s14<13><42>[13][42]+((3(-1<43>[12]+<32>[14]+<14>[32]-1<12>[43])[42]+(-4<43>[12]+4<32>[14]+<14>[32]-1(<12>+3[12])[43])<42>)[13]+((-4<43>[12]+<32>[14]+4<14>[32]-1(<12>+3[12])[43])[42]+(4[13][42]+3<12>[43]-3(<43>+2[43])[12])<42>)<13>)MW^2 -3((-1[42][143>+[43][231>)<42>[13]+((-1[43][142>+[13][234>)[42]+(-1[42][143>+(-1[142>+[231>)[43]+[13][234>)<42>)<13>)MW)
      //34 out:
      //(-2s14<13><42>[13][42]+((-3(<43>[12]+<32>[14]+<14>[32]+<12>[43])[42]+(4<43>[12]+4<32>[14]+<14>[32]+(<12>+3[12])[43])<42>)[13]+((4<43>[12]+<32>[14]+4<14>[32]+(<12>+3[12])[43])[42]+(4[13][42]+3<12>[43]-3(<43>+2[43])[12])<42>)<13>)MW^2 -3(([42][143>+[43][231>)<42>[13]+(-1([43][142>+[13][234>)[42]+(-1[42][143>+([142>-1[231>)[43]+[13][234>)<42>)<13>)MW)
      amplitude += normFactor*preA*(-two*m14*a13a.v(ds1a,ds3a)*a42a.v(ds4a,ds2a)*s13s.v(ds1b,ds3b)*s42s.v(ds4b,ds2b)+((four*a42a.v(ds4a,ds2b)-three*s42s.v(ds4a,ds2b))*a32a.v(ds3a,ds2a)*s13s.v(ds1a,ds3b)*s14s.v(ds1b,ds4b)+(-three*(a14a.v(ds1a,ds4b)*s32s.v(ds3a,ds2a)*s42s.v(ds4a,ds2b)+a43a.v(ds4a,ds3a)*s12s.v(ds1a,ds2b)*s42s.v(ds4b,ds2a))+(four*a43a.v(ds4b,ds3a)*s12s.v(ds1a,ds2b)+a14a.v(ds1a,ds4b)*s32s.v(ds3a,ds2b)+four*a13a.v(ds1a,ds3a)*s42s.v(ds4b,ds2b))*a42a.v(ds4a,ds2a))*s13s.v(ds1b,ds3b)+(three*a42a.v(ds4a,ds2a)*s12s.v(ds1a,ds2b)*s13s.v(ds1b,ds3b)+((three*a13a.v(ds1b,ds3b)+s13s.v(ds1b,ds3b))*a42a.v(ds4a,ds2a)+(a13a.v(ds1b,ds3b)-three*s13s.v(ds1b,ds3b))*s42s.v(ds4a,ds2a))*a12a.v(ds1a,ds2b))*s43s.v(ds4b,ds3a)+((a32a.v(ds3a,ds2a)*s14s.v(ds1b,ds4b)+four*a14a.v(ds1b,ds4b)*s32s.v(ds3a,ds2a))*s42s.v(ds4a,ds2b)+(four*a43a.v(ds4a,ds3a)*s42s.v(ds4b,ds2a)+three*s42s.v(ds4a,ds2a)*s43s.v(ds4b,ds3a)-three*(a43a.v(ds4b,ds3a)+two*s43s.v(ds4b,ds3a))*a42a.v(ds4a,ds2a))*s12s.v(ds1b,ds2b))*a13a.v(ds1a,ds3b))*MW*MW-three*(((a42a.v(ds4a,ds2a)-s42s.v(ds4a,ds2a))*s142a.v(ds1b,ds2b)-a42a.v(ds4a,ds2a)*s231a.v(ds2b,ds1b))*a13a.v(ds1a,ds3b)*s43s.v(ds4b,ds3a)+((-a13a.v(ds1a,ds3a)+s13s.v(ds1a,ds3a))*s42s.v(ds4b,ds2b)*s143a.v(ds1b,ds3b)+s13s.v(ds1a,ds3b)*s43s.v(ds4b,ds3a)*s231a.v(ds2b,ds1b))*a42a.v(ds4a,ds2a)+(a42a.v(ds4a,ds2a)-s42s.v(ds4a,ds2a))*a13a.v(ds1a,ds3a)*s13s.v(ds1b,ds3b)*s234a.v(ds2b,ds4b))*MW)/pDenTA;

      //U-Channel A
      //preA = e*e/(2*MW*MW*MW*MW);
      //W+W+W-W- all in:
      //((3(<34>[12]-<23>[14]-<14>[23]+<12>[34])MW[13]+((<34>[12]-4<23>[14]-(<14>+3[14])[23]+4<12>[34])MW+3[23][124>-3[13][432>)<13>)MW[24]+(4MW^2 <34>[12][13]-((4<23>[14]+(<14>+3[14])[23]-<12>[34])MW-3[24][123>+3[23][431>)MW[13]+(-3MW^2 <23>[14]+3MW^2 <14>[23]-6MW^2 [14][23]+4MW^2 [13][24]-2s12[13][24]+3MW[24][123>+3MW[23][124>-3MW[23][431>-3MW[13][432>)<13>)<24>)
      amplitude += normFactor*preA*(-two*m14*a12a.v(ds1a,ds2a)*a43a.v(ds4a,ds3a)*s12s.v(ds1b,ds2b)*s43s.v(ds4b,ds3b)+((four*a32a.v(ds3a,ds2a)*a43a.v(ds4a,ds3b)*s14s.v(ds1b,ds4b)+(a43a.v(ds4b,ds3a)-three*s43s.v(ds4b,ds3a))*a42a.v(ds4a,ds2a)*s13s.v(ds1b,ds3b))*s12s.v(ds1a,ds2b)+(four*a43a.v(ds4a,ds3a)*s42s.v(ds4b,ds2a)+three*s42s.v(ds4a,ds2a)*s43s.v(ds4b,ds3a)-three*(a43a.v(ds4b,ds3a)+two*s43s.v(ds4b,ds3a))*a42a.v(ds4a,ds2a))*a13a.v(ds1a,ds3b)*s12s.v(ds1b,ds2b)+(four*a14a.v(ds1b,ds4b)*s32s.v(ds3a,ds2a)*s43s.v(ds4a,ds3b)+(three*a43a.v(ds4a,ds3b)+s43s.v(ds4a,ds3b))*a32a.v(ds3a,ds2a)*s14s.v(ds1b,ds4b)+three*(a14a.v(ds1b,ds4b)*s32s.v(ds3b,ds2a)+a13a.v(ds1b,ds3b)*s42s.v(ds4b,ds2a))*a43a.v(ds4a,ds3a)+four*a13a.v(ds1b,ds3b)*s42s.v(ds4a,ds2a)*s43s.v(ds4b,ds3a)+(three*a43a.v(ds4b,ds3a)*s13s.v(ds1b,ds3b)+(-three*a13a.v(ds1b,ds3b)+s13s.v(ds1b,ds3b))*s43s.v(ds4b,ds3a))*a42a.v(ds4a,ds2a))*a12a.v(ds1a,ds2b)+(a14a.v(ds1a,ds4b)*s32s.v(ds3b,ds2a)+four*a12a.v(ds1a,ds2a)*s43s.v(ds4b,ds3b))*a43a.v(ds4a,ds3a)*s12s.v(ds1b,ds2b))*MW*MW+three*(-(a12a.v(ds1a,ds2a)+s12s.v(ds1a,ds2a))*a43a.v(ds4a,ds3a)*s43s.v(ds4b,ds3b)*s241a.v(ds2b,ds1b)+(-(a12a.v(ds1a,ds2b)+s12s.v(ds1a,ds2b))*s43s.v(ds4b,ds3a)*s123a.v(ds1b,ds3b)+(a43a.v(ds4b,ds3a)+s43s.v(ds4b,ds3a))*s12s.v(ds1a,ds2b)*s341a.v(ds3b,ds1b))*a42a.v(ds4a,ds2a)+(a43a.v(ds4a,ds3a)+s43s.v(ds4a,ds3a))*a12a.v(ds1a,ds2a)*s12s.v(ds1b,ds2b)*s423a.v(ds4b,ds3b))*MW)/pDenSA;
      
      //T-Channel Z
      //preZ = e*e/(2.0*MW*MW*MZ*MZ*SW*SW);
      //W+W+W-W- all in:
      //- preZ 3MW( (-[24][123>+[23][431>)<24>[13] + ( (-[23][124>+[13][432>)[24] + (-[24][123>+(-[124>+[431>)[23]+[13][432>)<24> )<13> )/(t-MZ^2)
      //- preZ MZ^2(<13><24>[14][23]+<14><23>[13][24]+(((-<34>[12]+4<23>[14]+(<14>+3[14])[23]-4<12>[34])<13>+3(-<34>[12]+<23>[14]+<14>[23]-<12>[34])[13])[24]+(-4<34>[12][13]+(-3<14>[23]+3(<23>+2[23])[14]-4[13][24])<13>+(4<23>[14]+(<14>+3[14])[23]-<12>[34])[13])<24>)CW^2 )/(t-MZ^2)
      //- preZ 2s12<13><24>[13][24]/(t-MZ^2)
      amplitude += - normFactor*preZ*(((-((a32a.v(ds3a,ds2a)*s14s.v(ds1b,ds4b)+four*a14a.v(ds1b,ds4b)*s32s.v(ds3a,ds2a))*s42s.v(ds4a,ds2b)+(four*a43a.v(ds4a,ds3a)*s42s.v(ds4b,ds2a)+three*s42s.v(ds4a,ds2a)*s43s.v(ds4b,ds3a))*s12s.v(ds1b,ds2b))*CW*MZ-three*s42s.v(ds4a,ds2a)*s43s.v(ds4b,ds3a)*s142a.v(ds1b,ds2b))*CW+(three*CW*CW*MZ*a43a.v(ds4b,ds3a)*s12s.v(ds1b,ds2b)+((MZ+6*CW*CW*MZ)*s12s.v(ds1b,ds2b)+three*(s142a.v(ds1b,ds2b)-s231a.v(ds2b,ds1b))*CW)*s43s.v(ds4b,ds3a))*a42a.v(ds4a,ds2a))*MZ*a13a.v(ds1a,ds3b)-(four*CW*CW*MZ*MZ*a43a.v(ds4b,ds3a)*s12s.v(ds1a,ds2b)*s13s.v(ds1b,ds3b)+CW*CW*MZ*MZ*a14a.v(ds1a,ds4b)*s13s.v(ds1b,ds3b)*s32s.v(ds3a,ds2b)+four*CW*CW*MZ*MZ*a13a.v(ds1a,ds3a)*s13s.v(ds1b,ds3b)*s42s.v(ds4b,ds2b)-two*m14*a13a.v(ds1a,ds3a)*s13s.v(ds1b,ds3b)*s42s.v(ds4b,ds2b)+three*CW*CW*MZ*MZ*a12a.v(ds1a,ds2b)*a13a.v(ds1b,ds3b)*s43s.v(ds4b,ds3a)+CW*CW*MZ*MZ*a12a.v(ds1a,ds2b)*s13s.v(ds1b,ds3b)*s43s.v(ds4b,ds3a)+three*CW*CW*MZ*MZ*s12s.v(ds1a,ds2b)*s13s.v(ds1b,ds3b)*s43s.v(ds4b,ds3a)+three*CW*MZ*a13a.v(ds1a,ds3a)*s42s.v(ds4b,ds2b)*s143a.v(ds1b,ds3b)-three*CW*MZ*s13s.v(ds1a,ds3a)*s42s.v(ds4b,ds2b)*s143a.v(ds1b,ds3b)-three*CW*MZ*s13s.v(ds1a,ds3b)*s43s.v(ds4b,ds3a)*s231a.v(ds2b,ds1b)-three*CW*MZ*a13a.v(ds1a,ds3a)*s13s.v(ds1b,ds3b)*s234a.v(ds2b,ds4b))*a42a.v(ds4a,ds2a)+((-four*a42a.v(ds4a,ds2b)+three*s42s.v(ds4a,ds2b))*CW*CW*MZ*a32a.v(ds3a,ds2a)*s13s.v(ds1a,ds3b)*s14s.v(ds1b,ds4b)+(three*CW*CW*a14a.v(ds1a,ds4b)*s32s.v(ds3a,ds2a)*s42s.v(ds4a,ds2b)+(a12a.v(ds1a,ds2b)+three*CW*CW*s12s.v(ds1a,ds2b))*a43a.v(ds4a,ds3a)*s42s.v(ds4b,ds2a))*MZ*s13s.v(ds1b,ds3b)-(a13a.v(ds1b,ds3b)-three*s13s.v(ds1b,ds3b))*CW*CW*MZ*a12a.v(ds1a,ds2b)*s42s.v(ds4a,ds2a)*s43s.v(ds4b,ds3a)-three*CW*a13a.v(ds1a,ds3a)*s13s.v(ds1b,ds3b)*s42s.v(ds4a,ds2a)*s234a.v(ds2b,ds4b))*MZ)/pDenTZ;

      //U-Channel Z
      //preZ = e*e/(2.0*MW*MW*MZ*MZ*SW*SW);
      //W+W+W-W- all in:
      //(2s12<13><24>[13][24]+(<13><24>[14][23]+<14><23>[13][24]+(((-<34>[12]+4<23>[14]+(<14>+3[14])[23]-4<12>[34])<13>+3(-<34>[12]+<23>[14]+<14>[23]-<12>[34])[13])[24]+(-4<34>[12][13]+(-3<14>[23]+3(<23>+2[23])[14]-4[13][24])<13>+(4<23>[14]+(<14>+3[14])[23]-<12>[34])[13])<24>)CW^2 )MZ^2 +3((-[24][123>+[23][431>)<24>[13]+((-[23][124>+[13][432>)[24]+(-[24][123>+(-[124>+[431>)[23]+[13][432>)<24>)<13>)CWMZ)
      amplitude += - normFactor*preZ*(-CW*CW*MZ*MZ*a42a.v(ds4a,ds2a)*a43a.v(ds4b,ds3a)*s12s.v(ds1a,ds2b)*s13s.v(ds1b,ds3b)-four*CW*CW*MZ*MZ*a32a.v(ds3a,ds2a)*a43a.v(ds4a,ds3b)*s12s.v(ds1a,ds2b)*s14s.v(ds1b,ds4b)-CW*CW*MZ*MZ*a14a.v(ds1a,ds4b)*a43a.v(ds4a,ds3a)*s12s.v(ds1b,ds2b)*s32s.v(ds3b,ds2a)+three*CW*CW*MZ*MZ*a42a.v(ds4a,ds2a)*s12s.v(ds1a,ds2b)*s13s.v(ds1b,ds3b)*s43s.v(ds4b,ds3a)+((three*CW*CW*a43a.v(ds4b,ds3a)+(1+6*CW*CW)*s43s.v(ds4b,ds3a))*a42a.v(ds4a,ds2a)-(four*a43a.v(ds4a,ds3a)*s42s.v(ds4b,ds2a)+three*s42s.v(ds4a,ds2a)*s43s.v(ds4b,ds3a))*CW*CW)*MZ*MZ*a13a.v(ds1a,ds3b)*s12s.v(ds1b,ds2b)-four*CW*CW*MZ*MZ*a12a.v(ds1a,ds2a)*a43a.v(ds4a,ds3a)*s12s.v(ds1b,ds2b)*s43s.v(ds4b,ds3b)+two*m14*a12a.v(ds1a,ds2a)*a43a.v(ds4a,ds3a)*s12s.v(ds1b,ds2b)*s43s.v(ds4b,ds3b)-(((three*a43a.v(ds4a,ds3b)+s43s.v(ds4a,ds3b))*CW*CW*a32a.v(ds3a,ds2a)*s14s.v(ds1b,ds4b)+(three*a43a.v(ds4a,ds3a)*s32s.v(ds3b,ds2a)+four*s32s.v(ds3a,ds2a)*s43s.v(ds4a,ds3b))*CW*CW*a14a.v(ds1b,ds4b)+(three*CW*CW*a13a.v(ds1b,ds3b)-s13s.v(ds1b,ds3b))*a43a.v(ds4a,ds3a)*s42s.v(ds4b,ds2a)+four*CW*CW*a13a.v(ds1b,ds3b)*s42s.v(ds4a,ds2a)*s43s.v(ds4b,ds3a))*MZ+(three*CW*MZ*a43a.v(ds4b,ds3a)*s13s.v(ds1b,ds3b)+((-three*a13a.v(ds1b,ds3b)+s13s.v(ds1b,ds3b))*CW*MZ-three*s123a.v(ds1b,ds3b))*s43s.v(ds4b,ds3a))*CW*a42a.v(ds4a,ds2a))*MZ*a12a.v(ds1a,ds2b)+three*CW*MZ*a42a.v(ds4a,ds2a)*s12s.v(ds1a,ds2b)*s43s.v(ds4b,ds3a)*s123a.v(ds1b,ds3b)+three*CW*MZ*a12a.v(ds1a,ds2a)*a43a.v(ds4a,ds3a)*s43s.v(ds4b,ds3b)*s241a.v(ds2b,ds1b)+three*CW*MZ*a43a.v(ds4a,ds3a)*s12s.v(ds1a,ds2a)*s43s.v(ds4b,ds3b)*s241a.v(ds2b,ds1b)-three*CW*MZ*a42a.v(ds4a,ds2a)*a43a.v(ds4b,ds3a)*s12s.v(ds1a,ds2b)*s341a.v(ds3b,ds1b)-three*CW*MZ*a42a.v(ds4a,ds2a)*s12s.v(ds1a,ds2b)*s43s.v(ds4b,ds3a)*s341a.v(ds3b,ds1b)-three*(a43a.v(ds4a,ds3a)+s43s.v(ds4a,ds3a))*CW*MZ*a12a.v(ds1a,ds2a)*s12s.v(ds1b,ds2b)*s423a.v(ds4b,ds3b))/pDenSZ;

      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble WWWW2::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over spins 1/3^2=1/9
    return amp2/9.0;
  }
  



  //  Tests
  int test_WWWW2(){
    int n=0;//Number of fails
    std::cout<<"\t* W+, W- -> W+, W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      WWWW2 WWWW4Amp = WWWW2(EE,mh,0,MW,0,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.563378238252062E+02,4.503988221084452E+01,1.755119910189961E+01,8.915661477818531E+00,5.219004741877505E+00,3.347470178235117E+00,2.293365715406586E+00,1.654170970525979E+00,1.245114790363838E+00,9.725084477578637E-01,7.851501906466927E-01,6.533715872407764E-01,5.591556154018003E-01,4.911325790096854E-01,4.418948367984821E-01,4.064846857352055E-01,3.815074735996752E-01,3.645929661079332E-01,3.540581967811323E-01,3.486907941525331E-01};
      i += WWWW4Amp.test_2to2_amp2([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH);
      i += WWWW4Amp.test_2to2_amp2_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH);
      i += WWWW4Amp.test_2to2_amp2_boosts([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH);
      i += WWWW4Amp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {1.588451390262576E+02,2.956881807562652E+01,1.320101564037441E+01,7.430827347748311E+00,4.688393385177053E+00,3.178167545875126E+00,2.268897377151777E+00,1.687447549737127E+00,1.299160083131472E+00,1.031353327605228E+00,8.420790234037233E-01,7.058402344789360E-01,6.064861109381752E-01,5.334351737305003E-01,4.795577370215638E-01,4.399350677083765E-01,4.111048037639216E-01,3.905873024383555E-01,3.765800690020982E-01,3.677558448235396E-01};
      i += WWWW4Amp.test_2to2_amp2([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH2);
      i += WWWW4Amp.test_2to2_amp2_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH2);
      i += WWWW4Amp.test_2to2_amp2_boosts([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH2);
      i += WWWW4Amp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {6.443442341732017E+08,7.160078876565878E+07,2.577879855938749E+07,1.315373121334773E+07,7.957971546895633E+06,5.327756706758454E+06,3.814919692966164E+06,2.865707978127272E+06,2.231305060274738E+06,1.786454377100097E+06,1.462523578728815E+06,1.219349275383327E+06,1.032157816812508E+06,8.849952053622548E+05,7.672109087766274E+05,6.714747598580248E+05,5.926080023286534E+05,5.268677430348845E+05,4.714944165585224E+05,4.244173094596451E+05};
      i += WWWW4Amp.test_2to2_amp2([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH4);
      i += WWWW4Amp.test_2to2_amp2_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH4);
      i += WWWW4Amp.test_2to2_amp2_boosts([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH4);
      i += WWWW4Amp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH4);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=250\n";
      mh=1;
      WWWW4Amp.set_masses(mh,MW);
      pspatial=250;
      ldouble dataCH5[20] = {2.649629077745788E+02,4.568845194160681E+01,1.773702479829053E+01,9.004912884984970E+00,5.274010269199242E+00,3.385702051164875E+00,2.321551713785867E+00,1.675578477994178E+00,1.261627685978568E+00,9.853545309854732E-01,7.952005333177272E-01,6.612812405954435E-01,5.654375427577201E-01,4.962002778497900E-01,4.460888480502931E-01,4.100901198165482E-01,3.847661742309341E-01,3.677122840233300E-01,3.572175152007832E-01,3.520465032086718E-01};
      i += WWWW4Amp.test_2to2_amp2([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH5);
      i += WWWW4Amp.test_2to2_amp2_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH5);
      i += WWWW4Amp.test_2to2_amp2_boosts([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH5);
      i += WWWW4Amp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH5);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH6[20] = {2.012835480187636E+02,3.437661937927619E+01,1.481323817451041E+01,8.181263748406998E+00,5.102345252914466E+00,3.432530890964385E+00,2.437621256119159E+00,1.806061792876518E+00,1.386516153826833E+00,1.098239060461784E+00,8.950521919768725E-01,7.490914810670239E-01,6.428095057929626E-01,5.647634071346190E-01,5.072715849785221E-01,4.650497255931578E-01,4.343860852123939E-01,4.126263695769401E-01,3.978428281586990E-01,3.886164097259885E-01};
      i += WWWW4Amp.test_2to2_amp2([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH6);
      i += WWWW4Amp.test_2to2_amp2_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH6);
      i += WWWW4Amp.test_2to2_amp2_boosts([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH6);
      i += WWWW4Amp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH6);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH7[20] = {7.812415807475721E+08,1.130735360168684E+08,4.843576207924464E+07,2.795935761682009E+07,1.855039308823812E+07,1.334220116129114E+07,1.011732765659646E+07,7.964966229538589E+06,6.448986787100058E+06,5.337002605081021E+06,4.495028784536018E+06,3.840983079006357E+06,3.322098488358643E+06,2.903097120439591E+06,2.559613046898231E+06,2.274348001681478E+06,2.034728818464021E+06,1.831428460906906E+06,1.657403322339229E+06,1.507249686540676E+06};
      i += WWWW4Amp.test_2to2_amp2([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH7);
      i += WWWW4Amp.test_2to2_amp2_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH7);
      i += WWWW4Amp.test_2to2_amp2_boosts([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH7);
      i += WWWW4Amp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH7);
      //std::cout<<"\n# mh=125, MW=50, pspatial=125\n";
      mh=125;
      MW=50;
      MZ=MW/CW;
      WWWW4Amp.set_masses(mh,MW);
      pspatial = 125;
      ldouble dataCH9[20] = {2.061911754646830E+02,3.869018431590226E+01,1.587694955102091E+01,8.345356635030772E+00,5.014481712222136E+00,3.291853584542936E+00,2.308158603309467E+00,1.706576971827064E+00,1.319832343574932E+00,1.061696643155913E+00,8.844277458516132E-01,7.600755087629659E-01,6.715285469827464E-01,6.079210900756936E-01,5.621408671671849E-01,5.294124805067258E-01,5.064598048201634E-01,4.909937692258121E-01,4.813898484648840E-01,4.764797604255753E-01};
      i += WWWW4Amp.test_2to2_amp2([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH9);
      i += WWWW4Amp.test_2to2_amp2_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH9);
      i += WWWW4Amp.test_2to2_amp2_boosts([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH9);
      i += WWWW4Amp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH9);
      //std::cout<<"\n# mh=125, MW=50, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH10[20] = {9.655201466940087E+07,1.073090803902375E+07,3.864173224960032E+06,1.972050805855260E+06,1.193291962210881E+06,7.990314791708251E+05,5.722424237149262E+05,4.299339147675121E+05,3.348141613267291E+05,2.681092936542360E+05,2.195320135580685E+05,1.830619688702597E+05,1.549855130689819E+05,1.329109881831349E+05,1.152417071149843E+05,1.008787113026679E+05,8.904555373659628E+04,7.918104086100706E+04,7.087139712766193E+04,6.380612153676905E+04};
      i += WWWW4Amp.test_2to2_amp2([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH10);
      i += WWWW4Amp.test_2to2_amp2_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH10);
      i += WWWW4Amp.test_2to2_amp2_boosts([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH10);
      i += WWWW4Amp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH10);
      //std::cout<<"\n# mh=1, MW=50, pspatial=125\n";
      mh=1;
      MW=50;
      MZ=MW/CW;
      WWWW4Amp.set_masses(mh,MW);
      pspatial = 125;
      ldouble dataCH11[20] = {2.236047446507091E+02,4.066384206325570E+01,1.650753080558468E+01,8.612911771602885E+00,5.134199849052453E+00,3.335419730634330E+00,2.306033615753990E+00,1.674089493900330E+00,1.265740718102132E+00,9.915217871854202E-01,8.019265290611673E-01,6.679610287002321E-01,5.718526301130803E-01,5.022921174368614E-01,4.518554707842158E-01,4.155393842766569E-01,3.898956347709435E-01,3.725016643518827E-01,3.616270752604485E-01,3.560180362264847E-01};
      i += WWWW4Amp.test_2to2_amp2([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH11);
      i += WWWW4Amp.test_2to2_amp2_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH11);
      i += WWWW4Amp.test_2to2_amp2_boosts([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH11);
      i += WWWW4Amp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH11);
      //std::cout<<"\n# mh=1, MW=50, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH12[20] = {1.170515013437128E+08,1.694135184097774E+07,7.257115261590748E+06,4.189318332019399E+06,2.779660757915390E+06,1.999361265535111E+06,1.516199803629081E+06,1.193720315654501E+06,9.665827065234523E+05,7.999714305610713E+05,6.738132677279257E+05,5.758108913754754E+05,4.980588801698995E+05,4.352719646708896E+05,3.837996646391776E+05,3.410503518031692E+05,3.051402973341285E+05,2.746720990827283E+05,2.485904539269491E+05,2.260857453149866E+05};
      i += WWWW4Amp.test_2to2_amp2([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH12);
      i += WWWW4Amp.test_2to2_amp2_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH12);
      i += WWWW4Amp.test_2to2_amp2_boosts([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH12);
      i += WWWW4Amp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWW4Amp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH12);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
  

}
