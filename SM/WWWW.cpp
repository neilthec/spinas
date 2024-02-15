
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

//File:  SPINAS/SM/WWWW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/WWWW.h"

namespace spinas {

  WWWW::WWWW(const ldouble& echarge, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& widthZ, const ldouble& sinW):
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
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //
    s431a = sproduct(SQUARE,&p4,&p3,&p1);
    s123a = sproduct(SQUARE,&p1,&p2,&p3);
    s321a = sproduct(SQUARE,&p3,&p2,&p1);
    s143a = sproduct(SQUARE,&p1,&p4,&p3);
    s124a = sproduct(SQUARE,&p1,&p2,&p4);
    s421a = sproduct(SQUARE,&p4,&p2,&p1);
    s243a = sproduct(SQUARE,&p2,&p4,&p3);
    s432a = sproduct(SQUARE,&p4,&p3,&p2);
    //
    s412a = sproduct(SQUARE,&p4,&p1,&p2);
    s214a = sproduct(SQUARE,&p2,&p1,&p4);
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
    //Couplings
    preh = e*e/(MW*MW*SW*SW);
    preA = e*e/(2.0*MW*MW*MW*MW);
    preZ = e*e/(2.0*MW*MW*MZ*MZ*SW*SW);
    pre4 = e*e/(MW*MW*MW*MW*SW*SW);
  }
  void WWWW::set_masses(const ldouble& massh, const ldouble& massW){
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
  void WWWW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    s12s.update();
    a12a.update();
    s34s.update();
    a34a.update();
    s23s.update();
    a23a.update();
    s24s.update();
    a24a.update();
    s14s.update();
    a14a.update();
    s13s.update();
    a13a.update();
    //
    s431a.update();
    s123a.update();
    s321a.update();
    s143a.update();
    s124a.update();
    s421a.update();
    s243a.update();
    s432a.update();
    //
    s412a.update();
    s214a.update();
    s312a.update();
    s213a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenTh=proph.denominator(propTP);
    pDenUh=proph.denominator(propUP);
    pDenTZ=propZ.denominator(propTP);
    pDenUZ=propZ.denominator(propUP);
    pDenTA=propA.denominator(propTP);
    pDenUA=propA.denominator(propUP);
    //Mandelstahms
    m12 = +2.*p1.dot(p2)+2.*MW*MW;//(p1+p2)^2
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble WWWW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
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
      //all ingoing = 34 outgoing:
      //preh [13]<13>[24]<24>/(t-Mh^2)
      amplitude += normFactor*preh*s13s.v(ds1a,ds3a)*a13a.v(ds1b,ds3b)*s24s.v(ds2a,ds4a)*a24a.v(ds2b,ds4b)/pDenTh;
      
      //U-Channel h
      //preh = e*e/(MW*MW*SW*SW);
      //all ingoing = 34 outgoing:
      //preh [14]<14>[23]<23>/(u-Mh^2)
      amplitude += normFactor*preh*s14s.v(ds1a,ds4a)*a14a.v(ds1b,ds4b)*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)/pDenUh;

      //4-Point
      //pre4 = e*e/(MW*MW*MW*MW*SW*SW);
      //all ingoing = 34 outgoing:
      //- pre4 [12]<12>[34]<34>
      amplitude += - normFactor*pre4*s12s.v(ds1a,ds2a)*a12a.v(ds1b,ds2b)*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b);

      //T-Channel A
      //preA = e*e/(2*MW*MW*MW*MW);
      //all in:
      //((3(<34>[12]-<23>[14]-<14>[23]+<12>[34])MW[13]+((<34>[12]-4<23>[14]-(<14>+3[14])[23]+4<12>[34])MW+3[23][124>-3[13][432>)<13>)MW[24]+(4MW^2 <34>[12][13]-((4<23>[14]+(<14>+3[14])[23]-<12>[34])MW-3[24][123>+3[23][431>)MW[13]+(-3MW^2 <23>[14]+3MW^2 <14>[23]-6MW^2 [14][23]+4MW^2 [13][24]-2s12[13][24]+3MW[24][123>+3MW[23][124>-3MW[23][431>-3MW[13][432>)<13>)<24>)
      //34 out:
      //((3(<34>[12]+<23>[14]+<14>[23]+<12>[34])MW[13]-((<34>[12]+4<23>[14]+(<14>-3[14])[23]+4<12>[34])MW-3[23][124>+3[13][432>)<13>)MW[24]-(4MW^2 <34>[12][13]+((4<23>[14]+(<14>-3[14])[23]+<12>[34])MW-3[24][123>+3[23][431>)MW[13]+(-3MW^2 <23>[14]+3MW^2 <14>[23]+6MW^2 [14][23]-4MW^2 [13][24]+2s12[13][24]+3MW[24][123>+3MW[23][124>-3MW[23][431>-3MW[13][432>)<13>)<24>)
      amplitude += normFactor*preA*((((-a13a.v(ds1a,ds3b)*s12s.v(ds1b,ds2b)+three*s12s.v(ds1a,ds2b)*s13s.v(ds1b,ds3b))*a34a.v(ds3a,ds4a)*s24s.v(ds2a,ds4b)+(-four*a24a.v(ds2b,ds4a)*s14s.v(ds1b,ds4b)+three*s14s.v(ds1b,ds4a)*s24s.v(ds2b,ds4b))*a23a.v(ds2a,ds3a)*s13s.v(ds1a,ds3b)+(-a24a.v(ds2a,ds4b)*s23s.v(ds2b,ds3b)+three*s23s.v(ds2a,ds3b)*s24s.v(ds2b,ds4b))*a14a.v(ds1a,ds4a)*s13s.v(ds1b,ds3a)-four*a12a.v(ds1a,ds2b)*a13a.v(ds1b,ds3b)*s24s.v(ds2a,ds4a)*s34s.v(ds3a,ds4b)+three*a12a.v(ds1a,ds2b)*s13s.v(ds1b,ds3b)*s24s.v(ds2a,ds4a)*s34s.v(ds3a,ds4b))*MW+((three*a23a.v(ds2a,ds3b)*a24a.v(ds2b,ds4a)*s14s.v(ds1b,ds4b)-three*a14a.v(ds1b,ds4a)*a24a.v(ds2a,ds4b)*s23s.v(ds2b,ds3b)-(four*a23a.v(ds2a,ds3b)*s14s.v(ds1b,ds4a)+(a14a.v(ds1b,ds4a)-three*s14s.v(ds1b,ds4a))*s23s.v(ds2a,ds3b))*s24s.v(ds2b,ds4b))*MW+three*s23s.v(ds2a,ds3b)*s24s.v(ds2b,ds4a)*s124a.v(ds1b,ds4b)-three*s13s.v(ds1b,ds3b)*s24s.v(ds2a,ds4a)*s432a.v(ds4b,ds2b))*a13a.v(ds1a,ds3a))*MW-(four*MW*MW*a34a.v(ds3a,ds4b)*s12s.v(ds1a,ds2b)*s13s.v(ds1b,ds3b)+(MW*a12a.v(ds1a,ds2b)*s13s.v(ds1b,ds3b)*s34s.v(ds3a,ds4b)-three*(s24s.v(ds2b,ds4b)*s123a.v(ds1b,ds3b)+(MW*s14s.v(ds1b,ds4b)-s431a.v(ds4b,ds1b))*s23s.v(ds2b,ds3b))*s13s.v(ds1a,ds3a))*MW+(6*MW*MW*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)+three*(s24s.v(ds2b,ds4b)*s123a.v(ds1b,ds3b)+(s124a.v(ds1b,ds4b)-s431a.v(ds4b,ds1b))*s23s.v(ds2b,ds3b))*MW+(two*(-two*MW*MW+m12)*s24s.v(ds2b,ds4b)-three*MW*s432a.v(ds4b,ds2b))*s13s.v(ds1b,ds3b))*a13a.v(ds1a,ds3a))*a24a.v(ds2a,ds4a))/pDenTA;

      //U-Channel A
      //preA = e*e/(2*MW*MW*MW*MW);
      amplitude += normFactor*preA*((-MW*a23a.v(ds2a,ds3a)*a24a.v(ds2b,ds4a)*s13s.v(ds1a,ds3b)*s14s.v(ds1b,ds4b)-three*MW*a24a.v(ds2a,ds4a)*s13s.v(ds1a,ds3a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)+(three*(a14a.v(ds1b,ds4a)*a24a.v(ds2a,ds4b)-two*a24a.v(ds2a,ds4a)*s14s.v(ds1b,ds4b))*s23s.v(ds2b,ds3b)+(-four*a14a.v(ds1b,ds4a)+three*s14s.v(ds1b,ds4a))*s23s.v(ds2a,ds3b)*s24s.v(ds2b,ds4b)+(three*a24a.v(ds2b,ds4a)*s14s.v(ds1b,ds4b)+(three*a14a.v(ds1b,ds4a)-four*s14s.v(ds1b,ds4a))*s24s.v(ds2b,ds4b))*a23a.v(ds2a,ds3b))*MW*a13a.v(ds1a,ds3a)+four*MW*a12a.v(ds1a,ds2b)*a14a.v(ds1b,ds4a)*s23s.v(ds2a,ds3b)*s34s.v(ds3a,ds4b)+(four*a34a.v(ds3a,ds4a)*s12s.v(ds1a,ds2b)*s14s.v(ds1b,ds4b)+(-three*a14a.v(ds1b,ds4a)+s14s.v(ds1b,ds4a))*a12a.v(ds1a,ds2b)*s34s.v(ds3a,ds4b))*MW*a23a.v(ds2a,ds3b)+three*a24a.v(ds2a,ds4a)*s14s.v(ds1a,ds4b)*s23s.v(ds2b,ds3a)*s143a.v(ds1b,ds3b)+three*a23a.v(ds2a,ds3a)*a24a.v(ds2b,ds4a)*s14s.v(ds1a,ds4b)*s321a.v(ds3b,ds1b)-three*a24a.v(ds2a,ds4a)*s14s.v(ds1a,ds4b)*s23s.v(ds2b,ds3a)*s321a.v(ds3b,ds1b)-three*a23a.v(ds2a,ds3a)*s14s.v(ds1a,ds4a)*s23s.v(ds2b,ds3b)*s421a.v(ds4b,ds1b))*MW+(-three*MW*MW*a23a.v(ds2a,ds3b)*a34a.v(ds3a,ds4b)*s12s.v(ds1b,ds2b)+(MW*a34a.v(ds3a,ds4b)*s12s.v(ds1b,ds2b)*s23s.v(ds2a,ds3b)-(MW*s13s.v(ds1b,ds3a)*s23s.v(ds2b,ds3b)+three*s23s.v(ds2b,ds3a)*s143a.v(ds1b,ds3b))*a24a.v(ds2a,ds4b)+three*s14s.v(ds1b,ds4b)*s23s.v(ds2a,ds3a)*s243a.v(ds2b,ds3b))*MW+(three*MW*MW*a24a.v(ds2b,ds4b)*s13s.v(ds1b,ds3b)+((four*MW*MW-two*m12)*s23s.v(ds2b,ds3b)-three*MW*s243a.v(ds2b,ds3b))*s14s.v(ds1b,ds4b)+three*MW*s23s.v(ds2b,ds3b)*s421a.v(ds4b,ds1b))*a23a.v(ds2a,ds3a))*a14a.v(ds1a,ds4a))/pDenUA;
      
      //T-Channel Z
      //preZ = e*e/(2.0*MW*MW*MZ*MZ*SW*SW);
      //all in:
      //- preZ 3MW( (-[24][123>+[23][431>)<24>[13] + ( (-[23][124>+[13][432>)[24] + (-[24][123>+(-[124>+[431>)[23]+[13][432>)<24> )<13> )/(t-MZ^2)
      //- preZ MZ^2(<13><24>[14][23]+<14><23>[13][24]+(((-<34>[12]+4<23>[14]+(<14>+3[14])[23]-4<12>[34])<13>+3(-<34>[12]+<23>[14]+<14>[23]-<12>[34])[13])[24]+(-4<34>[12][13]+(-3<14>[23]+3(<23>+2[23])[14]-4[13][24])<13>+(4<23>[14]+(<14>+3[14])[23]-<12>[34])[13])<24>)CW^2 )/(t-MZ^2)
      //- preZ 2s12<13><24>[13][24]/(t-MZ^2)
      //34 out:
      //+ preZ 3MW(([24][123>-[23][431>)<24>[13]+(([23][124>-[13][432>)[24]+(-[24][123>+(-[124>+[431>)[23]+[13][432>)<24>)<13>)/(t-MZ^2)
      //- preZ MZ^2(<13><24>[14][23]+<14><23>[13][24]+((4<34>[12][13]+(-3<23>[14]+3(<14>+2[14])[23]-4[13][24])<13>+(4<23>[14]+(<14>-3[14])[23]+<12>[34])[13])<24>+(-3(<34>[12]+<23>[14]+<14>[23]+<12>[34])[13]+(<34>[12]+4<23>[14]+(<14>-3[14])[23]+4<12>[34])<13>)[24])CW^2 )/(t-MZ^2)
      //- preZ 2s12<13><24>[13][24]/(t-MZ^2)
      amplitude += - normFactor*preZ*(two*m12*a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s13s.v(ds1b,ds3b)*s24s.v(ds2b,ds4b)+(a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)+a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s13s.v(ds1b,ds3b)*s24s.v(ds2b,ds4b)+((-three*a13a.v(ds1a,ds3a)*a23a.v(ds2a,ds3b)+four*a23a.v(ds2a,ds3a)*s13s.v(ds1a,ds3b))*a24a.v(ds2b,ds4a)*s14s.v(ds1b,ds4b)+a13a.v(ds1a,ds3b)*a34a.v(ds3a,ds4a)*s12s.v(ds1b,ds2b)*s24s.v(ds2a,ds4b)-three*a34a.v(ds3a,ds4a)*s12s.v(ds1a,ds2b)*s13s.v(ds1b,ds3b)*s24s.v(ds2a,ds4b)+(three*a13a.v(ds1a,ds3a)*a14a.v(ds1b,ds4a)+a14a.v(ds1a,ds4a)*s13s.v(ds1b,ds3a))*a24a.v(ds2a,ds4b)*s23s.v(ds2b,ds3b)+((four*a13a.v(ds1a,ds3a)*a23a.v(ds2a,ds3b)-three*a23a.v(ds2a,ds3a)*s13s.v(ds1a,ds3b))*s14s.v(ds1b,ds4a)+(-three*a14a.v(ds1a,ds4a)*s13s.v(ds1b,ds3a)+(a14a.v(ds1b,ds4a)-three*s14s.v(ds1b,ds4a))*a13a.v(ds1a,ds3a))*s23s.v(ds2a,ds3b))*s24s.v(ds2b,ds4b)+(four*a13a.v(ds1b,ds3b)-three*s13s.v(ds1b,ds3b))*a12a.v(ds1a,ds2b)*s24s.v(ds2a,ds4a)*s34s.v(ds3a,ds4b)+(three*(two*a13a.v(ds1a,ds3a)-s13s.v(ds1a,ds3a))*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)+(four*a34a.v(ds3a,ds4b)*s12s.v(ds1a,ds2b)-four*a13a.v(ds1a,ds3a)*s24s.v(ds2b,ds4b)+a12a.v(ds1a,ds2b)*s34s.v(ds3a,ds4b))*s13s.v(ds1b,ds3b))*a24a.v(ds2a,ds4a))*CW*CW)*MZ*MZ-three*((s24s.v(ds2b,ds4b)*s123a.v(ds1b,ds3b)-s23s.v(ds2b,ds3b)*s431a.v(ds4b,ds1b))*a24a.v(ds2a,ds4a)*s13s.v(ds1a,ds3a)+(s23s.v(ds2a,ds3b)*s24s.v(ds2b,ds4a)*s124a.v(ds1b,ds4b)-s13s.v(ds1b,ds3b)*s24s.v(ds2a,ds4a)*s432a.v(ds4b,ds2b)+(-s24s.v(ds2b,ds4b)*s123a.v(ds1b,ds3b)+(-s124a.v(ds1b,ds4b)+s431a.v(ds4b,ds1b))*s23s.v(ds2b,ds3b)+s13s.v(ds1b,ds3b)*s432a.v(ds4b,ds2b))*a24a.v(ds2a,ds4a))*a13a.v(ds1a,ds3a))*CW*MZ)/pDenTZ;

      //U-Channel Z
      //preZ = e*e/(2.0*MW*MW*MZ*MZ*SW*SW);
      amplitude += - normFactor*preZ*(two*m12*a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)+(a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)+a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s13s.v(ds1b,ds3b)*s24s.v(ds2b,ds4b)+((-four*a23a.v(ds2a,ds3b)*a34a.v(ds3a,ds4a)*s12s.v(ds1a,ds2b)+a23a.v(ds2a,ds3a)*a24a.v(ds2b,ds4a)*s13s.v(ds1a,ds3b)+three*a24a.v(ds2a,ds4a)*s13s.v(ds1a,ds3a)*s23s.v(ds2b,ds3b))*s14s.v(ds1b,ds4b)+((three*a23a.v(ds2a,ds3b)-s23s.v(ds2a,ds3b))*a34a.v(ds3a,ds4b)*s12s.v(ds1b,ds2b)+a24a.v(ds2a,ds4b)*s13s.v(ds1b,ds3a)*s23s.v(ds2b,ds3b)-(three*a24a.v(ds2b,ds4b)*s13s.v(ds1b,ds3b)+four*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b))*a23a.v(ds2a,ds3a))*a14a.v(ds1a,ds4a)+(-three*a23a.v(ds2a,ds3b)*a24a.v(ds2b,ds4a)*s14s.v(ds1b,ds4b)-three*(a14a.v(ds1b,ds4a)*a24a.v(ds2a,ds4b)-two*a24a.v(ds2a,ds4a)*s14s.v(ds1b,ds4b))*s23s.v(ds2b,ds3b)+((-three*a14a.v(ds1b,ds4a)+four*s14s.v(ds1b,ds4a))*a23a.v(ds2a,ds3b)+(four*a14a.v(ds1b,ds4a)-three*s14s.v(ds1b,ds4a))*s23s.v(ds2a,ds3b))*s24s.v(ds2b,ds4b))*a13a.v(ds1a,ds3a)+((three*a14a.v(ds1b,ds4a)-s14s.v(ds1b,ds4a))*a23a.v(ds2a,ds3b)-four*a14a.v(ds1b,ds4a)*s23s.v(ds2a,ds3b))*a12a.v(ds1a,ds2b)*s34s.v(ds3a,ds4b))*CW*CW)*MZ*MZ+three*((-s143a.v(ds1b,ds3b)+s321a.v(ds3b,ds1b))*a24a.v(ds2a,ds4a)*s14s.v(ds1a,ds4b)*s23s.v(ds2b,ds3a)+(a24a.v(ds2a,ds4b)*s23s.v(ds2b,ds3a)*s143a.v(ds1b,ds3b)+(a23a.v(ds2a,ds3a)-s23s.v(ds2a,ds3a))*s14s.v(ds1b,ds4b)*s243a.v(ds2b,ds3b)-a23a.v(ds2a,ds3a)*s23s.v(ds2b,ds3b)*s421a.v(ds4b,ds1b))*a14a.v(ds1a,ds4a)+(-a24a.v(ds2b,ds4a)*s14s.v(ds1a,ds4b)*s321a.v(ds3b,ds1b)+s14s.v(ds1a,ds4a)*s23s.v(ds2b,ds3b)*s421a.v(ds4b,ds1b))*a23a.v(ds2a,ds3a))*CW*MZ)/pDenUZ;

      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble WWWW::amp2(){
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
    //Symmetry factor 1/2
    return amp2/18.0;
  }
  



  //  Tests
  int test_WWWW(){
    int n=0;//Number of fails
    std::cout<<"\t* W+, W+ -> W+, W+      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      WWWW WWWWAmp = WWWW(EE,mh,0,MW,0,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.330498476822454E+02,2.574400305762363E+01,1.117817047479175E+01,6.389083613347302E+00,4.249873849952222E+00,3.130785827467675E+00,2.492964262988100E+00,2.117842664358412E+00,1.905654805317800E+00,1.809199057352321E+00,1.809199057352527E+00,1.905654805317751E+00,2.117842664358410E+00,2.492964262988184E+00,3.130785827467681E+00,4.249873849952687E+00,6.389083613346880E+00,1.117817047479180E+01,2.574400305762331E+01,1.330498476822448E+02};
      i += WWWWAmp.test_2to2_amp2([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH);
      i += WWWWAmp.test_2to2_amp2_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH);
      i += WWWWAmp.test_2to2_amp2_boosts([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH);
      i += WWWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {8.130631379749040E+01,1.625408821644292E+01,7.872366706414579E+00,4.853128067966389E+00,3.389299273693746E+00,2.574183518277992E+00,2.087897391405524E+00,1.792505256033727E+00,1.621788480938466E+00,1.543225271210162E+00,1.543225271210159E+00,1.621788480938461E+00,1.792505256033726E+00,2.087897391405526E+00,2.574183518277990E+00,3.389299273693748E+00,4.853128067966391E+00,7.872366706414574E+00,1.625408821644291E+01,8.130631379749011E+01};
      i += WWWWAmp.test_2to2_amp2([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH2);
      i += WWWWAmp.test_2to2_amp2_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH2);
      i += WWWWAmp.test_2to2_amp2_boosts([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH2);
      i += WWWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {3.279004391475152E+08,3.797423564810603E+07,1.438213861332197E+07,7.804461683462548E+06,5.085792978202032E+06,3.721867791593252E+06,2.962859154610754E+06,2.522608909310526E+06,2.275527956386337E+06,2.163658176582976E+06,2.163658176582976E+06,2.275527956386336E+06,2.522608909310525E+06,2.962859154610754E+06,3.721867791593251E+06,5.085792978202032E+06,7.804461683462544E+06,1.438213861332196E+07,3.797423564810599E+07,3.279004391475137E+08};
      i += WWWWAmp.test_2to2_amp2([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH4);
      i += WWWWAmp.test_2to2_amp2_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH4);
      i += WWWWAmp.test_2to2_amp2_boosts([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH4);
      i += WWWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH4);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=250\n";
      mh=1;
      WWWWAmp.set_masses(mh,MW);
      pspatial=250;
      ldouble dataCH5[20] = {1.289354271472530E+02,2.543386665305536E+01,1.107701607612621E+01,6.326853267894585E+00,4.199402895661335E+00,3.085154978598395E+00,2.449788658833762E+00,2.076059777061563E+00,1.864661540138805E+00,1.768568952283996E+00,1.768568952284205E+00,1.864661540138756E+00,2.076059777061562E+00,2.449788658833852E+00,3.085154978598395E+00,4.199402895661808E+00,6.326853267894166E+00,1.107701607612626E+01,2.543386665305505E+01,1.289354271472524E+02};
      i += WWWWAmp.test_2to2_amp2([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH5);
      i += WWWWAmp.test_2to2_amp2_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH5);
      i += WWWWAmp.test_2to2_amp2_boosts([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH5);
      i += WWWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH5);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH6[20] = {6.322075388133022E+01,1.417987528933380E+01,7.176682746487629E+00,4.526824780745282E+00,3.204494441286307E+00,2.454368906813108E+00,2.001227603501230E+00,1.723632291151748E+00,1.562318626091595E+00,1.487852428745605E+00,1.487852428745602E+00,1.562318626091590E+00,1.723632291151748E+00,2.001227603501232E+00,2.454368906813107E+00,3.204494441286310E+00,4.526824780745284E+00,7.176682746487625E+00,1.417987528933380E+01,6.322075388133001E+01};
      i += WWWWAmp.test_2to2_amp2([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH6);
      i += WWWWAmp.test_2to2_amp2_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH6);
      i += WWWWAmp.test_2to2_amp2_boosts([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH6);
      i += WWWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH6);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH7[20] = {2.608918746427962E+08,1.996628642120627E+07,5.214336344714892E+06,2.013519198280014E+06,9.581250189006138E+05,5.249405587342475E+05,3.221793458499209E+05,2.200997783803048E+05,1.685599275956673E+05,1.467065159549060E+05,1.467065159549062E+05,1.685599275956670E+05,2.200997783803041E+05,3.221793458499206E+05,5.249405587342468E+05,9.581250189006139E+05,2.013519198280013E+06,5.214336344714886E+06,1.996628642120624E+07,2.608918746427949E+08};
      i += WWWWAmp.test_2to2_amp2([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH7);
      i += WWWWAmp.test_2to2_amp2_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH7);
      i += WWWWAmp.test_2to2_amp2_boosts([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH7);
      i += WWWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH7);
      //std::cout<<"\n# mh=125, MW=50, pspatial=125\n";
      mh=125;
      MW=50;
      MZ=MW/CW;
      WWWWAmp.set_masses(mh,MW);
      pspatial = 125;
      ldouble dataCH9[20] = {1.095017954580556E+02,2.261852152671634E+01,1.027443279910154E+01,6.031164829907299E+00,4.082537678952108E+00,3.045044715189483E+00,2.446904275344808E+00,2.092473906756835E+00,1.891042727666908E+00,1.799238102159719E+00,1.799238102159731E+00,1.891042727666844E+00,2.092473906756904E+00,2.446904275344853E+00,3.045044715189508E+00,4.082537678952161E+00,6.031164829907317E+00,1.027443279910160E+01,2.261852152671633E+01,1.095017954580552E+02};
      i += WWWWAmp.test_2to2_amp2([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH9);
      i += WWWWAmp.test_2to2_amp2_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH9);
      i += WWWWAmp.test_2to2_amp2_boosts([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH9);
      i += WWWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH9);
      //std::cout<<"\n# mh=125, MW=50, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH10[20] = {4.913759765767320E+07,5.692401739240093E+06,2.156579201564024E+06,1.170625609325476E+06,7.630634979400714E+05,5.585720923930962E+05,4.447645637648232E+05,3.787475764988856E+05,3.416949842391653E+05,3.249182948145064E+05,3.249182948145064E+05,3.416949842391651E+05,3.787475764988855E+05,4.447645637648231E+05,5.585720923930962E+05,7.630634979400713E+05,1.170625609325475E+06,2.156579201564022E+06,5.692401739240088E+06,4.913759765767299E+07};
      i += WWWWAmp.test_2to2_amp2([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH10);
      i += WWWWAmp.test_2to2_amp2_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH10);
      i += WWWWAmp.test_2to2_amp2_boosts([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH10);
      i += WWWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH10);
      //std::cout<<"\n# mh=1, MW=50, pspatial=125\n";
      mh=1;
      MW=50;
      MZ=MW/CW;
      WWWWAmp.set_masses(mh,MW);
      pspatial = 125;
      ldouble dataCH11[20] = {1.013914231929271E+02,2.169090775285456E+01,9.941653645477583E+00,5.843129748126032E+00,3.944770146758831E+00,2.928324671792071E+00,2.340051571020195E+00,1.990548099383739E+00,1.791576100165265E+00,1.700804453787553E+00,1.700804453787564E+00,1.791576100165202E+00,1.990548099383809E+00,2.340051571020236E+00,2.928324671792093E+00,3.944770146758884E+00,5.843129748126052E+00,9.941653645477649E+00,2.169090775285454E+01,1.013914231929267E+02};
      i += WWWWAmp.test_2to2_amp2([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH11);
      i += WWWWAmp.test_2to2_amp2_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH11);
      i += WWWWAmp.test_2to2_amp2_boosts([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH11);
      i += WWWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH11);
      //std::cout<<"\n# mh=1, MW=50, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH12[20] = {3.910098778749280E+07,2.994425252220048E+06,7.826557892667830E+05,3.025216758387694E+05,1.441224074864304E+05,7.906947324540273E+04,4.860173953117031E+04,3.325412067510127E+04,2.550143925438343E+04,2.221317498939151E+04,2.221317498939151E+04,2.550143925438344E+04,3.325412067510116E+04,4.860173953117021E+04,7.906947324540267E+04,1.441224074864303E+05,3.025216758387692E+05,7.826557892667827E+05,2.994425252220043E+06,3.910098778749262E+07};
      i += WWWWAmp.test_2to2_amp2([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH12);
      i += WWWWAmp.test_2to2_amp2_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH12);
      i += WWWWAmp.test_2to2_amp2_boosts([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH12);
      i += WWWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WWWWAmp.amp2(); }, MW,MW,MW,MW,pspatial,dataCH12);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
  

}
