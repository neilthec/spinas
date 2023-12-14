
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

//File:  SPINAS/SM/ZWZW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ZWZW.h"

namespace spinas {

  ZWZW::ZWZW(const ldouble& echarge, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), mh(massh), wh(widthh), MW(massW), WW(widthW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW) {
    CW = std::sqrt(1.0-sinW*sinW);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    proph = propagator(mh,wh);
    propW = propagator(MW,WW);
    p1=particle(MZ);
    p2=particle(MW);
    p3=particle(MZ);
    p4=particle(MW);
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s42s = sproduct(SQUARE,&p4,&p2);
    a42a = sproduct(ANGLE,&p4,&p2);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s32s = sproduct(SQUARE,&p3,&p2);
    a32a = sproduct(ANGLE,&p3,&p2);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    //
    s241a = sproduct(SQUARE,&p2,&p4,&p1);
    s134a = sproduct(SQUARE,&p1,&p3,&p4);
    s431a = sproduct(SQUARE,&p4,&p3,&p1);
    s124a = sproduct(SQUARE,&p1,&p2,&p4);
    s132a = sproduct(SQUARE,&p1,&p3,&p2);
    s231a = sproduct(SQUARE,&p2,&p3,&p1);
    s324a = sproduct(SQUARE,&p3,&p2,&p4);
    s243a = sproduct(SQUARE,&p2,&p4,&p3);
    //Couplings
    preh = e*e/(MW*MW*SW*SW);
    preW = e*e/(2.0*MW*MW*MZ*MZ*SW*SW);
    pre4 = e*e/(MZ*MZ*MZ*MZ*SW*SW);
  }
  void ZWZW::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(MZ);
    p2.set_mass(MW);
    p3.set_mass(MZ);
    p4.set_mass(MW);
    proph.set_mass(mh);
    propW.set_mass(MW);
    //Couplings
    preh = e*e/(MW*MW*SW*SW);
    preW = e*e/(2.0*MW*MW*MZ*MZ*SW*SW);
    pre4 = e*e/(MZ*MZ*MZ*MZ*SW*SW);
  }
  void ZWZW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    s13s.update();
    a13a.update();
    s42s.update();
    a42a.update();
    s34s.update();
    a34a.update();
    s32s.update();
    a32a.update();
    s12s.update();
    a12a.update();
    s14s.update();
    a14a.update();
    //
    s241a.update();
    s134a.update();
    s431a.update();
    s124a.update();
    s132a.update();
    s231a.update();
    s324a.update();
    s243a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propW.den(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=proph.den(propTP);
    pDenU=propW.den(propUP);
    //Mandelstahms
    m13 = -2.*p1.dot(p3)+2.*MZ*MZ;//(p1-p3)^2
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ZWZW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1, two=2, three=3, four=4, six=6;
    int ds1a, ds1b, ds2a, ds2b, ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds1,ds2,ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds1,ds2,ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds1,ds1a,ds1b, ds2,ds2a,ds2b, ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);

      
      //S-Channel h
      //preh = e*e/(MW*MW*SW*SW);
      //ZZW-W+ all ingoing
      //preh [12]<12>[34]<34>/(s-Mh^2)
      //ZW+ZW-: 4->2->3->4
      //preh [13]<13>[42]<42>/(t-Mh^2)
      //34 out:
      //preh [13]<13>[42]<42>/(t-Mh^2)
      amplitude += normFactor*preh*s13s.v(ds1a,ds3a)*a13a.v(ds1b,ds3b)*s42s.v(ds4a,ds2a)*a42a.v(ds4b,ds2b)/pDenT;

      //4-Point
      //pre4 = e*e/(MZ*MZ*MZ*MZ*SW*SW);
      //ZZW-W+ all in:
      //pre4 [12]<12>[34]<34>
      //ZW+ZW-: 4->2->3->4
      //pre4 [13]<13>[42]<42>
      //34 out:
      //pre4 [13]<13>[42]<42>
      amplitude += normFactor*pre4*s13s.v(ds1a,ds3a)*a13a.v(ds1b,ds3b)*s42s.v(ds4a,ds2a)*a42a.v(ds4b,ds2b);

      
      //T-Channel W
      //preW = e*e/(2.0*MW*MW*MZ*MZ*SW*SW);
      //ZZW-W+ all in:
      //-(-2CW^2 s12 <13><24>[13][24]
      //  +MZ^2(-2<13><24>[13][24]
      //   +(-4<23><24>[13][14]-3<13><24>[14][23]+6<13><24>[13][24]+3<34>[12][13][24]-<13><14>[23][24])CW^2
      //   +(4<24><34>[12][13]-3<13><23><24>[14]+3<13><14><24>[23]-3<24>[13][14][23]+<13><34>[12][24]
      //     -3<23>[13][14][24]-3<14>[13][23][24]-3<13>[14][23][24]+<12><24>[13][34]+4<12><13>[24][34])CW^3
      //   +(-<14><24>[13][23]-4<13><24>[14][23]-<14><23>[13][24]-4<13><23>[14][24]+3<12>[13][24][34])CW^4
      //   )
      //  +MZ(
      //    (3<24>[13][24][123>+3<13>[23][24][124>-3<13><24>[23][431>-3<13><24>[13][432>)CW^2
      //    +(3<13><24>[24][123>+3<13><24>[23][124>-3<24>[13][23][431>-3<13>[13][24][432>)CW^3
      //   )
      // )/(t-MW^2)
      //ZW+ZW- 4->2->3->4:
      //-(-2CW^2 s13 <14><32>[14][32]
      //  +MZ^2(-2<14><32>[14][32]
      //   +(-4<34><32>[14][12]-3<14><32>[12][34]+6<14><32>[14][32]+3<42>[13][14][32]-<14><12>[34][32])CW^2
      //   +(4<32><42>[13][14]-3<14><34><32>[12]+3<14><12><32>[34]-3<32>[14][12][34]+<14><42>[13][32]
      //     -3<34>[14][12][32]-3<12>[14][34][32]-3<14>[12][34][32]+<13><32>[14][42]+4<13><14>[32][42])CW^3
      //   +(-<12><32>[14][34]-4<14><32>[12][34]-<12><34>[14][32]-4<14><34>[12][32]+3<13>[14][32][42])CW^4
      //   )
      //  +MZ(
      //    (3<32>[14][32][134>+3<14>[34][32][132>-3<14><32>[34][241>-3<14><32>[14][243>)CW^2
      //    +(3<14><32>[32][134>+3<14><32>[34][132>-3<32>[14][34][241>-3<14>[14][32][243>)CW^3
      //   )
      // )/(t-MW^2)
      //34 out:
      //-(-2CW^2 s13 <14><32>[14][32]
      //  +MZ^2(-2<14><32>[14][32]
      //   +(+4<34><32>[14][12]-3<14><32>[12][34]+6<14><32>[14][32]-3<42>[13][14][32]+<14><12>[34][32])CW^2
      //   +(4<32><42>[13][14]-3<14><34><32>[12]+3<14><12><32>[34]+3<32>[14][12][34]+<14><42>[13][32]
      //     -3<34>[14][12][32]-3<12>[14][34][32]+3<14>[12][34][32]+<13><32>[14][42]+4<13><14>[32][42])CW^3
      //   +(+<12><32>[14][34]-4<14><32>[12][34]-<12><34>[14][32]+4<14><34>[12][32]-3<13>[14][32][42])CW^4
      //   )
      //  +MZ(
      //    (-3<32>[14][32][134>+3<14>[34][32][132>+3<14><32>[34][241>-3<14><32>[14][243>)CW^2
      //    +(3<14><32>[32][134>-3<14><32>[34][132>-3<32>[14][34][241>+3<14>[14][32][243>)CW^3
      //   )
      // )/(u-MW^2)
      amplitude += -normFactor*preW*(
				     -two*CW*CW*m13*a14a.v(ds1a,ds4a)*a32a.v(ds3a,ds2a)*s14s.v(ds1b,ds4b)*s32s.v(ds3b,ds2b)
				     +MZ*MZ*(-two*a14a.v(ds1a,ds4a)*a32a.v(ds3a,ds2a)*s14s.v(ds1b,ds4b)*s32s.v(ds3b,ds2b)
				       +(+four*a34a.v(ds3a,ds4a)*a32a.v(ds3b,ds2a)*s14s.v(ds1a,ds4b)*s12s.v(ds1b,ds2b)
					 -three*a14a.v(ds1a,ds4a)*a32a.v(ds3a,ds2a)*s12s.v(ds1b,ds2b)*s34s.v(ds3b,ds4b)
					 +six*a14a.v(ds1a,ds4a)*a32a.v(ds3a,ds2a)*s14s.v(ds1b,ds4b)*s32s.v(ds3b,ds2b)
					 -three*a42a.v(ds4a,ds2a)*s13s.v(ds1a,ds3b)*s14s.v(ds1b,ds4b)*s32s.v(ds3a,ds2b)
					 +a14a.v(ds1a,ds4a)*a12a.v(ds1b,ds2a)*s34s.v(ds3a,ds4b)*s32s.v(ds3b,ds2b))*CW*CW
				       +(+four*a32a.v(ds3a,ds2a)*a42a.v(ds4a,ds2b)*s13s.v(ds1a,ds3b)*s14s.v(ds1b,ds4b)
					 -three*a14a.v(ds1a,ds4a)*a34a.v(ds3a,ds4b)*a32a.v(ds3b,ds2a)*s12s.v(ds1b,ds2b)
					 +three*a14a.v(ds1a,ds4a)*a12a.v(ds1b,ds2a)*a32a.v(ds3a,ds2b)*s34s.v(ds3b,ds4b)
					 +three*a32a.v(ds3a,ds2a)*s14s.v(ds1a,ds4a)*s12s.v(ds1b,ds2b)*s34s.v(ds3b,ds4b)
					 +a14a.v(ds1a,ds4b)*a42a.v(ds4a,ds2a)*s13s.v(ds1b,ds3b)*s32s.v(ds3a,ds2b)
					 -three*a34a.v(ds3a,ds4a)*s14s.v(ds1a,ds4b)*s12s.v(ds1b,ds2a)*s32s.v(ds3b,ds2b)
					 -three*a12a.v(ds1a,ds2a)*s14s.v(ds1b,ds4a)*s34s.v(ds3a,ds4b)*s32s.v(ds3b,ds2b)
					 +three*a14a.v(ds1a,ds4a)*s12s.v(ds1b,ds2a)*s34s.v(ds3a,ds4b)*s32s.v(ds3b,ds2b)
					 +a13a.v(ds1a,ds3b)*a32a.v(ds3a,ds2a)*s14s.v(ds1b,ds4b)*s42s.v(ds4a,ds2b)
					 +four*a13a.v(ds1a,ds3b)*a14a.v(ds1b,ds4b)*s32s.v(ds3a,ds2a)*s42s.v(ds4a,ds2b))*CW*CW*CW
				       +(+a12a.v(ds1a,ds2a)*a32a.v(ds3a,ds2b)*s14s.v(ds1b,ds4a)*s34s.v(ds3b,ds4b)
					 -four*a14a.v(ds1a,ds4a)*a32a.v(ds3a,ds2a)*s12s.v(ds1b,ds2b)*s34s.v(ds3b,ds4b)
					 -a12a.v(ds1a,ds2a)*a34a.v(ds3a,ds4a)*s14s.v(ds1b,ds4b)*s32s.v(ds3b,ds2b)
					 +four*a14a.v(ds1a,ds4a)*a34a.v(ds3a,ds4b)*s12s.v(ds1b,ds2a)*s32s.v(ds3b,ds2b)
					 -three*a13a.v(ds1a,ds3b)*s14s.v(ds1b,ds4b)*s32s.v(ds3a,ds2a)*s42s.v(ds4a,ds2b))*CW*CW*CW*CW
				       )
				     +MZ*(
				       (-three*a32a.v(ds3a,ds2a)*s14s.v(ds1a,ds4a)*s32s.v(ds3b,ds2b)*s134a.v(ds1b,ds4b)
					+three*a14a.v(ds1a,ds4a)*s34s.v(ds3a,ds4b)*s32s.v(ds3b,ds2a)*s132a.v(ds1b,ds2b)
					+three*a14a.v(ds1a,ds4a)*a32a.v(ds3a,ds2a)*s34s.v(ds3b,ds4b)*s241a.v(ds2b,ds1b)
					-three*a14a.v(ds1a,ds4a)*a32a.v(ds3a,ds2a)*s14s.v(ds1b,ds4b)*s243a.v(ds2b,ds3b))*CW*CW
				       +(+three*a14a.v(ds1a,ds4a)*a32a.v(ds3a,ds2a)*s32s.v(ds3b,ds2b)*s134a.v(ds1b,ds4b)
					 -three*a14a.v(ds1a,ds4a)*a32a.v(ds3a,ds2a)*s34s.v(ds3b,ds4b)*s132a.v(ds1b,ds2b)
					 -three*a32a.v(ds3a,ds2a)*s14s.v(ds1a,ds4a)*s34s.v(ds3b,ds4b)*s241a.v(ds2b,ds1b)
					 +three*a14a.v(ds1a,ds4a)*s14s.v(ds1b,ds4b)*s32s.v(ds3a,ds2a)*s243a.v(ds2b,ds3b))*CW*CW*CW
				       )
				     )/pDenU;

      //U-Channel W
      //preW = e*e/(2.0*MW*MW*MZ*MZ*SW*SW);
      //ZZW-W+ all in:
      //-(-2CW^2 s12<14><23>[14][23]
      // +(-2<14><23>[14][23]
      //  +(-3<14><23><34>[12]-1<14><24>[13][23]-4<13><24>[14][23]-1<14><23>[13][24]-4<13><23>[14][24])CW^4
      //  +(-1<23><24>[13][14]+6<14><23>[14][23]-3<13><24>[14][23]-4<13><14>[23][24]-3<12><14><23>[34])CW^2
      //  +(-3<14><23><24>[13]-3<13><23><24>[14]-4<23><34>[12][14]-3<13><14><24>[23]-1<14><34>[12][23]
      //   +3<24>[13][14][23]-3<13><14><23>[24]-3<13>[14][23][24]-1<12><23>[14][34]-4<12><14>[23][34])CW^3
      //  )MZ^2
      //  +(
      //    (-3<24>[14][23][143>-3<14>[14][23][243>+3<23><24>[14][321>+3<14><23>[23][421>)CW^2
      //   +(-3<14><24>[23][143>-3<14><23>[14][243>+3<24>[14][23][321>+3<23>[14][23][421>)CW^3
      //  )MZ
      // )/(u-MW^2)
      //ZW+ZW- 4->2->3->4:
      //-(-2CW^2 s13<12><34>[12][34]
      // +(-2<12><34>[12][34]
      //  +(-3<12><34><42>[13]-1<12><32>[14][34]-4<14><32>[12][34]-1<12><34>[14][32]-4<14><34>[12][32])CW^4
      //  +(-1<34><32>[14][12]+6<12><34>[12][34]-3<14><32>[12][34]-4<14><12>[34][32]-3<13><12><34>[42])CW^2
      //  +(-3<12><34><32>[14]-3<14><34><32>[12]-4<34><42>[13][12]-3<14><12><32>[34]-1<12><42>[13][34]
      //   +3<32>[14][12][34]-3<14><12><34>[32]-3<14>[12][34][32]-1<13><34>[12][42]-4<13><12>[34][42])CW^3
      //  )MZ^2
      //  +(
      //    (-3<32>[12][34][124>-3<12>[12][34][324>+3<34><32>[12][431>+3<12><34>[34][231>)CW^2
      //   +(-3<12><32>[34][124>-3<12><34>[12][324>+3<32>[12][34][431>+3<34>[12][34][231>)CW^3
      //  )MZ
      // )/(u-MW^2)
      //34 out:
      //-(-2CW^2 s13<12><34>[12][34]
      // +(-2<12><34>[12][34]
      //  +(+3<12><34><42>[13]+1<12><32>[14][34]-4<14><32>[12][34]-1<12><34>[14][32]+4<14><34>[12][32])CW^4
      //  +(+1<34><32>[14][12]+6<12><34>[12][34]-3<14><32>[12][34]+4<14><12>[34][32]+3<13><12><34>[42])CW^2
      //  +(+3<12><34><32>[14]-3<14><34><32>[12]+4<34><42>[13][12]-3<14><12><32>[34]+1<12><42>[13][34]
      //   -3<32>[14][12][34]+3<14><12><34>[32]+3<14>[12][34][32]+1<13><34>[12][42]+4<13><12>[34][42])CW^3
      //  )MZ^2
      //  +(
      //    (-3<32>[12][34][124>+3<12>[12][34][324>+3<34><32>[12][431>-3<12><34>[34][231>)CW^2
      //   +(-3<12><32>[34][124>+3<12><34>[12][324>+3<32>[12][34][431>-3<34>[12][34][231>)CW^3
      //  )MZ
      // )/(s-MW^2)
      amplitude += -normFactor*preW*(
				     -two*CW*CW*m13*a12a.v(ds1a,ds2a)*a34a.v(ds3a,ds4a)*s12s.v(ds1b,ds2b)*s34s.v(ds3b,ds4b)
				     +MZ*MZ*(-two*a12a.v(ds1a,ds2a)*a34a.v(ds3a,ds4a)*s12s.v(ds1b,ds2b)*s34s.v(ds3b,ds4b)
				       +(+three*a12a.v(ds1a,ds2a)*a34a.v(ds3a,ds4b)*a42a.v(ds4a,ds2b)*s13s.v(ds1b,ds3b)
					 +a12a.v(ds1a,ds2a)*a32a.v(ds3a,ds2b)*s14s.v(ds1b,ds4a)*s34s.v(ds3b,ds4b)
					 -four*a14a.v(ds1a,ds4a)*a32a.v(ds3a,ds2a)*s12s.v(ds1b,ds2b)*s34s.v(ds3b,ds4b)
					 -a12a.v(ds1a,ds2a)*a34a.v(ds3a,ds4a)*s14s.v(ds1b,ds4b)*s32s.v(ds3b,ds2b)
					 +four*a14a.v(ds1a,ds4a)*a34a.v(ds3a,ds4b)*s12s.v(ds1b,ds2a)*s32s.v(ds3b,ds2b))*CW*CW*CW*CW
				       +(+a34a.v(ds3a,ds4a)*a32a.v(ds3b,ds2a)*s14s.v(ds1a,ds4b)*s12s.v(ds1b,ds2b)
					 +six*a12a.v(ds1a,ds2a)*a34a.v(ds3a,ds4a)*s12s.v(ds1b,ds2b)*s34s.v(ds3b,ds4b)
					 -three*a14a.v(ds1a,ds4a)*a32a.v(ds3a,ds2a)*s12s.v(ds1b,ds2b)*s34s.v(ds3b,ds4b)
					 +four*a14a.v(ds1a,ds4a)*a12a.v(ds1b,ds2a)*s34s.v(ds3a,ds4b)*s32s.v(ds3b,ds2b)
					 +three*a13a.v(ds1a,ds3b)*a12a.v(ds1b,ds2a)*a34a.v(ds3a,ds4b)*s42s.v(ds4a,ds2b))*CW*CW
				       +(+three*a12a.v(ds1a,ds2a)*a34a.v(ds3a,ds4a)*a32a.v(ds3b,ds2b)*s14s.v(ds1b,ds4b)
					 -three*a14a.v(ds1a,ds4a)*a34a.v(ds3a,ds4b)*a32a.v(ds3b,ds2a)*s12s.v(ds1b,ds2b)
					 +four*a34a.v(ds3a,ds4b)*a42a.v(ds4a,ds2a)*s13s.v(ds1a,ds3b)*s12s.v(ds1b,ds2b)
					 -three*a14a.v(ds1a,ds4a)*a12a.v(ds1b,ds2a)*a32a.v(ds3a,ds2b)*s34s.v(ds3b,ds4b)
					 +a12a.v(ds1a,ds2a)*a42a.v(ds4a,ds2b)*s13s.v(ds1b,ds3b)*s34s.v(ds3a,ds4b)
					 -three*a32a.v(ds3a,ds2a)*s14s.v(ds1a,ds4a)*s12s.v(ds1b,ds2b)*s34s.v(ds3b,ds4b)
					 +three*a14a.v(ds1a,ds4a)*a12a.v(ds1b,ds2a)*a34a.v(ds3a,ds4b)*s32s.v(ds3b,ds2b)
					 +three*a14a.v(ds1a,ds4a)*s12s.v(ds1b,ds2a)*s34s.v(ds3a,ds4b)*s32s.v(ds3b,ds2b)
					 +a13a.v(ds1a,ds3b)*a34a.v(ds3a,ds4b)*s12s.v(ds1b,ds2a)*s42s.v(ds4a,ds2b)
					 +four*a13a.v(ds1a,ds3b)*a12a.v(ds1b,ds2a)*s34s.v(ds3a,ds4b)*s42s.v(ds4a,ds2b))*CW*CW*CW
				       )
				     +MZ*((-three*a32a.v(ds3a,ds2a)*s12s.v(ds1a,ds2b)*s34s.v(ds3b,ds4a)*s124a.v(ds1b,ds4b)
					   +three*a12a.v(ds1a,ds2a)*s12s.v(ds1b,ds2b)*s34s.v(ds3a,ds4a)*s324a.v(ds3b,ds4b)
					   +three*a34a.v(ds3a,ds4a)*a32a.v(ds3b,ds2a)*s12s.v(ds1a,ds2b)*s431a.v(ds4b,ds1b)
					   -three*a12a.v(ds1a,ds2a)*a34a.v(ds3a,ds4a)*s34s.v(ds3b,ds4b)*s231a.v(ds2b,ds1b))*CW*CW
					  +(-three*a12a.v(ds1a,ds2a)*a32a.v(ds3a,ds2b)*s34s.v(ds3b,ds4a)*s124a.v(ds1b,ds4b)
					    +three*a12a.v(ds1a,ds2a)*a34a.v(ds3a,ds4a)*s12s.v(ds1b,ds2b)*s324a.v(ds3b,ds4b)
					    +three*a32a.v(ds3a,ds2a)*s12s.v(ds1a,ds2b)*s34s.v(ds3b,ds4a)*s431a.v(ds4b,ds1b)
					    -three*a34a.v(ds3a,ds4a)*s12s.v(ds1a,ds2a)*s34s.v(ds3b,ds4b)*s231a.v(ds2b,ds1b))*CW*CW*CW
					  )
				     )/pDenS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ZWZW::amp2(){
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
  int test_ZWZW(){
    int n=0;//Number of fails
    std::cout<<"\t* Z , W+ -> Z , W+      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      ZWZW ZWZWAmp = ZWZW(EE,mh,0,MW,0,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.400180690614037E-01,2.607957753715975E-01,2.771556087075341E-01,2.945821429947157E-01,3.159801849391463E-01,3.437536326113061E-01,3.806001924649887E-01,4.300632148163758E-01,4.971828638613576E-01,5.894936937753744E-01,7.187333871400338E-01,9.039757012553801E-01,1.177748047249108E+00,1.598862119394270E+00,2.281793607715884E+00,3.471998166153202E+00,5.770739441057290E+00,1.097076703430870E+01,2.647364006591239E+01,1.129303506750818E+02};
      i += ZWZWAmp.test_2to2_amp2([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH);
      i += ZWZWAmp.test_2to2_amp2_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH);
      i += ZWZWAmp.test_2to2_amp2_boosts([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH);
      i += ZWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {2.205826058530009E-01,2.313490869869667E-01,2.440150543748788E-01,2.597088532111242E-01,2.797751400967960E-01,3.058872099354161E-01,3.402312196777793E-01,3.857792915093713E-01,4.467060714264587E-01,5.290499569329832E-01,6.418047694946938E-01,7.987968527972351E-01,1.022061475485640E+00,1.348242406927561E+00,1.841506512539880E+00,2.621696730434066E+00,3.932008899773368E+00,6.324014313565417E+00,1.126221473432283E+01,2.370622644626177E+01};
      i += ZWZWAmp.test_2to2_amp2([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH2);
      i += ZWZWAmp.test_2to2_amp2_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH2);
      i += ZWZWAmp.test_2to2_amp2_boosts([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH2);
      i += ZWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {1.752549513231766E+00,1.752723295484089E+00,1.752897098491987E+00,1.753070922258355E+00,1.753244766786090E+00,1.753418632078086E+00,1.753592518137243E+00,1.753766424966457E+00,1.753940352568625E+00,1.754114300946647E+00,1.754288270103421E+00,1.754462260041846E+00,1.754636270764821E+00,1.754810302275248E+00,1.754984354576027E+00,1.755158427670059E+00,1.755332521560246E+00,1.755506636249491E+00,1.755680771740695E+00,1.755854928036761E+00};
      i += ZWZWAmp.test_2to2_amp2([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH4);
      i += ZWZWAmp.test_2to2_amp2_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH4);
      i += ZWZWAmp.test_2to2_amp2_boosts([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH4);
      i += ZWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH4);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=250\n";
      mh=1;
      ZWZWAmp.set_masses(mh,MW);
      pspatial=250;
      ldouble dataCH5[20] = {5.434202095897271E-01,2.921334755350321E-01,2.841980532519765E-01,2.929011851812670E-01,3.094281322492493E-01,3.338010801802480E-01,3.679025530856882E-01,4.148944667010394E-01,4.795881213194446E-01,5.693440741068656E-01,6.957342866448265E-01,8.776417090137557E-01,1.147340419846507E+00,1.563267974370828E+00,2.239295343648624E+00,3.419806129963957E+00,5.703937229470132E+00,1.087950339079743E+01,2.633335106354514E+01,1.126438528535888E+02};
      i += ZWZWAmp.test_2to2_amp2([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH5);
      i += ZWZWAmp.test_2to2_amp2_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH5);
      i += ZWZWAmp.test_2to2_amp2_boosts([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH5);
      i += ZWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH5);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH6[20] = {5.320672150791189E+00,7.429423353377017E-01,4.106984197380950E-01,3.345080214918089E-01,3.167738149211194E-01,3.229862502307154E-01,3.448772477310487E-01,3.814959875725019E-01,4.352504658088600E-01,5.112504907694031E-01,6.179168365168304E-01,7.686489751599049E-01,9.850892267778447E-01,1.303438732780967E+00,1.787277745221810E+00,2.555569509346454E+00,3.850046354536690E+00,6.219571531336211E+00,1.112292804385914E+01,2.350509769426082E+01};
      i += ZWZWAmp.test_2to2_amp2([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH6);
      i += ZWZWAmp.test_2to2_amp2_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH6);
      i += ZWZWAmp.test_2to2_amp2_boosts([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH6);
      i += ZWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH6);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH7[20] = {1.095450439185843E+07,7.842475336927310E+06,5.890051724101557E+06,4.585276641399700E+06,3.670438410717359E+06,3.004334668116741E+06,2.504338216407079E+06,2.119485199666940E+06,1.816956911632716E+06,1.574845342268578E+06,1.378073134857445E+06,1.215988648148701E+06,1.080893876374665E+06,9.671134877840674E+05,8.703890443717999E+05,7.874747850693767E+05,7.158618152601587E+05,6.535861158696512E+05,5.990924754934162E+05,5.511364754710854E+05};
      i += ZWZWAmp.test_2to2_amp2([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH7);
      i += ZWZWAmp.test_2to2_amp2_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH7);
      i += ZWZWAmp.test_2to2_amp2_boosts([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH7);
      i += ZWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH7);
      //std::cout<<"\n# mh=125, MW=50, pspatial=125\n";
      mh=125;
      MW=50;
      MZ=MW/CW;
      ZWZWAmp.set_masses(mh,MW);
      pspatial = 125;
      ldouble dataCH8[20] = {2.517648103497350E-01,2.693678265918049E-01,2.875276754038855E-01,3.080619404764429E-01,3.327647698328290E-01,3.636728267698193E-01,4.033508033781278E-01,4.552663152050266E-01,5.243515910478490E-01,6.179153121946469E-01,7.472109782127575E-01,9.302811010805885E-01,1.197417904216073E+00,1.602372659884774E+00,2.247334677810866E+00,3.344736753689353E+00,5.393016064319590E+00,9.784774498243628E+00,2.165703988965396E+01,7.332396924537993E+01};
      i += ZWZWAmp.test_2to2_amp2([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH8);
      i += ZWZWAmp.test_2to2_amp2_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH8);
      i += ZWZWAmp.test_2to2_amp2_boosts([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH8);
      i += ZWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH8);
      //std::cout<<"\n# mh=125, MW=50, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH9[20] = {1.978520426979848E+00,1.978987904631867E+00,1.979455523736125E+00,1.979923284343201E+00,1.980391186503697E+00,1.980859230268238E+00,1.981327415687467E+00,1.981795742812052E+00,1.982264211692682E+00,1.982732822380068E+00,1.983201574924940E+00,1.983670469378054E+00,1.984139505790186E+00,1.984608684212134E+00,1.985078004694715E+00,1.985547467288772E+00,1.986017072045169E+00,1.986486819014790E+00,1.986956708248542E+00,1.987426739797354E+00};
      i += ZWZWAmp.test_2to2_amp2([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH9);
      i += ZWZWAmp.test_2to2_amp2_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH9);
      i += ZWZWAmp.test_2to2_amp2_boosts([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH9);
      i += ZWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH9);
      //std::cout<<"\n# mh=1, MW=50, pspatial=125\n";
      mh=1;
      MW=50;
      MZ=MW/CW;
      ZWZWAmp.set_masses(mh,MW);
      pspatial = 125;
      ldouble dataCH10[20] = {9.715784318652729E-01,3.283632851733490E-01,2.920713098032484E-01,2.933919080357928E-01,3.067255258652148E-01,3.290709903565099E-01,3.614073731250699E-01,4.064833692332597E-01,4.687885565077961E-01,5.552750643907035E-01,6.768382236335739E-01,8.510959729018056E-01,1.107773246563899E+00,1.499790402443351E+00,2.128016670069567E+00,3.202601500495524E+00,5.217604464474328E+00,9.555970852127952E+00,2.132829574511760E+01,7.274087790039592E+01};
      i += ZWZWAmp.test_2to2_amp2([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH10);
      i += ZWZWAmp.test_2to2_amp2_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH10);
      i += ZWZWAmp.test_2to2_amp2_boosts([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH10);
      i += ZWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH10);
      //std::cout<<"\n# mh=125, MW=50, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH11[20] = {1.638238403457604E+06,1.172673297266403E+06,8.806077027768949E+05,6.854383546270001E+05,5.486056546640652E+05,4.489832282831772E+05,3.742090198736779E+05,3.166584717810056E+05,2.714218432527616E+05,2.352217694721867E+05,2.058028075506777E+05,1.815716194031966E+05,1.613767302500372E+05,1.443692484612108E+05,1.299122375410561E+05,1.175202445178503E+05,1.068180391323142E+05,9.751189297380274E+04,8.936922489451976E+04,8.220393916050329E+04};
      i += ZWZWAmp.test_2to2_amp2([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH11);
      i += ZWZWAmp.test_2to2_amp2_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH11);
      i += ZWZWAmp.test_2to2_amp2_boosts([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH11);
      i += ZWZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWZWAmp.amp2(); }, MZ,MW,MZ,MW,pspatial,dataCH11);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
  

}
