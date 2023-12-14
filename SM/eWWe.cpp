
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

//File:  SPINAS/SM/eWWe.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eWWe.h"

namespace spinas {

  eWWe::eWWe(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propne(0,0), propA(0,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);
    proph = propagator(mh,wh);  
    p1=particle(me);
    p2=particle(MW);
    p3=particle(MW);
    p4=particle(me);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    //[413>,[231>,[132>
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
    s431a = sproduct(SQUARE,&p4,&p3,&p1);
    s134a = sproduct(SQUARE,&p1,&p3,&p4);
    //Couplings
    prene = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preh = 2.0*e*e*me/(4.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
    preA = 2*e*e/MW;
  }
  void eWWe::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(me);
    p2.set_mass(MW);
    p3.set_mass(MW);
    p4.set_mass(me);
    proph.set_mass(mh);
    propZ.set_mass(MZ);
    //Couplings
    prene = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preh = 2.0*e*e*me/(4.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
    preA = 2*e*e/MW;
  }
  void eWWe::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s34s.update();
    a34a.update();
    s24s.update();
    a24a.update();
    s14s.update();
    a14a.update();
    //[413>,[231>,[132>
    s213a.update();
    s431a.update();
    s134a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenhU=proph.den(propUP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenneS=propne.den(propSP);
    pDenZU=propZ.den(propUP);
    pDenAU=propA.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eWWe::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds2a, ds2b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds2,ds3);
    ldouble normFactor=get_spin_normalization(ds2,ds3);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, ds3,ds3a,ds3b, i);
      
      //U-Channel h
      //preh = e*e*me/(4.0*MW*MW*SW*SW);
      //eEW-W+ all ingoing: 
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //eW+W-E:  2<->4
      //preh [23] <23> ([14]+<14>)/(u-Mh^2)
      //34 out:
      //- preh [23] <23> ([14]-<14>)/(u-Mh^2)
      amplitude += - normFactor*preh*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(s14s.v(ds1,ds4)-a14a.v(ds1,ds4))/pDenhU;
      
      //S-Channel ne
      //prene = e*e/(2.0*MW*MW*SW*SW);
      //eEW-W+ all ingoing:
      //- prene [23] <14> ([413>-MW <34>))/u
      //eW+W-E: 2<->4
      //+ prene [34] <12> ([213>+MW <23>))/s
      //34 out:
      //- prene [34] <12> ([213>+MW <23>))/s
      amplitude += - normFactor*prene*s34s.v(ds3a,ds4)*a12a.v(ds1,ds2a)*(MW*a23a.v(ds2b,ds3b)+s213a.v(ds2b,ds3b))/pDenneS;

      //U-Channel A
      //preA = 2*e*e/MW;
      //eEW-W+ all in:
      //+ preA ( <13>[24] + [13]<24> + <14>[23] + [14]<23> )( <34> + [34] )/s
      //+ preA/MW [34]<34>([231>+[132>)/s
      //eW+W-E 2<->4:
      //+ preA ( <13>[24] + [13]<24> + <12>[34] + [12]<34> )( <23> + [23] )/u
      //+ preA/MW [23]<23>([431>+[134>)/u
      //34 out:
      //+ preA (- <13>[24] - [13]<24> + <12>[34] + [12]<34> )(- <23> + [23] )/u
      //+ preA/MW [23]<23>([431>-[134>)/u
      amplitude +=
	+ normFactor*preA*(
			   a12a.v(ds1,ds2a)*s34s.v(ds3a,ds4) + s12s.v(ds1,ds2a)*a34a.v(ds3a,ds4)
			   - a13a.v(ds1,ds3a)*s24s.v(ds2a,ds4) - s13s.v(ds1,ds3a)*a24a.v(ds2a,ds4)
			   )*(-a23a.v(ds2b,ds3b)+s23s.v(ds2b,ds3b))/pDenAU
	+ normFactor*preA/MW*(
			      s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(s431a.v(ds4,ds1)-s134a.v(ds1,ds4))
			      )/pDenAU;
      
      //U-Channel Z
      //preZ = e*e/(2.0*MW*MW*SW*SW); //=prene
      //eEW-W+ all ingoing:
      //- preZ (  + (gL-gR)Me[34]<34>([12]-<12>)
      //          + 2[34]<34>(gL[231>+gR[132>)
      //          + 2MW([34]+<34>)( gL([23]<14>+[24]<13>) + gR([13]<24>+[14]<23>) )
      //        )/(s-MZ^2)
      //eW+W-E: 2<->4
      //- preZ (  + (gL-gR)Me[23]<23>([14]-<14>)
      //          + 2[23]<23>(gL[431>+gR[134>)
      //          + 2MW([23]+<23>)( gL([34]<12>+[24]<13>) + gR([13]<24>+[12]<34>) )
      //        )/(u-MZ^2)
      //34 out:
      //- preZ (  - (gL-gR)Me[23]<23>([14]+<14>)
      //          + 2[23]<23>(gL[431>-gR[134>)
      //          + 2MW([23]-<23>)( gL([34]<12>-[24]<13>) + gR(-[13]<24>+[12]<34>) )
      //        )/(u-MZ^2)
      amplitude += - normFactor*preZ*(
				      - (gL-gR)*me*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(s14s.v(ds1,ds4)+a14a.v(ds1,ds4))
				      + two*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(gL*s431a.v(ds4,ds1)-gR*s134a.v(ds1,ds4))
				      + two*MW*(s23s.v(ds2a,ds3a)-a23a.v(ds2a,ds3a))*(
										      + gL*(+s34s.v(ds3b,ds4)*a12a.v(ds1,ds2b)-s24s.v(ds2b,ds4)*a13a.v(ds1,ds3b))
										      + gR*(-s13s.v(ds1,ds3b)*a24a.v(ds2b,ds4)+s12s.v(ds1,ds2b)*a34a.v(ds3b,ds4)) )
				      )/pDenZU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eWWe::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/3=1/6
    return amp2/6.0;
  }

  
  

  //  Tests
  int test_eWWe(){
    int n=0;//Number of fails
    std::cout<<"\t* e , W+ -> W+, e       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      eWWe eWWeAmp = eWWe(EE,me,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.774986712185618E-03,5.707323688938847E-03,9.376137602131287E-03,1.402224657738295E-02,1.997692853006095E-02,2.770316077799020E-02,3.785979816597120E-02,5.140421080475444E-02,6.976163800329019E-02,9.511490661198045E-02,1.309214644254101E-01,1.828835510150093E-01,2.608822684789542E-01,3.831320260011523E-01,5.859896108123273E-01,9.501808443176387E-01,1.684111200562499E+00,3.466933655556732E+00,9.678676814455841E+00,7.661696232521470E+01};
      i += eWWeAmp.test_2to2_amp2([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH);
      i += eWWeAmp.test_2to2_amp2_rotations([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH);
      i += eWWeAmp.test_2to2_amp2_boosts([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH);
      i += eWWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {6.321531803921920E-03,9.295747301137573E-03,1.303551227068612E-02,1.777960466941266E-02,2.385306554610571E-02,3.170520345386247E-02,4.196811449277572E-02,5.554924855239859E-02,7.378242307454089E-02,9.868318574421701E-02,1.333992280934185E-01,1.830458315616170E-01,2.563531875328298E-01,3.691709485654878E-01,5.526845257956349E-01,8.753305148402371E-01,1.512830348424923E+00,3.042788235522622E+00,8.459964464439661E+00,7.432588159683374E+01};
      i += eWWeAmp.test_2to2_amp2([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH2);
      i += eWWeAmp.test_2to2_amp2_rotations([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH2);
      i += eWWeAmp.test_2to2_amp2_boosts([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH2);
      i += eWWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, mh=125, MW=80.385, pspatial=95\n";
      me = 125.1;
      mh = 125;
      pspatial = 95;
      eWWeAmp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {1.263688156260893E-01,1.446136485231535E-01,1.658787564318502E-01,1.909222182527958E-01,2.207460344800112E-01,2.566958681279893E-01,3.006119406748662E-01,3.550636400093655E-01,4.237259211361128E-01,5.120056545780800E-01,6.281296576613171E-01,7.851341910051753E-01,1.004736684143462E+00,1.325476155987591E+00,1.821602474381698E+00,2.653013368653133E+00,4.223569613695125E+00,7.841543463799127E+00,2.018510330504544E+01,1.639932957973720E+02};
      i += eWWeAmp.test_2to2_amp2([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH4);
      i += eWWeAmp.test_2to2_amp2_rotations([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH4);
      i += eWWeAmp.test_2to2_amp2_boosts([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH4);
      i += eWWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH4);
      //std::cout<<"\n# me=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      me = 125;
      mh = 0.0005;
      pspatial = 125.1;
      eWWeAmp.set_masses(me,mh,MW);
      ldouble dataCH3[20] = {1.127431045021816E-01,1.291145474915429E-01,1.481452645467890E-01,1.705567901238186E-01,1.973187826410847E-01,2.297563558528123E-01,2.697151008925533E-01,3.198219547101344E-01,3.839109273357225E-01,4.677439824554086E-01,5.802856715427175E-01,7.360759812063647E-01,9.599319049943974E-01,1.297011239987770E+00,1.836581379568256E+00,2.775956335890260E+00,4.626695206908823E+00,9.089896863748725E+00,2.508966297104032E+01,2.224850780960082E+02};
      i += eWWeAmp.test_2to2_amp2([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH3);
      i += eWWeAmp.test_2to2_amp2_rotations([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH3);
      i += eWWeAmp.test_2to2_amp2_boosts([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH3);
      i += eWWeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eWWeAmp.amp2(); }, me,MW,MW,me,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}
