
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

//File:  SPINAS/SM/hhZZ.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/hhZZ.h"

namespace spinas {

  hhZZ::hhZZ(const ldouble& echarge, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);
    proph = propagator(mh,wh);  
    p1=particle(mh);
    p2=particle(mh);
    p3=particle(MZ);
    p4=particle(MZ);
    //<34>,[34]
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    //[314>,[413>
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
    //Couplings
    pre = e*e/(2.0*MW*MW*SW*SW);
    preS = 3.0*pre*mh*mh;
  }
  void hhZZ::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(mh);
    p2.set_mass(mh);
    p3.set_mass(MZ);
    p4.set_mass(MZ);
    propZ.set_mass(MZ);
    proph.set_mass(mh);
    //Couplings
    pre = e*e/(2.0*MW*MW*SW*SW);
    preS = 3.0*pre*mh*mh;
  }
  void hhZZ::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<34>,[34]
    s34s.update();
    a34a.update();
    //[314>,[413>
    s314a.update();
    s413a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=proph.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=propZ.denominator(propTP);
    pDenU=propZ.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble hhZZ::amp(const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1;
    int ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);

      //pre = e*e/(2.0*MW*MW*SW*SW);
      
      //4-Point      
      amplitude += - normFactor*pre*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b);
          
      //S-Channel h
      //preS = 3.0*pre*mh*mh;
      //all ingoing: 
      //preS [34] <34> /(s-Mh^2)
      //34 outgoing:
      //preS [34] <34> /(s-Mh^2)
      amplitude += normFactor*preS*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)/pDenS;
      
      //T-Channel Z
      //all ingoing:
      //+pre ( 2MZ^2[34]<34> + MZ([34][314>+<34>[413>) + [314>[413> )/(t-MZ^2)
      //34 outgoing:
      //+pre ( 2MZ^2[34]<34> - MZ([34][314>+<34>[413>) + [314>[413> )/(t-MZ^2)
      amplitude += normFactor*pre*2.0*MZ*MZ*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)/pDenT
      	-          normFactor*pre*MZ*(s34s.v(ds3a,ds4a)*s314a.v(ds3b,ds4b)+a34a.v(ds3a,ds4a)*s413a.v(ds4b,ds3b))/pDenT
      	+          normFactor*pre*s314a.v(ds3a,ds4a)*s413a.v(ds4b,ds3b)/pDenT;
      
      //U-Channel Z
      //all ingoing:
      //+pre ( 2MZ^2[34]<34> - MZ([34][413>+<34>[314>) + [314>[413> )/(u-MZ^2)
      //34 outgoing:
      //+pre ( 2MZ^2[34]<34> + MZ([34][413>+<34>[314>) + [314>[413> )/(u-MZ^2)
      amplitude += normFactor*pre*2.0*MZ*MZ*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)/pDenU
      	+          normFactor*pre*MZ*(s34s.v(ds3a,ds4a)*s413a.v(ds4b,ds3b)+a34a.v(ds3a,ds4a)*s314a.v(ds3b,ds4b))/pDenU
      	+          normFactor*pre*s314a.v(ds3a,ds4a)*s413a.v(ds4b,ds3b)/pDenU;
      


    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble hhZZ::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=2)
      for(int j4=-2;j4<=2;j4+=2){
	M = amp(j3,j4);
	amp2 += std::pow(std::abs(M),2);
      }
    //Symmetry Factor 1/2
    return amp2/2.0;
  }

   


  //  Tests
  int test_hhZZ(){
    int n=0;//Number of fails
    std::cout<<"\t* h , h  -> Z , Z       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      hhZZ hhZZAmp = hhZZ(EE,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.393489509525684E+01,3.748747753024061E+00,1.680617706172940E+00,9.556604827919933E-01,6.296050982616538E-01,4.610506302842268E-01,3.667247593592984E-01,3.122816535253189E-01,2.819664320572433E-01,2.683276751906571E-01,2.683276751906592E-01,2.819664320572405E-01,3.122816535253176E-01,3.667247593592987E-01,4.610506302842268E-01,6.296050982616561E-01,9.556604827919951E-01,1.680617706172940E+00,3.748747753024051E+00,1.393489509525681E+01};
      i += hhZZAmp.test_2to2_amp2([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH);
      i += hhZZAmp.test_2to2_amp2_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH);
      i += hhZZAmp.test_2to2_amp2_boosts([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH);
      i += hhZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {3.442907027433093E+00,1.846457347110825E+00,1.155021255182787E+00,8.017863275666668E-01,6.020032123716252E-01,4.816076940083718E-01,4.066065529046288E-01,3.599365309079385E-01,3.326084811332237E-01,3.199509683595942E-01,3.199509683595944E-01,3.326084811332232E-01,3.599365309079381E-01,4.066065529046294E-01,4.816076940083717E-01,6.020032123716252E-01,8.017863275666669E-01,1.155021255182787E+00,1.846457347110824E+00,3.442907027433089E+00};
      i += hhZZAmp.test_2to2_amp2([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH2);
      i += hhZZAmp.test_2to2_amp2_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH2);
      i += hhZZAmp.test_2to2_amp2_boosts([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH2);
      i += hhZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {5.950942712701554E-01,5.950549365246995E-01,5.950199749315097E-01,5.949893855639474E-01,5.949631676113112E-01,5.949413203787924E-01,5.949238432874399E-01,5.949107358741272E-01,5.949019977915340E-01,5.948976288081197E-01,5.948976288081197E-01,5.949019977915335E-01,5.949107358741278E-01,5.949238432874399E-01,5.949413203787925E-01,5.949631676113113E-01,5.949893855639475E-01,5.950199749315095E-01,5.950549365247001E-01,5.950942712701552E-01};
      i += hhZZAmp.test_2to2_amp2([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH4);
      i += hhZZAmp.test_2to2_amp2_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH4);
      i += hhZZAmp.test_2to2_amp2_boosts([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH4);
      i += hhZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH4);
      //std::cout<<"\n# mh=5, MW=80.385, pspatial=95\n";
      mh = 5;
      pspatial = 95;
      hhZZAmp.set_masses(mh,MW);
      ldouble dataCH6[20] = {5.093548809218734E-01,4.855439406413581E-01,4.654168186218827E-01,4.485623621152399E-01,4.346506064535563E-01,4.234190843403733E-01,4.146624442968827E-01,4.082246568396835E-01,4.039932929068483E-01,4.018955116627572E-01,4.018955116627572E-01,4.039932929068482E-01,4.082246568396835E-01,4.146624442968826E-01,4.234190843403733E-01,4.346506064535564E-01,4.485623621152399E-01,4.654168186218827E-01,4.855439406413580E-01,5.093548809218734E-01};
      i += hhZZAmp.test_2to2_amp2([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH6);
      i += hhZZAmp.test_2to2_amp2_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH6);
      i += hhZZAmp.test_2to2_amp2_boosts([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH6);
      i += hhZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH6);
      //std::cout<<"\n# mh=5, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH7[20] = {1.193768412751491E+01,3.705446901862067E+00,1.825793253720890E+00,1.119703469075197E+00,7.846486032157109E-01,6.034717066197560E-01,4.981550271148607E-01,4.354616811118717E-01,3.997403965086873E-01,3.834389792361504E-01,3.834389792361498E-01,3.997403965086879E-01,4.354616811118718E-01,4.981550271148620E-01,6.034717066197570E-01,7.846486032157105E-01,1.119703469075198E+00,1.825793253720889E+00,3.705446901862067E+00,1.193768412751489E+01};
      i += hhZZAmp.test_2to2_amp2([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH7);
      i += hhZZAmp.test_2to2_amp2_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH7);
      i += hhZZAmp.test_2to2_amp2_boosts([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH7);
      i += hhZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH7);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=95\n";
      mh = 0.0005;
      pspatial = 95;
      hhZZAmp.set_masses(mh,MW);
      ldouble dataCH3[20] = {5.071620540594148E-01,4.842315882628248E-01,4.648153627526531E-01,4.485320177637508E-01,4.350746267661532E-01,4.241984494399119E-01,4.157116085270992E-01,4.094680671948373E-01,4.053624593513622E-01,4.033264565875229E-01,4.033264565875229E-01,4.053624593513622E-01,4.094680671948372E-01,4.157116085270992E-01,4.241984494399119E-01,4.350746267661532E-01,4.485320177637508E-01,4.648153627526531E-01,4.842315882628248E-01,5.071620540594148E-01};
      i += hhZZAmp.test_2to2_amp2([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH3);
      i += hhZZAmp.test_2to2_amp2_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH3);
      i += hhZZAmp.test_2to2_amp2_boosts([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH3);
      i += hhZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH3);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH5[20] = {1.193409290785243E+01,3.705593394913418E+00,1.826192415661682E+00,1.120082043433994E+00,7.849799487574084E-01,6.037621176673770E-01,4.984146424642518E-01,4.356997699436336E-01,3.999649213496139E-01,3.836569573826608E-01,3.836569573826592E-01,3.999649213496142E-01,4.356997699436334E-01,4.984146424642531E-01,6.037621176673761E-01,7.849799487574086E-01,1.120082043433991E+00,1.826192415661682E+00,3.705593394913415E+00,1.193409290785241E+01};
      i += hhZZAmp.test_2to2_amp2([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH5);
      i += hhZZAmp.test_2to2_amp2_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH5);
      i += hhZZAmp.test_2to2_amp2_boosts([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH5);
      i += hhZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH5);
      //std::cout<<"\n# mh=0, MW=80.385, pspatial=95\n";
      mh = 0;
      pspatial = 95;
      hhZZAmp.set_masses(mh,MW);
      ldouble dataCH8[20] = {5.071620540375101E-01,4.842315882497074E-01,4.648153627466421E-01,4.485320177634573E-01,4.350746267704141E-01,4.241984494477386E-01,4.157116085376356E-01,4.094680672073261E-01,4.053624593651157E-01,4.033264566018978E-01,4.033264566018978E-01,4.053624593651157E-01,4.094680672073260E-01,4.157116085376356E-01,4.241984494477385E-01,4.350746267704140E-01,4.485320177634572E-01,4.648153627466421E-01,4.842315882497073E-01,5.071620540375101E-01};
      i += hhZZAmp.test_2to2_amp2([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH8);
      i += hhZZAmp.test_2to2_amp2_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH8);
      i += hhZZAmp.test_2to2_amp2_boosts([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH8);
      i += hhZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH8);
      //std::cout<<"\n# mh=0, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH9[20] = {1.193409290781652E+01,3.705593394914884E+00,1.826192415665676E+00,1.120082043437780E+00,7.849799487607214E-01,6.037621176702803E-01,4.984146424668490E-01,4.356997699460129E-01,3.999649213518584E-01,3.836569573848385E-01,3.836569573848381E-01,3.999649213518577E-01,4.356997699460120E-01,4.984146424668488E-01,6.037621176702790E-01,7.849799487607213E-01,1.120082043437780E+00,1.826192415665674E+00,3.705593394914882E+00,1.193409290781650E+01};
      i += hhZZAmp.test_2to2_amp2([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH9);
      i += hhZZAmp.test_2to2_amp2_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH9);
      i += hhZZAmp.test_2to2_amp2_boosts([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH9);
      i += hhZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhZZAmp.amp2(); }, mh,mh,MZ,MZ,pspatial,dataCH9);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  

}
