
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

//File:  SPINAS/SM/neZneh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/neZneh.h"

namespace spinas {

  neZneh::neZneh(const ldouble& echarge, const ldouble& massh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), mh(massh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);  
    p1=particle(0);
    p2=particle(MZ);
    p3=particle(0);
    p4=particle(mh);
    //[23], <13>
    s23s = sproduct(SQUARE,&p2,&p3);
    a12a = sproduct(ANGLE,&p1,&p2);
    //Couplings
    preZ = sqrt2*two*e*e/(4.0*CW*CW*SW*SW);
  }
  void neZneh::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p2.set_mass(MZ);
    p4.set_mass(mh);
    propZ.set_mass(MZ);
  }
  void neZneh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);    
    //[23], <13>
    s23s.update();
    a12a.update();
    //Propagator Momentum
    ldouble propTP[4];
    for(int j=0;j<4;j++)
      propTP[j] = mom1[j]-mom3[j];
    pDenT=propZ.denominator(propTP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble neZneh::amp(const int& ds2){//Double Spin
    cdouble amplitude(0,0);
    int ds2a, ds2b;
    
    //Symmetrize the Z-Boson Spin indices
    int nCombs=get_num_spin_loops(ds2);
    ldouble normFactor=get_spin_normalization(ds2);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, i);
      
      //S-Channel Z
      //preZ = 2.0*e*e*MZ*MZ/(4.0*MW*MW*SW*SW);
      //nNZh all ingoing: 
      //preZ [23] <13> /(s-MZ^2)
      //nZNh: 2<->3
      //- preZ [23] <12> /(t-MZ^2)
      //34 out:
      //- preZ [23] <12> /(t-MZ^2)
      amplitude += - normFactor*preZ*s23s.v(ds2a)*a12a.v(ds2b)/pDenT;

    }

    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble neZneh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2){
      M = amp(j2);
      amp2 += std::pow(std::abs(M),2);
    }
    //Average over initial spins 1/3
    return amp2/3.0;
  }
  



  //  Tests
  int test_neZneh(){
    int n=0;//Number of fails
    std::cout<<"\t* ne, Z  -> ne , h      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      neZneh neZnehAmp = neZneh(EE,mh,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {8.059353005451516E+00,2.258651915986207E+00,1.009733579169442E+00,5.538794937568469E-01,3.406383099614448E-01,2.251887949692742E-01,1.563377477717746E-01,1.123538527437846E-01,8.277446791600084E-02,6.207476903552580E-02,4.712311527851940E-02,3.604184349203651E-02,2.765256524651663E-02,2.118721545363439E-02,1.612894610150310E-02,1.212034846005574E-02,8.908377516802466E-03,6.310156266380950E-03,4.191104969597635E-03,2.450583822212670E-03};
      i += neZnehAmp.test_2to2_amp2([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH);
      i += neZnehAmp.test_2to2_amp2_rotations([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH);
      i += neZnehAmp.test_2to2_amp2_boosts([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH);
      i += neZnehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {1.203318325024187E+00,6.962067652514329E-01,4.455433350491694E-01,3.045904033179917E-01,2.180841365103223E-01,1.615232807441949E-01,1.227255266120302E-01,9.509396278352465E-02,7.481187568589291E-02,5.955116555947516E-02,4.782897465099734E-02,3.866613076139407E-02,3.139601146411148E-02,2.555273433164944E-02,2.080322665059452E-02,1.690458100721226E-02,1.367650825209593E-02,1.098308907148645E-02,8.720411996893886E-03,6.808029931832895E-03};
      i += neZnehAmp.test_2to2_amp2([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH2);
      i += neZnehAmp.test_2to2_amp2_rotations([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH2);
      i += neZnehAmp.test_2to2_amp2_boosts([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH2);
      i += neZnehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH2);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=125.1\n";
      mh = 0.0005;
      pspatial = 125.1;
      neZnehAmp.set_masses(mh,MW);
      ldouble dataCH3[20] = {1.400673229879967E+00,7.406714371979310E-01,4.486102438820463E-01,2.954960673322746E-01,2.059796901205218E-01,1.495033241086341E-01,1.118106108309811E-01,8.554300971463970E-02,6.659990106539156E-02,5.255440310954875E-02,4.189869012006953E-02,3.365763888158330E-02,2.717905691634746E-02,2.201406478088723E-02,1.784599556081821E-02,1.444659976144764E-02,1.164822168320758E-02,9.325636782747034E-03,7.383907862541557E-03,5.750090348538588E-03};
      i += neZnehAmp.test_2to2_amp2([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH3);
      i += neZnehAmp.test_2to2_amp2_rotations([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH3);
      i += neZnehAmp.test_2to2_amp2_boosts([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH3);
      i += neZnehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neZnehAmp.amp2(); }, 0,MZ,0,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
